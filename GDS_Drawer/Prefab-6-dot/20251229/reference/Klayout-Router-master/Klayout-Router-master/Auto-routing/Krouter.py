from Map import Map
from PinMatcher import PinMatcher

import numpy as np
import klayout.db as kdb

from skimage.graph import route_through_array

from pathlib import Path
from typing import TypedDict, Generator

class ObsRule(TypedDict):
    """Specified type hint for obstacle rules."""
    layers: str | list[str]  # Layer(s) can be a string or a list of strings
    bbx: bool               # Whether to use bounding box
    pad: float              # Padding size in microns


class Krouter(object):
    """
    Krouter is a layout-oriented (or user-oriented) class that performs auto-routing process.
        - This is a single-area routing class, considering only one total area bounding box where all obstacles and all pins are inside.
        - This class handles all the interaction with the layout (layer).
        - It extracts shapes or regions from the layout (layer) for underlying classes (Map, PinMatcher...) who
          only manage abstract concept to perform routing.
        - It takes length quantities in micron unit, and converts them to layout's dbu unit.
        - It feeds the Map class with
            - layout
            - total area bounding box, resolution
            - obstacle-regions
        - It feeds the PinMatcher class with 
            - Shapes collected from the layout (layer).
        - Cost of i-th pair of pins is 1000 + i, e.g. 0-th pair of pins is 1000. Compatible with less that 1000 pairs of pins.
        - Cost of path connecting i-th pair of pins is 2000 + i, e.g. 0-th path is 2000.
        - This cost does not have the original meaning of cost, but is used to distinguish different pairs of pins and paths from obstacles and plain pixels.
    """

    def __init__(self,
                 file_path: str,
                 cell_name: str,
                 total_area_bbox_layer: str, 
                 Pin_A_layer: str, Pin_B_layer: str, 
                 map_resolution: float, 
                 obs_rules: list[ObsRule],
                 ) -> None:
        """
        Args:
            total_area_bbox_layer (str): Layer containing the total area bounding box.
            Pin_A_layer (str): Layer containing pins of set A.
            Pin_B_layer (str): Layer containing pins of set B.
            map_resolution (float): Resolution of the map, in microns.
            obs_rules (list[ObsRule]): List of obstacle rules, each rule is a dictionary with keys:
                - "layers": Layer(s) containing obstacles, can be a string or a list of strings.
                - "bbx": Whether to use bounding box for the shapes.
                - "pad": Padding size for the shapes in microns.
        """
        # Empty layout
        self.__layout = kdb.Layout()

        # Get the layout and cell from the file
        self.__file_path = Path(file_path)
        self.__layout.read(self.__file_path)
        self.__cell = self.__layout.cell(cell_name)

        # Back up source file
        self.__layout.write(self.__file_path.with_stem(self.__file_path.stem + "_backup"))

        # --- Useful variables ---
        self.__obs_rules = [obs_rules] if isinstance(obs_rules, dict) else obs_rules  # Ensure obs_rules is a list

        # --- Map ---
        # Initialize the map
        self.__Map = Map(self.__layout, 
                         kdb.Region(self.__cell.shapes(self.try_find_layer(total_area_bbox_layer))).bbox(), 
                         self.to_dbu(map_resolution))
        # Generate total obs region based on obs_rules
        self.__obs_region = self.generate_obs_region_by_obs_rules()

        # --- Pin matcher ---
        # Initialize the pin matcher
        self.__PinMatcher = PinMatcher()
        # Get pins and feed them to the pin matcher
        self.__PinMatcher.get_pins_from_shapes(self.__cell.shapes(self.try_find_layer(Pin_A_layer)))
        self.__PinMatcher.get_pins_from_shapes(self.__cell.shapes(self.try_find_layer(Pin_B_layer)))

    # --- General methods ---
    def to_dbu(self, value: float) -> int:
        """Convert a value in micron to dbu unit.
        Args:
            value (float): The value in micron.
        Returns:
            int: The value in dbu unit, rounded to the nearest integer.
        """
        return int(value / self.__layout.dbu + 0.5)
    
    def try_find_layer(self, layer_str: str) -> int:
        """Find layer index of a layer by its string representation. 
        
        Raises ValueError to warn, compared to kdb.Layout.find_layer()
        Args:
            layer_str (str): Layer string like "layer_number/datatype", e.g. "1/0"; name, e.g. "test"; or both, e.g. "1/0, test".
        Returns:
            int: The index of the layer in the layout.
        Raises:
            ValueError: If the layer is not found in the layout.
        """
        layer_index = self.layout.find_layer(kdb.LayerInfo.from_string(layer_str))
        if layer_index is None:
            raise ValueError(f"Layer {layer_str} not found in the layout.")
        return layer_index

    # --- Map related methods ---
    def bbx_and_pad_region_on_layers(self, cbpr_layers: str | list[str], bbx: bool, pad: int) -> kdb.Region:
        """Collect all shapes from the cbpr_layer into a single region.
        
        You can choose whether to bounding box them or to pad them or both (pad after bbx).

        Args:
            cbpr_layer (list[str] or str): A list of layers of shapes like ["1/0", "2/0", ...].
                A single layer str is also accepted, e.g. "1/0".
            bbx (bool, optional): Whether to turn shapes into their bounding boxes before padding.
            pad (int, optional): Padding size of shapes, the unit is dbu.

        Returns:
            kdb.Region: The region collected.
        """
        # Check if cbpr_layers is a string or a list
        if isinstance(cbpr_layers, str): cbpr_layers = [cbpr_layers]
        
        # Collect all shapes from the layer into a region
        cbpr_region = kdb.Region()
        for clct_layer in cbpr_layers:
            cbpr_region += kdb.Region(self.__cell.shapes(self.try_find_layer(clct_layer)))
            print(f"Krouter: Collecting region from layer {clct_layer}.")

        # Pad the every polygon or pad each of their bounding boxes
        # Region.inside(box) ensures the final polygons are within the total area bounding box without merging.
        if bbx:
            # region.extents(pad) pads the bounding box of each polygon
            return cbpr_region.extents(pad)
        else:
            # region.sized(pad) pads each polygon
            return cbpr_region.sized(pad)

    def generate_obs_region_by_obs_rules(self) -> kdb.Region:
        """Generate the total obstacle region following the obs_rules."""
        # Empty region to contain obstacle-region
        obs_region = kdb.Region()

        # Collect all obstacle regions from the rules
        for obs_rule in self.__obs_rules:
            # Convert padding size to dbu
            obs_rule["pad"] = self.to_dbu(obs_rule["pad"])
            # Collect obs region
            obs_region.join_with(self.bbx_and_pad_region_on_layers(*obs_rule.values()))
        
        return obs_region

    # --- PinMatcher related methods ---
    def distance_matrix(self) -> np.ndarray:
        # initialize distance matrix
        distance_matrix = np.zeros((len(self.__PinMatcher.pin_sets["A"]), len(self.__PinMatcher.pin_sets["B"])), dtype=int)
        # Generate MCP_Geometric
        mcpg = self.__Map.MCPG
        # Turn pins to map coordinates (map index)
        pins_A = [self.__Map.layout_coor_to_map_coor(pin.x, pin.y) for pin in self.__PinMatcher.pin_sets["A"]]
        pins_B = [self.__Map.layout_coor_to_map_coor(pin.x, pin.y) for pin in self.__PinMatcher.pin_sets["B"]]

        # Use MCP_Geometric to find costs between every pair of pins. Read skimage.graph.MCP_Geometric for more details.
        for i, pin_a in enumerate(pins_A):
            # First, get costs from pin_a to all pins_B
            costs = mcpg.find_costs([pin_a], pins_B)[0]
            # Then, record cost between every pair
            for j, pin_b in enumerate(pins_B):
                distance_matrix[i, j] = costs[pin_b]

        return distance_matrix
    
    # --- routing animation methods ---
    def prepare_map_for_routing(self, obs_safe_distance: float, obs_hardness: int, obs_damping_step: int,
                                Pin_A_safe_distance: float, Pin_B_safe_distance: float, pin_hardness: int, 
                                Pin_A_damping_step: int, Pin_B_damping_step: int,
                                ) -> Generator[np.ndarray, None, None]:
        """Prepare the map with soft safe distance on obs and pins for routing, and generate frames for animation.

        If you do not need animation, use run_prepare_map_for_routing().
        
        Args:
            obs_safe_distance (float): Safe distance around obstacles in microns. 
                Not necessarily too large, this is just to avoid path-not-found issue.
            obs_hardness (int): Cost of the pixels around obstacles, the higher the harder for paths to get close, typically 10-100.
            obs_damping_step (int): Number of steps to gradually increase cost around obstacles within the obs_safe_distance region. Recommended 10.
            Pin_A_safe_distance (float): Safe distance for pins in set A in microns.
            Pin_B_safe_distance (float): Safe distance for pins in set B in microns.
            pin_hardness (int): Cost of the pixels around pins, the higher the harder for paths to get close, typically 10-100.
            pin_damping_step (int): Number of steps to gradually increase cost around pins within the pin_safe_distance region. Recommended 10.
        
        Yields:
            np.ndarray: The map at different stages of the preparation process for animation.
        """
        self.__obs_safe_distance = self.to_dbu(obs_safe_distance)
        self.__obs_hardness = obs_hardness
        self.__Pin_A_safe_distance = self.to_dbu(Pin_A_safe_distance)
        self.__Pin_B_safe_distance = self.to_dbu(Pin_B_safe_distance)
        self.__pin_hardness = pin_hardness
        self.__Pin_A_damping_step = Pin_A_damping_step
        self.__Pin_B_damping_step = Pin_B_damping_step

        # Yield the initial empty map twice to compensate the first frame. A bug of FuncAnimation?
        yield self.__Map.map
        yield self.__Map.map

        # Rasterize the obstacles to -1 in the map
        self.__Map.rasterize_region(self.__obs_region, -1)
        yield self.__Map.map

        # Rasterize the padded obstacles to obs_hardness in the map on non-negative pixels
        m0, n0, m1, n1, obs_raster = self.__get_damping_raster(self.__obs_region.merged(), self.__obs_safe_distance,
                                                                self.__obs_hardness, obs_damping_step)
        self.__Map.overwrite(obs_raster, obs_raster > 0, lambda x: x >= 0, m0, n0, m1, n1)
        yield self.__Map.map

        # Rasterize the pins to -2 in the map
        for pin_region_a, pin_region_b in self.__PinMatcher.matched_polygons:
            self.__Map.rasterize_region(pin_region_a, -2, lambda x: True, merge=True)
            self.__Map.rasterize_region(pin_region_b, -2, lambda x: True, merge=True)
        yield self.__Map.map

        #  Safely-size pins region for later recover use in the routing
        self.__safely_sized_PINs_region = self.__PinMatcher.matched_region_A.sized(self.__Pin_A_safe_distance) + \
                                          self.__PinMatcher.matched_region_B.sized(self.__Pin_B_safe_distance)
        self.__safely_sized_PINs_region.merged_semantics = False
        
        # Rasterize damping safely sized Pin A and B
        m0, n0, m1, n1, Pin_A_raster = self.__get_damping_raster(self.__PinMatcher.matched_region_A, self.__Pin_A_safe_distance,
                                                                 self.__pin_hardness, Pin_A_damping_step)
        self.__Map.overwrite(Pin_A_raster, Pin_A_raster > 0, lambda x: (0 < x) & (x < Pin_A_raster), m0, n0, m1, n1)  # This is called supplemental overwriting

        m0, n0, m1, n1, pin_B_raster = self.__get_damping_raster(self.__PinMatcher.matched_region_B, self.__Pin_B_safe_distance,
                                                                 self.__pin_hardness, Pin_B_damping_step)
        self.__Map.overwrite(pin_B_raster, pin_B_raster > 0, lambda x: (0 < x) & (x < pin_B_raster), m0, n0, m1, n1)
        yield self.__Map.map

        self.__prepared_map = self.__Map.map.copy()  # Backup the prepared map

        print("Krouter: Map prepared for routing.")

    def route(self, routing_path_width: float, path_safe_distance: float, path_hardness: int, 
              path_damping_step: int = 10) -> Generator[np.ndarray, None, None]:
        """Route on the prepared map and generate frames for animation.

        If you do not need animation, use run_route().

        Args:
            routing_path_width (float): Width of the routing paths, in microns.
            path_safe_distance (float): soft safe distance between the routing paths, in microns.
            path_hardness (int): Cost of the pixels around paths, the higher the harder for paths to get close, typically 10-100.
            path_damping_step (int): Number of steps to gradually increase cost around path within the path_safe_distance region. Default is 10.
                The larger the smoother the cost increase, but you also need to think about the resolution limit.
        Yields:
            np.ndarray: The map at different stages of the routing process for animation.
        """
        self.__routing_path_width = self.to_dbu(routing_path_width)  # Convert routing path width to dbu
        self.__path_safe_distance = self.to_dbu(path_safe_distance)
        self.__path_hardness = path_hardness

        # Empty region to collect paths as polygons
        Path_collection_region = kdb.Region()
        # Empty list to collect paths as kdb.Path
        Path_list = []

        # Take the i-th pair of pins
        for i, (pin_a, pin_b) in enumerate(self.__PinMatcher.matched_pins):
            # 1. Recover i-th pins themselves to 1 in the map
            ith_pins_region = self.__PinMatcher.matched_polygons[i][0] + self.__PinMatcher.matched_polygons[i][1]
            self.__Map.rasterize_region(ith_pins_region, 1, lambda x: True, merge=True)
            yield self.__Map.map

            # # 2. Clear hard pixels around routing pins. This is to make sure high performance of Dijkstra's algorithm.
            # # prepare region to be cleared, which should avoid all other sized pins and sized obstacles
            # sized_ith_pins_region = self.__PinMatcher.matched_polygons[i][0].sized(self.__Pin_A_safe_distance) + \
            #                         self.__PinMatcher.matched_polygons[i][1].sized(self.__Pin_B_safe_distance)
            # temp_region = sized_ith_pins_region - self.__safely_sized_PINs_region.not_in(sized_ith_pins_region)
            # # temp_region -= self.__obs_region.sized(self.__obs_safe_distance)
            # # clear hard pixels around the i-th pair of pins
            # self.__Map.rasterize_region(temp_region, 1, lambda x: x > 1, merge=True)
            # yield self.__Map.map

            # 3.1 routing between pin_a and pin_b, get the path nodes (list of (m, n) tuples)
            try:
                path_nodes = route_through_array(self.__Map.map,
                                                start=self.__Map.layout_coor_to_map_coor(pin_a.x, pin_a.y),
                                                end=self.__Map.layout_coor_to_map_coor(pin_b.x, pin_b.y),
                                                fully_connected=True,
                                                geometric=True,
                                                )[0]
                print(f"Krouter: Routing path {i + 1} between pin {pin_a} and pin {pin_b}.")
                found_path = True
            except Exception as e:
                print(f"Krouter: Failed to route path {i + 1} between pin {pin_a} and pin {pin_b}. Exception: {e}")
                found_path = False

            if found_path:
                # 3.2 Insert kdb.Path into the layout only if the path is found
                Path = self.generate_path_from_nodes(path_nodes, self.__routing_path_width, compress=True) # Generate kdb.Path from path nodes
                Path_collection_region += Path  # Add the path to the paths region
                Path_list.append(Path)  # Add the path to the paths list

                # 4. set path to -3 on pixels
                self.__Map.rasterize_region(kdb.Region(Path), -3, lambda x: True, merge=True)
                yield self.__Map.map

                # 5. set safely sized path to path_hardness on 1-pixels
                # self.__Map.rasterize_region(kdb.Region(Path).sized(self.__path_safe_distance), 
                #                             self.__path_hardness, lambda x: x==1, merge=True)
                m0, n0, m1, n1, path_raster = self.__get_damping_raster(kdb.Region(Path), self.__path_safe_distance, 
                                                                        self.__path_hardness, path_damping_step)
                self.__Map.overwrite(path_raster, path_raster > 0, lambda x: (0 < x) & (x < path_raster), m0, n0, m1, n1)

                yield self.__Map.map

            # 6. set i-th pair of pins themselves back to -2 on >=0-pixels
            self.__Map.rasterize_region(ith_pins_region, -2, lambda x: x >= 0, merge=True)
            yield self.__Map.map

            # # 7. set safely sized i-th pair of pins back to pin_hardness on 1-pixels
            # m0, n0, m1, n1, ith_pin_a_raster = self.__get_damping_raster(self.__PinMatcher.matched_polygon_A[i], self.__Pin_A_safe_distance,
            #                                                          self.__pin_hardness, self.__Pin_A_damping_step)
            # self.__Map.overwrite(ith_pin_a_raster, ith_pin_a_raster > 0, lambda x: (0 < x) & (x < ith_pin_a_raster), m0, n0, m1, n1)
            # m0, n0, m1, n1, ith_pin_b_raster = self.__get_damping_raster(self.__PinMatcher.matched_polygon_B[i], self.__Pin_B_safe_distance,
            #                                                          self.__pin_hardness, self.__Pin_B_damping_step)
            # self.__Map.overwrite(ith_pin_b_raster, ith_pin_b_raster > 0, lambda x: (0 < x) & (x < ith_pin_b_raster), m0, n0, m1, n1)

            yield self.__Map.map

        self.__Path_collection_region = Path_collection_region
        self.__Path_list = Path_list

        print(f"Krouter: Found {Path_collection_region.count()} paths.")

    def adapt_path_density(self, path_density_hardness: int) -> Generator[np.ndarray, None, None]:
        """Adapt the path density according to Path_collection_region. Region with more paths overlapping will be harder on the map.

        Args:
            path_density_hardness (int): Cost of the pixels of overlapped paths.
        """
        # Make sure Path_collection_region is not empty
        if self.__Path_collection_region.is_empty():
            # yield self.__Map.map
            return

        # Rasterize the safely sized path collection
        path_density_map = np.array((self.__Path_collection_region.sized(self.__path_safe_distance))\
                                    .rasterize(self.__Map.total_area_bbox.p1, 
                                               kdb.Vector(self.__Map.resolution, self.__Map.resolution), 
                                               self.__Map.N, self.__Map.M))
        
        # Normalize the path density map, where the cost means how many paths overlapping on the pixel, e.g. 1 means 2 overlapping paths, 2 means 3.
        path_density_map = path_density_map // (self.__Map.resolution * self.__Map.resolution) - 1

        # Turn density to hardness on the map
        path_density_map[path_density_map > 0] *= path_density_hardness

        # Overwrite 1-pixels on the map with path density map where pixels > 0
        self.__Map.overwrite(path_density_map, (path_density_map > 0) & (self.__Map.map == 1))
        print("Krouter: Adapted path density on the map.")
        yield self.__Map.map

    def self_adaptive_route(self, obs_safe_distance: float, obs_hardness: int, obs_damping_step: int,
                            Pin_A_safe_distance: float, Pin_B_safe_distance: float, pin_hardness: int,
                            Pin_A_damping_step: int, Pin_B_damping_step: int,
                            routing_path_width: float, path_safe_distance: float, path_hardness: int, path_damping_step: int,
                            path_density_hardness: int,
                            round: int = 1,
                            ) -> Generator[np.ndarray, None, None]:
        """A generator that yields frames of the self-adaptive routing process for animation.
        Each yield is a 2D numpy array representing the map at a certain stage of the routing process.
        """
        for r in range(round):
            print(f"Krouter: Self-adaptive routing round {r+1}.")
            if r == 0:
                # Prepare the map for the first time
                yield from self.prepare_map_for_routing(obs_safe_distance, obs_hardness, obs_damping_step,
                                                        Pin_A_safe_distance, Pin_B_safe_distance, pin_hardness,
                                                        Pin_A_damping_step, Pin_B_damping_step)
                self.__Path_collection_region = kdb.Region()  # Empty region for 1st round adapt_path_density()
            else:
                self.__Map.overwrite(self.__prepared_map, self.__Map.map != self.__prepared_map)  # Recover the prepared map
                yield self.__Map.map

            # Adapt path density and yield frames during the process
            yield from self.adapt_path_density(path_density_hardness)
            # Perform routing and yield frames during the process
            yield from self.route(routing_path_width, path_safe_distance, path_hardness, path_damping_step)
        
    # --- sub-methods for routing --- 
    def generate_path_from_nodes(self, path_nodes: list[tuple[int, int]], path_width: int,
                                 compress: bool = True) -> kdb.Path:
        """Generate a kdb.Path from a list of (m, n) tuples representing path nodes in map coordinates.

        Args:
            path_nodes (list[tuple[int, int]]): List of (m, n) tuples representing path nodes in map coordinates.
            compress (bool, optional): Whether to compress the path to keep only essential points. Defaults to True.

        Returns:
            kdb.Path: The generated kdb.Path object.
        """
        # First convert nodes to kdb.Points
        path_points = [kdb.Point(*self.__Map.map_coor_to_layout_coor(m, n)) for m, n in path_nodes]
        # Then create kdb.Path
        Path = kdb.Path(path_points, 
                        width=path_width, 
                        bgn_ext=path_width/2, end_ext=path_width/2, 
                        round=True,
                        )
        # Compress the path to keep only essential points
        if compress: self.compress_path(Path)

        return Path

    def compress_path(self, path: kdb.Path) -> kdb.Path:
        """Keep only the essential points (start, end, inflection) in a path, removing redundant points.

        This operation is in-situ, meaning it modifies the path directly, and returns its reference.

        Args:
            path (kdb.Path): The path to be compressed.

        Returns:
            kdb.Path: The compressed path with only essential points.
        """
        # If the path is too short to have any inflections, return it as is.
        if path.num_points() < 3:
            return path

        # The compressed path always starts with the first point.
        points = list(path.each_point())
        essential_points = [points[0]]

        # Iterate through the path, looking at three points at a time:
        # prev_point, current_point, and next_point.
        for i in range(1, len(points) - 1):
            prev_point = points[i-1]
            current_point = points[i]
            next_point = points[i+1]

            # Calculate the direction vector from the previous point to the current one.
            # For grid coordinates, the vector is just the difference in coordinates.
            vector1 = current_point - prev_point
            
            # Calculate the direction vector from the current point to the next one.
            vector2 = next_point - current_point

            # If the direction changes, the current point is a corner and must be kept.
            if vector1/vector1.abs() != vector2/vector2.abs():
                essential_points.append(current_point)

        # The essential_points path always ends with the last point.
        essential_points.append(points[-1])

        # Reassign the points of the path to the essential points.
        path.points = essential_points

        return path

    def __get_damping_raster(self, region: kdb.Region, safe_distance: int, hardness: int, step: int) -> tuple[int, int, int, int, np.ndarray]:
        """Get a raster with gradually damping hardness from the boundary of the region to the outside.

        Args:
            region (kdb.Region): The region to be rasterized.
            safe_distance (int): The safe distance around the region, in dbu.
            hardness (int): The maximum hardness around the region.
            step (int): The number of steps to gradually decrease hardness.
        Returns:
            tuple[int, int, int, int, np.ndarray]: The boundary (m0, n0, m1, n1) and the raster of the region with gradually damping hardness.
        """
        # Copy the region to avoid modifying the original one
        safely_sized_region = region.dup()

        # Safely size the region
        for i in range(step):
            safely_sized_region += region.sized(safe_distance * (i+1) / step)

        # Get the raster of the safely sized region
        m0, n0, m1, n1, damping_raster = self.__Map.get_boundary_and_raster_of_region(safely_sized_region)

        # Normalize the raster to gradually damping hardness
        damping_raster = (damping_raster // (self.__Map.resolution * self.__Map.resolution)) * (hardness // step)

        return m0, n0, m1, n1, damping_raster


    def run_prepare_map_for_routing(self, obs_safe_distance: float, obs_hardness: int, obs_damping_step: int,
                                    Pin_A_safe_distance: float, Pin_B_safe_distance: float, pin_hardness: int,
                                    Pin_A_damping_step: int, Pin_B_damping_step: int,
                                    ) -> None:
        """Deplete the generator to prepare the map for routing."""
        for _ in self.prepare_map_for_routing(obs_safe_distance, obs_hardness, obs_damping_step,
                                                Pin_A_safe_distance, Pin_B_safe_distance, pin_hardness,
                                                Pin_A_damping_step, Pin_B_damping_step):
            pass
    
    def run_route(self, routing_path_width: float, 
                        path_safe_distance: float, path_hardness: int, path_damping_step: int = 10
                    ) -> None:
        """Deplete the generator to perform routing."""
        for _ in self.route(routing_path_width, path_safe_distance, path_hardness, path_damping_step):
            pass

    def run_self_adaptive_route(self, obs_safe_distance: float, obs_hardness: int, obs_damping_step: int,
                                Pin_A_safe_distance: float, Pin_B_safe_distance: float, pin_hardness: int,
                                Pin_A_damping_step: int, Pin_B_damping_step: int,
                                routing_path_width: float, path_safe_distance: float, path_hardness: int, path_damping_step: int,
                                path_density_hardness: int,
                                round: int = 1,
                                ) -> None:
        """Deplete the generator to perform self-adaptive routing."""
        for _ in self.self_adaptive_route(obs_safe_distance, obs_hardness, obs_damping_step,
                                         Pin_A_safe_distance, Pin_B_safe_distance, pin_hardness,
                                         Pin_A_damping_step, Pin_B_damping_step,
                                         routing_path_width, path_safe_distance, path_hardness, path_damping_step,
                                         path_density_hardness,
                                         round):
            pass

    # --- Utility methods ---
    def visualize(self) -> None:
        """Visualize Map and matched pins in the layout."""
        # Visualize the map
        self.__Map.visualize(self.__layout, self.__cell)
        # Visualize the matched pins
        # self.__PinMatcher.visualize(self.__layout, self.__cell)
        return

    def save_layout(self, file_path: str = None) -> None:
        """Save the current layout."""
        file_path = self.__file_path if file_path is None else file_path
        self.__layout.write(file_path)
        print(f"Krouter: Layout saved to {file_path}.")

    def insert_paths_to_layout(self, routing_path_layer: str) -> None:
        """Insert the generated paths to the layout on the specified layer."""
        # Make sure there are paths to insert
        if len(self.__Path_list) == 0:
            print("Krouter: No paths to insert.")
            return
        
        # layer to draw the paths
        self.__routing_path_layer = routing_path_layer  
        # Insert the path on the routing path layer
        for path in self.__Path_list:
            self.__cell.shapes(self.__layout.layer(kdb.LayerInfo.from_string(
                                self.__routing_path_layer))).insert(path)
        print(f"Krouter: Inserted {len(self.__Path_list)} paths to layer {self.__routing_path_layer}.")

    @property
    def layout(self) -> kdb.Layout:
        """Get the current layout."""
        return self.__layout
    @property
    def cell(self) -> kdb.Cell:
        """Get the current cell."""
        return self.__cell
    @property
    def Map(self) -> Map:
        """Get the current map."""
        return self.__Map
    @property
    def PinMatcher(self) -> PinMatcher:
        """Get the current pin matcher."""
        return self.__PinMatcher
    @property
    def Path_collection_region(self) -> kdb.Region:
        """Get the current path collection."""
        return self.__Path_collection_region
    @property
    def Path_list(self) -> list[kdb.Path]:
        """Get the current path list."""
        return self.__Path_list