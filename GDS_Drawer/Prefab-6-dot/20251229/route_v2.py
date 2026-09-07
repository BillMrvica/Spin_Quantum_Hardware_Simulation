import klayout.db as kdb
import numpy as np
from typing import List, Tuple, Callable, TypedDict
from Map import Map
from PinMatcher import PinMatcher
from skimage.graph import route_through_array

class ObsRule(TypedDict):
    """Specified type hint for obstacle rules."""
    layers: str | list[str]  # Layer(s) can be a string or a list of strings
    bbx: bool               # Whether to use bounding box
    pad: float              # Padding size in microns

class RouteV2(object):
    """
    RouteV2 is a specialized routing class for connecting quantum dot electrodes to pad leads.
    It implements a multi-width routing strategy (40nm, 100nm, 500nm) with specific geometric constraints.
    """

    def __init__(self,
                 file_path: str,
                 cell_name: str,
                 total_area_bbox_layer: str,
                 quantum_dot_layer: str,
                 pad_lead_layer: str,
                 map_resolution: float,
                 obs_rules: list[ObsRule],
                 routing_layer: str,
                 initial_extension_distance: float = 0.5, # in microns
                 ) -> None:
        """
        Args:
            file_path (str): Path to the GDS file.
            cell_name (str): Name of the cell to route within.
            total_area_bbox_layer (str): Layer containing the total area bounding box for routing.
            quantum_dot_layer (str): Layer containing the quantum dot electrodes (Pin A).
            pad_lead_layer (str): Layer containing the pad leads (Pin B).
            map_resolution (float): Resolution of the routing map, in microns.
            obs_rules (list[ObsRule]): List of obstacle rules.
            routing_layer (str): Layer on which the routes will be drawn.
            initial_extension_distance (float): Distance to extend horizontally/vertically from quantum dots before S-bend, in microns.
        """
        self.__layout = kdb.Layout()
        self.__layout.read(file_path)
        self.__cell = self.__layout.cell(cell_name)
        self.__dbu_per_unit = int(1 / self.__layout.dbu + 0.5)

        self.__map_resolution = map_resolution
        self.__routing_layer = routing_layer
        self.__initial_extension_distance = self.to_dbu(initial_extension_distance)

        # Initialize Map
        self.__Map = Map(self.__layout,
                         kdb.Region(self.__cell.shapes(self.try_find_layer(total_area_bbox_layer))).bbox(),
                         self.to_dbu(map_resolution))

        # Generate total obstacle region based on obs_rules
        self.__obs_rules = [obs_rules] if isinstance(obs_rules, dict) else obs_rules
        self.__obs_region = self.generate_obs_region_by_obs_rules()

        # Initialize PinMatcher
        self.__PinMatcher = PinMatcher()
        self.__PinMatcher.get_pins_from_shapes(self.__cell.shapes(self.try_find_layer(quantum_dot_layer)))
        self.__PinMatcher.get_pins_from_shapes(self.__cell.shapes(self.try_find_layer(pad_lead_layer)))

        # Perform initial matching (can be refined later)
        self.__PinMatcher.hungarian_match()
        self.__PinMatcher.sort_matched_pins() # Sort by distance for consistent routing order

    def to_dbu(self, value: float) -> int:
        """Convert a value in micron to dbu unit."""
        return int(value / self.__layout.dbu + 0.5)

    def try_find_layer(self, layer_str: str) -> int:
        """Find layer index of a layer by its string representation."""
        layer_index = self.__layout.find_layer(kdb.LayerInfo.from_string(layer_str))
        if layer_index is None:
            raise ValueError(f"Layer {layer_str} not found in the layout.")
        return layer_index

    def bbx_and_pad_region_on_layers(self, cbpr_layers: str | list[str], bbx: bool, pad: int) -> kdb.Region:
        """Collect all shapes from the cbpr_layer into a single region."""
        if isinstance(cbpr_layers, str): cbpr_layers = [cbpr_layers]
        
        cbpr_region = kdb.Region()
        for clct_layer in cbpr_layers:
            cbpr_region += kdb.Region(self.__cell.shapes(self.try_find_layer(clct_layer)))

        if bbx:
            return cbpr_region.extents(pad)
        else:
            return cbpr_region.sized(pad)

    def generate_obs_region_by_obs_rules(self) -> kdb.Region:
        """Generate the total obstacle region following the obs_rules."""
        obs_region = kdb.Region()
        for obs_rule in self.__obs_rules:
            obs_rule["pad"] = self.to_dbu(obs_rule["pad"])
            obs_region.join_with(self.bbx_and_pad_region_on_layers(*obs_rule.values()))
        return obs_region

    def _create_path_segments(self, points: List[kdb.Point], width: int) -> kdb.Path:
        """
        Creates a kdb.Path from a list of points with a given width.
        All segments will be orthogonal or 45-degree.
        """
        if len(points) < 2:
            raise ValueError("Path must have at least two points.")
        
        path = kdb.Path(points, width=width, bgn_ext=width//2, end_ext=width//2, round=True)
        return path

    def _get_pad_orientation(self, pad_pin_index: int) -> str:
        """
        Determines the orientation of the pad lead (horizontal or vertical) based on its bounding box.
        """
        pad_shape_region = self.__PinMatcher.Pin_B_polygons[pad_pin_index]
        bbox = pad_shape_region.bbox()
        
        if bbox.width() > bbox.height():
            return "horizontal"
        elif bbox.height() > bbox.width():
            return "vertical"
        else:
            # Square pad, assume vertical for now or handle based on context
            return "vertical" # Default to vertical for square pads

    def _get_damping_raster(self, region: kdb.Region, safe_distance: int, hardness: int, step: int) -> tuple[int, int, int, int, np.ndarray]:
        """Get a raster with gradually damping hardness from the boundary of the region to the outside."""
        # Copy the region to avoid modifying the original one
        safely_sized_region = region.dup()

        # Safely size the region
        for i in range(step):
            safely_sized_region += region.sized(safe_distance * (i+1) / step)

        # Get the boundary (m0, n0, m1, n1) and raster of the safely sized region
        m0, n0, m1, n1, damping_raster = self.__Map.get_boundary_and_raster_of_region(safely_sized_region)

        # Normalize the raster to gradually damping hardness
        damping_raster = (damping_raster // (self.__Map.resolution * self.__Map.resolution)) * (hardness // step)

        return m0, n0, m1, n1, damping_raster

    def _prepare_map_for_routing_internal(self, obs_safe_distance: float, obs_hardness: int, obs_damping_step: int,
                                          pin_A_safe_distance: float, pin_B_safe_distance: float, pin_hardness: int,
                                          pin_A_damping_step: int, pin_B_damping_step: int) -> None:
        """
        Prepares the map with soft safe distance on obstacles and pins for routing.
        This replicates the logic of Krouter.prepare_map_for_routing.
        """
        obs_safe_distance_dbu = self.to_dbu(obs_safe_distance)
        pin_A_safe_distance_dbu = self.to_dbu(pin_A_safe_distance)
        pin_B_safe_distance_dbu = self.to_dbu(pin_B_safe_distance)

        # Rasterize the obstacles to -1 in the map
        self.__Map.rasterize_region(self.__obs_region, -1)

        # Rasterize the padded obstacles to obs_hardness in the map on non-negative pixels
        m0, n0, m1, n1, obs_raster = self._get_damping_raster(self.__obs_region.merged(), obs_safe_distance_dbu,
                                                                obs_hardness, obs_damping_step)
        self.__Map.overwrite(obs_raster, obs_raster > 0, lambda x: x >= 0, m0, n0, m1, n1)

        # Rasterize the pins to -2 in the map
        for pin_region_a, pin_region_b in self.__PinMatcher.matched_polygons:
            self.__Map.rasterize_region(pin_region_a, -2, lambda x: True, merge=True)
            self.__Map.rasterize_region(pin_region_b, -2, lambda x: True, merge=True)

        # Safely-size pins region for later recover use in the routing
        self.__safely_sized_PINs_region = self.__PinMatcher.matched_region_A.sized(pin_A_safe_distance_dbu) + \
                                          self.__PinMatcher.matched_region_B.sized(pin_B_safe_distance_dbu)
        self.__safely_sized_PINs_region.merged_semantics = False
        
        # Rasterize damping safely sized Pin A and B
        m0, n0, m1, n1, Pin_A_raster = self._get_damping_raster(self.__PinMatcher.matched_region_A, pin_A_safe_distance_dbu,
                                                                 pin_hardness, pin_A_damping_step)
        self.__Map.overwrite(Pin_A_raster, Pin_A_raster > 0, lambda x: (0 < x) & (x < Pin_A_raster), m0, n0, m1, n1)  # This is called supplemental overwriting

        m0, n0, m1, n1, pin_B_raster = self._get_damping_raster(self.__PinMatcher.matched_region_B, pin_B_safe_distance_dbu,
                                                                 pin_hardness, pin_B_damping_step)
        self.__Map.overwrite(pin_B_raster, pin_B_raster > 0, lambda x: (0 < x) & (x < pin_B_raster), m0, n0, m1, n1)
        
        print("RouteV2: Map prepared for routing.")
        self.__Map.visualize(self.__layout, self.__cell) # Visualize the map after preparation

    def route_paths(self,
                    obs_safe_distance: float, obs_hardness: int, obs_damping_step: int,
                    pin_safe_distance: float, pin_hardness: int, pin_damping_step: int,
                    path_safe_distance: float, path_hardness: int, path_damping_step: int,
                    path_density_hardness: int,
                    round: int = 1,
                    ) -> None:
        """
        Performs the multi-width routing for all matched quantum dot to pad lead pairs.
        """
        # Convert all distances to dbu
        obs_safe_distance_dbu = self.to_dbu(obs_safe_distance)
        pin_safe_distance_dbu = self.to_dbu(pin_safe_distance)
        path_safe_distance_dbu = self.to_dbu(path_safe_distance)

        # Prepare the map with obstacles and pins using the internal method
        self._prepare_map_for_routing_internal(obs_safe_distance, obs_hardness, obs_damping_step,
                                               pin_safe_distance, pin_safe_distance, pin_hardness,
                                               pin_damping_step, pin_damping_step)

        # Define path widths in dbu
        width_40nm = self.to_dbu(0.04)
        width_100nm = self.to_dbu(0.1)
        width_500nm = self.to_dbu(0.5)

        # Define spacing requirements in dbu
        spacing_100nm_min = self.to_dbu(0.1) # 100nm line spacing > 100nm
        spacing_500nm_min = self.to_dbu(0.5) # 500nm line spacing >= 500nm

        all_paths = []

        # Backup the initial prepared map to reset for each path if needed, or for adaptive routing rounds
        initial_prepared_map = self.__Map.map.copy()

        for i, (qd_pin, pad_pin) in enumerate(self.__PinMatcher.matched_pins):
            print(f"Routing path {i+1} from quantum dot {qd_pin} to pad lead {pad_pin}")

            # Reset map for each new path to avoid previous paths influencing current path's initial routing
            self.__Map.overwrite(initial_prepared_map, self.__Map.map != initial_prepared_map)


            # Determine pad lead orientation for perpendicular connection
            pad_orientation = self._get_pad_orientation(self.__PinMatcher.match_idx_sets["B"][i])

            current_point = qd_pin
            
            # Determine overall direction from quantum dot to pad
            dx = pad_pin.x - qd_pin.x
            dy = pad_pin.y - qd_pin.y

            # --- 1. Initial extension from quantum dot (40nm width) ---
            # Extend horizontally or vertically from the quantum dot to avoid overlap.
            # The direction of initial extension should be away from the quantum dot, towards the pad.
            
            target_point_initial_extension = None
            if abs(dx) > abs(dy): # Dominant horizontal movement
                extension_x = current_point.x + (self.__initial_extension_distance if dx > 0 else -self.__initial_extension_distance)
                target_point_initial_extension = kdb.Point(extension_x, current_point.y)
            else: # Dominant vertical movement
                extension_y = current_point.y + (self.__initial_extension_distance if dy > 0 else -self.__initial_extension_distance)
                target_point_initial_extension = kdb.Point(current_point.x, extension_y)
            
            print(f"RouteV2: Attempting initial extension for path {i+1} from {qd_pin} to {target_point_initial_extension}")
            try:
                path_nodes_initial_extension = route_through_array(self.__Map.map,
                                                                   start=self.__Map.layout_coor_to_map_coor(qd_pin.x, qd_pin.y),
                                                                   end=self.__Map.layout_coor_to_map_coor(target_point_initial_extension.x, target_point_initial_extension.y),
                                                                   fully_connected=True,
                                                                   geometric=True,
                                                                   )[0]
                path_40nm_initial = self._create_path_segments([kdb.Point(*self.__Map.map_coor_to_layout_coor(m, n)) for m, n in path_nodes_initial_extension], width_40nm)
                all_paths.append(path_40nm_initial)
                current_point = target_point_initial_extension
                
                # Mark this path as an obstacle for subsequent routing segments
                self.__Map.rasterize_region(kdb.Region(path_40nm_initial).sized(path_safe_distance_dbu), path_hardness, lambda x: x >= 0, merge=True)
                print(f"RouteV2: Initial extension for path {i+1} successful. Current point: {current_point}")

            except Exception as e:
                print(f"RouteV2: Failed initial extension for path {i+1}: {e}")
                continue # Skip to next path if initial extension fails

            # --- 2. 40nm S-bend section (45-degree turns) ---
            # This section needs to be a 45-degree S-bend and meet 100nm spacing.
            # The end point of this section will be the start of the 100nm section.
            
            # Heuristic for 45-degree S-bend: two 45-degree segments
            s_bend_offset_40nm = max(width_40nm, spacing_100nm_min) 

            target_point_40nm_sbend_1 = kdb.Point(current_point.x + (s_bend_offset_40nm if dx > 0 else -s_bend_offset_40nm), current_point.y + (s_bend_offset_40nm if dy > 0 else -s_bend_offset_40nm))
            target_point_40nm_sbend_2 = kdb.Point(target_point_40nm_sbend_1.x + (s_bend_offset_40nm if dx > 0 else -s_bend_offset_40nm), target_point_40nm_sbend_1.y + (s_bend_offset_40nm if dy > 0 else -s_bend_offset_40nm))

            print(f"RouteV2: Attempting 40nm S-bend segment 1 for path {i+1} from {current_point} to {target_point_40nm_sbend_1}")
            try:
                path_nodes_40nm_sbend_1 = route_through_array(self.__Map.map,
                                                              start=self.__Map.layout_coor_to_map_coor(current_point.x, current_point.y),
                                                              end=self.__Map.layout_coor_to_map_coor(target_point_40nm_sbend_1.x, target_point_40nm_sbend_1.y),
                                                              fully_connected=True,
                                                              geometric=True,
                                                              )[0]
                path_40nm_sbend_segment1 = self._create_path_segments([kdb.Point(*self.__Map.map_coor_to_layout_coor(m, n)) for m, n in path_nodes_40nm_sbend_1], width_40nm)
                all_paths.append(path_40nm_sbend_segment1)
                self.__Map.rasterize_region(kdb.Region(path_40nm_sbend_segment1).sized(path_safe_distance_dbu), path_hardness, lambda x: x >= 0, merge=True)
                current_point = target_point_40nm_sbend_1
                print(f"RouteV2: 40nm S-bend segment 1 for path {i+1} successful. Current point: {current_point}")

                print(f"RouteV2: Attempting 40nm S-bend segment 2 for path {i+1} from {current_point} to {target_point_40nm_sbend_2}")
                path_nodes_40nm_sbend_2 = route_through_array(self.__Map.map,
                                                              start=self.__Map.layout_coor_to_map_coor(current_point.x, current_point.y),
                                                              end=self.__Map.layout_coor_to_map_coor(target_point_40nm_sbend_2.x, target_point_40nm_sbend_2.y),
                                                              fully_connected=True,
                                                              geometric=True,
                                                              )[0]
                path_40nm_sbend_segment2 = self._create_path_segments([kdb.Point(*self.__Map.map_coor_to_layout_coor(m, n)) for m, n in path_nodes_40nm_sbend_2], width_40nm)
                all_paths.append(path_40nm_sbend_segment2)
                self.__Map.rasterize_region(kdb.Region(path_40nm_sbend_segment2).sized(path_safe_distance_dbu), path_hardness, lambda x: x >= 0, merge=True)
                current_point = target_point_40nm_sbend_2
                print(f"RouteV2: 40nm S-bend segment 2 for path {i+1} successful. Current point: {current_point}")

            except Exception as e:
                print(f"RouteV2: Failed 40nm S-bend for path {i+1}: {e}")
                continue

            # --- 3. 100nm width section (S-bend) ---
            # This section connects the 40nm section to the 500nm section and meets 100nm spacing.
            
            s_bend_offset_100nm = max(width_100nm, spacing_100nm_min) * 2 

            target_point_100nm_sbend_1 = kdb.Point(current_point.x + (s_bend_offset_100nm if dx > 0 else -s_bend_offset_100nm), current_point.y + (s_bend_offset_100nm if dy > 0 else -s_bend_offset_100nm))
            target_point_100nm_sbend_2 = kdb.Point(target_point_100nm_sbend_1.x + (s_bend_offset_100nm if dx > 0 else -s_bend_offset_100nm), target_point_100nm_sbend_1.y + (s_bend_offset_100nm if dy > 0 else -s_bend_offset_100nm))

            print(f"RouteV2: Attempting 100nm S-bend segment 1 for path {i+1} from {current_point} to {target_point_100nm_sbend_1}")
            try:
                path_nodes_100nm_sbend_1 = route_through_array(self.__Map.map,
                                                               start=self.__Map.layout_coor_to_map_coor(current_point.x, current_point.y),
                                                               end=self.__Map.layout_coor_to_map_coor(target_point_100nm_sbend_1.x, target_point_100nm_sbend_1.y),
                                                               fully_connected=True,
                                                               geometric=True,
                                                               )[0]
                path_100nm_sbend_segment1 = self._create_path_segments([kdb.Point(*self.__Map.map_coor_to_layout_coor(m, n)) for m, n in path_nodes_100nm_sbend_1], width_100nm)
                all_paths.append(path_100nm_sbend_segment1)
                self.__Map.rasterize_region(kdb.Region(path_100nm_sbend_segment1).sized(path_safe_distance_dbu), path_hardness, lambda x: x >= 0, merge=True)
                current_point = target_point_100nm_sbend_1
                print(f"RouteV2: 100nm S-bend segment 1 for path {i+1} successful. Current point: {current_point}")

                print(f"RouteV2: Attempting 100nm S-bend segment 2 for path {i+1} from {current_point} to {target_point_100nm_sbend_2}")
                path_nodes_100nm_sbend_2 = route_through_array(self.__Map.map,
                                                               start=self.__Map.layout_coor_to_map_coor(current_point.x, current_point.y),
                                                               end=self.__Map.layout_coor_to_map_coor(target_point_100nm_sbend_2.x, target_point_100nm_sbend_2.y),
                                                               fully_connected=True,
                                                               geometric=True,
                                                               )[0]
                path_100nm_sbend_segment2 = self._create_path_segments([kdb.Point(*self.__Map.map_coor_to_layout_coor(m, n)) for m, n in path_nodes_100nm_sbend_2], width_100nm)
                all_paths.append(path_100nm_sbend_segment2)
                self.__Map.rasterize_region(kdb.Region(path_100nm_sbend_segment2).sized(path_safe_distance_dbu), path_hardness, lambda x: x >= 0, merge=True)
                current_point = target_point_100nm_sbend_2
                print(f"RouteV2: 100nm S-bend segment 2 for path {i+1} successful. Current point: {current_point}")

            except Exception as e:
                print(f"RouteV2: Failed 100nm S-bend for path {i+1}: {e}")
                continue

            # --- 4. 500nm width section ---
            # This section connects the 100nm section to the pad lead, perpendicularly, and meets 500nm spacing.
            
            target_intermediate_point_500nm = None
            if pad_orientation == "horizontal":
                target_intermediate_point_500nm = kdb.Point(current_point.x, pad_pin.y + (spacing_500nm_min if dy > 0 else -spacing_500nm_min))
            else: # vertical
                target_intermediate_point_500nm = kdb.Point(pad_pin.x + (spacing_500nm_min if dx > 0 else -spacing_500nm_min), current_point.y)

            print(f"RouteV2: Attempting 500nm segment 1 for path {i+1} from {current_point} to {target_intermediate_point_500nm}")
            try:
                path_nodes_500nm_segment1 = route_through_array(self.__Map.map,
                                                                start=self.__Map.layout_coor_to_map_coor(current_point.x, current_point.y),
                                                                end=self.__Map.layout_coor_to_map_coor(target_intermediate_point_500nm.x, target_intermediate_point_500nm.y),
                                                                fully_connected=True,
                                                                geometric=True,
                                                                )[0]
                path_500nm_segment1 = self._create_path_segments([kdb.Point(*self.__Map.map_coor_to_layout_coor(m, n)) for m, n in path_nodes_500nm_segment1], width_500nm)
                all_paths.append(path_500nm_segment1)
                self.__Map.rasterize_region(kdb.Region(path_500nm_segment1).sized(path_safe_distance_dbu), path_hardness, lambda x: x >= 0, merge=True)
                current_point = target_intermediate_point_500nm
                print(f"RouteV2: 500nm segment 1 for path {i+1} successful. Current point: {current_point}")
                
                print(f"RouteV2: Attempting 500nm segment 2 for path {i+1} from {current_point} to {pad_pin}")
                path_nodes_500nm_segment2 = route_through_array(self.__Map.map,
                                                                start=self.__Map.layout_coor_to_map_coor(current_point.x, current_point.y),
                                                                end=self.__Map.layout_coor_to_map_coor(pad_pin.x, pad_pin.y),
                                                                fully_connected=True,
                                                                geometric=True,
                                                                )[0]
                path_500nm_segment2 = self._create_path_segments([kdb.Point(*self.__Map.map_coor_to_layout_coor(m, n)) for m, n in path_nodes_500nm_segment2], width_500nm)
                all_paths.append(path_500nm_segment2)
                self.__Map.rasterize_region(kdb.Region(path_500nm_segment2).sized(path_safe_distance_dbu), path_hardness, lambda x: x >= 0, merge=True)
                current_point = pad_pin
                print(f"RouteV2: 500nm segment 2 for path {i+1} successful. Current point: {current_point}")

            except Exception as e:
                print(f"RouteV2: Failed 500nm section for path {i+1}: {e}")
                continue

        # Insert all generated paths into the layout
        for path in all_paths:
            self.__cell.shapes(self.__layout.layer(kdb.LayerInfo.from_string(self.__routing_layer))).insert(path)
        
        print(f"RouteV2: Inserted {len(all_paths)} paths to layer {self.__routing_layer}.")

    def save_layout(self, file_path: str = None) -> None:
        """Save the current layout."""
        file_path = file_path if file_path else "routed_design_v2.gds"
        self.__layout.write(file_path)
        print(f"RouteV2: Layout saved to {file_path}.")

    @property
    def layout(self) -> kdb.Layout:
        return self.__layout

    @property
    def cell(self) -> kdb.Cell:
        return self.__cell

    @property
    def Map(self) -> Map:
        return self.__Map

    @property
    def PinMatcher(self) -> PinMatcher:
        return self.__PinMatcher
