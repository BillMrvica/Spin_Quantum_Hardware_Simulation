import numpy as np
import klayout.db as kdb
from typing import Callable
from functools import cached_property
from scipy.optimize import linear_sum_assignment

class PinMatcher(object):
    """
    Class to match start-end pins. 
        - A pin is kdb.Point object.
        - Always get matched pins/shapes/polygons through __matched_idx_sets, so that you only have to sort the idx to get all other things sorted.

    Vars:
        self.__pin_sets (dict): Dictionary to hold pin sets.
            - "A" (list[kdb.Point]): List of pins for set A.
            - "B" (list[kdb.Point]): List of pins for set B.
            - List of pins must contain kdb.Point objects.
        self.__shape_sets (dict): Dictionary to hold shapes where pins are located.
            - "A" (list[kdb.Shape]): Shapes for set A.
            - "B" (list[kdb.Shape]): Shapes for set B.
        self.__match_idx_sets (dict): Dictionary to hold matched indices.
            - "A" (list[int]): Indices of matched pins in set A.
            - "B" (list[int]): Indices of matched pins in set B.
    """

    def __init__(self) -> None:
        """
        """
        # Initialize pin sets
        self.__pin_sets = {
            "A": None,
            "B": None,
        }

        # Initialize shape sets, which contains list of kdb.Shape instead of kdb.Shapes.
        self.__shape_sets = {
            "A": None,
            "B": None,
        }

        self.__match_idx_sets = {
            "A": None,
            "B": None,
        }

    def add_pins(self, pins: list[kdb.Point]) -> None:
        """
        Add pins to self.__pin_sets.
        Args:
            pins [kdb.Point,]: List of kdb.Point objects representing pin positions.
        """
        if self.__pin_sets["A"] is None:
            self.__pin_sets["A"] = pins
            print(f"PinMatcher: Added {len(pins)} pins to Pin set A.")
        elif self.__pin_sets["B"] is None:
            self.__pin_sets["B"] = pins
            print(f"PinMatcher: Added {len(pins)} pins to Pin set B.")
        else:
            raise ValueError("Pin sets already contain both A and B pins.")
    
    def add_shape_list(self, shape_list: list[kdb.Shape]) -> None:
        """
        Add shapes to self.__shape_sets.
        Args:
            shapes (list[kdb.Shape]): list of shape to be added.
        Raises:
            TypeError: If shape_list is not a list.
            ValueError: If both shape sets A and B are already populated.
        """
        # Check if shape_list is a list of kdb.Shape objects instead of kdb.Shapes
        if not isinstance(shape_list, list):
            raise TypeError("shape_list must be a list of kdb.Shape objects.")
        
        if self.__shape_sets["A"] is None:
            self.__shape_sets["A"] = shape_list
            print(f"PinMatcher: Added {len(shape_list)} shapes to shape set A.")
        elif self.__shape_sets["B"] is None:
            self.__shape_sets["B"] = shape_list
            print(f"PinMatcher: Added {len(shape_list)} shapes to shape set B.")
        else:
            raise ValueError("Shape sets already contain both A and B shapes.")

    def get_pins_from_shapes(self, gpfs_shapes: kdb.Shapes) -> list[kdb.Point]:
        """
        Get pins from shapes list, and store the shapes. 
            - It takes the last point of the path, or the position of the text, of the center of other shape as pin.
            - Use `cell.shapes(layer_idx)` as input.
            - The shapes had best to be the same type, e.g. all paths or all texts.
        Args:
            gpfs_shapes (kdb.Shapes): Shapes list to extract pins from.
        """
        # Empty list to collect pins and shapes
        gpfs_pins = []

        # For different type of shapes, use different method to extract pins
        for gpfs_shape in gpfs_shapes.each():
            if gpfs_shape.is_text():
                # Get the position of the text
                gpfs_pins.append(gpfs_shape.text.position())
            elif gpfs_shape.is_path():
                # Get the last point of the path
                gpfs_pins.append(list(gpfs_shape.each_point())[-1])
            else:
                # For other shapes, use the center of the bounding box
                gpfs_pins.append(gpfs_shape.bbox().center())
                
        self.add_pins(gpfs_pins)
        self.add_shape_list(list(gpfs_shapes.each()))

        return

    def hungarian_match(self, distance_func: Callable[[kdb.Point, kdb.Point], float] = None, 
                        additional_cost_matrix: list[list[int]] = None, weight_ratio: float = 1) -> list[tuple[kdb.Point, kdb.Point]]:
        """
        Use the Hungarian algorithm to find the optimal matching between pin sets A and B.
        
        The default distance function is Euclidean distance, you can provide your own function.
        Additional cost matrix can be provided to adjust the cost, in case cost matrix is hard to calculate by distance_func.
        The final cost matrix is a weighted average of the two matrices.

        Args:
            distance_func (Callable[[kdb.Point, kdb.Point], float]): Function to calculate distance between two pins.
                - Should take two kdb.Point objects and return a float distance.
            additional_cost_matrix (list[list[int]]): Precomputed cost matrix. If provided, distance_func will be ignored.
            weight_ratio (float): For weighted average between the matrix calculated by distance_func and provided cost_matrix.
        
        Returns:
            List of tuples: [(pin_A, pin_B), ...] representing matched pins.

        This method finds the optimal match (minimum total distance) and guarantees no-crossing of pin-pin paths.
        """
        # Check if both pin sets A and B are populated
        if self.__pin_sets["A"] is None or self.__pin_sets["B"] is None:
            raise ValueError("Both pin sets A and B must be populated before matching.")

        # Calculate cost matrix
        if distance_func is not None:
            # Use the provided distance function to calculate cost matrix
            self.__cost_matrix = [[distance_func(pin_a, pin_b) for pin_b in self.__pin_sets["B"]] for pin_a in self.__pin_sets["A"]]
        else:
            # Default distance function is Euclidean distance
            self.__cost_matrix = [[pin_a.distance(pin_b) for pin_b in self.__pin_sets["B"]] for pin_a in self.__pin_sets["A"]]
        
        # Weighted average previous matrix and additional_cost_matrix if the latter are provided
        if additional_cost_matrix is not None: 
            self.__cost_matrix = weight_ratio * np.array(self.__cost_matrix) + (1 - weight_ratio) * np.array(additional_cost_matrix)
      
        # Apply the Hungarian algorithm
        row_indices, col_indices = linear_sum_assignment(self.__cost_matrix)
        
        # Create matched pairs based on the optimal assignment, and store it
        self.__match_idx_sets["A"] = list(row_indices)
        self.__match_idx_sets["B"] = list(col_indices)

        print(f"PinMatcher: Found {len(row_indices)} matched pins using Hungarian algorithm.")

    def visualize(self, layout: kdb.Layout, cell: kdb.Cell, vm_path_width: int = 1000) -> None:
        """
        Connect the matched pins by path on the layer "Match".
        Args:
            layout (kdb.Layout): The layout object to insert shapes into.
            cell (kdb.Cell): The cell object to insert shapes into.
            vm_path_width (int): Width of the path to draw between matched pins, in dbu.
        """
        # Firstly check if pins have been matched
        if self.__match_idx_sets["A"] is None:
            raise ValueError("Pins haven't been matched!")

        # Draw paths connecting matched pins on the layer
        for pin_a, pin_b in self.matched_pins:
            # Create a path between the matched pins
            path = kdb.Path([pin_a, pin_b], 
                            width=int(vm_path_width + 0.5),
                            bgn_ext=int(vm_path_width/2 + 0.5),
                            end_ext=int(vm_path_width/2 + 0.5),
                            round=True)
            # Insert the path
            cell.shapes(layout.layer("Match")).insert(path)

        num = len(self.__match_idx_sets["A"])
        print(f"PinMatcher: Visualized {num} matched pins on layer \"Match\".")

    def read_match(self, connections: kdb.Region) -> None:
        """Read matched pins from regions.
        This method allows you to assign/adjust the match mannually by drawing paths, and then read them.
        Args:
            connections (kdb.Region): Multiple polygons connecting matched pin-shapes.
        """
        # Disable merging semantics to avoid crossing paths merged
        connections.merged_semantics = False
        self.__match_idx_sets["A"] = []
        self.__match_idx_sets["B"] = []

        # Take a polygon from Pin A
        for i, poly_a in enumerate(self.Pin_A_polygons):
            # Find the connection that overlaps with poly_a
            connection = poly_a.pull_overlapping(connections)

            # Find which polygon in Pin B overlaps with the connection
            for j, poly_b in enumerate(self.Pin_B_polygons):
                    # If poly_b overlaps with the connection, it is a matched polygon
                    if not poly_b.overlapping(connection).is_empty():
                        self.__match_idx_sets["A"].append(i)
                        self.__match_idx_sets["B"].append(j)
                        break  # Break the inner loop as soon as we found a match for poly_a
        
        print(f"PinMatcher: Read {len(self.__match_idx_sets['A'])} matched pins.")

    def sort_matched_pins(self) -> None:
        """Sort the matched pins in a ascending order of cost (distance between them)."""
        if self.__match_idx_sets["A"] is None or self.__match_idx_sets["B"] is None:
            raise ValueError("No matched pins found. Please run match first.")
        
        # costs = [self.__cost_matrix[a, b] for a, b in zip(self.__match_idx_sets["A"], self.__match_idx_sets["B"])]
        costs = [self.__pin_sets["A"][i].distance(self.__pin_sets["B"][j]) for i, j in zip(self.__match_idx_sets["A"], self.__match_idx_sets["B"])]
        # Sort the indices based on costs
        sorted_indices = np.argsort(costs)
        # Sort the matched pins and polygons based on the sorted indices
        self.__match_idx_sets["A"] = [self.__match_idx_sets["A"][i] for i in sorted_indices]
        self.__match_idx_sets["B"] = [self.__match_idx_sets["B"][i] for i in sorted_indices]

        print(f"PinMatcher: Sorted {len(self.__match_idx_sets['A'])} matched pins by cost.")

    def clockwise_sort(self, reference: str = "B", counter: bool = False) -> None:
        """Take the closest pair of pins as start, and then sort the rest by clockwise or counter-clockwise order.
        
        Args:
            reference (str): Which pin set to use as reference for nearest neighbor search.
                - "A": Use pin set A as reference.
                - "B": Use pin set B as reference.
                - Default is "B".
            counter (bool): If True, sort in counter-clockwise order. Default is False (clockwise).
        """
        if reference not in ["A", "B"]:
            raise ValueError("reference must be either 'A' or 'B'.")

        if self.__match_idx_sets["A"] is None or self.__match_idx_sets["B"] is None:
            raise ValueError("No matched pins found. Please run match first.")
        
        # Take the closest pair of pins as start
        costs = [self.__pin_sets["A"][i].distance(self.__pin_sets["B"][j]) for i, j in zip(self.__match_idx_sets["A"], self.__match_idx_sets["B"])]
        start_pin_idx = self.__match_idx_sets[reference][np.argmin(costs)]  # Starting pin index in the reference set

        # Get the center of all pins in the reference set
        refer_pins_region = kdb.Region()
        for pin_region in (self.Pin_A_polygons if reference == "A" else self.Pin_B_polygons):
            refer_pins_region += pin_region
        center = refer_pins_region.bbox().center()

        # Calculate rotation angles around center between each other pins and the start pin
        # clockwise angle \theta = sgn(a \times b) * arccos(a . b / |a||b|)
        start_vector = self.__pin_sets[reference][start_pin_idx] - center
        refer_vectors = [self.__pin_sets[reference][i] - center for i in self.__match_idx_sets[reference]]
        angles = np.array([self.clockwise_angle_from_v1_to_v2(start_vector, v) for v in refer_vectors])
        angles = (angles + 2 * np.pi) % (2 * np.pi)  # Shift to [0, 2pi]

        # If counter is True, sort in counter-clockwise order
        if counter: angles *= -1

        # Get the sorted indices based on angles
        sorted_idx = sorted(range(len(angles)), key=lambda i: angles[i])
        
        # Reorder the matched idx sets based on sorted_idx
        self.__match_idx_sets["A"] = [self.__match_idx_sets["A"][i] for i in sorted_idx]
        self.__match_idx_sets["B"] = [self.__match_idx_sets["B"][i] for i in sorted_idx]
        
        print(f"PinMatcher: Sorted {len(self.__match_idx_sets['A'])} matched pins by {counter}-counter clockwise order on Pin set {reference}.")


    def text_to_box(self, text: kdb.Text) -> kdb.Box:
        """Convert kdb.Text to kdb.Box of small size (10 dbu) at the text position."""
        pos = text.position()
        return kdb.Box(pos.x - 5, pos.y - 5, pos.x + 5, pos.y + 5)
    
    def shape_to_polygon(self, shape: kdb.Shape) -> kdb.Region:
        """Convert a shape to a region."""
        if shape.is_text():
            # Specially take care of text shapes as they will be ignored by kdb.Region()
            return kdb.Region(self.text_to_box(shape.text))
        elif shape.is_box():
            return kdb.Region(shape.box)
        elif shape.is_polygon():
            return kdb.Region(shape.polygon)
        elif shape.is_path():
            return kdb.Region(shape.path)
        else:
            raise TypeError(f"PinMatcher: Unsupported shape type to polygon: {shape.type}.")

    def clockwise_angle_from_v1_to_v2(self, v1: kdb.Vector, v2: kdb.Vector) -> float:
        """Calculate the clockwise angle from vector v1 to vector v2 in radians. 
        
        Range: [-pi, pi], use (v + 2pi) % (2pi) to shift to [0, 2pi] if needed.

        Args:
            v1 (kdb.Vector): The first vector.
            v2 (kdb.Vector): The second vector.
        Returns:
            float: The clockwise angle from v1 to v2 in radians.
        """
        # Calculate cosine of the angle between v1 and v2, and clip to avoid numerical issues
        cos_theta = np.clip(v1.sprod(v2) / (v1.abs() * v2.abs()), -1, 1)

        # Determine the quadrant using the cross product
        if v1.vprod(v2) == 0:
            return np.arccos(cos_theta)
        else:
            return v1.vprod_sign(v2) * np.arccos(cos_theta)

    @property
    def pin_sets(self) -> dict[str, list[kdb.Point]]:
        """Get the pin sets."""
        return self.__pin_sets
    @property
    def shape_sets(self) -> dict[str, list[kdb.Shape]]:
        """Get the shape sets."""
        return self.__shape_sets
    @property
    def match_idx_sets(self) -> dict[str, list[int]]:
        """Get the matched indices."""
        return self.__match_idx_sets

    @cached_property
    def Pin_A_polygons(self) -> list[kdb.Region]:
        """Get the regions of Pin A shapes."""
        if self.__shape_sets["A"] is None:
            raise ValueError("Pin A shapes not found. Please add shapes to Pin A first.")
        return [self.shape_to_polygon(shape) for shape in self.__shape_sets["A"]]
    @cached_property
    def Pin_B_polygons(self) -> list[kdb.Region]:
        """Get the regions of Pin B shapes."""
        if self.__shape_sets["B"] is None:
            raise ValueError("Pin B shapes not found. Please add shapes to Pin B first.")
        return [self.shape_to_polygon(shape) for shape in self.__shape_sets["B"]]
    
    @property
    def matched_Pin_A(self) -> list[kdb.Point]:
        """Get the matched pins in Pin A set.
        Returns:
            List of kdb.Point: Matched pins in Pin A set.
        """
        if self.__match_idx_sets["A"] is None:
            raise ValueError("No matched pins found. Please run match first.")
        return [self.__pin_sets["A"][i] for i in self.__match_idx_sets["A"]]
    @property
    def matched_Pin_B(self) -> list[kdb.Point]:
        """Get the matched pins in Pin B set.
        Returns:
            List of kdb.Point: Matched pins in Pin B set.
        """
        if self.__match_idx_sets["B"] is None:
            raise ValueError("No matched pins found. Please run match first.")
        return [self.__pin_sets["B"][i] for i in self.__match_idx_sets["B"]]
    @property
    def matched_pins(self) -> list[tuple[kdb.Point, kdb.Point]]:
        """Get the matched pins as a list of tuples (pin_a, pin_b).
        Returns:
            List of tuples: [(pin_a, pin_b), ...] representing matched pins.
        """
        return list(zip(self.matched_Pin_A, self.matched_Pin_B))

    @property
    def matched_polygon_A(self) -> list[kdb.Region]:
        """
        Get the matched shapes as polygons (regions) for Pin A.
        Returns:
            List of kdb.Region: Matched polygons (regions) of Pin A.
        Raises:
            ValueError: If no matched shapes are found.
        """
        if self.__match_idx_sets["A"] is None:
            raise ValueError("No matched shapes found.")
        
        return [self.Pin_A_polygons[i] for i in self.__match_idx_sets["A"]]
    @property
    def matched_polygon_B(self) -> list[kdb.Region]:
        """
        Get the matched shapes as polygons (regions) for Pin B.
        Returns:
            List of kdb.Region: Matched polygons (regions) of Pin B.
        Raises:
            ValueError: If no matched shapes are found.
        """
        if self.__match_idx_sets["B"] is None:
            raise ValueError("No matched shapes found.")
        
        return [self.Pin_B_polygons[i] for i in self.__match_idx_sets["B"]]
    @property
    def matched_polygons(self) -> list[tuple[kdb.Region, kdb.Region]]:
        """ Get the matched polygons as a list of tuples (polygon_A, polygon_B).
        Returns:
            List of tuples: [(polygon_A, polygon_B), ...] representing matched polygons.
        """
        return list(zip(self.matched_polygon_A, self.matched_polygon_B))
    
    @property
    def matched_region_A(self) -> kdb.Region:
        """Get the union of matched polygons in Pin A set."""
        region = kdb.Region()
        for poly in self.matched_polygon_A:
            region += poly
        return region
    @property
    def matched_region_B(self) -> kdb.Region:
        """Get the union of matched polygons in Pin B set."""
        region = kdb.Region()
        for poly in self.matched_polygon_B:
            region += poly
        return region