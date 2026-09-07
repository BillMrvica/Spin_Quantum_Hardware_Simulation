import klayout.db as kdb
import numpy as np

from skimage.graph import MCP_Geometric
from typing import List, Tuple, Callable

# use int(float_value + 0.5) to round to the nearest integer.

class Map(object):
    """
    A class to create and manage a discrete routing map (a matrix) from a layout (bounding box).
        - All numbers are integers inside this class.
        - All caculations are done in database units (dbu).

    About self.__map:
        - map dimentsion is (M, N) where M is the number of rows and N is the number of columns.
        - map idx is (m, n) where m is the row index and n is the column index.
        - To access the cost of a pixel, use map[m, n].
    
    About layout coordinate system:
        - The default unit is dbu (database unit).
        - The origin is (X0, Y0), the bottom left corner of the total area bounding box.
        - The dimension is (width, height) of the total area bounding box.
        - layout coordinate is (x, y) where x is the horizontal coordinate and y is the vertical coordinate.
    
    About conversion between layout and map coordinates:
        - self.layout_coor_to_map_coor(x, y) converts (x, y) to (m, n).
        - self.map_coor_to_layout_coor(m, n) converts (m, n) to (x, y).
        - x <-> n, y <-> m
    """
     
    def __init__(self, layout: kdb.Layout, total_area_bbox: kdb.Box, resolution: int) -> None:
        """Initializes the map object with the given parameters.

        Args:
            layout (kdb.Layout): The layout is only used to get the dbu (database unit) value.
            total_area_bbox (kdb.Box): Bounding box of the total area for routing.
            resolution (int): The size of one map pixel in dbu.

        Note:
            self.__map (numpy.ndarray): A 2D array where walkable pixels have positive values, otherwise negative values.
        """
        # Make sure total_area_bbox is a kdb.Box
        if not isinstance(total_area_bbox, kdb.Box):
            raise TypeError("total_area_bbox must be a kdb.Box object!")

        # Make sure box is not divisible by resolution. This is for algorithm stability.
        if total_area_bbox.width() % resolution == 0:
            total_area_bbox.right += resolution // 2
            print("Map: Warning! The width of the total area bounding box is divisible by the resolution. The box has been expanded to the right by 1/2 resolution.")
        if total_area_bbox.height() % resolution == 0:
            total_area_bbox.top += resolution // 2
            print("Map: Warning! The height of the total area bounding box is divisible by the resolution. The box has been expanded to the top by 1/2 resolution.")

        # !calculation involving float numbers! but the result is an integer
        self.__dbu_per_unit = int(1 / layout.dbu + 0.5)  # dbu/um, usually 1000. This is an integer.
        
        # Dimensions of the map
        self.__resol = resolution
        self.__N = total_area_bbox.width() // self.__resol + 1  # Number of columns
        self.__M = total_area_bbox.height() // self.__resol + 1  # Number of rows

        # Dimensions of the total area bounding box
        self.__X0 = total_area_bbox.left
        self.__Y0 = total_area_bbox.bottom
        self.__X1 = total_area_bbox.right
        self.__Y1 = total_area_bbox.top
        self.__width = total_area_bbox.width()
        self.__height = total_area_bbox.height()
        self.__origin = (self.__X0, self.__Y0)

        # object attributes
        self.__total_area_bbox = total_area_bbox
        self.__resol = resolution # Resolution in database units

        print(f"Map: Initializing map with total_area_bbox (dbu): {total_area_bbox.bbox()}, resolution (dbu): {resolution}")
        print(f"Map: Calculated map dimensions: N={self.__N}, M={self.__M}")

        self.__initialize_map()
    
    def __initialize_map(self) -> None:
        """
        Mesh the layout geometry.
            - Create an initialized numpy matrix (the map) from the layout geometry, with cost of each pixel set to 1.
            - Create X and Y meshgrid for the map coordinates.
        """
        # Initialize the map with a default cost of 0 for walkable terrain (graph nodes)
        self.__map = np.ones((self.__M, self.__N), dtype=int)
        # Generate pixel array

        print(f"Map: Created a {self.__N}x{self.__M} map matrix with a resolution of {self.__resol / self.__dbu_per_unit} µm.")

    def print_info(self) -> None:
        """Print map attributes information."""
        print("--- Map Attributes ---")
        print(f"Map resolution: {self.__resol / self.__dbu_per_unit} µm")
        print(f"Database unit: {1 / self.__dbu_per_unit} µm/dbu")
        print(f"Map dimensions: {self.__N} columns, {self.__M} rows")
        print(f"Map origin: ({self.__X0 / self.__dbu_per_unit}, {self.__Y0 / self.__dbu_per_unit}) µm")
        print(f"Map width: {self.__width / self.__dbu_per_unit} µm, height: {self.__height / self.__dbu_per_unit} µm")

    def layout_coor_to_map_coor(self, x: int, y: int) -> Tuple[int, int]:
        """Identify which map pixel contains the layout coordinate (x, y).
        
        Uses the map resolution and bounding box to calculate the conversion.

        Args:
            x (int): X coordinate in dbu units.
            y (int): Y coordinate in dbu units.

        Returns:
            tuple[int, int]: 
                - (m, n): indices of the map matrix
                - The corresponding matrix element is M[m, n].
                
        Raises:
            ValueError: If coordinates are outside the total area bounding box.
        """
        if x < self.__X0 or x > self.__X1 or y < self.__Y0 or y > self.__Y1:
            raise ValueError(f"Coordinates ({x}, {y}) are out of the total area bounding box with corners ({self.__X0}, {self.__Y0}) and ({self.__X1}, {self.__Y1}).")

        # Convert layout coordinates to map indices
        m = (y - self.__Y0) // self.__resol
        n = (x - self.__X0) // self.__resol

        return m, n
    
    def map_coor_to_layout_coor(self, m: int, n: int) -> Tuple[int, int]:
        """Convert map index (m, n) to pixel center coordinate (x, y).

        Args:
            m (int): Row index of the map matrix.
            n (int): Column index of the map matrix.

        Returns:
            tuple[int, int]: 
                - (x, y): coordinates in dbu units
                - The center of the pixel in the layout
                
        Raises:
            ValueError: If indices are out of bounds for the map dimensions.
        """
        if n < 0 or n >= self.__N or m < 0 or m >= self.__M:
            raise ValueError(f"Indices ({n}, {m}) are out of bounds for the map with dimensions ({self.__N}, {self.__M}).")

        # Convert map indices to layout coordinates
        # self.__resol//2 garantees integer calculation
        x = self.__X0 + n * self.__resol + self.__resol // 2
        y = self.__Y0 + m * self.__resol + self.__resol // 2

        return x, y

    def pixel(self, m: int, n: int) -> kdb.Region:
        """Create a pixel (kdb.Region) at coordinate (m, n).
        
        Args:
            m (int): Row coordinate of the map.
            n (int): Column coordinate of the map.
            
        Returns:
            kdb.Region: A rectangular region at coordinate (m, n) with size of self.__resol.
        """
        x, y = self.map_coor_to_layout_coor(m, n)
        return kdb.Region(kdb.Box(x - self.__resol//2, y - self.__resol//2, x + self.__resol//2, y + self.__resol//2))

    def __rasterize_polygon(self, poly_region: kdb.Region, cost: int, 
                          overwrite_condition: Callable[[np.ndarray], np.ndarray] = lambda x: x >= 0) -> Tuple[List[int], List[int]]:
        """Convert a polygon (single-connected kdb.Region) of any shape to map pixels, and set the cost of these pixels.
        
        This method uses kdb.Region.rasterize() as the core.
        
        Args:
            poly_region (kdb.Region): Region to be rasterized.
            cost (int): Cost to be set. You can set cost to 1 if you only want to get indices of the region.
            overwrite_condition (callable): A function, used to determine which pixels to overwrite.
                The function takes a numpy array as input and returns a boolean array of the same shape.
                Default is lambda x: x >= 0, which means all walkable pixels (cost >= 0) will be overwritten.
            
        Returns:
            tuple[list[int], list[int]]: (list of m indices, list of n indices) of the rasterized pixels.
        """
        if poly_region.is_empty():
            print("Map: The polygon to rasterize is empty.")
            return

        # Get the boundary pixels of the polygon
        m0, n0, m1, n1, poly_raster = self.get_boundary_and_raster_of_region(poly_region)

        # Create a cost array for later overwriting
        cost_array = np.full(poly_raster.shape, cost, dtype=int)

        # Overwrite the map with the rasterized polygon
        self.overwrite(cost_array, poly_raster > 0, overwrite_condition, m0, n0, m1, n1)

        return
    
    def rasterize_region(self, rast_region: kdb.Region, cost: int, 
                        overwrite_condition: Callable[[np.ndarray], np.ndarray] = lambda x: x >= 0, 
                        merge: bool = True) -> None:
        """Rasterize a kdb.Region (specifically, a multi-polygon region) and set the cost of the pixels in the map.

        This method handles multiple polygons by rasterizing each one individually. This avoids too large memory usage.
        
        Args:
            rast_region (kdb.Region): Region (containing multiple polygons) to be rasterized.
            cost (int): Cost to be set for the rasterized pixels.
            overwrite_condition (callable): A function, used to determine which pixels to overwrite.
                The function takes a numpy array as input and returns a boolean array of the same shape.
                Default is lambda x: x >= 0, which means all walkable pixels (cost >= 0) will be overwritten.
            merge (bool): Whether to merge the region after clipped to total area box.
                It is recommended True if each merged polygon is still small; False if polygons become huge after merged.
        """
        # Clip the region to the total area box; this will not merge the region.
        rast_region = rast_region & self.__total_area_bbox
        
        if rast_region.is_empty():
            print("Map: The region to rasterize is empty after clipping with the total area box.")
            return

        # Merge the region to reduce the number of polygons if specified
        if merge: rast_region = rast_region.merged()
        
        # Rasterize each polygon in the region
        for rast_polygon in rast_region.each():
            self.__rasterize_polygon(kdb.Region(rast_polygon), cost, overwrite_condition)

        return

    def get_boundary_and_raster_of_region(self, gbrr_region: kdb.Region) -> Tuple[int, int, int, int, np.ndarray]:
        """Get the boundary and klayout-raster of a region on the map, so that you can process the raster and then overwrite it.

        - How to use these returned values: 
            - self.overwrite(raster, mask, overwrite_condition, m0, n0, m1, n1)

        Args:
            gbrr_region (kdb.Region): The polygon region to get the boundary pixels from.

        Returns:
            tuple[int, int, int, int, numpy.ndarray]: 
                - (m0, n0, m1, n1): boundary pixel indices of the region on the map.
                - raster (numpy.ndarray): rasterized array of the region.
        """
        # Check if the region is empty
        if gbrr_region.is_empty():
            raise ValueError("Map: The region to get boundary and raster is empty!")

        # Clip the region to the total area box
        gbrr_region = gbrr_region & self.__total_area_bbox

        # Get the bottom-left and top-right pixel coor of the polygon's bounding box
        m0, n0 = self.layout_coor_to_map_coor(gbrr_region.bbox().left, gbrr_region.bbox().bottom)
        m1, n1 = self.layout_coor_to_map_coor(gbrr_region.bbox().right, gbrr_region.bbox().top)

        # Get the bottom-left corner of the bottom-left pixel of the polygon's bounding box
        x0, y0 = self.map_coor_to_layout_coor(m0, n0)
        x0 -= self.__resol // 2
        y0 -= self.__resol // 2

        # Get the top-right corner of the top-right pixel of the polygon's bounding box
        x1, y1 = self.map_coor_to_layout_coor(m1, n1)
        x1 += self.__resol // 2
        y1 += self.__resol // 2

        # Get relative coordinates of the polygon pixels
        nx = (x1 - x0) // self.__resol  # Number of columns in the rasterized area
        ny = (y1 - y0) // self.__resol  # Number of rows in the rasterized area

        # Get raster of the region
        raster = np.array(gbrr_region.rasterize(kdb.Point(x0, y0), kdb.Vector(self.__resol, self.__resol), nx, ny))

        return m0, n0, m1, n1, raster

    def overwrite(self, content: np.ndarray, content_mask: np.ndarray, 
                  overwrite_condition: Callable[[np.ndarray], np.ndarray] = lambda x: True, 
                  m0: int = 0, n0: int = 0, m1: int = None, n1: int = None) -> None:
        """Overwrite a part of the map with new content based on a condition and a mask.
        
        - Use self.get_boundary_and_raster_of_region() to get (m0, n0, m1, n1) of a polygon region.
        - A trick: overwrite(content, content != self.__map) will directly set the map to content.

        Args:
            content (numpy.ndarray): New content to overwrite the map with.
            content_mask (numpy.ndarray): A boolean mask indicating which pixels in the content should be applied.
            overwrite_condition (callable): A function, used to determine which pixels to overwrite.
                The function takes a numpy array as input and returns a boolean array of the same shape.
            m0 (int): Starting row index of the region to overwrite. Default is 0.
            n0 (int): Starting column index of the region to overwrite. Default is 0
            m1 (int): Ending row index of the region to overwrite. Default is None, which means the last row.
            n1 (int): Ending column index of the region to overwrite. Default is None, which means the last column.
        Raises:
            ValueError: If the dimensions of the content do not match the specified region dimensions.
        """
        if m1 is None: m1 = self.__M - 1
        if n1 is None: n1 = self.__N - 1

        if content.shape != (m1 - m0 + 1, n1 - n0 + 1):
            raise ValueError(f"New map dimensions {content.shape} do not match the specified region dimensions {(m1 - m0 + 1, n1 - n0 + 1)}.")

        # Identify the pixels to overwrite based on the condition and the content mask
        overwrite_mask = overwrite_condition(self.__map[m0:m1 + 1, n0:n1 + 1]) & content_mask

        # Set the cost of the identified pixels
        self.__map[m0:m1 + 1, n0:n1 + 1][overwrite_mask] = content[overwrite_mask]

        return

    def initialize(self, cost: int = 1) -> None:
        """Re-initialize the map to a uniform cost value.
        
        Args:
            cost (int): Cost value to initialize the map with. Default is 1 (walkable terrain).
        """
        self.__map.fill(cost)

    def __create_pixel_cell(self) -> None:
        """Create a mesh of the map with all edges and nodes."""
        # --- Setup: Get KLayout objects ---
        lv = kdb.Application.instance().main_window().current_view()
        layout = lv.active_cellview().layout()
        top_cell = lv.active_cellview().cell # Our main cell
        dbu = layout.dbu

        # --- 1. Create a "Unit Cell" to hold the base shape ---
        unit_cell_name = "BOX_10x20_CELL"
        unit_cell = layout.create_cell(unit_cell_name)

        # Define the geometry layer inside the unit cell
        geo_layer = layout.insert_layer(kdb.LayerInfo(1, 0))
        base_shape = kdb.Box(0, 0, int(10 / dbu), int(20 / dbu))
        unit_cell.shapes(geo_layer).insert(base_shape)
        print(f"Map: Created a unit cell '{unit_cell_name}' with the base shape.")

        # --- 2. Define Array Parameters ---
        num_rows = 5
        num_cols = 10
        col_vector_um = kdb.DVector(15.0, 0)
        row_vector_um = kdb.DVector(0, 25.0)

        # The origin point of the entire array in the top cell
        array_origin = kdb.Point(0, 0)
        origin_trans = kdb.Trans(array_origin)

        # --- 3. Create and insert a single Cell Instance Array object ---
        print("Map: Creating a hierarchical cell instance array...")

        # Create the CellInstArray object
        # It needs: the index of the cell to instance, the origin transformation,
        # the column and row vectors (in dbu), and the number of columns and rows.
        inst_array = kdb.CellInstArray(unit_cell.cell_index(), 
                                    origin_trans, 
                                    col_vector_um / dbu, 
                                    row_vector_um / dbu, 
                                    num_cols, 
                                    num_rows)

        # Insert this single object into our top cell
        top_cell.insert(inst_array)

        print(f"Map: Successfully created an array of {num_rows * num_cols} instances of '{unit_cell_name}'.")

    def visualize(self, layout: kdb.Layout, cell: kdb.Cell, cost: int = 1) -> None:
        """Draw pixels of certain cost on the layer "Mesh/cost" in the layout.
        Args:
            layout (kdb.Layout): The layout to visualize the map on.
            cell (kdb.Cell): The cell in which to draw the map.
            cost (int): The cost of the pixels to visualize. Default is 1, which means walkable pixels.
        """
        for m, n in list(zip(*np.where(self.__map == cost))):
            # Get the layout coordinates of the pixel center
            x, y = self.map_coor_to_layout_coor(m, n)
            # Insert a box representing the pixel
            cell.shapes(layout.layer(f"Mesh_{cost}")).insert(
                kdb.Box(int(x - self.__resol // 2), int(y - self.__resol // 2), 
                        int(x + self.__resol // 2), int(y + self.__resol // 2)))

        print(f"Map: Visualized map on layer \"Mesh_{cost}\".")


    # Properties to get private attributes
    @property
    def map(self) -> np.ndarray:
        """numpy.ndarray: The map matrix containing the cost of each pixel."""
        return self.__map

    @property
    def M(self) -> int:
        """int: Number of rows in the map."""
        return self.__M
        
    @property
    def N(self) -> int:
        """int: Number of columns in the map."""
        return self.__N
        
    @property
    def width(self) -> int:
        """int: Width of the total area bounding box in dbu."""
        return self.__width
        
    @property
    def height(self) -> int:
        """int: Height of the total area bounding box in dbu."""
        return self.__height
        
    @property
    def origin(self) -> Tuple[int, int]:
        """tuple[int, int]: Origin (X0, Y0) of the total area bounding box."""
        return self.__origin
        
    @property
    def resolution(self) -> int:
        """int: Resolution of the map in dbu."""
        return self.__resol
    
    @property
    def total_area_bbox(self) -> kdb.Box:
        """kdb.Box: Bounding box of the total area for routing."""
        return self.__total_area_bbox

    @property
    def dbu_per_unit(self) -> int:
        """int: Database units per micron."""
        return self.__dbu_per_unit
    
    @property
    def MCPG(self, offsets=None, fully_connected=True) -> MCP_Geometric:
        """skimage.graph.MCP_Geometric: Minimum Cost Path (MCP) object for pathfinding on the map."""
        return MCP_Geometric(self.__map, offsets=offsets, fully_connected=fully_connected)
