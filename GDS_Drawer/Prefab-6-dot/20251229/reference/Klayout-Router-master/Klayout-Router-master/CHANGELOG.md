# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 

### Changed

- Auto-patching
  - Intersection of the pattern with the grid is added to patches, together with boxes.
  - No grid layer is generated anymore.





## [2.0.0] - 9/29/2025

I decide to give up the `NetworkX` scheme, it is too slow even if I have adopted lots of optimization, e.g. invisible node sets, subgraph...

### Changed

- **Breaking change**: 
  - The project has been completely refactored from a KLayout macro to a standalone Python application. It no longer runs inside the KLayout GUI but instead operates on GDS files directly from the command line. This requires a dedicated Python environment.
  - The core geometry engine has been migrated from the GUI-dependent `pya` module to the high-performance `klayout.db` (from the `klayout` PyPI package).
  - The program now read and write GDS files directly.
- Refactored the procedure of `Krouter.route()` to fit new `Map.rasterize` method and `soft safe distance` feature.
- Refactored rasterization method.
  - Re-written rasterization method base on the **`klayout` native engine** `klayout.db.Region.rasterize()`. New method is faster, clearer, more lightweight and versatile. Compress time consumption **from 1s/path down to 0.1s/path**
  - Support rasterizing region and **overwriting pixels by condition**.
  - There are now **two ways to do rasterization**: use `Map.rasterize_region()` to set cost conveniently, or use `Map.get_boundary_and_raster_of_region()` and `Map.overwrite()` to do more flexible and personalized overwriting on the map.

### Added

With all new features below, and appropriate parameters set, the routing now can be finished nearly perfect within seconds.

- Utilize **`scikit-image`** PyPI package to do routing, much much more compatible and faster than `NetworkX`.
- Adopt **damping soft safe distance** between paths and around obstacles & pins, avoiding *no-path-found* problem.
- **Clockwise sort matched pins**.
- **Generate animation** of routing process.
- Hungarian match with weighted-average of Euclidean and real-distance cost matrix.
- Two new widgets (klayout macro): `merge_layers.py` and `change_path_width.py`.

## [1.6.0] - 9/8/2025

### Added

- Keep away from all other Pins with a safe distance when routing.
- `PinMatcher` can read manual match from the layout.
- `PinMatcher` can choose to use shortest real path distance as the cost matrix element.
- `Krouter`: `compress_path()`

### Changed

- `GraphM`:  
  - New feature: Keep track of sets of labelled nodes and invisible nodes, instead of deleting them.
    - `labelled_nodes`: dictionary of nodes with different labels
    - `invisible_nodes`: set of invisible nodes.
    - `add_labelled_nodes()`: add nodes into `labelled_nodes`.
    - `visible_graph()`: Create a subgraph avoiding invisible nodes.
- `PinMatcher`:
  - Has three intrinsic sets: `pin_sets`, `shape_sets` and `match_idx_sets`, all other things comes from these three sets
  - All match methods change `match_idx_sets`
  - New methods:
    - `read_match()`
    - `sort_matched_pins()`
    - `shape_to_polygon()`

### Deprecated Method

- `GraphM`:
  - `bidirectional_BFS_routing()`: replaced by `bidirectional_dijkstra()`



## [1.5.0] - 7/25/2025

New feature implemented: routing with path padding, no more overlapping between paths.

### Big ideas:

- All the interactions with the layout (read from layer, write to layer) are now handled by `Krouter`.
- `Map`, `PinMatcher`, and `GraphM` are more compact, and they do not handle layout directly.
- Layers to put reference items (mesh, matches) are now hard-coded named, like f"Mesh_{cost}", "Matches", etc.

### Details:

- `Krouter` (new)
    - `Krouter` creates a `PinMatcher`, a `Map`, and a `GraphM` instance.
    - `Krouter` gets shapes/regions from layers, feed them to `Map` and `PinMatcher`.
    - `Krouter` takes floats in micron (μm) and turn them into ints in dbu for `Map`, `PinMatcher`, and `GraphM`.      
    - Some great methods:
        - `to_dbu`: to convert float in micron to int in dbu.
        - `__find_layer`: find a layer and raise an error if not found.
        - `routing`: routing which garantees padding between paths.
- `Map` (modified)
    - No `self.__layout`, `self.__cell` anymore.
    - Cost of walkable pixels is now 0, and cost of obstacle pixels is 1.
    - Other value can be assigned to pixels of different meaning, like padded routing paths.
    - New methods:
        - `mesh_X`, `mesh_Y`: to get the coordinates of the mesh.
        - `where`: to get the coordinates of pixels with a specific cost.
        - `where_not`/`where_less_than`/`where_greater_than`: to get coordinates of pixels with specific conditions.   
        - `get_coor_list`: return list of all coordinates in the map.
        - `create_pixel_cell`: to be implemented.
    - Modified methods:
        - `self.map` -> `self.__map`
        - `rasterize_rect_region` name-changed to `rasterize_rectangle`.
        - `rasterize_region` name-changed to `rasterize_polygon`.
        - `rasterize_region` now takes a real `pya.Region` as input, containing multiple polygons.
        - `rasterize_region`: use `pya.Region.rectangles()/non_rectangles()` to extract rectangles from a region.      
        - `visualize`: Draw pixels with certain cost on layer f"Mesh{cost}", and use `self.where()` to get coordinates.
- `GraphM` (modified)
    - Use `self.__Map.where()` to get coordinates of pixels. More efficient.
- `PinMatcher` (modified)
    - No `self.__layout`, `self.__cell` anymore.
    - Modified methods:
        - all four different `get_pins_from_xxx` methods are unified into one `get_pins_from_shapes`.
- Method_test.py (new): A script to test methods in Klayout.