# Klayout-Router Technical Manual

## 1. Introduction

**Klayout-Router** is a specialized automation toolkit designed for **laboratory-level nano-device design** in Condensed Matter Physics (CMP). Unlike industrial EDA tools that focus on standard cell routing, this project addresses the highly customized, flexible, and often "imperfect" nature of experimental chip design.

This document serves as a comprehensive guide to the internal logic, key algorithms, and practical usage of the router.

## 2. Design Philosophy & Logic

The core of the router is based on a **Cost-Map** approach, transforming the continuous vector geometry of a GDS layout into a discrete grid (rasterization) where pathfinding algorithms can operate.

### 2.1 The Map System (`Map.py`)
The layout is discretized into a 2D numpy matrix (`self.__map`).
- **Resolution**: Defined by the user (e.g., 5µm). A smaller resolution offers finer control but increases memory usage and computation time.
- **Coordinate System**: The router maintains a strict conversion between:
    - **Layout Coordinates**: Continuous (x, y) in database units (dbu).
    - **Map Coordinates**: Discrete (row, col) indices in the matrix.

### 2.2 Rasterization
Instead of using slow, pixel-by-pixel checks, the router utilizes KLayout's native C++ engine (`kdb.Region.rasterize`) to convert complex polygons into grid masks. This allows for extremely fast "conditional overwriting" of the map—e.g., "mark these pixels as obstacles only if they are currently empty."

### 2.3 The Cost Function
The map is not just binary (0 for obstacle, 1 for free). It uses a **Weighted Cost System**:
- **1**: Free space (lowest cost).
- **>1**: "Soft" obstacles or buffer zones (high cost).
- **<0**: Absolute obstacles (impassable).

### 2.4 Routing algorithm

`Krouter` will do *Dijkstra algorithm* between every pair of pins on the `Map`, which is a sparse graph. There are also a lot of designs to arrange all paths on the map, so that they:

- do not cross each other
- keep away from each other
- keep away from obstacles and other pins

See [key features](#4.-Key-Features) and the code for more details.

## 3. Architecture Overview

### `Krouter` (The Controller)

- **Role**: The boss. Handles all interaction with the layout files.
- **Responsibilities**:
  - Reads `.gds` files.
  - Extracts shapes from specified layers.
  - Coordinates `Map` and `PinMatcher`.
  - Converts physical units (µm) to database units (dbu).
  - Executes the routing loop and saves the result.

### `Map` (The Terrain)

- **Role**: The grid manager.
- **Responsibilities**:
  - Stores the `numpy` cost matrix.
  - Handles `Layout <-> Map` coordinate conversion.
  - Performs `rasterize` operations to "paint" obstacles and costs onto the grid.

### `PinMatcher` (The Strategist)

- **Role**: Manage pins and match them.
- **Responsibilities**:
  - Identifies start and end points from shapes.
  - Uses the **Hungarian Algorithm** (Linear Sum Assignment) to find the optimal pairing of pins that minimizes total path length.
  - Performs Clockwise Sorting or Distance-ascending sorting.

## 4. Key Features

### Soft Constraints & Damping
A common issue in routing is the "No Path Found" error when constraints are too strict. Klayout-Router solves this with **Soft Constraints**:
- Instead of strictly forbidding paths near obstacles, we assign a **High Cost** to these areas.
- **Damping**: The cost doesn't drop instantly from "High" to "Low". It decays gradually (like a halo) around obstacles and pins.
- **Result**: The router *will* find a path if one exists, even if it has to squeeze through a "dangerous" narrow channel, but it will naturally prefer safer, wider routes.

### Distance-ascending Sorting Pins
sort all pairs of pins by their distances.

This will define the order by which the routing is done. Firstly finding the path that connects pins with shorter distance is beneficial, because this path is less like to affect other paths.

### Clock-wise Sorting Pins

Sort all pairs of pins clockwise around the center.

Sometimes this order is better than distance-ascending sorting.

### Self-Adaptive Routing (Experimental)
To prevent paths from crowding together:
1.  Route all paths.
2.  Calculate a "Path Density Map" based on the result.
3.  Feed this density back into the cost map as a penalty.
4.  Re-route. Paths will now "repel" each other, spreading out evenly across the available space.

## Troubleshooting & Tuning

**Problem: "No path found" or paths are crossing.**

- **Cause**: `resolution` might be too coarse, blocking narrow passages.
- **Fix**: Decrease `resolution` (e.g., from 5 to 2).
- **Fix**: Decrease `obs_safe_distance` to allow paths to squeeze through tighter spots.

**Problem: Paths are hugging obstacles too closely.**
- **Cause**: `obs_hardness` is too low.
- **Fix**: Increase `obs_hardness` or `obs_damping_step` to create a stronger "repulsion" field.

**Problem: Paths are crossing other pins.**
- **Cause**: The router doesn't see other pins as obstacles.
- **Fix**: Ensure `pin_hardness` is high (e.g., >20).

**Problem: Routing is too slow.**
- **Cause**: Map is too large or resolution is too fine.
- **Fix**: Increase `resolution` or reduce the `BBox` area.
