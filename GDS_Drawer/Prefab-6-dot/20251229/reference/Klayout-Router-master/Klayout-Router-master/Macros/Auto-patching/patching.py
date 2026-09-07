import pya

def create_grid_from_shapes(layout, cell, shape_layer, grid_line_width, 
                            shape_layer_datatype=0,):
    """
    Create a grid with `grid_line_width` along the edges of all shapes on `shape_layer`.
    """
    # Get the database unit (dbu) for unit conversion
    dbu = layout.dbu
    grid_line_width_dbu = int(grid_line_width / dbu + 0.5)

    # Get layer indices
    shape_layer_index = layout.layer(pya.LayerInfo(shape_layer, shape_layer_datatype))

    # For every shape, find its edges and extend them to create grid lines, 
    # and finally combine all grid regions.
    final_grid_region = pya.Region()
    for shape in cell.shapes(shape_layer_index):
        edges = pya.Region(shape.polygon).edges() # standard way to get edges from any shapes
        grid_region = edges.extended(0, 0, grid_line_width_dbu//2, grid_line_width_dbu//2, True)
        final_grid_region += grid_region
    
    return final_grid_region


def create_patch(layout, cell, grid_region, electrode_layer, patch_layer, patch_size,
                 electrode_layer_datatype=0, patch_layer_datatype=0):
    """
    Create patches at the center of intersections between electrode and grid layers.
    Datatype for all layers must be 0.
    """
    # Get the database unit (dbu) for unit conversion
    dbu = layout.dbu
    patch_size_dbu = int(patch_size / dbu + 0.5)

    # Get layer indices
    electrode_layer_index = layout.layer(pya.LayerInfo(electrode_layer, electrode_layer_datatype))
    patch_layer_index = layout.layer(pya.LayerInfo(patch_layer, patch_layer_datatype))

    # Create Regions for boolean operations
    electrode_regioin = pya.Region(cell.shapes(electrode_layer_index))

    # Find intersections
    intersection_region = electrode_regioin & grid_region
    intersection_region.merge()

    # Check if there are any intersections
    if intersection_region.is_empty():
        print("Warning: No intersections found between the electrode layer and the grid layer.")
        return
    
    # Intersection region is actually the patches, 
    # but for safety, we additionally create boxes at the center of each intersection shape
    for intersection_poly in intersection_region.each():
        center = intersection_poly.bbox().center()
        p1 = pya.Point(center.x - patch_size_dbu//2, center.y - patch_size_dbu//2)
        p2 = pya.Point(center.x + patch_size_dbu//2, center.y + patch_size_dbu//2)
        cell.shapes(patch_layer_index).insert(pya.Box(p1, p2))

    # Finally, insert the patches into the layout
    cell.shapes(patch_layer_index).insert(intersection_region)

    print(f"Patches created on layer {pya.LayerInfo(patch_layer, patch_layer_datatype)}.")
    
    return

