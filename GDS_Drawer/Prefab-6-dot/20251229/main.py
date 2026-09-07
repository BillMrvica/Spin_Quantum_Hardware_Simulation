import gdstk
import numpy as np
import klayout.db as kdb # Import klayout.db
from typing import List # Import List
from quantum_dot_generator import QuantumDotGenerator
from pad_lead_generator import PadLeadGenerator
from router import Router # Keep for now, will remove after RouteV2 is integrated
from route_v2 import RouteV2, ObsRule # Import the new router class and ObsRule TypedDict
from utils import plot_gds, visualize_and_save_pad_frame, visualize_and_save_assembly

if __name__ == "__main__":
    ACTIVE_SIZE = 180
    ROUTING_LAYER = 5
    PAD_CONNECTION_EXTENSION_LEN = 10.0 # User specified 10um

    # Initialize generators
    qd_gen = QuantumDotGenerator()
    pad_gen = PadLeadGenerator()

    # Step 1: Generate Quantum Dot Device with Leads
    lib_device, cell_device, device_points = qd_gen.create_6qd_layout_with_labels()
    plot_gds(lib_device, cell_device, title="步骤 1: 仅带长引线的量子点器件", plot_flag=False, save_gds_flag=True, gds_filename="step1_device_with_leads.gds") # Pass lib_device and gds_filename
    # lib_device.write_gds("step1_device_with_leads.gds") # Removed as it's now handled by plot_gds

    # Step 2: Generate Pad Frame with Wires
    lib_pads, cell_pads, pad_points, active_center, all_pads_info = pad_gen.create_rect_wire_layout(
        active_size=ACTIVE_SIZE, edge_margin=50, active_entry_len=40
    )
    visualize_and_save_pad_frame(lib_pads, cell_pads, "步骤 2: 仅带引线的焊盘框架", "step2_pad_frame.gds", plot_flag=False, save_gds_flag=True)

    # Prepare label to side info for routing
    label_to_side_info = {}
    pads_by_edge = {'top': [], 'bottom': [], 'left': [], 'right': []}

    for pad_info in all_pads_info:
        pads_by_edge[pad_info['edge']].append(pad_info)

    for side, pads in pads_by_edge.items():
        if side in ['top', 'bottom']:
            pads.sort(key=lambda p: p['center'][0])
        else: # 'left', 'right'
            pads.sort(key=lambda p: p['center'][1], reverse=True)
        
        for index, pad_info in enumerate(pads):
            label_to_side_info[pad_info['label']] = (side, index, len(pads))

    # Translate device to active_center
    bbox = cell_device.bounding_box()
    if bbox:
        (min_x, min_y), (max_x, max_y) = bbox
        cx, cy = (min_x + max_x) / 2, (min_y + max_y) / 2
        for poly in cell_device.polygons:
            poly.translate(-cx, -cy)
        for label in cell_device.labels:
            label.origin = (label.origin[0] - cx, label.origin[1] - cy)
        for label, pos in device_points.items():
            device_points[label] = (pos[0] - cx, pos[1] - cy)
    
    # Assemble device and pads
    assembly_lib = gdstk.Library()
    assembly_cell = assembly_lib.new_cell('ASSEMBLY_CELL')
    assembly_cell.add(gdstk.Reference(cell_pads))
    device_ref = gdstk.Reference(cell_device, active_center)
    assembly_cell.add(device_ref)
    
    visualize_and_save_assembly(assembly_lib, assembly_cell, "步骤 3: 器件已放置 (未布线)", "step3_placed_unrouted.gds", plot_flag=False, save_gds_flag=True)
    
    print("--- 步骤 4 & 5: 正在使用全新 Router 连接所有线... ---")

    # Extract active region bounding box
    active_region_bbox = (
        (active_center[0] - ACTIVE_SIZE / 2, active_center[1] - ACTIVE_SIZE / 2),
        (active_center[0] + ACTIVE_SIZE / 2, active_center[1] + ACTIVE_SIZE / 2)
    )

    # Extract quantum dot polygons from cell_device
    quantum_dot_polygons = []
    for layer in [0, 2, 3]: # S/D, BG, PG layers
        quantum_dot_polygons.extend(cell_device.get_polygons(layer=layer, datatype=0))
    
    # Translate quantum dot polygons to their global position
    translated_qd_polygons = []
    for poly in quantum_dot_polygons:
        translated_poly = gdstk.Polygon(poly.points) # Create a copy
        translated_poly.translate(active_center[0], active_center[1])
        translated_qd_polygons.append(translated_poly)

    # Add a dummy rectangle for total_area_bbox_layer to assembly_cell
    # This ensures the RouteV2 class has a valid bounding box to work with.
    bbox_rect = gdstk.rectangle(
        corner1=active_region_bbox[0],
        corner2=active_region_bbox[1],
        layer=100, # Use layer 100 for the total area bounding box
        datatype=0
    )
    assembly_cell.add(bbox_rect)
    
    # Re-save step3_placed_unrouted.gds after adding the bbox_rect
    visualize_and_save_assembly(assembly_lib, assembly_cell, "步骤 3: 器件已放置 (未布线)", "step3_placed_unrouted.gds", plot_flag=False, save_gds_flag=True)

    # Initialize RouteV2
    # Define obstacle rules (example, adjust as needed)
    obs_rules: List[ObsRule] = [
        {"layers": "0/0", "bbx": False, "pad": 0.01}, # Reduced obstacle padding to 10nm
        {"layers": "2/0", "bbx": False, "pad": 0.01}, # Reduced obstacle padding to 10nm
        {"layers": "3/0", "bbx": False, "pad": 0.01}, # Reduced obstacle padding to 10nm
    ]

    router_v2 = RouteV2(
        file_path="step3_placed_unrouted.gds", # Input GDS file
        cell_name='ASSEMBLY_CELL', # The cell containing the device and pads
        total_area_bbox_layer="100/0", # A layer that defines the total routing area (e.g., a dummy layer)
        quantum_dot_layer="0/0", # Layer for quantum dot electrodes (Pin A)
        pad_lead_layer="1/0", # Layer for pad leads (Pin B)
        map_resolution=0.1, # Adjusted resolution to 0.1um for better performance
        obs_rules=obs_rules,
        routing_layer=str(ROUTING_LAYER) + "/0", # Layer for the new routes
        initial_extension_distance=0.5 # Reverted initial extension distance
    )

    # Perform routing
    router_v2.route_paths(
        obs_safe_distance=0.2, # Increased safe distance around obstacles
        obs_hardness=50, # Reduced obstacle hardness
        obs_damping_step=10,
        pin_safe_distance=0.1, # Increased safe distance around pins
        pin_hardness=20, # Reduced pin hardness
        pin_damping_step=5,
        path_safe_distance=0.1, # Increased safe distance between paths
        path_hardness=100, # Reduced path hardness
        path_damping_step=10,
        path_density_hardness=50, # Reduced path density hardness
        round=1
    )

    # Save the routed layout
    router_v2.save_layout("step5_fully_routed_v2.gds")
    
    # Add the generated paths to the assembly cell for visualization
    # Note: RouteV2 directly inserts paths into its internal layout.
    # To visualize with gdstk, we need to re-read the GDS or extract paths.
    # For simplicity, let's assume the save_layout creates the final GDS.
    # If direct addition to assembly_cell is needed, RouteV2 would need a method to return its generated paths.
    
    # For now, we'll just rely on the saved GDS for visualization.
    # visualize_and_save_assembly(assembly_lib, assembly_cell, "步骤 5: 最终完整布线", "step5_fully_routed_v2.gds", plot_flag=False, save_gds_flag=True)
    print("--- 步骤 4 & 5: 使用 RouteV2 完成布线 ---")
    print("布线结果已保存到 step5_fully_routed_v2.gds")

    print("全部分步流程执行完毕！")
