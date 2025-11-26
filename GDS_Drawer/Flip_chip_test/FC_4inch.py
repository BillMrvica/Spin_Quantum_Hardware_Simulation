import gdstk
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch, Circle

# --- 详细版图可视化函数 (无变化) ---
def plot_cell_details(cell, layer_map, title="Layout Preview"):
    fig, ax = plt.subplots(figsize=(14, 14))
    color_map = {0: 'pink', 1: 'lightcyan', 2: 'darkblue', 3: 'lightgrey', 4: '#A7B5CC', 5: '#E5E4E2', 6: '#FFD700', 7: '#708090'}
    z_order = {0: -1, 1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7:6}
    polygons = cell.get_polygons()
    for poly in polygons:
        color, z = color_map.get(poly.layer, 'black'), z_order.get(poly.layer, 0)
        patch = Polygon(poly.points, closed=True, facecolor=color, edgecolor='dimgray', linewidth=0.5, zorder=z)
        ax.add_patch(patch)
    legend_elements = [Patch(facecolor=color_map[layer], edgecolor='dimgray', label=name)
                       for layer, name in sorted(layer_map.items()) if layer in color_map]
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_aspect('equal')
    bbox = cell.bounding_box()
    if bbox:
        x_margin, y_margin = (bbox[1][0] - bbox[0][0]) * 0.05, (bbox[1][1] - bbox[0][1]) * 0.05
        ax.set_xlim(bbox[0][0] - x_margin, bbox[1][0] + x_margin)
        ax.set_ylim(bbox[0][1] - y_margin, bbox[1][1] + y_margin)
    ax.set_xlabel("X coordinate (μm)"); ax.set_ylabel("Y coordinate (μm)")
    ax.set_title(title); ax.grid(True, linestyle='--', alpha=0.6)
    print(f"\n正在显示详细版图预览: {title}")
    plt.show()

# --- 晶圆宏观布局示意图函数 (已更新，可显示安全区) ---
def plot_wafer_schematic(wafer_cell, wafer_diameter_um, exclusion_um):
    fig, ax = plt.subplots(figsize=(10, 10))
    # 绘制晶圆物理轮廓
    wafer_outline_poly = [p for p in wafer_cell.get_polygons() if p.layer == 0][0]
    ax.add_patch(Polygon(wafer_outline_poly.points, facecolor='pink', edgecolor='red', alpha=0.3))
    # **【修改处】**: 绘制安全区边界
    safe_zone = Circle((0, 0), wafer_diameter_um / 2 - exclusion_um, facecolor='none', edgecolor='green', linestyle='--', linewidth=1.5, label='Safe Zone')
    ax.add_patch(safe_zone)
    
    for ref in wafer_cell.references:
        bbox = ref.bounding_box()
        width, height = bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1]
        rect = plt.Rectangle(bbox[0], width, height, facecolor='lightblue', edgecolor='blue', alpha=0.6)
        ax.add_patch(rect)
        ax.text(bbox[0][0] + width/2, bbox[0][1] + height/2, ref.cell.name.replace("Unit", ""), ha='center', va='center', fontsize=6)
    
    ax.legend()
    ax.set_aspect('equal')
    ax.set_xlim(-wafer_diameter_um/2 * 1.05, wafer_diameter_um/2 * 1.05)
    ax.set_ylim(-wafer_diameter_um/2 * 1.05, wafer_diameter_um/2 * 1.05)
    ax.set_xlabel("X coordinate (μm)"); ax.set_ylabel("Y coordinate (μm)")
    ax.set_title("Wafer-Level Schematic Layout"); ax.grid(True, linestyle='--', alpha=0.6)
    print(f"\n正在显示晶圆宏观布局示意图...")
    plt.show()

# --- Marker创建函数 (无变化) ---
def create_alignment_marker(lib, layers, size=400, width=80, marker_type='positive'):
    marker_cell = lib.new_cell(f"AlignmentMarker_{marker_type}")
    l, w = size / 2, width / 2
    cross = gdstk.cross((0, 0), l, w, layer=layers['marker'])
    marker_cell.add(cross)
    tick_width = 2; tick_lengths = {5: 15, 10: 20, 50: 25}; start_offset = w + tick_lengths[50] + 10
    for i in range(int(start_offset), int(l) + 1):
        tick_len = 0
        if i % 50 == 0: tick_len = tick_lengths[50]
        elif i % 10 == 0: tick_len = tick_lengths[10]
        elif i % 5 == 0: tick_len = tick_lengths[5]
        if tick_len > 0:
            if marker_type == 'positive':
                marker_cell.add(gdstk.rectangle((i - tick_width / 2, w), (i + tick_width / 2, w + tick_len), layer=layers['marker']), gdstk.rectangle((i - tick_width / 2, -w - tick_len), (i + tick_width / 2, -w), layer=layers['marker']), gdstk.rectangle((-i - tick_width / 2, w), (-i + tick_width / 2, w + tick_len), layer=layers['marker']), gdstk.rectangle((-i - tick_width / 2, -w - tick_len), (-i + tick_width / 2, -w), layer=layers['marker']), gdstk.rectangle((w, i - tick_width / 2), (w + tick_len, i + tick_width / 2), layer=layers['marker']), gdstk.rectangle((-w - tick_len, i - tick_width / 2), (-w, i + tick_width / 2), layer=layers['marker']), gdstk.rectangle((w, -i - tick_width / 2), (w + tick_len, -i + tick_width / 2), layer=layers['marker']), gdstk.rectangle((-w - tick_len, -i - tick_width / 2), (-w, -i + tick_width / 2), layer=layers['marker']))
            else:
                marker_cell.add(gdstk.rectangle((i - tick_width / 2, w - tick_len), (i + tick_width / 2, w), layer=layers['marker']), gdstk.rectangle((i - tick_width / 2, -w), (i + tick_width / 2, -w + tick_len), layer=layers['marker']), gdstk.rectangle((-i - tick_width / 2, w - tick_len), (-i + tick_width / 2, w), layer=layers['marker']), gdstk.rectangle((-i - tick_width / 2, -w), (-i + tick_width / 2, -w + tick_len), layer=layers['marker']), gdstk.rectangle((w - tick_len, i - tick_width / 2), (w, i + tick_width / 2), layer=layers['marker']), gdstk.rectangle((-w, i - tick_width / 2), (-w + tick_len, i + tick_width / 2), layer=layers['marker']), gdstk.rectangle((w - tick_len, -i - tick_width / 2), (w, -i + tick_width / 2), layer=layers['marker']), gdstk.rectangle((-w, -i - tick_width / 2), (-w + tick_len, -i + tick_width / 2), layer=layers['marker']))
    return marker_cell

# --- 独立焊盘单元生成函数 (无变化) ---
def create_carrier_chip_with_individual_pads(input_gds, lib, layers, marker_cell, marker_positions):
    print("正在生成单元: 独立焊盘载板")
    wb_pad_size = (100, 100); clearance = 700; trace_width = 50; bump_diameter = 20.0; sio2_margin = 50.0; INPUT_LAYER = 86
    carrier_cell = lib.new_cell('IndividualPadUnit')
    original_cell = gdstk.read_gds(input_gds).top_level()[0]
    polygons = [p for p in original_cell.get_polygons() if p.layer == INPUT_LAYER]
    if not polygons: raise ValueError("在输入层上未找到多边形。")
    overall_bbox = original_cell.bounding_box()
    center_x, center_y = (overall_bbox[0][0] + overall_bbox[1][0]) / 2, (overall_bbox[0][1] + overall_bbox[1][1]) / 2
    for poly in polygons:
        pad_bbox = poly.bounding_box(); pad_center = ((pad_bbox[0][0] + pad_bbox[1][0]) / 2, (pad_bbox[0][1] + pad_bbox[1][1]) / 2)
        metal_pad = poly.copy(); metal_pad.layer = layers['metal']; carrier_cell.add(metal_pad)
        dx, dy = pad_center[0] - center_x, pad_center[1] - center_y
        if abs(dy) > abs(dx): wb_pad_center = (pad_center[0], pad_center[1] + np.sign(dy) * (clearance + wb_pad_size[1] / 2 + (pad_bbox[1][1]-pad_bbox[0][1]) / 2))
        else: wb_pad_center = (pad_center[0] + np.sign(dx) * (clearance + wb_pad_size[0] / 2 + (pad_bbox[1][0]-pad_bbox[0][0]) / 2), pad_center[1])
        wb_pad = gdstk.rectangle((wb_pad_center[0] - wb_pad_size[0] / 2, wb_pad_center[1] - wb_pad_size[1] / 2), (wb_pad_center[0] + wb_pad_size[0] / 2, wb_pad_center[1] + wb_pad_size[1] / 2), layer=layers['metal'])
        carrier_cell.add(wb_pad, gdstk.FlexPath([pad_center, wb_pad_center], trace_width, layer=layers['metal']))
        ubm_ti = poly.copy(); ubm_ti.layer = layers['ubm_ti']; carrier_cell.add(ubm_ti); ubm_pt = poly.copy(); ubm_pt.layer = layers['ubm_pt']; carrier_cell.add(ubm_pt); ubm_au = poly.copy(); ubm_au.layer = layers['ubm_au']; carrier_cell.add(ubm_au)
        carrier_cell.add(gdstk.regular_polygon(pad_center, bump_diameter / 2.0, 16, layer=layers['bump']))
    for pos in marker_positions: carrier_cell.add(gdstk.Reference(marker_cell, pos))
    final_bbox = carrier_cell.bounding_box()
    carrier_cell.add(gdstk.rectangle((final_bbox[0][0] - sio2_margin, final_bbox[0][1] - sio2_margin), (final_bbox[1][0] + sio2_margin, final_bbox[1][1] + sio2_margin), layer=layers['sio2']))
    return carrier_cell

# --- 两两相连焊盘单元生成函数 (无变化) ---
def create_carrier_chip_with_paired_pads(input_gds, lib, layers, marker_cell, marker_positions):
    print("正在生成单元: 两两相连焊盘载板")
    trace_width = 20.0; bump_diameter = 20.0; sio2_margin = 50.0; INPUT_LAYER = 86
    paired_pad_cell = lib.new_cell('PairedPadUnit')
    original_cell = gdstk.read_gds(input_gds).top_level()[0]
    polygons = [p for p in original_cell.get_polygons() if p.layer == INPUT_LAYER]
    if not polygons: raise ValueError("在输入层上未找到多边形。")
    center_x, center_y = np.vstack([p.points for p in polygons]).mean(axis=0)
    pads_with_info = sorted([{'polygon': p, 'center': ((p.bounding_box()[0][0] + p.bounding_box()[1][0]) / 2, (p.bounding_box()[0][1] + p.bounding_box()[1][1]) / 2), 'angle': np.arctan2(((p.bounding_box()[0][1] + p.bounding_box()[1][1]) / 2) - center_y, ((p.bounding_box()[0][0] + p.bounding_box()[1][0]) / 2) - center_x)} for p in polygons], key=lambda item: item['angle'])
    for i in range(0, len(pads_with_info), 2):
        if i + 1 >= len(pads_with_info): continue
        pad1, pad2 = pads_with_info[i], pads_with_info[i+1]
        pair_metal = gdstk.boolean([pad1['polygon'], pad2['polygon']], gdstk.FlexPath([pad1['center'], pad2['center']], trace_width, layer=layers['metal']), 'or', layer=layers['metal'])
        paired_pad_cell.add(*pair_metal)
    for pad_info in pads_with_info:
        ubm_ti = pad_info['polygon'].copy(); ubm_ti.layer = layers['ubm_ti']; paired_pad_cell.add(ubm_ti); ubm_pt = pad_info['polygon'].copy(); ubm_pt.layer = layers['ubm_pt']; paired_pad_cell.add(ubm_pt); ubm_au = pad_info['polygon'].copy(); ubm_au.layer = layers['ubm_au']; paired_pad_cell.add(ubm_au)
        paired_pad_cell.add(gdstk.regular_polygon(pad_info['center'], bump_diameter / 2.0, 16, layer=layers['bump']))
    for pos in marker_positions: paired_pad_cell.add(gdstk.Reference(marker_cell, pos))
    final_bbox = paired_pad_cell.bounding_box()
    paired_pad_cell.add(gdstk.rectangle((final_bbox[0][0] - sio2_margin, final_bbox[0][1] - sio2_margin), (final_bbox[1][0] + sio2_margin, final_bbox[1][1] + sio2_margin), layer=layers['sio2']))
    return paired_pad_cell

# --- 主程序 (已更新晶圆生成和芯片裁剪逻辑) ---
if __name__ == "__main__":
    input_filename = "C:\\Users\\Administrator\\Desktop\\Spin_Quantum_Hardware_Simulation\\GDS_Drawer\\Flip_chip_test\\PAD_ring_Xiling.gds"
    if not os.path.exists(input_filename):
         print(f"错误: 输入文件 '{input_filename}' 不存在。")
    else:
        WAFER_DIAMETER_UM = 100000
        WAFER_RADIUS_UM = WAFER_DIAMETER_UM / 2
        EMPTY_WIDTH_UM = 15000
        CHIP_GAP_UM = 1000
        # **【新增参数】**: 定义晶圆边缘的安全距离
        WAFER_EDGE_EXCLUSION_UM = 5000 # 5mm

        layer_map = {0: 'Wafer Outline', 1: 'Passivation - SiO2', 2: 'Marker - Ti/Au', 3: 'Carrier Metal', 4: 'UBM - Ti', 5: 'UBM - Pt', 6: 'UBM - Au', 7: 'Bump'}
        layers = {'wafer_outline': 0, 'sio2': 1, 'marker': 2, 'metal': 3, 'ubm_ti': 4, 'ubm_pt': 5, 'ubm_au': 6, 'bump': 7}

        main_lib = gdstk.Library(unit=1e-6, precision=1e-9)
        
        marker_positive = create_alignment_marker(main_lib, layers, marker_type='positive')
        marker_negative = create_alignment_marker(main_lib, layers, marker_type='negative')

        original_cell = gdstk.read_gds(input_filename).top_level()[0]
        original_bbox = original_cell.bounding_box()
        marker_offset = 70
        marker_positions = [
            (original_bbox[0][0] - marker_offset, original_bbox[0][1] - marker_offset), (original_bbox[1][0] + marker_offset, original_bbox[0][1] - marker_offset),
            (original_bbox[0][0] - marker_offset, original_bbox[1][1] + marker_offset), (original_bbox[1][0] + marker_offset, original_bbox[1][1] + marker_offset)
        ]
        
        individual_pad_cell = create_carrier_chip_with_individual_pads(input_filename, main_lib, layers, marker_positive, marker_positions)
        paired_pad_cell = create_carrier_chip_with_paired_pads(input_filename, main_lib, layers, marker_negative, marker_positions)
            
        wafer_cell = main_lib.new_cell('Wafer_Layout_Clipped_With_Flat')
        
        flat_length_um = 32.5 * 1000
        flat_y_pos = -math.sqrt(WAFER_RADIUS_UM**2 - (flat_length_um / 2)**2)
        print(f"晶圆平边Y轴位置: {flat_y_pos:.2f} um")
        wafer_circle = gdstk.ellipse((0, 0), WAFER_RADIUS_UM, layer=layers['wafer_outline'])
        cutter = gdstk.rectangle((-WAFER_RADIUS_UM, -WAFER_RADIUS_UM), (WAFER_RADIUS_UM, flat_y_pos))
        wafer_shape = gdstk.boolean(wafer_circle, cutter, "not", layer=layers['wafer_outline'])
        wafer_cell.add(*wafer_shape)

        # **【修改处】**: 定义一个更严格的检查函数
        def is_chip_fully_inside(chip_origin, chip_width, chip_height):
            safe_radius_sq = (WAFER_RADIUS_UM - WAFER_EDGE_EXCLUSION_UM)**2
            safe_flat_y = flat_y_pos + WAFER_EDGE_EXCLUSION_UM
            corners = [
                chip_origin, (chip_origin[0] + chip_width, chip_origin[1]),
                (chip_origin[0], chip_origin[1] + chip_height), (chip_origin[0] + chip_width, chip_origin[1] + chip_height)
            ]
            for x, y in corners:
                if not (x**2 + y**2 < safe_radius_sq and y > safe_flat_y):
                    return False
            return True

        # 上半部分
        bbox_ind = individual_pad_cell.bounding_box()
        chip_w_ind, chip_h_ind = bbox_ind[1][0] - bbox_ind[0][0], bbox_ind[1][1] - bbox_ind[0][1]
        pitch_x_ind, pitch_y_ind = chip_w_ind + CHIP_GAP_UM, chip_h_ind + CHIP_GAP_UM
        num_cols_ind = int(WAFER_DIAMETER_UM / pitch_x_ind)
        num_rows_ind = int((WAFER_RADIUS_UM - EMPTY_WIDTH_UM / 2) / pitch_y_ind)
        for r in range(num_rows_ind):
            for c in range(-num_cols_ind // 2, num_cols_ind // 2 + 1):
                origin = (c * pitch_x_ind, EMPTY_WIDTH_UM / 2 + r * pitch_y_ind)
                if is_chip_fully_inside(origin, chip_w_ind, chip_h_ind):
                    wafer_cell.add(gdstk.Reference(individual_pad_cell, origin))

        # 下半部分
        bbox_pair = paired_pad_cell.bounding_box()
        chip_w_pair, chip_h_pair = bbox_pair[1][0] - bbox_pair[0][0], bbox_pair[1][1] - bbox_pair[0][1]
        pitch_x_pair, pitch_y_pair = chip_w_pair + CHIP_GAP_UM, chip_h_pair + CHIP_GAP_UM
        num_cols_pair = int(WAFER_DIAMETER_UM / pitch_x_pair)
        num_rows_pair = int((WAFER_RADIUS_UM - EMPTY_WIDTH_UM / 2) / pitch_y_pair)
        for r in range(num_rows_pair):
            for c in range(-num_cols_pair // 2, num_cols_pair // 2 + 1):
                origin = (c * pitch_x_pair, -EMPTY_WIDTH_UM / 2 - (r + 1) * pitch_y_pair)
                if is_chip_fully_inside(origin, chip_w_pair, chip_h_pair):
                    wafer_cell.add(gdstk.Reference(paired_pad_cell, origin))

        output_filename = r"\\10.20.42.18\SEQSI_main\Public_Resources\GDS Design\Flip-Chip Test Structure\\FC_test_stru_4_inch_wafer_20251125.gds"
        main_lib.write_gds(output_filename)
        print(f"\n成功生成GDS文件: {output_filename}")
        
        # 可视化
        plot_wafer_schematic(wafer_cell, WAFER_DIAMETER_UM, WAFER_EDGE_EXCLUSION_UM)
        plot_cell_details(individual_pad_cell.flatten(), layer_map, title="Detailed View: Individual Pad Unit")
        plot_cell_details(paired_pad_cell.flatten(), layer_map, title="Detailed View: Paired Pad Unit")