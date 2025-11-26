import gdstk
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch

# --- 可视化函数 (无变化) ---
def plot_cell(cell, layer_map, title="Layout Preview"):
    fig, ax = plt.subplots(figsize=(12, 12)) # 稍微增大画布以适应更大的版图
    color_map = {1: 'lightcyan', 2: 'lightgrey', 3: '#A7B5CC', 4: '#E5E4E2', 5: '#FFD700', 6: '#708090'}
    z_order = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5}
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
        x_margin, y_margin = (bbox[1][0] - bbox[0][0]) * 0.1, (bbox[1][1] - bbox[0][1]) * 0.1
        ax.set_xlim(bbox[0][0] - x_margin, bbox[1][0] + x_margin)
        ax.set_ylim(bbox[0][1] - y_margin, bbox[1][1] + y_margin)
    ax.set_xlabel("X coordinate (μm)"); ax.set_ylabel("Y coordinate (μm)")
    ax.set_title(title); ax.grid(True, linestyle='--', alpha=0.6)
    print(f"\n正在显示版图预览: {title}")
    plt.show()

# --- Marker创建函数 (无变化) ---
def create_alignment_marker(lib, layers, size=400, width=80):
    marker_cell = lib.new_cell("AlignmentMarker")
    cross = gdstk.cross((0, 0), size / 2, width / 2, layer=layers['metal'])
    marker_cell.add(cross)
    return marker_cell

# --- 独立焊盘单元生成函数 (无变化) ---
def create_carrier_chip_with_individual_pads(input_gds, lib, layers, marker_cell, marker_positions):
    print("正在生成版图: 独立焊盘载板 + Marker")
    wb_pad_size = (100, 100); clearance = 300; trace_width = 50; bump_diameter = 20.0; sio2_margin = 50.0; INPUT_LAYER = 86
    carrier_cell = lib.new_cell('IndividualPadUnit')
    original_cell = gdstk.read_gds(input_gds).top_level()[0]
    polygons = [p for p in original_cell.get_polygons() if p.layer == INPUT_LAYER]
    if not polygons: raise ValueError("在输入层上未找到多边形。")
    
    overall_bbox = original_cell.bounding_box()
    center_x, center_y = (overall_bbox[0][0] + overall_bbox[1][0]) / 2, (overall_bbox[0][1] + overall_bbox[1][1]) / 2

    for poly in polygons:
        pad_bbox = poly.bounding_box()
        pad_center = ((pad_bbox[0][0] + pad_bbox[1][0]) / 2, (pad_bbox[0][1] + pad_bbox[1][1]) / 2)
        metal_pad = poly.copy(); metal_pad.layer = layers['metal']; carrier_cell.add(metal_pad)
        dx, dy = pad_center[0] - center_x, pad_center[1] - center_y
        if abs(dy) > abs(dx):
            offset = np.sign(dy) * (clearance + wb_pad_size[1] / 2 + (pad_bbox[1][1]-pad_bbox[0][1]) / 2)
            wb_pad_center = (pad_center[0], pad_center[1] + offset)
        else:
            offset = np.sign(dx) * (clearance + wb_pad_size[0] / 2 + (pad_bbox[1][0]-pad_bbox[0][0]) / 2)
            wb_pad_center = (pad_center[0] + offset, pad_center[1])
        wb_pad = gdstk.rectangle((wb_pad_center[0] - wb_pad_size[0] / 2, wb_pad_center[1] - wb_pad_size[1] / 2),
                                 (wb_pad_center[0] + wb_pad_size[0] / 2, wb_pad_center[1] + wb_pad_size[1] / 2), layer=layers['metal'])
        carrier_cell.add(wb_pad, gdstk.FlexPath([pad_center, wb_pad_center], trace_width, layer=layers['metal']))
        ubm_ti = poly.copy(); ubm_ti.layer = layers['ubm_ti']; carrier_cell.add(ubm_ti)
        ubm_pt = poly.copy(); ubm_pt.layer = layers['ubm_pt']; carrier_cell.add(ubm_pt)
        ubm_au = poly.copy(); ubm_au.layer = layers['ubm_au']; carrier_cell.add(ubm_au)
        carrier_cell.add(gdstk.regular_polygon(pad_center, bump_diameter / 2.0, 16, layer=layers['bump']))

    for pos in marker_positions:
        carrier_cell.add(gdstk.Reference(marker_cell, pos))

    final_bbox = carrier_cell.bounding_box()
    sio2_layer = gdstk.rectangle((final_bbox[0][0] - sio2_margin, final_bbox[0][1] - sio2_margin),
                                 (final_bbox[1][0] + sio2_margin, final_bbox[1][1] + sio2_margin), layer=layers['sio2'])
    carrier_cell.add(sio2_layer)
    return carrier_cell

# --- 两两相连焊盘单元生成函数 (无变化) ---
def create_carrier_chip_with_paired_pads(input_gds, lib, layers, marker_cell, marker_positions):
    print("正在生成版图: 两两相连焊盘载板 + Marker")
    trace_width = 20.0; bump_diameter = 20.0; sio2_margin = 50.0; INPUT_LAYER = 86
    paired_pad_cell = lib.new_cell('PairedPadUnit')
    original_cell = gdstk.read_gds(input_gds).top_level()[0]
    polygons = [p for p in original_cell.get_polygons() if p.layer == INPUT_LAYER]
    if not polygons: raise ValueError("在输入层上未找到多边形。")
    
    center_x, center_y = np.vstack([p.points for p in polygons]).mean(axis=0)
    pads_with_info = []
    for p in polygons:
        bbox = p.bounding_box()
        center = ((bbox[0][0] + bbox[1][0]) / 2, (bbox[0][1] + bbox[1][1]) / 2)
        angle = np.arctan2(center[1] - center_y, center[0] - center_x)
        pads_with_info.append({'polygon': p, 'center': center, 'angle': angle})
    pads_with_info.sort(key=lambda item: item['angle'])
    
    for i in range(0, len(pads_with_info), 2):
        if i + 1 >= len(pads_with_info): continue
        pad1, pad2 = pads_with_info[i], pads_with_info[i+1]
        path = gdstk.FlexPath([pad1['center'], pad2['center']], trace_width, layer=layers['metal'])
        pair_metal = gdstk.boolean([pad1['polygon'], pad2['polygon']], path, 'or', layer=layers['metal'])
        paired_pad_cell.add(*pair_metal)
    for pad_info in pads_with_info:
        ubm_ti = pad_info['polygon'].copy(); ubm_ti.layer = layers['ubm_ti']; paired_pad_cell.add(ubm_ti)
        ubm_pt = pad_info['polygon'].copy(); ubm_pt.layer = layers['ubm_pt']; paired_pad_cell.add(ubm_pt)
        ubm_au = pad_info['polygon'].copy(); ubm_au.layer = layers['ubm_au']; paired_pad_cell.add(ubm_au)
        paired_pad_cell.add(gdstk.regular_polygon(pad_info['center'], bump_diameter / 2.0, 16, layer=layers['bump']))

    for pos in marker_positions:
        paired_pad_cell.add(gdstk.Reference(marker_cell, pos))

    final_bbox = paired_pad_cell.bounding_box()
    sio2_layer = gdstk.rectangle((final_bbox[0][0] - sio2_margin, final_bbox[0][1] - sio2_margin),
                                 (final_bbox[1][0] + sio2_margin, final_bbox[1][1] + sio2_margin), layer=layers['sio2'])
    paired_pad_cell.add(sio2_layer)
    return paired_pad_cell

# --- 主程序 (已修改) ---
if __name__ == "__main__":
    input_filename = "PAD_ring_Xiling.gds"
    if not os.path.exists(input_filename):
         print(f"错误: 输入文件 '{input_filename}' 不存在。")
    else:
        layer_map = {
            1: 'Passivation - SiO2', 2: 'Carrier Metal', 3: 'UBM - Ti (20nm)',
            4: 'UBM - Pt (50nm)', 5: 'UBM - Au (100nm)', 6: 'Bump - Indium Pillar'
        }
        layers = {'sio2': 1, 'metal': 2, 'ubm_ti': 3, 'ubm_pt': 4, 'ubm_au': 5, 'bump': 6}

        # 1. 创建一个GDS库和通用的组件
        main_lib = gdstk.Library(unit=1e-6, precision=1e-9)
        marker_cell = create_alignment_marker(main_lib, layers)

        # 2. 预先计算统一的Marker位置
        original_cell = gdstk.read_gds(input_filename).top_level()[0]
        original_bbox = original_cell.bounding_box()
        marker_offset = 250 # 增加间距以适应更大的marker
        marker_positions = [
            (original_bbox[0][0] - marker_offset, original_bbox[0][1] - marker_offset), # 左下
            (original_bbox[1][0] + marker_offset, original_bbox[0][1] - marker_offset), # 右下
            (original_bbox[0][0] - marker_offset, original_bbox[1][1] + marker_offset), # 左上
            (original_bbox[1][0] + marker_offset, original_bbox[1][1] + marker_offset)  # 右上
        ]
        
        # 3. 生成两种不同的单元
        individual_pad_cell = create_carrier_chip_with_individual_pads(input_filename, main_lib, layers, marker_cell, marker_positions)
        paired_pad_cell = create_carrier_chip_with_paired_pads(input_filename, main_lib, layers, marker_cell, marker_positions)
            
        # 4. 创建顶层单元，并计算统一间距
        top_cell = main_lib.new_cell('Combined_2x2_Reticle')
        bbox_ind = individual_pad_cell.bounding_box()
        bbox_pair = paired_pad_cell.bounding_box()
        max_width = max(bbox_ind[1][0] - bbox_ind[0][0], bbox_pair[1][0] - bbox_pair[0][0])
        max_height = max(bbox_ind[1][1] - bbox_ind[0][1], bbox_pair[1][1] - bbox_pair[0][1])
        gap = 500 # 单元间距
        spacing_x = max_width + gap
        spacing_y = max_height + gap
        
        # 5. 手动放置四种单元，形成棋盘状混合版图
        top_cell.add(gdstk.Reference(paired_pad_cell, (0, 0)))           # 左下: 两两相连
        top_cell.add(gdstk.Reference(individual_pad_cell, (spacing_x, 0))) # 右下: 独立
        top_cell.add(gdstk.Reference(individual_pad_cell, (0, spacing_y)))   # 左上: 独立
        top_cell.add(gdstk.Reference(paired_pad_cell, (spacing_x, spacing_y)))# 右上: 两两相连

        # 6. 扁平化，保存并绘图
        top_cell.flatten()
        output_filename = f"combined_carrier_chip_2x2_reticle.gds"
        main_lib.write_gds(output_filename)
        print(f"\n成功生成GDS文件: {output_filename}")
        
        plot_cell(top_cell, layer_map, title="Combined Carrier Chip 2x2 Reticle")