import gdstk
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch

# --- 可视化函数 (已更新以适应新的Marker图层) ---
def plot_cell(cell, layer_map, title="Layout Preview"):
    fig, ax = plt.subplots(figsize=(14, 14))
    # 为新的Marker图层添加颜色
    color_map = {1: 'lightcyan', 2: 'darkblue', 3: 'lightgrey', 4: '#A7B5CC', 5: '#E5E4E2', 6: '#FFD700', 7: '#708090'}
    z_order = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7:6}
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
    print(f"\n正在显示版图预览: {title}")
    plt.show()

# --- Marker创建函数 (已更新，使用新的Marker图层) ---
def create_alignment_marker(lib, layers, size=400, width=80, marker_type='positive'):
    marker_cell = lib.new_cell(f"AlignmentMarker_{marker_type}")
    l, w = size / 2, width / 2
    # 使用新的Marker图层
    cross = gdstk.cross((0, 0), l, w, layer=layers['marker'])
    marker_cell.add(cross)
    tick_width = 2; tick_lengths = {5: 15, 10: 20, 50: 25}; start_offset = w + tick_lengths[50] + 10
    for i in range(int(start_offset), int(l) + 1):
        tick_len = 0
        if i % 50 == 0: tick_len = tick_lengths[50]
        elif i % 10 == 0: tick_len = tick_lengths[10]
        elif i % 5 == 0: tick_len = tick_lengths[5]
        if tick_len > 0:
            # 使用新的Marker图层
            if marker_type == 'positive':
                marker_cell.add(gdstk.rectangle((i - tick_width / 2, w), (i + tick_width / 2, w + tick_len), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((i - tick_width / 2, -w - tick_len), (i + tick_width / 2, -w), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((-i - tick_width / 2, w), (-i + tick_width / 2, w + tick_len), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((-i - tick_width / 2, -w - tick_len), (-i + tick_width / 2, -w), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((w, i - tick_width / 2), (w + tick_len, i + tick_width / 2), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((-w - tick_len, i - tick_width / 2), (-w, i + tick_width / 2), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((w, -i - tick_width / 2), (w + tick_len, -i + tick_width / 2), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((-w - tick_len, -i - tick_width / 2), (-w, -i + tick_width / 2), layer=layers['marker']))
            else:
                marker_cell.add(gdstk.rectangle((i - tick_width / 2, w - tick_len), (i + tick_width / 2, w), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((i - tick_width / 2, -w), (i + tick_width / 2, -w + tick_len), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((-i - tick_width / 2, w - tick_len), (-i + tick_width / 2, w), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((-i - tick_width / 2, -w), (-i + tick_width / 2, -w + tick_len), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((w - tick_len, i - tick_width / 2), (w, i + tick_width / 2), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((-w, i - tick_width / 2), (-w + tick_len, i + tick_width / 2), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((w - tick_len, -i - tick_width / 2), (w, -i + tick_width / 2), layer=layers['marker']))
                marker_cell.add(gdstk.rectangle((-w, -i - tick_width / 2), (-w + tick_len, -i + tick_width / 2), layer=layers['marker']))
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

# --- 主程序 (已更新图层定义) ---
if __name__ == "__main__":
    input_filename = "C:\\Users\\Administrator\\Desktop\\Spin_Quantum_Hardware_Simulation\\GDS_Drawer\\Flip_chip_test\\PAD_ring_Xiling.gds"
    if not os.path.exists(input_filename):
         print(f"错误: 输入文件 '{input_filename}' 不存在。")
    else:
        ROWS = 8
        COLUMNS = 8
        
        # **【修改处】**: 更新图层字典和映射关系
        layer_map = {
            1: 'Passivation - SiO2', 
            2: 'Marker - Ti/Au',
            3: 'Carrier Metal', 
            4: 'UBM - Ti (20nm)',
            5: 'UBM - Pt (50nm)', 
            6: 'UBM - Au (100nm)', 
            7: 'Bump - Indium Pillar'
        }
        layers = {
            'sio2': 1, 
            'marker': 2, 
            'metal': 3, 
            'ubm_ti': 4, 
            'ubm_pt': 5, 
            'ubm_au': 6, 
            'bump': 7
        }

        main_lib = gdstk.Library(unit=1e-6, precision=1e-9)
        
        marker_positive = create_alignment_marker(main_lib, layers, marker_type='positive')
        marker_negative = create_alignment_marker(main_lib, layers, marker_type='negative')

        original_cell = gdstk.read_gds(input_filename).top_level()[0]
        original_bbox = original_cell.bounding_box()
        marker_offset = 100
        marker_positions = [
            (original_bbox[0][0] - marker_offset, original_bbox[0][1] - marker_offset),
            (original_bbox[1][0] + marker_offset, original_bbox[0][1] - marker_offset),
            (original_bbox[0][0] - marker_offset, original_bbox[1][1] + marker_offset),
            (original_bbox[1][0] + marker_offset, original_bbox[1][1] + marker_offset)
        ]
        
        individual_pad_cell = create_carrier_chip_with_individual_pads(input_filename, main_lib, layers, marker_positive, marker_positions)
        paired_pad_cell = create_carrier_chip_with_paired_pads(input_filename, main_lib, layers, marker_negative, marker_positions)
            
        top_cell = main_lib.new_cell(f'Combined_{ROWS}x{COLUMNS}_Reticle')
        bbox_ind = individual_pad_cell.bounding_box()
        bbox_pair = paired_pad_cell.bounding_box()
        max_width = max(bbox_ind[1][0] - bbox_ind[0][0], bbox_pair[1][0] - bbox_pair[0][0])
        max_height = max(bbox_ind[1][1] - bbox_ind[0][1], bbox_pair[1][1] - bbox_pair[0][1])
        gap = 2000
        spacing_x, spacing_y = max_width + gap, max_height + gap
        
        for r in range(ROWS):
            for c in range(COLUMNS):
                x_pos = c * spacing_x
                y_pos = r * spacing_y
                
                if r < 4:
                    top_cell.add(gdstk.Reference(paired_pad_cell, (x_pos, y_pos)))
                else:
                    top_cell.add(gdstk.Reference(individual_pad_cell, (x_pos, y_pos)))

        top_cell.flatten()
        output_filename = f"combined_carrier_chip_{ROWS}x{COLUMNS}_20251125.gds"
        main_lib.write_gds(output_filename)
        print(f"\n成功生成GDS文件: {output_filename}")
        
        plot_cell(top_cell, layer_map, title=f"Combined Carrier Chip {ROWS}x{COLUMNS} Reticle")