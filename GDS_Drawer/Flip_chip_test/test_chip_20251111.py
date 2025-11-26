import gdstk
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch

def plot_cell(cell, layer_map, title="Layout Preview"):
    """
    使用 Matplotlib 绘制 gdstk.Cell 的内容。
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    
    color_map = {
        1: 'lightcyan', 2: 'lightgrey', 3: '#A7B5CC', 4: '#E5E4E2', 5: '#FFD700', 6: '#708090'
    }
    z_order = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5}

    polygons = cell.get_polygons()
    for poly in polygons:
        color = color_map.get(poly.layer, 'black')
        z = z_order.get(poly.layer, 0)
        patch = Polygon(poly.points, closed=True, facecolor=color, edgecolor='dimgray', linewidth=0.5, zorder=z)
        ax.add_patch(patch)

    legend_elements = [Patch(facecolor=color_map[layer], edgecolor='dimgray', label=name)
                       for layer, name in sorted(layer_map.items()) if layer in color_map]
    ax.legend(handles=legend_elements, loc='upper right')

    ax.set_aspect('equal')
    bbox = cell.bounding_box()
    if bbox is not None:
        x_margin, y_margin = (bbox[1][0] - bbox[0][0]) * 0.1, (bbox[1][1] - bbox[0][1]) * 0.1
        ax.set_xlim(bbox[0][0] - x_margin, bbox[1][0] + x_margin)
        ax.set_ylim(bbox[0][1] - y_margin, bbox[1][1] + y_margin)

    ax.set_xlabel("X coordinate (μm)"); ax.set_ylabel("Y coordinate (μm)")
    ax.set_title(title); ax.grid(True, linestyle='--', alpha=0.6)
    
    print("\n正在显示版图预览...")
    plt.show()


def create_paired_pad_cell(input_gds, lib):
    """
    读取GDS文件，生成一个焊盘两两相连的单元 (cell) 并返回。
    """
    print(f"正在读取输入GDS文件: {input_gds} 以创建两两相连的焊盘单元...")
    
    # --- 参数定义 ---
    trace_width = 20.0
    bump_diameter = 30.0
    sio2_margin = 20.0
    INPUT_LAYER = 86

    # 图层定义 (SiO2 = 1, Metal = 2, UBM = 3-5, Bump = 6)
    LAYER_SIO2, LAYER_METAL, LAYER_UBM_TI, LAYER_UBM_PT, LAYER_UBM_AU, LAYER_BUMP = 1, 2, 3, 4, 5, 6
    
    original_lib = gdstk.read_gds(input_gds)
    original_cell = original_lib.top_level()[0]
    
    polygons_to_process = [p for p in original_cell.get_polygons() if p.layer == INPUT_LAYER]
    if not polygons_to_process:
        raise ValueError(f"在输入文件 '{input_gds}' 的第 {INPUT_LAYER} 层上没有找到任何多边形。")
        
    print(f"成功读取单元 '{original_cell.name}'，找到 {len(polygons_to_process)} 个多边形进行处理。")

    paired_pad_cell = lib.new_cell('PairedPadUnit')

    # --- 1. 排序焊盘以确定连接顺序 ---
    all_points = np.vstack([p.points for p in polygons_to_process])
    center_x, center_y = all_points.mean(axis=0)

    pads_with_info = []
    for poly in polygons_to_process:
        bbox = poly.bounding_box()
        pc = ((bbox[0][0] + bbox[1][0]) / 2, (bbox[0][1] + bbox[1][1]) / 2)
        angle = np.arctan2(pc[1] - center_y, pc[0] - center_x)
        pads_with_info.append({'polygon': poly, 'center': pc, 'angle': angle})
    
    sorted_pads = sorted(pads_with_info, key=lambda p: p['angle'])
    
    # --- 2. 遍历排序后的焊盘，两两相连 ---
    for i in range(0, len(sorted_pads), 2):
        # 确保我们有一个完整的对
        if i + 1 >= len(sorted_pads):
            # 如果焊盘总数是奇数，最后一个将不被连接
            print(f"警告: 焊盘总数为奇数。最后一个焊盘将不被连接。")
            continue

        pad1_info = sorted_pads[i]
        pad2_info = sorted_pads[i+1]

        # 创建连接两个焊盘中心的路径
        path = gdstk.FlexPath([pad1_info['center'], pad2_info['center']], trace_width, layer=LAYER_METAL)
        
        # 将两个焊盘多边形和路径合并成一个哑铃形状
        polygons_for_pair = [pad1_info['polygon'], pad2_info['polygon']]
        pair_metal_list = gdstk.boolean(polygons_for_pair, path, 'or', layer=LAYER_METAL)
        paired_pad_cell.add(*pair_metal_list)
        
        # 为这个哑铃形状创建对应的 SiO2 底层
        temp_cell = gdstk.Cell("temp_bbox").add(*pair_metal_list)
        pair_bbox = temp_cell.bounding_box()
        sio2_rect = gdstk.rectangle(
            (pair_bbox[0][0] - sio2_margin, pair_bbox[0][1] - sio2_margin),
            (pair_bbox[1][0] + sio2_margin, pair_bbox[1][1] + sio2_margin),
            layer=LAYER_SIO2
        )
        paired_pad_cell.add(sio2_rect)

    # --- 3. 在每个原始焊盘位置添加 UBM 和 Bumps ---
    for pad_info in sorted_pads:
        poly = pad_info['polygon']
        pad_center = pad_info['center']
        
        ubm_ti = poly.copy(); ubm_ti.layer = LAYER_UBM_TI; paired_pad_cell.add(ubm_ti)
        ubm_pt = poly.copy(); ubm_pt.layer = LAYER_UBM_PT; paired_pad_cell.add(ubm_pt)
        ubm_au = poly.copy(); ubm_au.layer = LAYER_UBM_AU; paired_pad_cell.add(ubm_au)
        indium_bump = gdstk.regular_polygon(pad_center, bump_diameter / 2.0, 16, layer=LAYER_BUMP)
        paired_pad_cell.add(indium_bump)
        
    return paired_pad_cell


if __name__ == "__main__":
    input_filename = "PAD_ring_Xiling.gds"
    if not os.path.exists(input_filename):
         print(f"错误: 输入文件 '{input_filename}' 不存在。")
    else:
        layer_map = {
            1: 'Passivation - SiO2',
            2: 'Paired-Pad Metal',
            3: 'UBM - Ti (20nm)',
            4: 'UBM - Pt (50nm)',
            5: 'UBM - Au (100nm)',
            6: 'Bump - Indium Pillar'
        }
        main_lib = gdstk.Library(unit=1e-6, precision=1e-9)
        
        # 1. 创建两两相连的焊盘单元
        paired_pad_unit_cell = create_paired_pad_cell(input_filename, main_lib)
        
        # 2. 将单元排列成 2x2 矩阵
        matrix_cell = main_lib.new_cell('PairedPad_2x2_Matrix')
        bbox = paired_pad_unit_cell.bounding_box()
        chip_width, chip_height = bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1]
        gap = 200 # 200 um gap
        
        chip_array = gdstk.Reference(
            paired_pad_unit_cell,
            columns=8, rows=8,
            spacing=(chip_width + gap, chip_height + gap)
        )
        matrix_cell.add(chip_array)
        
        # 3. 扁平化，保存并绘图
        matrix_cell.flatten()
        
        output_filename = "carrier_chip_paired_pad_2x2.gds"
        main_lib.write_gds(output_filename)
        
        print(f"\n2x2 两两相连焊盘矩阵GDS文件 '{output_filename}' 已成功生成。")
        print("版图包含以下图层:")
        for layer_num, name in sorted(layer_map.items()):
            print(f" - 第 {layer_num} 层: {name}")

        plot_cell(matrix_cell, layer_map, title="Paired-Pad Carrier Chip 2x2 Matrix")