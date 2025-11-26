import gdstk
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch

def plot_cell(cell, layer_map, title="Layout Preview"):
    """
    使用 Matplotlib 绘制 gdstk.Cell 的内容，并为不同图层使用不同颜色。
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # **【修改处】**: 更新颜色和绘图顺序以匹配新的图层编号
    color_map = {
        1: 'lightcyan', # SiO2
        2: 'lightgrey', # Carrier Metal
        3: '#A7B5CC',  # Ti
        4: '#E5E4E2',  # Pt
        5: '#FFD700',  # Au
        6: '#708090'   # Indium Bump
    }
    z_order = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5} # SiO2 (layer 1) is drawn first (bottom)

    # 扁平化单元以获取所有多边形
    polygons = cell.get_polygons()
    for poly in polygons:
        color = color_map.get(poly.layer, 'black')
        z = z_order.get(poly.layer, 0)
        patch = Polygon(poly.points, closed=True, facecolor=color, edgecolor='dimgray', linewidth=0.5, zorder=z)
        ax.add_patch(patch)

    # 创建图例
    legend_elements = [Patch(facecolor=color_map[layer], edgecolor='dimgray', label=name)
                       for layer, name in sorted(layer_map.items()) if layer in color_map]
    ax.legend(handles=legend_elements, loc='upper right')

    # 设置绘图区域
    ax.set_aspect('equal')
    bbox = cell.bounding_box()
    if bbox is not None:
        x_margin = (bbox[1][0] - bbox[0][0]) * 0.1
        y_margin = (bbox[1][1] - bbox[0][1]) * 0.1
        ax.set_xlim(bbox[0][0] - x_margin, bbox[1][0] + x_margin)
        ax.set_ylim(bbox[0][1] - y_margin, bbox[1][1] + y_margin)

    # 添加标签和标题
    ax.set_xlabel("X coordinate (μm)")
    ax.set_ylabel("Y coordinate (μm)")
    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.6)
    
    print("\n正在显示版图预览...")
    plt.show()


def create_single_carrier_chip_cell(input_gds, lib):
    """
    读取GDS文件，生成一个包含完整载板芯片设计的 gdstk.Cell 并返回。
    """
    print(f"正在读取输入GDS文件: {input_gds}")
    
    # --- 参数定义 ---
    wb_pad_size = (100, 100)
    clearance = 300
    trace_width = 50
    bump_diameter = 30.0
    sio2_margin = 20.0 # SiO2 layer overhang
    INPUT_LAYER = 86

    # **【修改处】**: 重新定义图层编号，SiO2 为第 1 层
    LAYER_SIO2, LAYER_CARRIER_METAL, LAYER_UBM_TI, LAYER_UBM_PT, LAYER_UBM_AU, LAYER_BUMP = 1, 2, 3, 4, 5, 6
    
    original_lib = gdstk.read_gds(input_gds)
    original_cell = original_lib.top_level()[0]
    
    polygons_to_process = [p for p in original_cell.get_polygons() if p.layer == INPUT_LAYER]
    if not polygons_to_process:
        raise ValueError(f"在输入文件 '{input_gds}' 的第 {INPUT_LAYER} 层上没有找到任何多边形。")
        
    print(f"成功读取单元 '{original_cell.name}'，找到 {len(polygons_to_process)} 个多边形进行处理。")

    carrier_cell = lib.new_cell('SingleChipUnit')

    bbox = original_cell.bounding_box()
    center_x, center_y = (bbox[0][0] + bbox[1][0]) / 2, (bbox[0][1] + bbox[1][1]) / 2
    
    for poly in polygons_to_process:
        pad_bbox = poly.bounding_box()
        pad_center = ((pad_bbox[0][0] + pad_bbox[1][0]) / 2, (pad_bbox[0][1] + pad_bbox[1][1]) / 2)
        
        # 1. 创建载板金属层 (Layer 2)
        flip_chip_pad_on_carrier = poly.copy(); flip_chip_pad_on_carrier.layer = LAYER_CARRIER_METAL
        carrier_cell.add(flip_chip_pad_on_carrier)
        
        dx, dy = pad_center[0] - center_x, pad_center[1] - center_y
        if abs(dy) > abs(dx):
            offset = np.sign(dy) * (clearance + wb_pad_size[1] / 2 + (pad_bbox[1][1]-pad_bbox[0][1]) / 2)
            wb_pad_center_y, wb_pad_center_x = pad_center[1] + offset, pad_center[0]
        else:
            offset = np.sign(dx) * (clearance + wb_pad_size[0] / 2 + (pad_bbox[1][0]-pad_bbox[0][0]) / 2)
            wb_pad_center_x, wb_pad_center_y = pad_center[0] + offset, pad_center[1]

        wb_pad = gdstk.rectangle(
            (wb_pad_center_x - wb_pad_size[0] / 2, wb_pad_center_y - wb_pad_size[1] / 2),
            (wb_pad_center_x + wb_pad_size[0] / 2, wb_pad_center_y + wb_pad_size[1] / 2),
            layer=LAYER_CARRIER_METAL)
        carrier_cell.add(wb_pad)
        
        trace = gdstk.FlexPath([pad_center, (wb_pad_center_x, wb_pad_center_y)],
                               trace_width, layer=LAYER_CARRIER_METAL)
        carrier_cell.add(trace)

        # 2. 创建 UBM 层 (Layers 3, 4, 5)
        ubm_ti = poly.copy(); ubm_ti.layer = LAYER_UBM_TI; carrier_cell.add(ubm_ti)
        ubm_pt = poly.copy(); ubm_pt.layer = LAYER_UBM_PT; carrier_cell.add(ubm_pt)
        ubm_au = poly.copy(); ubm_au.layer = LAYER_UBM_AU; carrier_cell.add(ubm_au)

        # 3. 创建 Indium Bump (Layer 6)
        indium_bump = gdstk.regular_polygon(pad_center, bump_diameter / 2.0, 8, layer=LAYER_BUMP)
        carrier_cell.add(indium_bump)
        
        # 4. 创建 SiO2 钝化层 (Layer 1)
        wb_pad_bbox = wb_pad.bounding_box()
        sio2_min_x = min(pad_bbox[0][0], wb_pad_bbox[0][0]) - sio2_margin
        sio2_min_y = min(pad_bbox[0][1], wb_pad_bbox[0][1]) - sio2_margin
        sio2_max_x = max(pad_bbox[1][0], wb_pad_bbox[1][0]) + sio2_margin
        sio2_max_y = max(pad_bbox[1][1], wb_pad_bbox[1][1]) + sio2_margin
        
        sio2_undermetal = gdstk.rectangle((sio2_min_x, sio2_min_y), (sio2_max_x, sio2_max_y), layer=LAYER_SIO2)
        carrier_cell.add(sio2_undermetal)
        
    return carrier_cell


if __name__ == "__main__":
    input_filename = "PAD_ring_Xiling.gds"
    if not os.path.exists(input_filename):
         print(f"错误: 输入文件 '{input_filename}' 不存在。请确保文件与脚本在同一目录下。")
    else:
        # **【修改处】**: 更新图层映射以反映新的编号顺序
        layer_map = {
            1: 'Passivation - SiO2',
            2: 'Carrier Metal (Pads + Traces)',
            3: 'UBM - Ti (20nm)',
            4: 'UBM - Pt (50nm)',
            5: 'UBM - Au (100nm)',
            6: 'Bump - Indium Pillar'
        }
        main_lib = gdstk.Library(unit=1e-6, precision=1e-9)
        single_chip_cell = create_single_carrier_chip_cell(input_filename, main_lib)
        
        matrix_cell = main_lib.new_cell('2x2_Matrix_Layout')
        bbox = single_chip_cell.bounding_box()
        chip_width = bbox[1][0] - bbox[0][0]
        chip_height = bbox[1][1] - bbox[0][1]
        gap = 200 # 200 um gap
        
        chip_array = gdstk.Reference(
            single_chip_cell,
            columns=2,
            rows=2,
            spacing=(chip_width + gap, chip_height + gap)
        )
        
        matrix_cell.add(chip_array)
        
        # 在保存和绘图之前，扁平化单元格
        matrix_cell.flatten()
        
        output_filename = "carrier_chip_2x2_matrix_reordered.gds"
        main_lib.write_gds(output_filename)
        
        print(f"\n2x2 矩阵GDS文件 '{output_filename}' 已成功生成。")
        print("版图包含以下图层:")
        for layer_num, name in sorted(layer_map.items()):
            print(f" - 第 {layer_num} 层: {name}")

        plot_cell(matrix_cell, layer_map, title="Carrier Chip 2x2 Matrix (SiO2 on Layer 1)")