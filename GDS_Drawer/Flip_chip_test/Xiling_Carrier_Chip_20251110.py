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
    
    # 定义颜色查找表
    color_map = {
        1: 'lightgrey',
        2:'#A7B5CC', # Titanium-like for Ti 
        3: '#E5E4E2', # Platinum-like for Pt
        4: '#FFD700', # Gold-like for Au
        5: '#708090'  # SlateGray for Indium Bump
    }
    
    # 遍历单元中的所有多边形并按图层排序以确保bump在最上层绘制
    polygons = sorted(cell.get_polygons(), key=lambda p: p.layer)
    for poly in polygons:
        color = color_map.get(poly.layer, 'black')
        patch = Polygon(poly.points, closed=True, facecolor=color, edgecolor='dimgray', linewidth=0.5)
        ax.add_patch(patch)

    # 创建图例
    legend_elements = [Patch(facecolor=color_map[layer], edgecolor='dimgray', label=name)
                       for layer, name in layer_map.items() if layer in color_map]
    ax.legend(handles=legend_elements, loc='upper right')

    # 设置绘图区域
    ax.set_aspect('equal')
    bbox = cell.bounding_box()
    if bbox is not None:
        x_margin = (bbox[1][0] - bbox[0][0]) * 0.05
        y_margin = (bbox[1][1] - bbox[0][1]) * 0.05
        ax.set_xlim(bbox[0][0] - x_margin, bbox[1][0] + x_margin)
        ax.set_ylim(bbox[0][1] - y_margin, bbox[1][1] + y_margin)

    # 添加标签和标题
    ax.set_xlabel("X coordinate (μm)")
    ax.set_ylabel("Y coordinate (μm)")
    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.6)
    
    print("\n正在显示版图预览...")
    plt.show()


def create_carrier_chip(input_gds, output_gds="carrier_chip.gds"):
    """
    读取GDS文件，生成带有UBM和Bump的载板芯片版图，并使用 Matplotlib 显示结果。
    """
    print(f"\n正在读取输入GDS文件: {input_gds}")
    
    # --- 参数定义 ---
    wb_pad_size = (100, 100)
    clearance = 300
    trace_width = 50
    bump_diameter = 30.0

    # --- 图层定义 ---
    INPUT_LAYER = 86
    layer_map = {
        1: 'Carrier Metal (Pads + Traces)',
        2: 'UBM - Ti (20nm)',
        3: 'UBM - Pt (50nm)',
        4: 'UBM - Au (100nm)',
        5: 'Bump - Indium Pillar'
    }
    LAYER_CARRIER_METAL = 1
    LAYER_UBM_AU = 2
    LAYER_UBM_PT = 3
    LAYER_UBM_TI = 4
    LAYER_BUMP = 5
    
    try:
        original_lib = gdstk.read_gds(input_gds)
    except FileNotFoundError:
        print(f"错误: 输入文件 '{input_gds}' 未找到。")
        return

    original_cell = original_lib.top_level()[0]
    
    polygons_to_process = [p for p in original_cell.get_polygons() if p.layer == INPUT_LAYER]
    
    if not polygons_to_process:
        print(f"警告: 在输入文件 '{input_gds}' 的第 {INPUT_LAYER} 层上没有找到任何多边形。")
        return
        
    print(f"成功读取单元 '{original_cell.name}'，在第 {INPUT_LAYER} 层上找到 {len(polygons_to_process)} 个多边形进行处理。")

    carrier_lib = gdstk.Library(unit=1e-6, precision=1e-9)
    carrier_cell = carrier_lib.new_cell('CarrierChip_Layout')

    bbox = original_cell.bounding_box()
    center_x = (bbox[0][0] + bbox[1][0]) / 2
    center_y = (bbox[0][1] + bbox[1][1]) / 2
    
    for poly in polygons_to_process:
        pad_bbox = poly.bounding_box()
        x_min, y_min, x_max, y_max = pad_bbox[0][0], pad_bbox[0][1], pad_bbox[1][0], pad_bbox[1][1]
        pad_center = ((x_min + x_max) / 2, (y_min + y_max) / 2)
        pad_width, pad_height = x_max - x_min, y_max - y_min

        # 1. 创建载板金属层
        flip_chip_pad_on_carrier = poly.copy()
        flip_chip_pad_on_carrier.layer = LAYER_CARRIER_METAL
        carrier_cell.add(flip_chip_pad_on_carrier)
        
        dx, dy = pad_center[0] - center_x, pad_center[1] - center_y
        if abs(dy) > abs(dx):
            offset = np.sign(dy) * (clearance + wb_pad_size[1] / 2 + pad_height / 2)
            wb_pad_center_y, wb_pad_center_x = pad_center[1] + offset, pad_center[0]
        else:
            offset = np.sign(dx) * (clearance + wb_pad_size[0] / 2 + pad_width / 2)
            wb_pad_center_x, wb_pad_center_y = pad_center[0] + offset, pad_center[1]

        wb_pad = gdstk.rectangle(
            (wb_pad_center_x - wb_pad_size[0] / 2, wb_pad_center_y - wb_pad_size[1] / 2),
            (wb_pad_center_x + wb_pad_size[0] / 2, wb_pad_center_y + wb_pad_size[1] / 2),
            layer=LAYER_CARRIER_METAL
        )
        carrier_cell.add(wb_pad)
        
        trace = gdstk.FlexPath([pad_center, (wb_pad_center_x, wb_pad_center_y)],
                               trace_width, layer=LAYER_CARRIER_METAL)
        carrier_cell.add(trace)

        # 2. 创建 UBM (Under-Bump Metallurgy) 层
        ubm_ti = poly.copy(); ubm_ti.layer = LAYER_UBM_TI
        carrier_cell.add(ubm_ti)
        ubm_pt = poly.copy(); ubm_pt.layer = LAYER_UBM_PT
        carrier_cell.add(ubm_pt)
        ubm_au = poly.copy(); ubm_au.layer = LAYER_UBM_AU
        carrier_cell.add(ubm_au)

        # --- **【错误修复处】** ---
        # 3. 创建 Indium Bump
        bump_radius = bump_diameter / 2.0 
        # 使用 gdstk.regular_polygon 来创建一个128边的多边形以近似圆形
        indium_bump = gdstk.regular_polygon(pad_center, bump_radius, 8, layer=LAYER_BUMP)
        carrier_cell.add(indium_bump)
        # --- **【修复结束】** ---

    # --- 文件保存和绘图 ---
    carrier_lib.write_gds(output_gds)
    print(f"\n载板芯片GDS文件 '{output_gds}' 已成功生成。")
    print("版图包含以下图层:")
    for layer_num, name in layer_map.items():
        print(f" - 第 {layer_num} 层: {name}")

    plot_cell(carrier_cell, layer_map, title="Carrier Chip Full Layout")


if __name__ == "__main__":
    input_filename = "PAD_ring_Xiling.gds"
    if not os.path.exists(input_filename):
         print(f"错误: 输入文件 '{input_filename}' 不存在。请确保文件与脚本在同一目录下。")
    else:
        output_filename = "carrier_chip_full.gds"
        create_carrier_chip(input_gds=input_filename, output_gds=output_filename)