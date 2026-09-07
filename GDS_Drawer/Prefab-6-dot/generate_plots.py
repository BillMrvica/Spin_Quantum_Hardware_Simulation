#!/usr/bin/env python3
"""
生成可视化图形
"""

import os
import sys

os.chdir(r"C:\Users\Administrator\Desktop\Spin_Quantum_Hardware_Simulation\GDS_Drawer\Prefab-6-dot")
sys.path.insert(0, os.getcwd())

import gdstk
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import numpy as np

# 设置输出目录
output_dir = r"C:\Users\Administrator\Desktop\Spin_Quantum_Hardware_Simulation\GDS_Drawer\Prefab-6-dot\output_gds"
os.makedirs(output_dir, exist_ok=True)

def plot_gds_to_file(cell, title, filename, figsize=(16, 12)):
    """保存GDS可视化到文件"""
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(title, fontsize=18)
    
    layer_config = {
        0: {'color': '#87CEEB', 'alpha': 0.6, 'label': 'S/D (L0)', 'zorder': 2},
        1: {'color': '#D3D3D3', 'alpha': 0.5, 'label': 'SG (L1)', 'hatch': '///', 'zorder': 1},
        2: {'color': '#FF1493', 'alpha': 0.9, 'label': 'BG (L2)', 'zorder': 3},
        3: {'color': '#8A2BE2', 'alpha': 0.9, 'label': 'PG (L3)', 'zorder': 4}
    }
    
    drawn_labels = set()
    
    for layer in [1, 0, 2, 3]:
        polygons = cell.get_polygons(layer=layer, datatype=0)
        if not polygons:
            continue
        
        cfg = layer_config.get(layer)
        label = cfg['label']
        
        for gds_poly in polygons:
            pts = gds_poly.points
            if pts.ndim != 2 or len(pts) < 3:
                continue
            
            current_label = label if label not in drawn_labels else None
            poly = MplPolygon(pts, closed=True, facecolor=cfg['color'],
                            edgecolor='black', linewidth=0.5, alpha=cfg['alpha'],
                            hatch=cfg.get('hatch', ''), label=current_label,
                            zorder=cfg['zorder'])
            ax.add_patch(poly)
            drawn_labels.add(label)
    
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    ax.set_aspect('equal')
    
    bbox = cell.bounding_box()
    if bbox is not None:
        min_x, min_y = bbox[0]
        max_x, max_y = bbox[1]
        width = max_x - min_x
        height = max_y - min_y
        margin_x = width * 0.1 if width > 0 else 1
        margin_y = height * 0.1 if height > 0 else 1
        ax.set_xlim(min_x - margin_x, max_x + margin_x)
        ax.set_ylim(min_y - margin_y, max_y + margin_y)
    
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')
    plt.grid(True, which='both', linestyle='--', alpha=0.3)
    plt.tight_layout()
    
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✓ 图形已保存: {filepath}")

def plot_assembly_to_file(cell, title, filename, figsize=(18, 18)):
    """保存组装后的可视化到文件"""
    plot_cell = cell.flatten()
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(title, fontsize=18)
    
    layer_config = {
        0: {'color': '#87CEEB', 'alpha': 0.9, 'label': 'S/D (L0)', 'zorder': 4},
        2: {'color': '#FF1493', 'alpha': 0.9, 'label': 'BG (L2)', 'zorder': 5},
        3: {'color': '#8A2BE2', 'alpha': 0.9, 'label': 'PG (L3)', 'zorder': 6},
        5: {'color': '#0000FF', 'alpha': 1.0, 'label': 'Final Routing (L5)', 'zorder': 3},
        10: {'color': 'pink', 'alpha': 0.5, 'label': 'Active Area (L10)', 'zorder': 2},
        11: {'color': 'gray', 'alpha': 1.0, 'label': 'SiO2 (L11)', 'zorder': 0}
    }
    
    drawn_labels = set()
    
    def draw_polygons(polygons, config):
        label = config['label']
        for gds_poly in polygons:
            vertices = gds_poly.points
            if len(vertices) > 0:
                current_label = label if label not in drawn_labels else None
                is_fill = True if config.get('layer', -1) != 11 else False
                edge_color = 'black' if config.get('layer', -1) != 11 else 'gray'
                line_style = 'solid' if config.get('layer', -1) != 11 else '--'
                
                mpl_poly = MplPolygon(vertices, closed=True, fill=is_fill,
                                      facecolor=config['color'], edgecolor=edge_color,
                                      linestyle=line_style, linewidth=0.2,
                                      alpha=config['alpha'], label=current_label,
                                      zorder=config.get('zorder', 1))
                ax.add_patch(mpl_poly)
                drawn_labels.add(label)
    
    for layer, cfg in layer_config.items():
        polys = plot_cell.get_polygons(layer=layer, datatype=0)
        if polys:
            draw_polygons(polys, {**cfg, 'layer': layer})
    
    # Pad层
    pad_polys = plot_cell.get_polygons(layer=1, datatype=0)
    if pad_polys:
        sg_config = {'color': '#D3D3D3', 'alpha': 0.5, 'label': 'SG (L1)',
                    'hatch': '///', 'zorder': 7}
        pad_config = {'color': 'gold', 'alpha': 0.8, 'label': 'Pads & Traces (L1)',
                     'zorder': 1}
        
        for poly in pad_polys:
            (min_x, min_y), (max_x, max_y) = poly.bounding_box()
            if (max_x - min_x) > 1 or (max_y - min_y) > 1:
                draw_polygons([poly], pad_config)
            else:
                draw_polygons([poly], sg_config)
    
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    bbox = plot_cell.bounding_box()
    if bbox:
        (min_x, min_y), (max_x, max_y) = bbox
        margin = (max_x - min_x) * 0.05
        ax.set_xlim(min_x - margin, max_x + margin)
        ax.set_ylim(min_y - margin, max_y + margin)
    
    ax.set_aspect('equal')
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✓ 图形已保存: {filepath}")

print("=" * 60)
print("生成可视化图形...")
print("=" * 60)
print()

# 步骤1: 量子点器件
print("步骤1: 量子点器件核心")
lib1 = gdstk.read_gds(os.path.join(output_dir, "step1_device_with_leads.gds"))
cell1 = lib1.top_level()[0]
plot_gds_to_file(cell1, "步骤1: 6量子点器件核心 (带长引线)", "step1_device_with_leads.png")

# 步骤2: Pad框架
print("\n步骤2: Pad框架")
lib2 = gdstk.read_gds(os.path.join(output_dir, "step2_pad_frame.gds"))
cell2 = lib2.top_level()[0]
plot_assembly_to_file(cell2, "步骤2: Pad框架 (带引线)", "step2_pad_frame.png")

# 步骤3: 组装未布线
print("\n步骤3: 器件已放置 (未布线)")
lib3 = gdstk.read_gds(os.path.join(output_dir, "step3_placed_unrouted.gds"))
cell3 = lib3.top_level()[0]
plot_assembly_to_file(cell3, "步骤3: 器件已放置 (未布线)", "step3_placed_unrouted.png")

# 步骤4: 单线布线
print("\n步骤4: 单线布线")
lib4 = gdstk.read_gds(os.path.join(output_dir, "step4_single_route.gds"))
cell4 = lib4.top_level()[0]
plot_assembly_to_file(cell4, "步骤4: 单线布线 (QD_PG4)", "step4_single_route.png")

# 步骤5: 完整布线
print("\n步骤5: 完整布线")
lib5 = gdstk.read_gds(os.path.join(output_dir, "step5_fully_routed.gds"))
cell5 = lib5.top_level()[0]
plot_assembly_to_file(cell5, "步骤5: 最终完整布线", "step5_fully_routed.png")

print()
print("=" * 60)
print("所有可视化图形生成完成!")
print("=" * 60)

# 列出所有生成的文件
print("\n生成的所有文件:")
for f in sorted(os.listdir(output_dir)):
    filepath = os.path.join(output_dir, f)
    size = os.path.getsize(filepath)
    print(f"  - {f} ({size:,} bytes)")
