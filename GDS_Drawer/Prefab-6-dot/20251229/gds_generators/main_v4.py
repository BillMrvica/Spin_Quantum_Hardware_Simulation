import sys
import os
# Add the parent directory to sys.path to allow absolute imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.patches as patches
import math
import re

from gds_generators.quantum_dot_pad_generator import QuantumDotPadGenerator
from gds_generators.route_v4 import RouteManager_v4

# ===================================================================
#                      解决中文显示问题的代码
# ===================================================================
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
# ===================================================================

def plot_gds(cell, title="量子点器件布局", show_plot=True):
    if not show_plot:
        return

    fig, ax = plt.subplots(figsize=(16, 12)); ax.set_title(title, fontsize=18)
    layer_config = {0:{'color':'#87CEEB','alpha':0.6,'label':'S/D (L0)','zorder':2}, 1:{'color':'#D3D3D3','alpha':0.5,'label':'SG (L1)','hatch':'///','zorder':1}, 2:{'color':'#FF1493','alpha':0.9,'label':'BG (L2)','zorder':3}, 3:{'color':'#8A2BE2','alpha':0.9,'label':'PG (L3)','zorder':4}}
    drawn_labels = set()
    for layer in [1, 0, 2, 3]: 
        polygons = cell.get_polygons(layer=layer, datatype=0)
        if not polygons: continue
        cfg = layer_config.get(layer); label = cfg['label']
        for gds_poly in polygons:
            pts = gds_poly.points
            if pts.ndim != 2 or len(pts) < 3: continue
            current_label = label if label not in drawn_labels else None
            poly = MplPolygon(pts, closed=True, facecolor=cfg['color'], edgecolor='black', linewidth=0.5, alpha=cfg['alpha'], hatch=cfg.get('hatch', ''), label=current_label, zorder=cfg['zorder'])
            ax.add_patch(poly); drawn_labels.add(label)
    handles, labels = ax.get_legend_handles_labels()
    if handles: by_label = dict(zip(labels, handles)); ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    ax.set_aspect('equal')
    bbox = cell.bounding_box()
    if bbox is not None:
        min_x, min_y = bbox[0]; max_x, max_y = bbox[1]; width = max_x - min_x; height = max_y - min_y; margin_x = width * 0.1 if width > 0 else 1; margin_y = height * 0.1 if height > 0 else 1
        ax.set_xlim(min_x - margin_x, max_x + margin_x); ax.set_ylim(min_y - margin_y, max_y + margin_y); print(f"Auto-adjusting View: X[{min_x:.2f}, {max_x:.2f}], Y[{min_y:.2f}, {max_y:.2f}]")
    else: ax.autoscale_view()
    ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); plt.grid(True, which='both', linestyle='--', alpha=0.3); plt.tight_layout(); plt.show()

def visualize_and_save_pad_frame(lib, cell_to_show, title, gds_filename, show_plot=True):
    print(f"--- {title} ---"); print(f"正在可视化: {title}...")
    
    if show_plot:
        fig, ax = plt.subplots(figsize=(18, 18)); plot_cell = cell_to_show.flatten()
        pads_and_traces = plot_cell.get_polygons(layer=1, datatype=0); active_area_poly = plot_cell.get_polygons(layer=10, datatype=0); sio2_poly = plot_cell.get_polygons(layer=11, datatype=0)
        for gds_poly in sio2_poly: ax.add_patch(MplPolygon(gds_poly.points, closed=True, fill=False, edgecolor='gray', linestyle='--'))
        for i, gds_poly in enumerate(pads_and_traces):
            label = "Pads & Traces (L1)" if i == 0 else None
            ax.add_patch(MplPolygon(gds_poly.points, closed=True, facecolor='gold', edgecolor='black', linewidth=0.2, label=label))
        if active_area_poly: ax.add_patch(MplPolygon(active_area_poly[0].points, closed=True, facecolor='pink', alpha=0.5, label="Active Area (L10)"))
        ax.legend(loc='upper right'); ax.set_aspect('equal')
        bbox = cell_to_show.bounding_box()
        if bbox: (min_x,min_y),(max_x,max_y)=bbox; margin=(max_x-min_x)*0.05; ax.set_xlim(min_x-margin,max_x+margin); ax.set_ylim(min_y-margin,max_y+margin)
        ax.set_title(title, fontsize=18); ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); plt.grid(True, linestyle='--', alpha=0.3); plt.tight_layout(); plt.show()
    
    lib.write_gds(gds_filename); print(f"GDS 文件已保存到: '{gds_filename}'\n")

def visualize_and_save_assembly(lib, cell_to_show, title, gds_filename, show_plot=True):
    print(f"--- {title} ---"); print(f"正在可视化: {title}...")
    
    if not show_plot:
        lib.write_gds(gds_filename); print(f"GDS 文件已保存到: '{gds_filename}'\n")
        return

    plot_cell = cell_to_show.flatten(); fig, ax = plt.subplots(figsize=(18, 18))
    layer_config = {0:{'color':'#87CEEB','alpha':0.9,'label':'S/D (L0)','zorder':4}, 2:{'color':'#FF1493','alpha':0.9,'label':'BG (L2)','zorder':5}, 3:{'color':'#8A2BE2','alpha':0.9,'label':'PG (L3)','zorder':6}, 5:{'color':'#0000FF','alpha':1.0,'label':'Final Routing (L5)','zorder':3}, 10:{'color':'pink','alpha':0.5,'label':'Active Area (L10)','zorder':2}, 11:{'color':'gray','alpha':1.0,'label':'SiO2 (L11)','zorder':0}}
    drawn_labels = set()
    def draw_polygons(polygons, config):
        label = config['label']
        for gds_poly in polygons:
            vertices = gds_poly.points
            if len(vertices) > 0:
                current_label = label if label not in drawn_labels else None; is_fill = True if config.get('layer', -1) != 11 else False; edge_color = 'black' if config.get('layer', -1) != 11 else 'gray'; line_style = 'solid' if config.get('layer', -1) != 11 else '--'
                mpl_poly = MplPolygon(vertices, closed=True, fill=is_fill, facecolor=config['color'], edgecolor=edge_color, linestyle=line_style, linewidth=0.2, alpha=config['alpha'], label=current_label, zorder=config.get('zorder', 1))
                ax.add_patch(mpl_poly); drawn_labels.add(label)
    for layer, cfg in layer_config.items():
        polys = plot_cell.get_polygons(layer=layer, datatype=0)
        if polys: draw_polygons(polys, {**cfg, 'layer': layer})
    pad_polys = plot_cell.get_polygons(layer=1, datatype=0)
    if pad_polys:
        sg_config = {'color':'#D3D3D3', 'alpha':0.5, 'label':'SG (L1)', 'hatch':'///', 'zorder':7}; pad_config = {'color':'gold','alpha':0.8,'label':'Pads & Traces (L1)','zorder':1}
        for poly in pad_polys:
            (min_x, min_y), (max_x, max_y) = poly.bounding_box()
            if (max_x - min_x) > 1 or (max_y - min_y) > 1: draw_polygons([poly], pad_config)
            else: draw_polygons([poly], sg_config)
    handles, labels = ax.get_legend_handles_labels(); by_label = dict(zip(labels, handles)); ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    bbox=plot_cell.bounding_box(); 
    if bbox: (min_x,min_y),(max_x,max_y)=bbox; margin=(max_x-min_x)*0.05; ax.set_xlim(min_x-margin,max_x+margin); ax.set_ylim(min_y-margin,max_y+margin)
    ax.set_aspect('equal'); ax.set_title(title, fontsize=18); ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); plt.grid(True, linestyle='--', alpha=0.3); plt.tight_layout(); plt.show()
    lib.write_gds(gds_filename); print(f"GDS 文件已保存到: '{gds_filename}'\n")


if __name__ == "__main__":

    ACTIVE_SIZE = 180
    ROUTING_LAYER = 5

    qd_pad_gen = QuantumDotPadGenerator()

    lib_device, cell_device, device_points = qd_pad_gen.generate_quantum_dot_layout()
    plot_gds(cell_device, title="步骤 1: 仅带长引线的量子点器件", show_plot=True)
    lib_device.write_gds("step1_device_with_leads.gds")

    lib_pads, cell_pads, pad_points, active_center, all_pads_info = qd_pad_gen.generate_pad_frame_layout(active_size=ACTIVE_SIZE, edge_margin=50, active_entry_len=40)
    visualize_and_save_pad_frame(lib_pads, cell_pads, "步骤 2: 仅带引线的焊盘框架", "step2_pad_frame.gds", show_plot=True)

    label_to_side_info = {}
    for pad_info in all_pads_info:
        label = pad_info['label']; side = pad_info['edge']
        side_pads = [p for p in all_pads_info if p['edge'] == side]
        if side in ['top', 'bottom']: side_pads.sort(key=lambda p: p['center'][0])
        else: side_pads.sort(key=lambda p: p['center'][1], reverse=True)
        side_labels = [p['label'] for p in side_pads]
        if label in side_labels:
            index = side_labels.index(label); total = len(side_labels)
            label_to_side_info[label] = (side, index, total)

    bbox = cell_device.bounding_box();
    if bbox:
        (min_x,min_y),(max_x,max_y)=bbox; cx,cy=(min_x+max_x)/2,(min_y+max_y)/2
        for poly in cell_device.polygons: poly.translate(-cx,-cy)
        for label in cell_device.labels: label.origin=(label.origin[0]-cx,label.origin[1]-cy)
        for label,pos in device_points.items(): device_points[label]=(pos[0]-cx,pos[1]-cy)
    
    assembly_lib = gdstk.Library(); assembly_cell = assembly_lib.new_cell('ASSEMBLY_CELL')
    assembly_cell.add(gdstk.Reference(cell_pads))
    device_ref = gdstk.Reference(cell_device, active_center)
    assembly_cell.add(device_ref)
    
    visualize_and_save_assembly(assembly_lib, assembly_cell, "步骤 3: 器件已放置 (未布线)", "step3_placed_unrouted.gds", show_plot=True)
    
    print("--- 步骤 4 & 5: 正在使用全新 RouteManager_v4 连接所有线... ---")
    manager = RouteManager_v4(chip_center=active_center, layer=ROUTING_LAYER)
    
    for label, local_dev_coord in device_points.items():
        if label in pad_points and label in label_to_side_info:
            pad_coord = pad_points[label]
            global_dev_coord = (local_dev_coord[0] + active_center[0], local_dev_coord[1] + active_center[1])
            side_info = label_to_side_info[label]
            manager.add_route_request(label=label, start_pt=global_dev_coord, end_pt=pad_coord, side_info=side_info)

    all_final_geometries = manager.plan_and_generate_routes()
    assembly_cell.add(*all_final_geometries)
            
    manager.check_drc()
    visualize_and_save_assembly(assembly_lib, assembly_cell, "步骤 5: 最终完整布线", "step5_fully_routed.gds", show_plot=True)

    print("全部分步流程执行完毕！")
