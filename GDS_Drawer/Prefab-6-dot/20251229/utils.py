import gdstk
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.patches as patches
import numpy as np

# ===================================================================
#                      解决中文显示问题的代码
# ===================================================================
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
# ===================================================================

def plot_gds(lib, cell, title="量子点器件布局", plot_flag=True, save_gds_flag=False, gds_filename=None): # Added lib and gds_filename
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
    ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); plt.grid(True, which='both', linestyle='--', alpha=0.3); plt.tight_layout();
    if plot_flag: # Conditional plotting
        plt.show()
    if save_gds_flag and gds_filename: # Conditional saving for plot_gds
        lib.write_gds(gds_filename); print(f"GDS 文件已保存到: '{gds_filename}'\n")

def visualize_and_save_pad_frame(lib, cell_to_show, title, gds_filename, plot_flag=True, save_gds_flag=True):
    print(f"--- {title} ---")
    if plot_flag:
        print(f"正在可视化: {title}...")
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
    if save_gds_flag:
        lib.write_gds(gds_filename); print(f"GDS 文件已保存到: '{gds_filename}'\n")

def visualize_and_save_assembly(lib, cell_to_show, title, gds_filename, plot_flag=True, save_gds_flag=True):
    print(f"--- {title} ---")

    # Always flatten the cell for consistent GDS saving, regardless of plot_flag
    flattened_cell = cell_to_show.flatten()

    if plot_flag:
        print(f"正在可视化: {title}...")
        fig, ax = plt.subplots(figsize=(18, 18))
        layer_config = {0:{'color':'#87CEEB','alpha':0.9,'label':'S/D (L0)','zorder':0}, # Changed zorder to 0
                        1:{'color':'#D3D3D3','alpha':0.5,'label':'SG (L1)','hatch':'///','zorder':7}, # SG layer
                        2:{'color':'#FF1493','alpha':0.9,'label':'BG (L2)','zorder':5},
                        3:{'color':'#8A2BE2','alpha':0.9,'label':'PG (L3)','zorder':6},
                        5:{'color':'#0000FF','alpha':1.0,'label':'Final Routing (L5)','zorder':3},
                        10:{'color':'pink','alpha':0.5,'label':'Active Area (L10)','zorder':2},
                        11:{'color':'gray','alpha':1.0,'label':'SiO2 (L11)','zorder':1}} # Changed zorder to 1
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
            polys = flattened_cell.get_polygons(layer=layer, datatype=0)
            if polys: draw_polygons(polys, {**cfg, 'layer': layer})
        pad_polys = flattened_cell.get_polygons(layer=1, datatype=0)
        if pad_polys:
            # The original logic here was trying to differentiate pads from SG based on size,
            # but layer 1 is used for both. Let's explicitly draw pads and SG if they are distinct.
            # Assuming pads are larger, and SG are smaller polygons on layer 1.
            # This might need further refinement if the distinction is not clear by size.
            sg_config = {'color':'#D3D3D3', 'alpha':0.5, 'label':'SG (L1)', 'hatch':'///', 'zorder':7};
            pad_config = {'color':'gold','alpha':0.8,'label':'Pads & Traces (L1)','zorder':1}
            for poly in pad_polys:
                (min_x, min_y), (max_x, max_y) = poly.bounding_box()
                # Heuristic: if a polygon on layer 1 is significantly large, treat it as a pad/trace
                # otherwise, treat it as an SG gate. This might need adjustment based on actual sizes.
                if (max_x - min_x) > 1 or (max_y - min_y) > 1: # Assuming pads/traces are generally larger than 1um
                    draw_polygons([poly], pad_config)
                else:
                    draw_polygons([poly], sg_config)
        handles, labels = ax.get_legend_handles_labels(); by_label = dict(zip(labels, handles)); ax.legend(by_label.values(), by_label.keys(), loc='upper right')
        bbox=flattened_cell.bounding_box(); 
        if bbox: (min_x,min_y),(max_x,max_y)=bbox; margin=(max_x-min_x)*0.05; ax.set_xlim(min_x-margin,max_x+margin); ax.set_ylim(min_y-margin,max_y+margin)
        ax.set_aspect('equal'); ax.set_title(title, fontsize=18); ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); plt.grid(True, linestyle='--', alpha=0.3); plt.tight_layout(); plt.show()
    if save_gds_flag:
        # Create a temporary library to save the flattened cell
        temp_lib = gdstk.Library()
        temp_lib.add(flattened_cell)
        temp_lib.write_gds(gds_filename)
        print(f"GDS 文件已保存到: '{gds_filename}'\n")
