import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.patches as patches
import math
import re

def plot_gds(cell):
    fig, ax = plt.subplots(figsize=(14, 10))
    
    layer_config = {
        0: {'color': '#87CEEB', 'alpha': 0.6, 'label': 'S/D (L0)', 'zorder': 2},
        1: {'color': '#D3D3D3', 'alpha': 0.5, 'label': 'SG (L1)', 'hatch': '///', 'zorder': 1},
        2: {'color': '#FF1493', 'alpha': 0.9, 'label': 'BG (L2)', 'zorder': 3}, 
        3: {'color': '#8A2BE2', 'alpha': 0.9, 'label': 'PG (L3)', 'zorder': 4}, 
    }
    
    # 绘制多边形
    for layer in [1, 0, 2, 3]: 
        polygons = cell.get_polygons(layer=layer, datatype=0)
        if not polygons: continue
        cfg = layer_config.get(layer)
        for pts in polygons:
            if hasattr(pts, 'points'): pts = pts.points
            pts = np.array(pts)
            if pts.ndim != 2 or len(pts) < 3: continue
            pts = np.concatenate((pts, [pts[0]]))
            poly = MplPolygon(pts, closed=True, 
                              facecolor=cfg['color'], edgecolor='black', linewidth=0.5,
                              alpha=cfg['alpha'], hatch=cfg.get('hatch', ''), 
                              label=cfg['label'], zorder=cfg['zorder'])
            ax.add_patch(poly)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='upper right')

    ax.set_aspect('equal')
    
    # --- [修改] 自动计算并设置视图范围 ---
    bbox = cell.bounding_box()
    if bbox is not None:
        min_x, min_y = bbox[0]
        max_x, max_y = bbox[1]
        
        # 计算宽和高
        width = max_x - min_x
        height = max_y - min_y
        
        # 添加 10% 的边距
        margin_x = width * 0.1
        margin_y = height * 0.1
        
        ax.set_xlim(min_x - margin_x, max_x + margin_x)
        ax.set_ylim(min_y - margin_y, max_y + margin_y)
        print(f"Auto-adjusting View: X[{min_x:.2f}, {max_x:.2f}], Y[{min_y:.2f}, {max_y:.2f}]")
    else:
        # Fallback if cell is empty
        ax.autoscale_view()
    
    ax.set_title('Final QD Layout V6: Fully Parametrized SGs')
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')
    plt.grid(True, which='both', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()

def create_6qd_layout_with_labels(
    # --- 几何尺寸输入 ---
    pg_max_width,
    pg_vert_side_len,
    pg_chamfer_h,
    pg_bot_flat_w,
    pg_top_flat_w,
    bg_max_width,
    bg_vert_side_len,
    bg_top_flat_w,
    bg_bot_flat_w,
    # --- 间距参数输入 ---
    gap_pg_bg,
    d1_gap,
    d2_gap,
    sd_gap_to_gate,
    gap_gate_outer_sg,
    # --- SG 参数输入 ---
    sg_mid_thick,   
    sg_top_thick,   
    sg_bot_thick,   
    sg_extension,   
    # --- 引线参数输入 ---
    lead_width,
    lead_length_bot,
    lead_length_top,
    lead_overlap,
    # --- S/D 参数输入 ---
    sd_height,
    sd_width,       # 用于中间 S/D 的宽度
    sd_lead_width,
    sd_outer_length
):
    # ==========================================
    # 1. 参数定义
    # ==========================================
    
    # 层定义
    LAYER_SD = 0
    LAYER_SG = 1
    LAYER_BG = 2
    LAYER_PG = 3
    LAYER_LABEL = 100  # 用于放置标签的专用层
    
    # --- 纵向布局 ---
    gate_half_height = pg_vert_side_len/2 + pg_chamfer_h 
    row_pitch = gate_half_height + d1_gap + sg_mid_thick + d2_gap + gate_half_height
    
    lib = gdstk.Library()
    cell = lib.new_cell('6QD_DEVICE_LABELED')

    # ==========================================
    # 2. 形状生成函数 (无变化)
    # ==========================================
    
    def get_pg_points(center):
        cx, cy = center
        y_v_top = pg_vert_side_len / 2; y_v_bot = -pg_vert_side_len / 2
        y_top = y_v_top + pg_chamfer_h; y_bot = y_v_bot - pg_chamfer_h
        x_max = pg_max_width / 2; x_top_flat = pg_top_flat_w / 2; x_bot_flat = pg_bot_flat_w / 2
        pts = [
            (cx + x_top_flat, cy + y_top), (cx + x_max, cy + y_v_top),
            (cx + x_max, cy + y_v_bot), (cx + x_bot_flat, cy + y_bot),
            (cx - x_bot_flat, cy + y_bot), (cx - x_max, cy + y_v_bot),
            (cx - x_max, cy + y_v_top), (cx - x_top_flat, cy + y_top)
        ]
        return pts, y_top, y_bot

    def get_bg_points(center):
        cx, cy = center
        y_v_top = bg_vert_side_len / 2; y_v_bot = -bg_vert_side_len / 2
        bg_chamfer_h_top = (bg_max_width - bg_top_flat_w) / 2
        bg_chamfer_h_bot = (bg_max_width - bg_bot_flat_w) / 2
        y_top = y_v_top + bg_chamfer_h_top; y_bot = y_v_bot - bg_chamfer_h_bot
        x_max = bg_max_width / 2; x_top_flat = bg_top_flat_w / 2; x_bot_flat = bg_bot_flat_w / 2
        pts = [
            (cx + x_top_flat, cy + y_top), (cx + x_max, cy + y_v_top),
            (cx + x_max, cy + y_v_bot), (cx + x_bot_flat, cy + y_bot),
            (cx - x_bot_flat, cy + y_bot), (cx - x_max, cy + y_v_bot),
            (cx - x_max, cy + y_v_top), (cx - x_top_flat, cy + y_top)
        ]
        return pts, y_top, y_bot

    # ==========================================
    # 计算引线的全局统一终点 Y 坐标 (无变化)
    # ==========================================
    pg_local_bot_y = -(pg_vert_side_len / 2 + pg_chamfer_h) 
    bg_local_bot_y = -(bg_vert_side_len / 2 + (bg_max_width - bg_bot_flat_w) / 2)
    global_min_gate_local_y = min(pg_local_bot_y, bg_local_bot_y)
    fixed_lead_bot_end_y = 0 + global_min_gate_local_y - lead_length_bot
    
    y_center_top = row_pitch
    pg_local_top_y = (pg_vert_side_len / 2 + pg_chamfer_h)
    bg_local_top_y = (bg_vert_side_len / 2 + (bg_max_width - bg_top_flat_w) / 2)
    global_max_gate_local_y = max(pg_local_top_y, bg_local_top_y)
    fixed_lead_top_end_y = y_center_top + global_max_gate_local_y + lead_length_top

    # ==========================================
    # 3. 布局生成：Bottom Row (无变化)
    # ==========================================
    
    all_shapes = []
    all_labels = []
    y_center_bot = 0
    current_x = 0
    
    bottom_gate_names = []
    for i in range(1, 7):
        bottom_gate_names.append(f'QD_B{i}')
        bottom_gate_names.append(f'QD_PG{i}')
    bottom_gate_names.append('QD_B7')

    bot_tips_y = [] 
    
    for i, name in enumerate(bottom_gate_names):
        g_type = 'PG' if 'PG' in name else 'BG'
        w = pg_max_width if g_type == 'PG' else bg_max_width
        layer = LAYER_PG if g_type == 'PG' else LAYER_BG

        if i == 0:
            current_x = 0
        else:
            prev_name = bottom_gate_names[i-1]
            prev_w = pg_max_width if 'PG' in prev_name else bg_max_width
            current_x += (prev_w / 2) + gap_pg_bg + (w / 2)
            
        center = (current_x, y_center_bot)
        
        if g_type == 'PG': pts, tip_y, base_y = get_pg_points(center)
        else: pts, tip_y, base_y = get_bg_points(center)
            
        bot_tips_y.append(tip_y)
        poly = gdstk.Polygon(pts, layer=layer)
        
        lead_y_start = base_y + lead_overlap 
        lead_y_end = fixed_lead_bot_end_y 
        
        lead = gdstk.rectangle((current_x - lead_width/2, lead_y_start),
                               (current_x + lead_width/2, lead_y_end),
                               layer=layer)
        
        all_shapes.extend(gdstk.boolean([poly, lead], [], 'or', layer=layer))
        label = gdstk.Label(name, (current_x, lead_y_end), layer=LAYER_LABEL)
        all_labels.append(label)

    bot_min_x = 0 - bg_max_width/2
    bot_max_x = current_x + bg_max_width/2
    
    # ==========================================
    # 4. 布局生成：Top Row (SETs) (无变化)
    # ==========================================
    
    center_x_global = (bot_min_x + bot_max_x) / 2
    inner_gap_required = 2 * sd_gap_to_gate + sd_width 
    offset_from_center = inner_gap_required/2 + bg_max_width/2 + gap_pg_bg + pg_max_width/2
    
    left_set_center_pg_x = center_x_global - offset_from_center
    right_set_center_pg_x = center_x_global + offset_from_center

    def create_set_at(center_pg_x, side):
        shapes = []
        labels = []
        offsets = [
            -(pg_max_width/2 + gap_pg_bg + bg_max_width/2), 
            0, 
            (pg_max_width/2 + gap_pg_bg + bg_max_width/2)
        ]
        
        if side == 'L': set_gate_names = ['SET1_B1', 'SET1_G', 'SET1_B2']
        else: set_gate_names = ['SET2_B1', 'SET2_G', 'SET2_B2']

        for offset, name in zip(offsets, set_gate_names):
            cx = center_pg_x + offset
            g_type = 'PG' if '_G' in name else 'BG'
            layer = LAYER_PG if g_type == 'PG' else LAYER_BG
            
            if g_type == 'PG': pts, _, _ = get_pg_points((cx, y_center_top))
            else: pts, _, _ = get_bg_points((cx, y_center_top))
            
            poly = gdstk.Polygon(pts, layer=layer)
            Wire_top_y = np.max([pt[1] for pt in pts])
            lead_y_start = Wire_top_y - lead_overlap
            lead_y_end = fixed_lead_top_end_y
            
            lead = gdstk.rectangle((cx - lead_width/2, lead_y_start), (cx + lead_width/2, lead_y_end), layer=layer)
            
            shapes.extend(gdstk.boolean([poly, lead], [], 'or', layer=layer))
            labels.append(gdstk.Label(name, (cx, lead_y_end), layer=LAYER_LABEL))
            
        return shapes, labels

    left_set_shapes, left_set_labels = create_set_at(left_set_center_pg_x, side='L')
    all_shapes.extend(left_set_shapes); all_labels.extend(left_set_labels)
    
    right_set_shapes, right_set_labels = create_set_at(right_set_center_pg_x, side='R')
    all_shapes.extend(right_set_shapes); all_labels.extend(right_set_labels)
    
    # ==========================================
    # 5. Screening Gates (SGs) (无变化)
    # ==========================================
    
    sg_x_start = bot_min_x - sg_extension; sg_x_end = bot_max_x + sg_extension
    y_sg_mid_bot = max(bot_tips_y) + d1_gap; y_sg_mid_top = y_sg_mid_bot + sg_mid_thick
    sg_mid = gdstk.rectangle((sg_x_start, y_sg_mid_bot), (sg_x_end, y_sg_mid_top), layer=LAYER_SG)
    
    y_sg_top_bot = y_center_top + gate_half_height + gap_gate_outer_sg; y_sg_top_top = y_sg_top_bot + sg_top_thick
    sg_top = gdstk.rectangle((sg_x_start, y_sg_top_bot), (sg_x_end, y_sg_top_top), layer=LAYER_SG)
                          
    y_sg_bot_top = y_center_bot - gate_half_height - gap_gate_outer_sg; y_sg_bot_bot = y_sg_bot_top - sg_bot_thick
    sg_bot = gdstk.rectangle((sg_x_start, y_sg_bot_bot), (sg_x_end, y_sg_bot_top), layer=LAYER_SG)
    
    all_shapes.extend([sg_top, sg_mid, sg_bot])
    all_labels.append(gdstk.Label("SG1", (sg_x_start, (y_sg_top_bot + y_sg_top_top) / 2), layer=LAYER_LABEL))
    all_labels.append(gdstk.Label("SG2", (sg_x_start, (y_sg_mid_bot + y_sg_mid_top) / 2), layer=LAYER_LABEL))
    all_labels.append(gdstk.Label("SG3", (sg_x_start, (y_sg_bot_bot + y_sg_bot_top) / 2), layer=LAYER_LABEL))

    # ==========================================
    # 6. S/D Electrodes
    # ==========================================
    
    sd_shapes = []
    sd_corner_r = 0.020 

    def create_rounded_sd(x_min, y_min, x_max, y_max, corners_to_round):
        poly = gdstk.Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max), (x_min, y_max)], layer=LAYER_SD)
        poly.fillet([sd_corner_r if c else 0 for c in corners_to_round], tolerance=0.001)
        return poly
    
    bot_left_bg_edge = bot_min_x; bot_left_sd_inner_x = bot_left_bg_edge - sd_gap_to_gate
    bot_right_bg_edge = bot_max_x; bot_right_sd_inner_x = bot_right_bg_edge + sd_gap_to_gate
    global_outer_left_x = bot_left_sd_inner_x - sd_outer_length
    global_outer_right_x = bot_right_sd_inner_x + sd_outer_length
    
    sd_shapes.append(create_rounded_sd(global_outer_left_x, y_center_bot - sd_height/2, bot_left_sd_inner_x, y_center_bot + sd_height/2, [False, True, True, False]))
    all_labels.append(gdstk.Label("QD_D", (global_outer_left_x, y_center_bot), layer=LAYER_LABEL))
        
    sd_shapes.append(create_rounded_sd(bot_right_sd_inner_x, y_center_bot - sd_height/2, global_outer_right_x, y_center_bot + sd_height/2, [True, False, False, True]))
    all_labels.append(gdstk.Label("QD_S", (global_outer_right_x, y_center_bot), layer=LAYER_LABEL))
        
    set_half_width = bg_max_width + gap_pg_bg + pg_max_width/2
    l_set_left_edge = left_set_center_pg_x - set_half_width; top_left_sd_inner_x = l_set_left_edge - sd_gap_to_gate
    r_set_right_edge = right_set_center_pg_x + set_half_width; top_right_sd_inner_x = r_set_right_edge + sd_gap_to_gate
    
    sd_shapes.append(create_rounded_sd(global_outer_left_x, y_center_top - sd_height/2, top_left_sd_inner_x, y_center_top + sd_height/2, [False, True, True, False]))
    all_labels.append(gdstk.Label("SET1_D", (global_outer_left_x, y_center_top), layer=LAYER_LABEL))
        
    sd_shapes.append(create_rounded_sd(top_right_sd_inner_x, y_center_top - sd_height/2, global_outer_right_x, y_center_top + sd_height/2, [True, False, False, True]))
    all_labels.append(gdstk.Label("SET2_D", (global_outer_right_x, y_center_top), layer=LAYER_LABEL))
    
    # --- Top Middle S/D ---
    l_set_right_edge = left_set_center_pg_x + set_half_width; s_mid_x_min = l_set_right_edge + sd_gap_to_gate
    r_set_left_edge = right_set_center_pg_x - set_half_width; s_mid_x_max = r_set_left_edge - sd_gap_to_gate
    
    mid_sd_poly = create_rounded_sd(s_mid_x_min, y_center_top - sd_height/2, s_mid_x_max, y_center_top + sd_height/2, [True, True, True, True])
    
    mid_sd_center_x = (s_mid_x_min + s_mid_x_max) / 2
    mid_lead_y_end = y_center_top + lead_length_top
    mid_sd_lead = gdstk.rectangle((mid_sd_center_x - sd_lead_width/2, y_center_top), (mid_sd_center_x + sd_lead_width/2, mid_lead_y_end), layer=LAYER_SD)
    
    sd_shapes.extend(gdstk.boolean([mid_sd_poly, mid_sd_lead], [], 'or', layer=LAYER_SD))
    
    ### 已修改 ###
    # 两个SET共享一个公共源，因此我们应用一个统一的标签
    all_labels.append(gdstk.Label("SET_S", (mid_sd_center_x, mid_lead_y_end), layer=LAYER_LABEL))

    all_shapes.extend(sd_shapes)
    cell.add(*all_shapes)
    cell.add(*all_labels)

    return lib, cell

def create_rect_wire_layout(
    N=6,
    layout_width=1400,
    layout_height=1400,
    pad_width=100,
    pad_height=100,      
    pad_spacing=30,
    edge_margin=100,
    active_size=200,     
    trace_width=10,
    trace_spacing=10,    
    taper_length=40,
    active_entry_len=100,
    gds_filename='qd_layout_rect_wire.gds'
):
    """
    Generates a Pad Layout where wires entering the active region are 
    rectangular polylines (constant width, flat ends), removing any arrow/taper effects at the center.
    """

    # ==========================================
    # 1. Pad Generation (Symmetric Logic)
    # ==========================================
    def natural_sort_key(s):
        return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

    num_sets = math.ceil(N / 3)
    pg_labels = sorted([f"QD_PG{i+1}" for i in range(N)], key=natural_sort_key)
    bg_labels = sorted([f"QD_B{i+1}" for i in range(N + 1)], key=natural_sort_key)

    gate_sequence = []
    for i in range(N):
        gate_sequence.append(bg_labels[i])
        gate_sequence.append(pg_labels[i])
    gate_sequence.append(bg_labels[N])

    top_pads = []
    if num_sets >= 1: top_pads.extend(["SET1_B1", "SET1_B2", "SET1_G"])
    top_pads.append("SET_S")
    if num_sets >= 2: top_pads.extend(["SET2_G", "SET2_B2", "SET2_B1"])

    left_fixed = ["SG1"]
    if num_sets >= 1: left_fixed.append("SET1_D") 
    left_fixed.extend(["SG2", "QD_D"])
    
    right_fixed = []
    if num_sets >= 2: right_fixed.append("SET2_D")
    right_fixed.append("QD_S")
    right_fixed.extend(["GND", "SG3", "GND"]) 

    # Symmetry Logic
    len_seq = len(gate_sequence)
    len_L_fix = len(left_fixed)
    len_R_fix = len(right_fixed)
    len_Top = len(top_pads)
    
    available_w = layout_width - 2 * edge_margin
    best_config = (0, len_seq, 0)
    min_diff = float('inf')
    
    start_side = max(len_L_fix, len_R_fix)
    for target_side in range(start_side, start_side + len_seq):
        add_L = target_side - len_L_fix
        add_R = target_side - len_R_fix
        if add_L + add_R > len_seq: break
        n_bot = len_seq - (add_L + add_R)
        if n_bot > 0 and (n_bot * pad_width + (n_bot-1)*pad_spacing > available_w): continue
        diff = abs(n_bot - len_Top)
        if diff < min_diff:
            min_diff = diff
            best_config = (add_L, n_bot, add_R)
            
    n_L, n_B, n_R = best_config
    pads_left_extra = gate_sequence[:n_L]
    pads_bottom = gate_sequence[n_L : n_L + n_B]
    pads_right_extra = gate_sequence[n_L + n_B :]

    left_pads = left_fixed + pads_left_extra
    right_pads = right_fixed + list(reversed(pads_right_extra)) # 右侧的顺序反转

    print(f"Pad Distribution: Left={left_pads},\n Bottom={pads_bottom},\n Right={right_pads},\n Top={top_pads}")

    all_pads = [] 
    def assign_coords(pads, edge):
        count = len(pads)
        if count == 0: return
        if edge in ['top', 'bottom']:
            total_span = count * pad_width + (count - 1) * pad_spacing
            start_x = (layout_width - total_span) / 2
            y = layout_height - edge_margin - pad_height if edge == 'top' else edge_margin
            for i, label in enumerate(pads):
                x = start_x + i * (pad_width + pad_spacing)
                cx, cy = x + pad_width/2, y + pad_height/2
                all_pads.append({'label': label, 'rect': (x, y, pad_width, pad_height), 'center': (cx, cy), 'edge': edge})
        else:
            total_span = count * pad_width + (count - 1) * pad_spacing
            top_y = (layout_height + total_span) / 2
            x = edge_margin if edge == 'left' else layout_width - edge_margin - pad_height
            for i, label in enumerate(pads):
                y = top_y - (i+1)*pad_width - i*pad_spacing
                cx, cy = x + pad_height/2, y + pad_width/2
                all_pads.append({'label': label, 'rect': (x, y, pad_height, pad_width), 'center': (cx, cy), 'edge': edge})

    assign_coords(top_pads, 'top')
    assign_coords(pads_bottom, 'bottom')
    assign_coords(left_pads, 'left')
    assign_coords(right_pads, 'right')

    # ==========================================
    # 2. Target Assignment
    # ==========================================
    center_x = layout_width / 2
    center_y = layout_height / 2
    aa_x1 = center_x - active_size/2
    aa_x2 = center_x + active_size/2
    aa_y1 = center_y - active_size/2
    aa_y2 = center_y + active_size/2

    pads_by_edge = {'top': [], 'bottom': [], 'left': [], 'right': []}
    for p in all_pads: pads_by_edge[p['edge']].append(p)

    pads_by_edge['top'].sort(key=lambda p: p['center'][0])
    pads_by_edge['bottom'].sort(key=lambda p: p['center'][0])
    pads_by_edge['left'].sort(key=lambda p: p['center'][1], reverse=True)
    pads_by_edge['right'].sort(key=lambda p: p['center'][1], reverse=True)

    pitch = trace_width + trace_spacing
    
    def assign_targets(pads, edge):
        n = len(pads)
        if n == 0: return
        indices = np.arange(n) - (n - 1) / 2
        offsets = indices * pitch
        
        if edge == 'top':
            for i, p in enumerate(pads):
                p['target'] = (center_x + offsets[i], aa_y2)
                p['target_edge'] = 'top'
        elif edge == 'bottom':
            for i, p in enumerate(pads):
                p['target'] = (center_x + offsets[i], aa_y1)
                p['target_edge'] = 'bottom'
        elif edge == 'left':
            for i, p in enumerate(pads):
                p['target'] = (aa_x1, center_y + offsets[-(i+1)])
                p['target_edge'] = 'left'
        elif edge == 'right':
            for i, p in enumerate(pads):
                p['target'] = (aa_x2, center_y + offsets[-(i+1)])
                p['target_edge'] = 'right'

    for edge in ['top', 'bottom', 'left', 'right']:
        assign_targets(pads_by_edge[edge], edge)

    # ==========================================
    # 3. GDS Generation (Updated Routing)
    # ==========================================
    lib = gdstk.Library()
    cell = lib.new_cell('QD_RECT_WIRE')
    
    # [修改点 1] 设置 Layer 定义
    LAYER_SUBSTRATE = 0
    LAYER_PAD = 1
    LAYER_TRACE = 1    # 将 Trace 设为与 Pad 同层 (原为3，现改为1)
    LAYER_TEXT = 2
    LAYER_ACTIVE = 10
    
    # [修改点 2] 定义 SiO2 层及可调参数
    LAYER_SIO2 = 11    # 新增 SiO2 层 ID
    sio2_margin = 10   # SiO2 层比 Pad 边缘大的距离 (um) - 可在此处调节

    cell.add(gdstk.rectangle((0, 0), (layout_width, layout_height), layer=LAYER_SUBSTRATE))
    cell.add(gdstk.rectangle((aa_x1, aa_y1), (aa_x2, aa_y2), layer=LAYER_ACTIVE))

    def create_taper(p_rect, edge, t_len, t_width):
        x, y, w, h = p_rect
        cx, cy = x + w/2, y + h/2
        if edge == 'top':
            pts = [(x, y), (x+w, y), (cx+t_width/2, y-t_len), (cx-t_width/2, y-t_len)]
            start = (cx, y - t_len)
        elif edge == 'bottom':
            pts = [(x, y+h), (x+w, y+h), (cx+t_width/2, y+h+t_len), (cx-t_width/2, y+h+t_len)]
            start = (cx, y + h + t_len)
        elif edge == 'left':
            pts = [(x+w, y), (x+w, y+h), (x+w+t_len, cy+t_width/2), (x+w+t_len, cy-t_width/2)]
            start = (x+w+t_len, cy)
        elif edge == 'right':
            pts = [(x, y), (x, y+h), (x-t_len, cy+t_width/2), (x-t_len, cy-t_width/2)]
            start = (x-t_len, cy)
        return gdstk.Polygon(pts, layer=LAYER_TRACE), start

    def route_safe_z(start_pt, target_pt, edge, width):
        """
        Calculates points for Safe Z-Routing.
        Returns a gdstk.FlexPath with 'flush' ends (rectangular).
        """
        x1, y1 = start_pt
        x2, y2 = target_pt
        
        safe_margin = 20 
        pts = [(x1, y1)]
        
        if edge in ['top', 'bottom']:
            is_top = (edge == 'top')
            sign_y = -1 if is_top else 1
            y_safe = y1 + sign_y * safe_margin
            pts.append((x1, y_safe))
            y_entry = y2 + sign_y * active_entry_len 
            
            avail_h = abs(y_entry - y_safe)
            needed_w = abs(x2 - x1)
            sign_x = 1 if x2 > x1 else -1
            
            if needed_w > avail_h:
                h_leg = avail_h / 2
                y_mid = y_safe + sign_y * h_leg
                x_mid_1 = x1 + sign_x * h_leg
                x_mid_2 = x2 - sign_x * h_leg
                pts.append((x_mid_1, y_mid))
                pts.append((x_mid_2, y_mid))
            else:
                dy_diag = needed_w
                y_turn_end = y_safe + sign_y * dy_diag
                pts.append((x2, y_turn_end))
            pts.append((x2, y2))
                
        else: # Left / Right
            is_left = (edge == 'left')
            sign_x = 1 if is_left else -1
            x_safe = x1 + sign_x * safe_margin
            pts.append((x_safe, y1))
            x_entry = x2 + sign_x * active_entry_len
            
            avail_w = abs(x_entry - x_safe)
            needed_h = abs(y2 - y1)
            sign_y = 1 if y2 > y1 else -1
            
            if needed_h > avail_w:
                w_leg = avail_w / 2
                x_mid = x_safe + sign_x * w_leg
                y_mid_1 = y1 + sign_y * w_leg
                y_mid_2 = y2 - sign_y * w_leg
                pts.append((x_mid, y_mid_1))
                pts.append((x_mid, y_mid_2))
            else:
                dx_diag = needed_h
                x_turn_end = x_safe + sign_x * dx_diag
                pts.append((x_turn_end, y2))
            pts.append((x2, y2))
            
        path = gdstk.FlexPath(pts, width, ends="flush", joins="round", layer=LAYER_TRACE)
        return path, pts

    # Loop
    for p in all_pads:
        x, y, w, h = p['rect']
        label = p['label']
        edge = p['edge']
        
        # [修改点 3] 绘制 SiO2 层 (比 Pad 大 sio2_margin)
        cell.add(gdstk.rectangle(
            (x - sio2_margin, y - sio2_margin), 
            (x + w + sio2_margin, y + h + sio2_margin), 
            layer=LAYER_SIO2
        ))
        
        cell.add(gdstk.rectangle((x, y), (x + w, y + h), layer=LAYER_PAD))
        
        cx, cy = x + w/2, y + h/2
        rot_rad = 0 if w > h else math.pi / 2
        cell.add(gdstk.Label(label, (cx, cy), layer=LAYER_TEXT, magnification=10, rotation=rot_rad))
        
        poly_taper, trace_start = create_taper(p['rect'], edge, taper_length, trace_width)
        cell.add(poly_taper)
        
        path_obj, path_pts = route_safe_z(trace_start, p['target'], edge, trace_width)
        cell.add(path_obj)
        p['trace_points'] = path_pts

    lib.write_gds(gds_filename)
    print(f"GDS file saved to: {gds_filename}")

    # Visualization
    fig, ax = plt.subplots(figsize=(15, 15))
    margin_view = 200
    ax.set_xlim(-margin_view, layout_width + margin_view)
    ax.set_ylim(-margin_view, layout_height + margin_view)
    ax.set_aspect('equal')
    
    chip = patches.Rectangle((0, 0), layout_width, layout_height, lw=2, ec='black', fc='#f8f8f8')
    ax.add_patch(chip)
    active_area = patches.Rectangle((aa_x1, aa_y1), active_size, active_size, lw=1, ec='red', fc='none', zorder=10)
    ax.add_patch(active_area)

    for p in all_pads:
        x, y, w, h = p['rect']
        label = p['label']
        fc = 'gold'
        if 'GND' in label: fc = 'lightgreen'
        if 'SG' in label: fc = 'orange'
        
        # 简单可视化 SiO2 层 (灰色虚线)
        sio2_rect = patches.Rectangle(
            (x - sio2_margin, y - sio2_margin), 
            w + 2*sio2_margin, h + 2*sio2_margin, 
            linewidth=1, edgecolor='gray', linestyle='--', facecolor='none', zorder=4
        )
        ax.add_patch(sio2_rect)

        rect = patches.Rectangle((x, y), w, h, linewidth=1, edgecolor='black', facecolor=fc, zorder=5)
        ax.add_patch(rect)
        t_rot = 90 if w < h else 0
        ax.text(x + w/2, y + h/2, label, ha='center', va='center', fontsize=6, rotation=t_rot, zorder=6)
        
        if 'trace_points' in p:
            pts = p['trace_points']
            line_x, line_y = zip(*pts)
            ax.plot(line_x, line_y, color='#1E90FF', linewidth=5, zorder=4)
            _, t_start = create_taper(p['rect'], p['edge'], taper_length, trace_width)
            cx, cy = x + w/2, y + h/2
            ax.plot([cx, t_start[0]], [cy, t_start[1]], color='#1E90FF', linewidth=10, zorder=4)

    plt.title(f"Rectangular Wire Ends (N={N})\nTrace & Pad on Layer {LAYER_PAD}, SiO2 on Layer {LAYER_SIO2} (+{sio2_margin}um)", fontsize=16, pad=30)
    plt.axis('off')
    plt.tight_layout()
    plt.show()