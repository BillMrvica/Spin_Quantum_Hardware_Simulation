import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

def create_final_layout_v4():
    # ==========================================
    # 1. 参数定义 (单位: um)
    # ==========================================
    
    # 层定义
    LAYER_SD = 0
    LAYER_SG = 1
    LAYER_BG = 2
    LAYER_PG = 3
    
    # --- 几何尺寸 ---
    pg_max_width = 0.120
    pg_vert_side_len = 0.040
    pg_chamfer_h = 0.040 
    pg_bot_flat_w = 0.042
    pg_top_flat_w = 0.040
    
    bg_max_width = 0.060
    bg_vert_side_len = 0.040
    bg_top_flat_w = 0.020
    bg_bot_flat_w = 0.042
    
    # --- 间距参数 ---
    gap_pg_bg = 0.002       # 栅极横向间距
    d1_gap = 0.030          # Bottom Gate Tip 到 Mid SG
    d2_gap = 0.030          # Mid SG 到 Top Gate Tip
    sd_gap_to_gate = 0.010  # [修改] S/D 到 BG 的间距，改为 100nm
    
    sg_mid_thick = 0.030    # 中间 SG 厚度
    
    # --- 纵向布局计算 ---
    # 半高 = 垂直边/2 + 倒角高度
    gate_half_height = pg_vert_side_len/2 + pg_chamfer_h # 0.06
    
    # Row Pitch 计算
    row_pitch = gate_half_height + d1_gap + sg_mid_thick + d2_gap + gate_half_height
    # Pitch ≈ 0.21 um
    
    # --- 引线参数 (优化连接) ---
    lead_width = 0.042
    lead_length_bot = 0.6  
    lead_length_top = 0.8   # [修改] 稍微加长顶部引线
    lead_overlap = 0.020    # [修改] 引线重叠量，确保"连在一起"无缝隙
    
    # --- S/D 参数 ---
    sd_height = 0.080      
    sd_width = 0.4         # 标准 S/D 宽度
    
    lib = gdstk.Library()
    cell = lib.new_cell('QD_DEVICE_V4')

    # ==========================================
    # 2. 形状生成函数
    # ==========================================
    
    def get_pg_points(center, flip_y=False):
        cx, cy = center
        y_v_top = pg_vert_side_len / 2
        y_v_bot = -pg_vert_side_len / 2
        y_top = y_v_top + pg_chamfer_h
        y_bot = y_v_bot - pg_chamfer_h
        
        x_max = pg_max_width / 2
        x_top_flat = pg_top_flat_w / 2
        x_bot_flat = pg_bot_flat_w / 2
        
        pts = [
            (cx + x_top_flat, cy + y_top), (cx + x_max, cy + y_v_top),
            (cx + x_max, cy + y_v_bot), (cx + x_bot_flat, cy + y_bot),
            (cx - x_bot_flat, cy + y_bot), (cx - x_max, cy + y_v_bot),
            (cx - x_max, cy + y_v_top), (cx - x_top_flat, cy + y_top)
        ]
        
        if flip_y:
            pts = [(x, 2*cy - y) for x, y in pts]
            visual_tip = 2*cy - y_bot 
            visual_base = 2*cy - y_top # Base connects lead
            return pts, visual_tip, visual_base 
        else:
            return pts, y_top, y_bot

    def get_bg_points(center, flip_y=False):
        cx, cy = center
        y_v_top = bg_vert_side_len / 2
        y_v_bot = -bg_vert_side_len / 2
        
        bg_chamfer_h_top = (bg_max_width - bg_top_flat_w) / 2
        bg_chamfer_h_bot = (bg_max_width - bg_bot_flat_w) / 2
        y_top = y_v_top + bg_chamfer_h_top
        y_bot = y_v_bot - bg_chamfer_h_bot
        
        x_max = bg_max_width / 2
        x_top_flat = bg_top_flat_w / 2
        x_bot_flat = bg_bot_flat_w / 2
        
        pts = [
            (cx + x_top_flat, cy + y_top), (cx + x_max, cy + y_v_top),
            (cx + x_max, cy + y_v_bot), (cx + x_bot_flat, cy + y_bot),
            (cx - x_bot_flat, cy + y_bot), (cx - x_max, cy + y_v_bot),
            (cx - x_max, cy + y_v_top), (cx - x_top_flat, cy + y_top)
        ]
        
        if flip_y:
            pts = [(x, 2*cy - y) for x, y in pts]
            return pts, 2*cy - y_bot, 2*cy - y_top
        else:
            return pts, y_top, y_bot

    # ==========================================
    # 3. 布局生成：Bottom Row
    # ==========================================
    
    all_shapes = []
    y_center_bot = 0
    current_x = 0
    
    elements = []
    for _ in range(6):
        elements.extend([('BG', bg_max_width, LAYER_BG), ('PG', pg_max_width, LAYER_PG)])
    elements.append(('BG', bg_max_width, LAYER_BG))
    
    bot_tips_y = [] 
    
    for i, (g_type, w, layer) in enumerate(elements):
        if i == 0: current_x = 0
        else:
            prev_w = elements[i-1][1]
            current_x += (prev_w / 2) + gap_pg_bg + (w / 2)
            
        center = (current_x, y_center_bot)
        if g_type == 'PG': pts, tip_y, base_y = get_pg_points(center, flip_y=False)
        else: pts, tip_y, base_y = get_bg_points(center, flip_y=False)
        
        bot_tips_y.append(tip_y)
        
        poly = gdstk.Polygon(pts, layer=layer, datatype=0)
        
        # [修改] 引线向下，增加 overlap 确保连接
        lead_y_start = base_y + lead_overlap 
        lead_y_end = base_y - lead_length_bot
        lead = gdstk.rectangle((current_x - lead_width/2, lead_y_start),
                               (current_x + lead_width/2, lead_y_end),
                               layer=layer, datatype=0)
        
        # Boolean 合并
        all_shapes.extend(gdstk.boolean([poly, lead], [], 'or', layer=layer, datatype=0))

    bot_min_x = 0 - bg_max_width/2
    bot_max_x = current_x + bg_max_width/2
    
    # ==========================================
    # 4. 布局生成：Top Row (SETs)
    # ==========================================
    
    y_center_top = row_pitch
    center_x_global = (bot_min_x + bot_max_x) / 2
    
    # [修改] 计算 SET 间距以容纳中间 S/D
    # 需求: LeftBG_edge --(100nm)-- SD(150nm) --(100nm)-- RightBG_edge
    # 总间隙 = 0.1 + 0.15 + 0.1 = 0.35 um (两个内侧BG边缘的距离)
    
    # 单个 SET 视觉宽度 (BG-PG-BG)
    set_visual_half_w = (bg_max_width + gap_pg_bg + pg_max_width + gap_pg_bg + bg_max_width) / 2
    
    # 计算两个 PG 中心的距离
    # Distance = (Visual_Half_W - BG/2) + Gap_Between_BGs + (Visual_Half_W - BG/2)
    # 实际上 Gap_Between_BGs = 0.35
    # Left_PG_x = Global_Center - (0.35/2 + BG/2 + Gap + PG/2)
    
    inner_gap_required = 2 * sd_gap_to_gate + sd_width # 0.1*2 + 0.15 = 0.35
    
    # 偏移量计算
    offset_from_center = inner_gap_required/2 + bg_max_width/2 + gap_pg_bg + pg_max_width/2
    
    left_set_center_pg_x = center_x_global - offset_from_center
    right_set_center_pg_x = center_x_global + offset_from_center

    top_tips_y = []

    def create_set_at(center_pg_x):
        shapes = []
        offsets = [
            -(pg_max_width/2 + gap_pg_bg + bg_max_width/2), 
            0, 
            (pg_max_width/2 + gap_pg_bg + bg_max_width/2)
        ]
        types = [('BG', LAYER_BG), ('PG', LAYER_PG), ('BG', LAYER_BG)]
        
        for offset, (gname, lay) in zip(offsets, types):
            cx = center_pg_x + offset
            if gname == 'PG': pts, tip_y, base_y = get_pg_points((cx, y_center_top), flip_y=True)
            else: pts, tip_y, base_y = get_bg_points((cx, y_center_top), flip_y=True)
            
            top_tips_y.append(tip_y)
            poly = gdstk.Polygon(pts, layer=lay, datatype=0)
            
            # [修改] 引线向上，增加 overlap
            lead_y_start = base_y - lead_overlap
            lead_y_end = base_y + lead_length_top
            lead = gdstk.rectangle((cx - lead_width/2, lead_y_start),
                                   (cx + lead_width/2, lead_y_end),
                                   layer=lay, datatype=0)
            
            shapes.extend(gdstk.boolean([poly, lead], [], 'or', layer=lay, datatype=0))
        return shapes

    all_shapes.extend(create_set_at(left_set_center_pg_x))
    all_shapes.extend(create_set_at(right_set_center_pg_x))
    
    cell.add(*all_shapes)

    # ==========================================
    # 5. Screening Gates (SGs)
    # ==========================================
    
    sg_ext = 0.6
    sg_x_start = bot_min_x - sg_ext
    sg_x_end = bot_max_x + sg_ext
    
    # Mid SG
    max_bot_tip = max(bot_tips_y)
    y_sg_mid_bot = max_bot_tip + d1_gap
    y_sg_mid_top = y_sg_mid_bot + sg_mid_thick
    sg2 = gdstk.rectangle((sg_x_start, y_sg_mid_bot), (sg_x_end, y_sg_mid_top), 
                          layer=LAYER_SG, datatype=0)
    
    # Top SG (SG1)
    top_gate_base_y = y_center_top + gate_half_height
    y_sg1_bot = top_gate_base_y + 0.04
    y_sg1_top = y_sg1_bot + 0.3
    sg1 = gdstk.rectangle((sg_x_start, y_sg1_bot), (sg_x_end, y_sg1_top),
                          layer=LAYER_SG, datatype=0)
                          
    # Bottom SG
    bot_gate_base_y = y_center_bot - gate_half_height
    y_sg_bot_top = bot_gate_base_y - 0.04
    y_sg_bot_bot = y_sg_bot_top - 0.3
    sg_bot = gdstk.rectangle((sg_x_start, y_sg_bot_bot), (sg_x_end, y_sg_bot_top),
                             layer=LAYER_SG, datatype=0)
                             
    cell.add(sg1, sg2, sg_bot)

    # ==========================================
    # 6. S/D Electrodes (修正间距)
    # ==========================================
    
    sd_shapes = []
    
    # --- Bottom Row S/D ---
    # Left S: 右边缘 = 最左BG左边缘 - 100nm
    bot_left_bg_edge = bot_min_x 
    s_bot_x_max = bot_left_bg_edge - sd_gap_to_gate
    sd_shapes.append(gdstk.rectangle(
        (s_bot_x_max - sd_width, y_center_bot - sd_height/2),
        (s_bot_x_max,            y_center_bot + sd_height/2), 
        layer=LAYER_SD, datatype=0))
        
    # Right D: 左边缘 = 最右BG右边缘 + 100nm
    bot_right_bg_edge = bot_max_x
    d_bot_x_min = bot_right_bg_edge + sd_gap_to_gate
    sd_shapes.append(gdstk.rectangle(
        (d_bot_x_min,            y_center_bot - sd_height/2),
        (d_bot_x_min + sd_width, y_center_bot + sd_height/2), 
        layer=LAYER_SD, datatype=0))
        
    # --- Top Row S/D ---
    # 辅助计算 SET 边缘
    set_half_width = bg_max_width + gap_pg_bg + pg_max_width/2
    
    # Left of Left SET
    # 边缘 = LeftCenter - HalfWidth - 100nm
    l_set_left_edge = left_set_center_pg_x - set_half_width
    d_top_l_x_max = l_set_left_edge - sd_gap_to_gate
    sd_shapes.append(gdstk.rectangle(
        (d_top_l_x_max - sd_width, y_center_top - sd_height/2),
        (d_top_l_x_max,            y_center_top + sd_height/2),
        layer=LAYER_SD, datatype=0))
        
    # Middle S (Between SETs)
    # 左边缘 = LeftSet_RightEdge + 100nm
    l_set_right_edge = left_set_center_pg_x + set_half_width
    s_mid_x_min = l_set_right_edge + sd_gap_to_gate
    
    # 右边缘 = RightSet_LeftEdge - 100nm
    r_set_left_edge = right_set_center_pg_x - set_half_width
    s_mid_x_max = r_set_left_edge - sd_gap_to_gate
    
    sd_shapes.append(gdstk.rectangle(
        (s_mid_x_min, y_center_top - sd_height/2),
        (s_mid_x_max, y_center_top + sd_height/2),
        layer=LAYER_SD, datatype=0))
        
    # Right of Right SET
    r_set_right_edge = right_set_center_pg_x + set_half_width
    d_top_r_x_min = r_set_right_edge + sd_gap_to_gate
    sd_shapes.append(gdstk.rectangle(
        (d_top_r_x_min,            y_center_top - sd_height/2),
        (d_top_r_x_min + sd_width, y_center_top + sd_height/2),
        layer=LAYER_SD, datatype=0))
        
    cell.add(*sd_shapes)

    return lib, cell

# ==========================================
# 7. 可视化
# ==========================================
def plot_gds(cell):
    fig, ax = plt.subplots(figsize=(14, 10))
    
    layer_config = {
        0: {'color': '#87CEEB', 'alpha': 0.6, 'label': 'S/D (L0)', 'zorder': 2},
        1: {'color': '#D3D3D3', 'alpha': 0.5, 'label': 'SG (L1)', 'hatch': '///', 'zorder': 1},
        2: {'color': '#FF1493', 'alpha': 0.9, 'label': 'BG (L2)', 'zorder': 3}, 
        3: {'color': '#8A2BE2', 'alpha': 0.9, 'label': 'PG (L3)', 'zorder': 4}, 
    }
    
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
    ax.autoscale_view()
    ax.set_ylim(-0.3, 0.7)
    
    ax.set_title('Final QD Layout V4: 100nm S/D Gaps, Continuous Leads')
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')
    plt.grid(True, which='both', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    lib, cell = create_final_layout_v4()
    lib.write_gds('qd_device_v4.gds')
    print("GDS saved as 'qd_device_v4.gds'")
    plot_gds(cell)
