# ===================================================================
#                      用户提供的代码基础 (不可删减)
# ===================================================================

import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.patches as patches
import math
import re

# ===================================================================
#                      [新增] 解决中文显示问题的代码
# ===================================================================
# 设置matplotlib支持中文的字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 'SimHei' 是 Windows 上常用的黑体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号'-'显示为方块的问题
# 如果您在 macOS 或 Linux 上运行，可以尝试其他字体，例如:
# plt.rcParams['font.sans-serif'] = ['PingFang SC'] # macOS
# plt.rcParams['font.sans-serif'] = ['WenQuanYi Micro Hei'] # Linux
# ===================================================================


def plot_gds(cell, title="量子点器件布局"):
    """
    用户提供的绘图函数，专门用于绘制量子点器件的核心层。
    """
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.set_title(title, fontsize=18)
    
    layer_config = {
        0: {'color': '#87CEEB', 'alpha': 0.6, 'label': 'S/D (L0)', 'zorder': 2},
        1: {'color': '#D3D3D3', 'alpha': 0.5, 'label': 'SG (L1)', 'hatch': '///', 'zorder': 1},
        2: {'color': '#FF1493', 'alpha': 0.9, 'label': 'BG (L2)', 'zorder': 3}, 
        3: {'color': '#8A2BE2', 'alpha': 0.9, 'label': 'PG (L3)', 'zorder': 4}, 
    }
    
    # 绘制多边形
    drawn_labels = set()
    for layer in [1, 0, 2, 3]: 
        polygons = cell.get_polygons(layer=layer, datatype=0)
        if not polygons: continue
        cfg = layer_config.get(layer)
        label = cfg['label']
        for gds_poly in polygons:
            pts = gds_poly.points
            if pts.ndim != 2 or len(pts) < 3: continue
            
            current_label = label if label not in drawn_labels else None
            poly = MplPolygon(pts, closed=True, 
                              facecolor=cfg['color'], edgecolor='black', linewidth=0.5,
                              alpha=cfg['alpha'], hatch=cfg.get('hatch', ''), 
                              label=current_label, zorder=cfg['zorder'])
            ax.add_patch(poly)
            drawn_labels.add(label)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='upper right')

    ax.set_aspect('equal')
    
    bbox = cell.bounding_box()
    if bbox is not None:
        min_x, min_y = bbox[0]; max_x, max_y = bbox[1]
        width = max_x - min_x; height = max_y - min_y
        margin_x = width * 0.1 if width > 0 else 1
        margin_y = height * 0.1 if height > 0 else 1
        ax.set_xlim(min_x - margin_x, max_x + margin_x)
        ax.set_ylim(min_y - margin_y, max_y + margin_y)
        print(f"Auto-adjusting View: X[{min_x:.2f}, {max_x:.2f}], Y[{min_y:.2f}, {max_y:.2f}]")
    else:
        ax.autoscale_view()
    
    ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)')
    plt.grid(True, which='both', linestyle='--', alpha=0.3); plt.tight_layout(); plt.show()

def create_6qd_layout_with_labels(
    pg_max_width=0.120, pg_vert_side_len=0.040, pg_chamfer_h=0.040, pg_bot_flat_w=0.042, pg_top_flat_w=0.040, bg_max_width=0.060, bg_vert_side_len=0.040, bg_top_flat_w=0.020, bg_bot_flat_w=0.042, gap_pg_bg=0.002, d1_gap=0.020, d2_gap=0.020, sd_gap_to_gate=0.01, gap_gate_outer_sg=0.030, sg_mid_thick=0.100,   sg_top_thick=0.300,   sg_bot_thick=0.300,   sg_extension=0.200, lead_width=0.042, lead_length_bot=0.4, lead_length_top=0.5, lead_overlap=0.025, sd_height=0.10, sd_width=0.35, sd_lead_width=0.060, sd_outer_length=0.300
):
    """
    用户提供的函数，用于创建带长引线和标签的6QD器件。
    [新增] 返回值中增加了 connection_points 字典。
    """
    LAYER_SD = 0; LAYER_SG = 1; LAYER_BG = 2; LAYER_PG = 3; LAYER_LABEL = 100
    gate_half_height = pg_vert_side_len/2 + pg_chamfer_h; row_pitch = gate_half_height + d1_gap + sg_mid_thick + d2_gap + gate_half_height
    lib = gdstk.Library(); cell = lib.new_cell('6QD_DEVICE_WITH_LEADS'); connection_points = {}
    def get_pg_points(center):
        cx, cy=center; y_v_top=pg_vert_side_len/2; y_v_bot=-pg_vert_side_len/2; y_top=y_v_top+pg_chamfer_h; y_bot=y_v_bot-pg_chamfer_h; x_max=pg_max_width/2; x_top_flat=pg_top_flat_w/2; x_bot_flat=pg_bot_flat_w/2
        return [(cx+x_top_flat,cy+y_top),(cx+x_max,cy+y_v_top),(cx+x_max,cy+y_v_bot),(cx+x_bot_flat,cy+y_bot),(cx-x_bot_flat,cy+y_bot),(cx-x_max,cy+y_v_bot),(cx-x_max,cy+y_v_top),(cx-x_top_flat,cy+y_top)], y_top, y_bot
    def get_bg_points(center):
        cx,cy=center; y_v_top=bg_vert_side_len/2; y_v_bot=-bg_vert_side_len/2; bg_chamfer_h_top=(bg_max_width-bg_top_flat_w)/2; bg_chamfer_h_bot=(bg_max_width-bg_bot_flat_w)/2; y_top=y_v_top+bg_chamfer_h_top; y_bot=y_v_bot-bg_chamfer_h_bot; x_max=bg_max_width/2; x_top_flat=bg_top_flat_w/2; x_bot_flat=bg_bot_flat_w/2
        return [(cx+x_top_flat,cy+y_top),(cx+x_max,cy+y_v_top),(cx+x_max,cy+y_v_bot),(cx+x_bot_flat,cy+y_bot),(cx-x_bot_flat,cy+y_bot),(cx-x_max,cy+y_v_bot),(cx-x_max,cy+y_v_top),(cx-x_top_flat,cy+y_top)], y_top, y_bot
    pg_local_bot_y=-(pg_vert_side_len/2+pg_chamfer_h); bg_local_bot_y=-(bg_vert_side_len/2+(bg_max_width-bg_bot_flat_w)/2); fixed_lead_bot_end_y=0+min(pg_local_bot_y,bg_local_bot_y)-lead_length_bot; y_center_top=row_pitch; pg_local_top_y=(pg_vert_side_len/2+pg_chamfer_h); bg_local_top_y=(bg_vert_side_len/2+(bg_max_width-bg_top_flat_w)/2); fixed_lead_top_end_y=y_center_top+max(pg_local_top_y,bg_local_top_y)+lead_length_top
    all_shapes=[]; y_center_bot=0; current_x=0; bottom_gate_names=[]; [bottom_gate_names.extend([f'QD_B{i}',f'QD_PG{i}']) for i in range(1,7)]; bottom_gate_names.append('QD_B7'); bot_tips_y=[]
    for i,name in enumerate(bottom_gate_names):
        g_type='PG' if 'PG' in name else 'BG'; w=pg_max_width if g_type=='PG' else bg_max_width; layer=LAYER_PG if g_type=='PG' else LAYER_BG;
        if i>0: prev_w=pg_max_width if 'PG' in bottom_gate_names[i-1] else bg_max_width; current_x+=(prev_w/2)+gap_pg_bg+(w/2)
        center=(current_x,y_center_bot); pts,tip_y,base_y=get_pg_points(center) if g_type=='PG' else get_bg_points(center); bot_tips_y.append(tip_y); poly=gdstk.Polygon(pts,layer=layer); lead_y_start=base_y+lead_overlap; lead_y_end=fixed_lead_bot_end_y; lead=gdstk.rectangle((current_x-lead_width/2,lead_y_start),(current_x+lead_width/2,lead_y_end),layer=layer)
        all_shapes.extend(gdstk.boolean([poly,lead],[],'or',layer=layer)); label=gdstk.Label(name,(current_x,lead_y_end),layer=LAYER_LABEL); cell.add(label); connection_points[name]=label.origin
    bot_min_x=0-bg_max_width/2; bot_max_x=current_x+bg_max_width/2; center_x_global=(bot_min_x+bot_max_x)/2; inner_gap_required=2*sd_gap_to_gate+sd_width; offset_from_center=inner_gap_required/2+bg_max_width/2+gap_pg_bg+pg_max_width/2; left_set_center_pg_x=center_x_global-offset_from_center; right_set_center_pg_x=center_x_global+offset_from_center
    def create_set_at(center_pg_x,side):
        shapes=[]; labels=[]; offsets=[-(pg_max_width/2+gap_pg_bg+bg_max_width/2),0,(pg_max_width/2+gap_pg_bg+bg_max_width/2)]; set_gate_names=['SET1_B1','SET1_G','SET1_B2'] if side=='L' else ['SET2_B2','SET2_G','SET2_B1']
        for offset,name in zip(offsets,set_gate_names):
            cx=center_pg_x+offset; g_type='PG' if '_G' in name else 'BG'; layer=LAYER_PG if g_type=='PG' else LAYER_BG; pts,_,_=get_pg_points((cx,y_center_top)) if g_type=='PG' else get_bg_points((cx,y_center_top)); poly=gdstk.Polygon(pts,layer=layer); Wire_top_y=np.max([pt[1] for pt in pts]); lead_y_start=Wire_top_y-lead_overlap; lead_y_end=fixed_lead_top_end_y; lead=gdstk.rectangle((cx-lead_width/2,lead_y_start),(cx+lead_width/2,lead_y_end),layer=layer)
            shapes.extend(gdstk.boolean([poly,lead],[],'or',layer=layer)); label=gdstk.Label(name,(cx,lead_y_end),layer=LAYER_LABEL); cell.add(label); connection_points[name]=label.origin
        return shapes
    all_shapes.extend(create_set_at(left_set_center_pg_x,side='L')); all_shapes.extend(create_set_at(right_set_center_pg_x,side='R'))
    sg_x_start=bot_min_x-sg_extension; sg_x_end=bot_max_x+sg_extension; y_sg_mid_bot=max(bot_tips_y)+d1_gap; y_sg_mid_top=y_sg_mid_bot+sg_mid_thick; sg_mid=gdstk.rectangle((sg_x_start,y_sg_mid_bot),(sg_x_end,y_sg_mid_top),layer=LAYER_SG); y_sg_top_bot=y_center_top+gate_half_height+gap_gate_outer_sg; y_sg_top_top=y_sg_top_bot+sg_top_thick; sg_top=gdstk.rectangle((sg_x_start,y_sg_top_bot),(sg_x_end,y_sg_top_top),layer=LAYER_SG); y_sg_bot_top=y_center_bot-gate_half_height-gap_gate_outer_sg; y_sg_bot_bot=y_sg_bot_top-sg_bot_thick; sg_bot=gdstk.rectangle((sg_x_start,y_sg_bot_bot),(sg_x_end,y_sg_bot_top),layer=LAYER_SG); all_shapes.extend([sg_top,sg_mid,sg_bot])
    label_sg1=gdstk.Label("SG1",(sg_x_start,(y_sg_top_bot+y_sg_top_top)/2),layer=LAYER_LABEL); cell.add(label_sg1); connection_points["SG1"]=label_sg1.origin
    label_sg2=gdstk.Label("SG2",(sg_x_start,(y_sg_mid_bot+y_sg_mid_top)/2),layer=LAYER_LABEL); cell.add(label_sg2); connection_points["SG2"]=label_sg2.origin
    label_sg3=gdstk.Label("SG3",(sg_x_end,(y_sg_bot_bot+y_sg_bot_top)/2),layer=LAYER_LABEL); cell.add(label_sg3); connection_points["SG3"]=label_sg3.origin
    sd_shapes=[]; sd_corner_r=0.020
    def create_rounded_sd(x_min,y_min,x_max,y_max,corners_to_round):
        poly=gdstk.Polygon([(x_min,y_min),(x_max,y_min),(x_max,y_max),(x_min,y_max)],layer=LAYER_SD); poly.fillet([sd_corner_r if c else 0 for c in corners_to_round],tolerance=0.001); return poly
    bot_left_bg_edge=bot_min_x; bot_left_sd_inner_x=bot_left_bg_edge-sd_gap_to_gate; bot_right_bg_edge=bot_max_x; bot_right_sd_inner_x=bot_right_bg_edge+sd_gap_to_gate; global_outer_left_x=bot_left_sd_inner_x-sd_outer_length; global_outer_right_x=bot_right_sd_inner_x+sd_outer_length
    sd_shapes.append(create_rounded_sd(global_outer_left_x,y_center_bot-sd_height/2,bot_left_sd_inner_x,y_center_bot+sd_height/2,[False,True,True,False])); label_qdd=gdstk.Label("QD_D",(global_outer_left_x,y_center_bot),layer=LAYER_LABEL); cell.add(label_qdd); connection_points["QD_D"]=label_qdd.origin
    sd_shapes.append(create_rounded_sd(bot_right_sd_inner_x,y_center_bot-sd_height/2,global_outer_right_x,y_center_bot+sd_height/2,[True,False,False,True])); label_qds=gdstk.Label("QD_S",(global_outer_right_x,y_center_bot),layer=LAYER_LABEL); cell.add(label_qds); connection_points["QD_S"]=label_qds.origin
    set_half_width=bg_max_width+gap_pg_bg+pg_max_width/2; l_set_left_edge=left_set_center_pg_x-set_half_width; top_left_sd_inner_x=l_set_left_edge-sd_gap_to_gate; r_set_right_edge=right_set_center_pg_x+set_half_width; top_right_sd_inner_x=r_set_right_edge+sd_gap_to_gate
    sd_shapes.append(create_rounded_sd(global_outer_left_x,y_center_top-sd_height/2,top_left_sd_inner_x,y_center_top+sd_height/2,[False,True,True,False])); label_s1d=gdstk.Label("SET1_D",(global_outer_left_x,y_center_top),layer=LAYER_LABEL); cell.add(label_s1d); connection_points["SET1_D"]=label_s1d.origin
    sd_shapes.append(create_rounded_sd(top_right_sd_inner_x,y_center_top-sd_height/2,global_outer_right_x,y_center_top+sd_height/2,[True,False,False,True])); label_s2d=gdstk.Label("SET2_D",(global_outer_right_x,y_center_top),layer=LAYER_LABEL); cell.add(label_s2d); connection_points["SET2_D"]=label_s2d.origin
    l_set_right_edge=left_set_center_pg_x+set_half_width; s_mid_x_min=l_set_right_edge+sd_gap_to_gate; r_set_left_edge=right_set_center_pg_x-set_half_width; s_mid_x_max=r_set_left_edge-sd_gap_to_gate
    mid_sd_poly=create_rounded_sd(s_mid_x_min,y_center_top-sd_height/2,s_mid_x_max,y_center_top+sd_height/2,[True,True,True,True]); mid_sd_center_x=(s_mid_x_min+s_mid_x_max)/2; mid_lead_y_end=y_center_top+lead_length_top; mid_sd_lead=gdstk.rectangle((mid_sd_center_x-sd_lead_width/2,y_center_top),(mid_sd_center_x+sd_lead_width/2,mid_lead_y_end),layer=LAYER_SD)
    sd_shapes.extend(gdstk.boolean([mid_sd_poly,mid_sd_lead],[],'or',layer=LAYER_SD)); label_sets=gdstk.Label("SET_S",(mid_sd_center_x,mid_lead_y_end),layer=LAYER_LABEL); cell.add(label_sets); connection_points["SET_S"]=label_sets.origin
    all_shapes.extend(sd_shapes); cell.add(*all_shapes)
    return lib, cell, connection_points

def create_rect_wire_layout(
    N=6, layout_width=1400, layout_height=1400, pad_width=100, pad_height=100, pad_spacing=30, edge_margin=50, active_size=180, trace_width=10, trace_spacing=10, taper_length=40, active_entry_len=30
):
    """
    用户提供的函数，用于创建带引线的焊盘框架。
    [新增] 返回值增加了 connection_points 和 active_center。
    """
    lib=gdstk.Library(); cell=lib.new_cell('PAD_FRAME_WITH_WIRES'); connection_points={}; all_pads=[]
    def natural_sort_key(s): return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)',s)]
    pg_labels=sorted([f"QD_PG{i+1}" for i in range(N)],key=natural_sort_key); bg_labels=sorted([f"QD_B{i+1}" for i in range(N+1)],key=natural_sort_key); gate_sequence=[]; [gate_sequence.extend([bg_labels[i],pg_labels[i]]) for i in range(N)]; gate_sequence.append(bg_labels[N]); top_pads=[]; num_sets=math.ceil(N/3)
    if num_sets>=1: top_pads.extend(["SET1_B1","SET1_G","SET1_B2"]); top_pads.append("SET_S");
    if num_sets>=2: top_pads.extend(["SET2_B2","SET2_G","SET2_B1"])
    left_fixed=["SG1"];
    if num_sets>=1: left_fixed.append("SET1_D"); left_fixed.extend(["SG2","QD_D"]); right_fixed=[]
    if num_sets>=2: right_fixed.append("SET2_D"); right_fixed.append("QD_S"); right_fixed.extend(["GND","SG3","GND"])
    len_seq=len(gate_sequence); len_L_fix=len(left_fixed); len_R_fix=len(right_fixed); len_Top=len(top_pads); available_w=layout_width-2*edge_margin; best_config=(0,len_seq,0); min_diff=float('inf'); start_side=max(len_L_fix,len_R_fix)
    for target_side in range(start_side,start_side+len_seq):
        add_L=target_side-len_L_fix; add_R=target_side-len_R_fix;
        if add_L+add_R > len_seq: break
        n_bot=len_seq-(add_L+add_R)
        if n_bot > 0 and (n_bot*pad_width+(n_bot-1)*pad_spacing>available_w): continue
        diff=abs(n_bot-len_Top);
        if diff < min_diff: min_diff=diff; best_config=(add_L,n_bot,add_R)
    n_L,n_B,n_R=best_config; pads_left_extra=gate_sequence[:n_L]; pads_bottom=gate_sequence[n_L:n_L+n_B]; pads_right_extra=gate_sequence[n_L+n_B:]; left_pads=left_fixed+pads_left_extra; right_pads=right_fixed+list(reversed(pads_right_extra))
    def assign_coords(pads,edge):
        count=len(pads)
        if count==0: return
        if edge in ['top','bottom']:
            total_span=count*pad_width+(count-1)*pad_spacing; start_x=(layout_width-total_span)/2; y=layout_height-edge_margin-pad_height if edge=='top' else edge_margin
            for i,label in enumerate(pads): x=start_x+i*(pad_width+pad_spacing); cx,cy=x+pad_width/2,y+pad_height/2; all_pads.append({'label':label,'rect':(x,y,pad_width,pad_height),'center':(cx,cy),'edge':edge})
        else:
            total_span=count*pad_width+(count-1)*pad_spacing; top_y=(layout_height+total_span)/2; x=edge_margin if edge=='left' else layout_width-edge_margin-pad_height
            for i,label in enumerate(pads): y=top_y-(i+1)*pad_width-i*pad_spacing; cx,cy=x+pad_height/2,y+pad_width/2; all_pads.append({'label':label,'rect':(x,y,pad_height,pad_width),'center':(cx,cy),'edge':edge})
    assign_coords(top_pads,'top'); assign_coords(pads_bottom,'bottom'); assign_coords(left_pads,'left'); assign_coords(right_pads,'right')
    center_x=layout_width/2; center_y=layout_height/2; active_center=(center_x,center_y); aa_x1=center_x-active_size/2; aa_x2=center_x+active_size/2; aa_y1=center_y-active_size/2; aa_y2=center_y+active_size/2; pads_by_edge={'top':[],'bottom':[],'left':[],'right':[]}; [pads_by_edge[p['edge']].append(p) for p in all_pads]
    pads_by_edge['top'].sort(key=lambda p: p['center'][0]); pads_by_edge['bottom'].sort(key=lambda p: p['center'][0]); pads_by_edge['left'].sort(key=lambda p: p['center'][1],reverse=True); pads_by_edge['right'].sort(key=lambda p: p['center'][1],reverse=True); pitch=trace_width+trace_spacing
    def assign_targets(pads,edge):
        n=len(pads)
        if n==0: return
        indices=np.arange(n)-(n-1)/2; offsets=indices*pitch
        if edge=='top': [p.update({'target':(center_x+offsets[i],aa_y2),'target_edge':'top'}) for i,p in enumerate(pads)]
        elif edge=='bottom': [p.update({'target':(center_x+offsets[i],aa_y1),'target_edge':'bottom'}) for i,p in enumerate(pads)]
        elif edge=='left': [p.update({'target':(aa_x1,center_y+offsets[-(i+1)]),'target_edge':'left'}) for i,p in enumerate(pads)]
        elif edge=='right': [p.update({'target':(aa_x2,center_y+offsets[-(i+1)]),'target_edge':'right'}) for i,p in enumerate(pads)]
    [assign_targets(pads_by_edge[edge],edge) for edge in ['top','bottom','left','right']]
    LAYER_SUBSTRATE=0; LAYER_PAD=1; LAYER_TRACE=1; LAYER_TEXT=2; LAYER_ACTIVE=10; LAYER_SIO2=11; sio2_margin=10
    cell.add(gdstk.rectangle((0,0),(layout_width,layout_height),layer=LAYER_SUBSTRATE)); cell.add(gdstk.rectangle((aa_x1,aa_y1),(aa_x2,aa_y2),layer=LAYER_ACTIVE))
    def create_taper(p_rect,edge,t_len,t_width):
        x,y,w,h=p_rect; cx,cy=x+w/2,y+h/2
        if edge=='top': pts=[(x,y),(x+w,y),(cx+t_width/2,y-t_len),(cx-t_width/2,y-t_len)]; start=(cx,y-t_len)
        elif edge=='bottom': pts=[(x,y+h),(x+w,y+h),(cx+t_width/2,y+h+t_len),(cx-t_width/2,y+h+t_len)]; start=(cx,y+h+t_len)
        elif edge=='left': pts=[(x+w,y),(x+w,y+h),(x+w+t_len,cy+t_width/2),(x+w+t_len,cy-t_width/2)]; start=(x+w+t_len,cy)
        else: pts=[(x,y),(x,y+h),(x-t_len,cy+t_width/2),(x-t_len,cy-t_width/2)]; start=(x-t_len,cy)
        return gdstk.Polygon(pts,layer=LAYER_TRACE),start
    def route_safe_z(start_pt,target_pt,edge,width):
        x1,y1=start_pt; x2,y2=target_pt; safe_margin=20; pts=[(x1,y1)]
        if edge in ['top','bottom']:
            sign_y=-1 if edge=='top' else 1; y_safe=y1+sign_y*safe_margin; pts.append((x1,y_safe)); y_entry=y2+sign_y*active_entry_len; avail_h=abs(y_entry-y_safe); needed_w=abs(x2-x1); sign_x=1 if x2>x1 else -1
            if needed_w>avail_h: h_leg=avail_h/2; y_mid=y_safe+sign_y*h_leg; x_mid_1=x1+sign_x*h_leg; x_mid_2=x2-sign_x*h_leg; pts.extend([(x_mid_1,y_mid),(x_mid_2,y_mid)])
            else: dy_diag=needed_w; y_turn_end=y_safe+sign_y*dy_diag; pts.append((x2,y_turn_end))
            pts.append((x2,y2))
        else:
            sign_x=1 if edge=='left' else -1; x_safe=x1+sign_x*safe_margin; pts.append((x_safe,y1)); x_entry=x2+sign_x*active_entry_len; avail_w=abs(x_entry-x_safe); needed_h=abs(y2-y1); sign_y=1 if y2>y1 else -1
            if needed_h>avail_w: w_leg=avail_w/2; x_mid=x_safe+sign_x*w_leg; y_mid_1=y1+sign_y*w_leg; y_mid_2=y2-sign_y*w_leg; pts.extend([(x_mid,y_mid_1),(x_mid,y_mid_2)])
            else: dx_diag=needed_h; x_turn_end=x_safe+sign_x*dx_diag; pts.append((x_turn_end,y2))
            pts.append((x2,y2))
        return gdstk.FlexPath(pts,width,ends="flush",joins="round",layer=LAYER_TRACE),pts
    for p in all_pads:
        x,y,w,h=p['rect']; label=p['label']; edge=p['edge']
        cell.add(gdstk.rectangle((x-sio2_margin,y-sio2_margin),(x+w+sio2_margin,y+h+sio2_margin),layer=LAYER_SIO2))
        cell.add(gdstk.rectangle((x,y),(x+w,y+h),layer=LAYER_PAD)); cx,cy=x+w/2,y+h/2; rot_rad=0 if w>h else math.pi/2
        cell.add(gdstk.Label(label,(cx,cy),layer=LAYER_TEXT,magnification=10,rotation=rot_rad))
        poly_taper,trace_start=create_taper(p['rect'],edge,taper_length,trace_width); cell.add(poly_taper)
        path_obj,path_pts=route_safe_z(trace_start,p['target'],edge,trace_width); cell.add(path_obj)
        p['trace_points']=path_pts; connection_points[p['label']] = p['target']
    
    return lib, cell, connection_points, active_center

# ===================================================================
#                      新增的分步执行与可视化逻辑
# ===================================================================

def visualize_and_save_pad_frame(lib, cell_to_show, title, gds_filename):
    """
    一个专门用于显示【焊盘框架】的函数，以匹配用户期望的样式。
    """
    print(f"--- {title} ---"); print(f"正在可视化: {title}...")
    fig, ax = plt.subplots(figsize=(18, 18))
    
    # 使用 flatten 以确保所有几何都在顶层，便于提取
    plot_cell = cell_to_show.flatten()

    pads_and_traces = plot_cell.get_polygons(layer=1, datatype=0)
    active_area_poly = plot_cell.get_polygons(layer=10, datatype=0)
    sio2_poly = plot_cell.get_polygons(layer=11, datatype=0)

    for gds_poly in sio2_poly: ax.add_patch(MplPolygon(gds_poly.points, closed=True, fill=False, edgecolor='gray', linestyle='--'))
    
    for i, gds_poly in enumerate(pads_and_traces):
        label = "Pads & Traces (L1)" if i == 0 else None
        ax.add_patch(MplPolygon(gds_poly.points, closed=True, facecolor='gold', edgecolor='black', linewidth=0.2, label=label))
    
    if active_area_poly:
        ax.add_patch(MplPolygon(active_area_poly[0].points, closed=True, facecolor='pink', alpha=0.5, label="Active Area (L10)"))

    ax.legend(loc='upper right'); ax.set_aspect('equal')
    bbox = cell_to_show.bounding_box()
    if bbox:
        (min_x,min_y),(max_x,max_y)=bbox; margin=(max_x-min_x)*0.05
        ax.set_xlim(min_x-margin,max_x+margin); ax.set_ylim(min_y-margin,max_y+margin)
    ax.set_title(title, fontsize=18); ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); plt.grid(True, linestyle='--', alpha=0.3); plt.tight_layout()
    plt.show()

    lib.write_gds(gds_filename); print(f"GDS 文件已保存到: '{gds_filename}'\n")

def visualize_and_save_assembly(lib, cell_to_show, title, gds_filename):
    """
    一个通用的函数，用于显示组装和布线后的单元。
    """
    print(f"--- {title} ---"); print(f"正在可视化: {title}...")

    plot_cell = cell_to_show.flatten()
    fig, ax = plt.subplots(figsize=(18, 18))
    
    layer_config = {
        0: {'color':'#87CEEB','alpha':0.9,'label':'S/D (L0)','zorder':4},
        # Layer 1 from pad frame and device are different, handle separately
        2: {'color':'#FF1493','alpha':0.9,'label':'BG (L2)','zorder':5},
        3: {'color':'#8A2BE2','alpha':0.9,'label':'PG (L3)','zorder':6},
        5: {'color':'#0000FF','alpha':1.0,'label':'Final Routing (L5)','zorder':3},
        10:{'color':'pink','alpha':0.5,'label':'Active Area (L10)','zorder':2},
        11:{'color':'gray', 'alpha':1.0, 'label':'SiO2 (L11)', 'zorder':0}
    }
    
    drawn_labels = set()
    def draw_polygons(polygons, config):
        label = config['label']
        for gds_poly in polygons:
            vertices = gds_poly.points
            if len(vertices) > 0:
                current_label = label if label not in drawn_labels else None
                is_fill = True if config.get('layer') != 11 else False
                edge_color = 'black' if config.get('layer') != 11 else 'gray'
                line_style = 'solid' if config.get('layer') != 11 else '--'

                mpl_poly = MplPolygon(vertices, closed=True, fill=is_fill, facecolor=config['color'], edgecolor=edge_color, 
                                      linestyle=line_style, linewidth=0.02, alpha=config['alpha'], 
                                      label=current_label, zorder=config.get('zorder', 1))
                ax.add_patch(mpl_poly)
                drawn_labels.add(label)

    # Draw specific layers
    for layer, cfg in layer_config.items():
        polys = plot_cell.get_polygons(layer=layer, datatype=0)
        if polys: draw_polygons(polys, cfg)

    # Handle Pad/Trace (L1) from pad_frame style
    pad_polys = plot_cell.get_polygons(layer=1, datatype=0)
    if pad_polys: draw_polygons(pad_polys, {'color':'gold','alpha':0.08,'label':'Pads & Traces (L1)','zorder':1})

    # Handle SG (L1) from device style
    sg_polys = plot_cell.get_polygons(layer=1, datatype=0)
    if sg_polys: draw_polygons(sg_polys, {'color':'#D3D3D3', 'alpha':0.05, 'label':'SG (L1)', 'hatch':'///', 'zorder':7})


    handles, labels = ax.get_legend_handles_labels(); by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    bbox=plot_cell.bounding_box(); 
    if bbox:
        (min_x,min_y),(max_x,max_y)=bbox; margin=(max_x-min_x)*0.05
        ax.set_xlim(min_x-margin,max_x+margin); ax.set_ylim(min_y-margin,max_y+margin)
    
    ax.set_aspect('equal'); ax.set_title(title, fontsize=18)
    ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); plt.grid(True, linestyle='--', alpha=0.3); plt.tight_layout()
    plt.show()

    lib.write_gds(gds_filename); print(f"GDS 文件已保存到: '{gds_filename}'\n")

# ===================================================================
#                           主程序入口 (全新分步流程)
# ===================================================================
if __name__ == "__main__":

    # --- 步骤 1: 生成并展示核心器件 ---
    lib_device, cell_device, device_points = create_6qd_layout_with_labels()
    # 使用用户提供的 plot_gds 函数进行可视化
    plot_gds(cell_device, title="步骤 1: 仅带长引线的量子点器件")
    lib_device.write_gds("step1_device_with_leads.gds")
    print("GDS 文件已保存到: 'step1_device_with_leads.gds'\n")


    # --- 步骤 2: 生成并展示焊盘框架 ---
    lib_pads, cell_pads, pad_points, active_center = create_rect_wire_layout()
    visualize_and_save_pad_frame(lib_pads, cell_pads, 
                       "步骤 2: 仅带引线的焊盘框架", 
                       "step2_pad_frame.gds")

    # --- 步骤 3: 组装 - 将器件放置在中心 (未布线) ---
    # 居中器件，以便精确放置
    bbox = cell_device.bounding_box()
    if bbox:
        (min_x, min_y), (max_x, max_y) = bbox
        cx, cy = (min_x+max_x)/2, (min_y+max_y)/2
        # 平移几何对象
        for poly in cell_device.polygons: poly.translate(-cx, -cy)
        # 更新连接点坐标
        for label, pos in device_points.items(): device_points[label] = (pos[0]-cx, pos[1]-cy)
    
    assembly_lib = gdstk.Library(); assembly_cell = assembly_lib.new_cell('ASSEMBLY_CELL')
    assembly_cell.add(gdstk.Reference(cell_pads))
    device_ref = gdstk.Reference(cell_device, active_center)
    assembly_cell.add(device_ref)
    
    visualize_and_save_assembly(assembly_lib, assembly_cell, 
                       "步骤 3: 器件已放置 (未布线)", 
                       "step3_placed_unrouted.gds")

    # --- 步骤 4: 连接单条线进行演示 ---
    # 这里的逻辑是连接两个长引线的末端，所以最终连线会非常短
    label_to_connect = 'QD_PG4'
    if label_to_connect in device_points and label_to_connect in pad_points:
        pad_coord = pad_points[label_to_connect]
        local_dev_coord = device_points[label_to_connect]
        global_dev_coord = (local_dev_coord[0] + active_center[0], local_dev_coord[1] + active_center[1])
        
        path = gdstk.FlexPath(points=[pad_coord, global_dev_coord], width=0.03, layer=5, ends="flush")
        assembly_cell.add(path)

    visualize_and_save_assembly(assembly_lib, assembly_cell, 
                       f"步骤 4: 连接单条线 ({label_to_connect})", 
                       "step4_single_route.gds")

    # --- 步骤 5: 完成所有剩余的连接 ---
    print("--- 步骤 5: 正在连接所有剩余的线... ---")
    for label, local_dev_coord in device_points.items():
        if label != label_to_connect and label in pad_points:
            pad_coord = pad_points[label]
            global_dev_coord = (local_dev_coord[0] + active_center[0], local_dev_coord[1] + active_center[1])
            path = gdstk.FlexPath(points=[pad_coord, global_dev_coord], width=0.03, layer=5, ends="flush")
            assembly_cell.add(path)
            
    visualize_and_save_assembly(assembly_lib, assembly_cell, 
                       "步骤 5: 最终完整布线", 
                       "step5_fully_routed.gds")

    print("全部分步流程执行完毕！")

