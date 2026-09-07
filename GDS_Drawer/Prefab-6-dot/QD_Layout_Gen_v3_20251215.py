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
#                      解决中文显示问题的代码
# ===================================================================
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
# ===================================================================


def plot_gds(cell, title="量子点器件布局"):
    # ... (此函数与上一版完全相同，保持不变)
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

def create_6qd_layout_with_labels(
    # ... (此函数与上一版完全相同，保持不变)
    pg_max_width=0.120, pg_vert_side_len=0.040, pg_chamfer_h=0.040, pg_bot_flat_w=0.042, pg_top_flat_w=0.040, bg_max_width=0.060, bg_vert_side_len=0.040, bg_top_flat_w=0.020, bg_bot_flat_w=0.042, gap_pg_bg=0.002, d1_gap=0.020, d2_gap=0.020, sd_gap_to_gate=0.01, gap_gate_outer_sg=0.030, sg_mid_thick=0.100,   sg_top_thick=0.300,   sg_bot_thick=0.300,   sg_extension=0.200, lead_width=0.042, lead_length_bot=0.4, lead_length_top=0.5, lead_overlap=0.025, sd_height=0.10, sd_width=0.35, sd_lead_width=0.060, sd_outer_length=0.300
):
    # (Code from previous version)
    LAYER_SD=0;LAYER_SG=1;LAYER_BG=2;LAYER_PG=3;LAYER_LABEL=100;gate_half_height=pg_vert_side_len/2+pg_chamfer_h;row_pitch=gate_half_height+d1_gap+sg_mid_thick+d2_gap+gate_half_height;lib=gdstk.Library();cell=lib.new_cell('6QD_DEVICE_WITH_LEADS');connection_points={};
    def get_pg_points(c): cx,cy=c;yvt=pg_vert_side_len/2;yvb=-pg_vert_side_len/2;yt=yvt+pg_chamfer_h;yb=yvb-pg_chamfer_h;xm=pg_max_width/2;xtf=pg_top_flat_w/2;xbf=pg_bot_flat_w/2;return [(cx+xtf,cy+yt),(cx+xm,cy+yvt),(cx+xm,cy+yvb),(cx+xbf,cy+yb),(cx-xbf,cy+yb),(cx-xm,cy+yvb),(cx-xm,cy+yvt),(cx-xtf,cy+yt)],yt,yb
    def get_bg_points(c): cx,cy=c;yvt=bg_vert_side_len/2;yvb=-bg_vert_side_len/2;bcht=(bg_max_width-bg_top_flat_w)/2;bchb=(bg_max_width-bg_bot_flat_w)/2;yt=yvt+bcht;yb=yvb-bchb;xm=bg_max_width/2;xtf=bg_top_flat_w/2;xbf=bg_bot_flat_w/2;return [(cx+xtf,cy+yt),(cx+xm,cy+yvt),(cx+xm,cy+yvb),(cx+xbf,cy+yb),(cx-xbf,cy+yb),(cx-xm,cy+yvb),(cx-xm,cy+yvt),(cx-xtf,cy+yt)],yt,yb
    pglb=-(pg_vert_side_len/2+pg_chamfer_h);bglb=-(bg_vert_side_len/2+(bg_max_width-bg_bot_flat_w)/2);fbe=min(pglb,bglb)-lead_length_bot;yct=row_pitch;pglt=(pg_vert_side_len/2+pg_chamfer_h);bglt=(bg_vert_side_len/2+(bg_max_width-bg_top_flat_w)/2);fte=yct+max(pglt,bglt)+lead_length_top
    all_shapes=[];ycb=0;curr_x=0;bgn=[];[bgn.extend([f'QD_B{i}',f'QD_PG{i}'])for i in range(1,7)];bgn.append('QD_B7');bot_tips_y=[]
    for i,n in enumerate(bgn):
        gt='PG' if 'PG' in n else 'BG';w=pg_max_width if gt=='PG' else bg_max_width;lyr=LAYER_PG if gt=='PG' else LAYER_BG;
        if i>0:pw=pg_max_width if 'PG' in bgn[i-1] else bg_max_width;curr_x+=(pw/2)+gap_pg_bg+(w/2)
        cen=(curr_x,ycb);pts,tip_y,by=get_pg_points(cen)if gt=='PG'else get_bg_points(cen);bot_tips_y.append(tip_y);poly=gdstk.Polygon(pts,layer=lyr);lys=by+lead_overlap;lye=fbe;lead=gdstk.rectangle((curr_x-lead_width/2,lys),(curr_x+lead_width/2,lye),layer=lyr)
        all_shapes.extend(gdstk.boolean([poly,lead],[],'or',layer=lyr));lbl=gdstk.Label(n,(curr_x,lye),layer=LAYER_LABEL);cell.add(lbl);connection_points[n]=lbl.origin
    bmx=-bg_max_width/2;bax=curr_x+bg_max_width/2;cxg=(bmx+bax)/2;igr=2*sd_gap_to_gate+sd_width;off_cen=igr/2+bg_max_width/2+gap_pg_bg+pg_max_width/2;lscpgx=cxg-off_cen;rscpgx=cxg+off_cen
    def create_set_at(cpgx,side):
        shps=[];offs=[-(pg_max_width/2+gap_pg_bg+bg_max_width/2),0,(pg_max_width/2+gap_pg_bg+bg_max_width/2)];sgn=['SET1_B1','SET1_G','SET1_B2']if side=='L'else['SET2_B1','SET2_G','SET2_B2']
        for off,n in zip(offs,sgn):
            cx=cpgx+off;gt='PG' if '_G' in n else 'BG';lyr=LAYER_PG if gt=='PG' else LAYER_BG;pts,_,_=get_pg_points((cx,yct))if gt=='PG'else get_bg_points((cx,yct));poly=gdstk.Polygon(pts,layer=lyr);Wty=np.max([pt[1]for pt in pts]);lys=Wty-lead_overlap;lye=fte;lead=gdstk.rectangle((cx-lead_width/2,lys),(cx+lead_width/2,lye),layer=lyr)
            shps.extend(gdstk.boolean([poly,lead],[],'or',layer=lyr));lbl=gdstk.Label(n,(cx,lye),layer=LAYER_LABEL);cell.add(lbl);connection_points[n]=lbl.origin
        return shps
    all_shapes.extend(create_set_at(lscpgx,side='L'));all_shapes.extend(create_set_at(rscpgx,side='R'))
    sgxs=bmx-sg_extension;sgxe=bax+sg_extension;ysgmb=max(bot_tips_y)+d1_gap;ysgmt=ysgmb+sg_mid_thick;sg_mid=gdstk.rectangle((sgxs,ysgmb),(sgxe,ysgmt),layer=LAYER_SG);ysgtb=yct+gate_half_height+gap_gate_outer_sg;ysgtt=ysgtb+sg_top_thick;sg_top=gdstk.rectangle((sgxs,ysgtb),(sgxe,ysgtt),layer=LAYER_SG);ysgbt=ycb-gate_half_height-gap_gate_outer_sg;ysgbb=ysgbt-sg_bot_thick;sg_bot=gdstk.rectangle((sgxs,ysgbb),(sgxe,ysgbt),layer=LAYER_SG);all_shapes.extend([sg_top,sg_mid,sg_bot])
    lsg1=gdstk.Label("SG1",(sgxs,(ysgtb+ysgtt)/2),layer=LAYER_LABEL);cell.add(lsg1);connection_points["SG1"]=lsg1.origin;lsg2=gdstk.Label("SG2",(sgxs,(ysgmb+ysgmt)/2),layer=LAYER_LABEL);cell.add(lsg2);connection_points["SG2"]=lsg2.origin;lsg3=gdstk.Label("SG3",(sgxs,(ysgbb+ysgbt)/2),layer=LAYER_LABEL);cell.add(lsg3);connection_points["SG3"]=lsg3.origin
    sd_shps=[];sdcr=0.020
    def cr(xmin,ymin,xmax,ymax,crnrs):poly=gdstk.Polygon([(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)],layer=LAYER_SD);poly.fillet([sdcr if c else 0 for c in crnrs],tolerance=0.001);return poly
    blbe=bmx;blsix=blbe-sd_gap_to_gate;brbe=bax;brsix=brbe+sd_gap_to_gate;golx=blsix-sd_outer_length;gorx=brsix+sd_outer_length;sd_shps.append(cr(golx,ycb-sd_height/2,blsix,ycb+sd_height/2,[False,True,True,False]));lqdd=gdstk.Label("QD_D",(golx,ycb),layer=LAYER_LABEL);cell.add(lqdd);connection_points["QD_D"]=lqdd.origin;sd_shps.append(cr(brsix,ycb-sd_height/2,gorx,ycb+sd_height/2,[True,False,False,True]));lqds=gdstk.Label("QD_S",(gorx,ycb),layer=LAYER_LABEL);cell.add(lqds);connection_points["QD_S"]=lqds.origin;shw=bg_max_width+gap_pg_bg+pg_max_width/2;lsle=lscpgx-shw;tlsix=lsle-sd_gap_to_gate;rsre=rscpgx+shw;trsix=rsre-sd_gap_to_gate;sd_shps.append(cr(golx,yct-sd_height/2,tlsix,yct+sd_height/2,[False,True,True,False]));ls1d=gdstk.Label("SET1_D",(golx,yct),layer=LAYER_LABEL);cell.add(ls1d);connection_points["SET1_D"]=ls1d.origin;sd_shps.append(cr(trsix,yct-sd_height/2,gorx,yct+sd_height/2,[True,False,False,True]));ls2d=gdstk.Label("SET2_D",(gorx,yct),layer=LAYER_LABEL);cell.add(ls2d);connection_points["SET2_D"]=ls2d.origin;lsre=lscpgx+shw;smxm=lsre+sd_gap_to_gate;rsle=rscpgx-shw;smxa=rsle-sd_gap_to_gate;msdp=cr(smxm,yct-sd_height/2,smxa,yct+sd_height/2,[True,True,True,True]);msdcx=(smxm+smxa)/2;mlye=yct+lead_length_top;msdl=gdstk.rectangle((msdcx-sd_lead_width/2,yct),(msdcx+sd_lead_width/2,mlye),layer=LAYER_SD);sd_shps.extend(gdstk.boolean([msdp,msdl],[],'or',layer=LAYER_SD));lsets=gdstk.Label("SET_S",(msdcx,mlye),layer=LAYER_LABEL);cell.add(lsets);connection_points["SET_S"]=lsets.origin;all_shapes.extend(sd_shps);cell.add(*all_shapes)
    return lib, cell, connection_points

def create_rect_wire_layout(
    N=6, layout_width=1400, layout_height=1400, pad_width=100, pad_height=100, pad_spacing=30, edge_margin=50, active_size=180, trace_width=10, trace_spacing=10, taper_length=20, active_entry_len=20
):
    """
    用户提供的函数，用于创建带引线的焊盘框架。
    [已修复] 返回值中增加了 all_pads，以匹配主程序的调用。
    """
    lib=gdstk.Library(); cell=lib.new_cell('PAD_FRAME_WITH_WIRES'); connection_points={}; all_pads=[]
    def natural_sort_key(s): return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)',s)]
    pg_labels=sorted([f"QD_PG{i+1}" for i in range(N)],key=natural_sort_key); bg_labels=sorted([f"QD_B{i+1}" for i in range(N+1)],key=natural_sort_key); gate_sequence=[]; [gate_sequence.extend([bg_labels[i],pg_labels[i]]) for i in range(N)]; gate_sequence.append(bg_labels[N]); top_pads=[]; num_sets=math.ceil(N/3)
    if num_sets>=1: top_pads.extend(["SET1_B1","SET1_B2","SET1_G"]); top_pads.append("SET_S");
    if num_sets>=2: top_pads.extend(["SET2_G","SET2_B2","SET2_B1"])
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
            total_span=count*pad_height+(count-1)*pad_spacing; top_y=(layout_height+total_span)/2; x=edge_margin if edge=='left' else layout_width-edge_margin-pad_width
            for i,label in enumerate(pads): y=top_y-(i+1)*pad_height-i*pad_spacing; cx,cy=x+pad_height/2,y+pad_height/2; all_pads.append({'label':label,'rect':(x,y,pad_height,pad_width),'center':(cx,cy),'edge':edge})
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
    
    # --- [已修复] ---
    # 添加 all_pads 到返回值，以匹配主程序的调用
    return lib, cell, connection_points, active_center, all_pads
    # --- 修复结束 ---

# ===================================================================
#                      分步执行与可视化逻辑
# ===================================================================

def visualize_and_save_pad_frame(lib, cell_to_show, title, gds_filename):
    # ... (此函数与上一版完全相同，保持不变)
    print(f"--- {title} ---"); print(f"正在可视化: {title}...")
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

def visualize_and_save_assembly(lib, cell_to_show, title, gds_filename):
    # ... (此函数与上一版完全相同，保持不变)
    print(f"--- {title} ---"); print(f"正在可视化: {title}...")
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

# ===================================================================
#           [核心修改] 全新：分层阶梯式 Z 型扇出布线算法
# ===================================================================
def route_tiered_z_fanout(start_pt, end_pt, chip_center, side_info, layer):
    """
    一个专业的扇出布线算法，它创建了一个分层的Z型走线结构以满足DRC。
    """
    # --- 1. 定义布线参数 ---
    W = {1: 0.04, 2: 0.1, 3: 0.5}  # 各段线宽 (40nm, 100nm, 500nm)
    S = {2: 0.1, 3: 0.5}           # 各段间距 (等于线宽)
    P = {2: W[2] + S[2], 3: W[3] + S[3]} # 中心间距 (Pitch)
    D = {1: 10, 2: 25, 3: 45}      # 三个阶梯层距中心的距离 (可调)

    p0 = np.array(start_pt); p_final = np.array(end_pt); center = np.array(chip_center)
    side, index, total = side_info

    # --- 2. 预计算阶梯落点 (解决DRC的关键) ---
    points = {}
    
    # 根据每条线在边上的索引，计算它在每个阶梯层上的偏移位置
    # 这是算法的核心，确保了间距
    offset_tier1 = (index - (total - 1) / 2.0) * P[2] # 最初用较宽的间距大致分开
    offset_tier2 = (index - (total - 1) / 2.0) * P[2]
    offset_tier3 = (index - (total - 1) / 2.0) * P[3]

    if side == 'bottom':
        points[1] = (center[0] + offset_tier1, center[1] - D[1])
        points[2] = (center[0] + offset_tier2, center[1] - D[2])
        points[3] = (center[0] + offset_tier3, center[1] - D[3])
    elif side == 'top':
        points[1] = (center[0] + offset_tier1, center[1] + D[1])
        points[2] = (center[0] + offset_tier2, center[1] + D[2])
        points[3] = (center[0] + offset_tier3, center[1] + D[3])
    elif side == 'left':
        points[1] = (center[0] - D[1], center[1] - offset_tier1)
        points[2] = (center[0] - D[2], center[1] - offset_tier2)
        points[3] = (center[0] - D[3], center[1] - offset_tier3)
    else: # right
        points[1] = (center[0] + D[1], center[1] - offset_tier1)
        points[2] = (center[0] + D[2], center[1] - offset_tier2)
        points[3] = (center[0] + D[3], center[1] - offset_tier3)

    # --- 3. 构建分段路径 ---
    p1_target = np.array(points[1]); p2_target = np.array(points[2]); p3_target = np.array(points[3])
    
    # 路径点列表
    path_points = [p0]

    # 段 1: 器件 -> Tier 1 (正交引出)
    if side in ['top','bottom']: path_points.append( (p1_target[0], p0[1]) )
    else: path_points.append( (p0[0], p1_target[1]) )
    path_points.append(p1_target)

    # 段 2: Tier 1 -> Tier 2 (Z-shape)
    if side in ['top','bottom']: path_points.append( (p2_target[0], p1_target[1]) )
    else: path_points.append( (p1_target[0], p2_target[1]) )
    path_points.append(p2_target)
    
    # 段 3: Tier 2 -> Tier 3 (Z-shape)
    if side in ['top','bottom']: path_points.append( (p3_target[0], p2_target[1]) )
    else: path_points.append( (p2_target[0], p3_target[1]) )
    path_points.append(p3_target)

    # 段 4: Tier 3 -> Pad 引线 (与Pad垂直对接)
    if side in ['top','bottom']: path_points.append( (p3_target[0], p_final[1]) )
    else: path_points.append( (p_final[0], p3_target[1]) )
    path_points.append(p_final)
    
    # --- 4. 定义宽度变化 ---
    # 宽度在每个顶点处定义
    widths = [W[1], W[1], W[2], W[2], W[3], W[3], W[3], W[3], W[3]]

    # 创建最终的FlexPath对象
    # joins='miter' 产生锐利的45度角
    path = gdstk.FlexPath(path_points, widths, layer=layer, joins="miter")

    return [path]


# ===================================================================
#                           主程序入口 (已更新)
# ===================================================================
if __name__ == "__main__":

    ACTIVE_SIZE = 200

    # --- 步骤 1: 生成并展示核心器件 ---
    lib_device, cell_device, device_points = create_6qd_layout_with_labels()
    plot_gds(cell_device, title="步骤 1: 仅带长引线的量子点器件")
    lib_device.write_gds("step1_device_with_leads.gds")

    # --- 步骤 2: 生成并展示焊盘框架 ---
    lib_pads, cell_pads, pad_points, active_center, all_pads_info = create_rect_wire_layout(active_size=ACTIVE_SIZE)
    visualize_and_save_pad_frame(lib_pads, cell_pads, "步骤 2: 仅带引线的焊盘框架", "step2_pad_frame.gds")

    # --- 创建一个查找表，用于快速获取每个标签的布线侧和索引 ---
    label_to_side_info = {}
    for pad_info in all_pads_info:
        label = pad_info['label']
        side = pad_info['edge']
        
        # 获取该侧的所有标签
        side_pads = [p for p in all_pads_info if p['edge'] == side]
        
        # 根据中心坐标对它们进行排序
        if side in ['top', 'bottom']:
            side_pads.sort(key=lambda p: p['center'][0])
        else: # left, right
            side_pads.sort(key=lambda p: p['center'][1], reverse=True)

        side_labels = [p['label'] for p in side_pads]

        if label in side_labels:
            index = side_labels.index(label)
            total = len(side_labels)
            label_to_side_info[label] = (side, index, total)

    # --- 步骤 3: 组装 ---
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
    
    visualize_and_save_assembly(assembly_lib, assembly_cell, "步骤 3: 器件已放置 (未布线)", "step3_placed_unrouted.gds")

    # --- 步骤 4 & 5: 使用全新 Fan-out 算法连接所有线 ---
    print("--- 步骤 4 & 5: 正在使用全新阶梯式 Fan-out 算法连接所有线... ---")
    ROUTING_LAYER = 5
    
    # 演示单条线
    label_to_connect = 'QD_PG4'
    if label_to_connect in device_points and label_to_connect in pad_points and label_to_connect in label_to_side_info:
        pad_coord = pad_points[label_to_connect]; local_dev_coord = device_points[label_to_connect]
        global_dev_coord = (local_dev_coord[0] + active_center[0], local_dev_coord[1] + active_center[1])
        side_info = label_to_side_info[label_to_connect]
        
        geometries = route_tiered_z_fanout(start_pt=global_dev_coord, end_pt=pad_coord, chip_center=active_center, side_info=side_info, layer=ROUTING_LAYER)
        assembly_cell.add(*geometries)

    visualize_and_save_assembly(assembly_lib, assembly_cell, f"步骤 4: 连接单条线 ({label_to_connect})", "step4_single_route.gds")

    # 连接所有剩余的线
    for label, local_dev_coord in device_points.items():
        if label != label_to_connect and label in pad_points and label in label_to_side_info:
            pad_coord = pad_points[label]; global_dev_coord = (local_dev_coord[0] + active_center[0], local_dev_coord[1] + active_center[1])
            side_info = label_to_side_info[label]
            geometries = route_tiered_z_fanout(start_pt=global_dev_coord, end_pt=pad_coord, chip_center=active_center, side_info=side_info, layer=ROUTING_LAYER)
            assembly_cell.add(*geometries)
            
    visualize_and_save_assembly(assembly_lib, assembly_cell, "步骤 5: 最终完整布线", "step5_fully_routed.gds")

    print("全部分步流程执行完毕！")