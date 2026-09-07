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

# --- User-Provided Functions (plot_gds, create_6qd_layout_with_labels, create_rect_wire_layout) ---
def plot_gds(cell, title="量子点器件布局"):
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
    pg_max_width=0.120, pg_vert_side_len=0.040, pg_chamfer_h=0.040, pg_bot_flat_w=0.042, pg_top_flat_w=0.040, bg_max_width=0.060, bg_vert_side_len=0.040, bg_top_flat_w=0.020, bg_bot_flat_w=0.042, gap_pg_bg=0.002, d1_gap=0.020, d2_gap=0.020, sd_gap_to_gate=0.01, gap_gate_outer_sg=0.030, sg_mid_thick=0.100,   sg_top_thick=0.300,   sg_bot_thick=0.300,   sg_extension=0.200, lead_width=0.042, lead_length_bot=0.4, lead_length_top=0.5, lead_overlap=0.025, sd_height=0.10, sd_width=0.35, sd_lead_width=0.060, sd_outer_length=0.300
):
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
    sgxs=bmx-sg_extension;sgxe=bax+sg_extension;ysgmb=max(bot_tips_y)+d1_gap;ysgmt=ysgmb+sg_mid_thick;
    ysgtb=yct+gate_half_height+gap_gate_outer_sg;ysgtt=ysgtb+sg_top_thick; # Define ysgtb and ysgtt here
    sg_mid=gdstk.rectangle((sgxs,ysgmb),(sgxe,ysgmt),layer=LAYER_SG);
    sg_top=gdstk.rectangle((sgxs,ysgtb),(sgxe,ysgtt),layer=LAYER_SG);
    ysgbt=ycb-gate_half_height-gap_gate_outer_sg;ysgbb=ysgbt-sg_bot_thick;
    sg_bot=gdstk.rectangle((sgxs,ysgbb),(sgxe,ysgbt),layer=LAYER_SG);all_shapes.extend([sg_top,sg_mid,sg_bot])
    lsg1=gdstk.Label("SG1",(sgxs,(ysgtb+ysgtt)/2),layer=LAYER_LABEL);cell.add(lsg1);connection_points["SG1"]=lsg1.origin;lsg2=gdstk.Label("SG2",(sgxs,(ysgmb+ysgmt)/2),layer=LAYER_LABEL);cell.add(lsg2);connection_points["SG2"]=lsg2.origin;lsg3=gdstk.Label("SG3",(sgxs,(ysgbb+ysgbt)/2),layer=LAYER_LABEL);cell.add(lsg3);connection_points["SG3"]=lsg3.origin
    sd_shps=[];sdcr=0.020
    def cr(xmin,ymin,xmax,ymax,crnrs):poly=gdstk.Polygon([(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)],layer=LAYER_SD);poly.fillet([sdcr if c else 0 for c in crnrs],tolerance=0.001);return poly
    blbe=bmx;blsix=blbe-sd_gap_to_gate;brbe=bax;brsix=brbe+sd_gap_to_gate;golx=blsix-sd_outer_length;gorx=brsix+sd_outer_length;sd_shps.append(cr(golx,ycb-sd_height/2,blsix,ycb+sd_height/2,[False,True,True,False]));lqdd=gdstk.Label("QD_D",(golx,ycb),layer=LAYER_LABEL);cell.add(lqdd);connection_points["QD_D"]=lqdd.origin;sd_shps.append(cr(brsix,ycb-sd_height/2,gorx,ycb+sd_height/2,[True,False,False,True]));lqds=gdstk.Label("QD_S",(gorx,ycb),layer=LAYER_LABEL);cell.add(lqds);connection_points["QD_S"]=lqds.origin;shw=bg_max_width+gap_pg_bg+pg_max_width/2;lsle=lscpgx-shw;tlsix=lsle-sd_gap_to_gate;rsre=rscpgx+shw;trsix=rsre-sd_gap_to_gate;sd_shps.append(cr(golx,yct-sd_height/2,tlsix,yct+sd_height/2,[False,True,True,False]));ls1d=gdstk.Label("SET1_D",(golx,yct),layer=LAYER_LABEL);cell.add(ls1d);connection_points["SET1_D"]=ls1d.origin;sd_shps.append(cr(trsix,yct-sd_height/2,gorx,yct+sd_height/2,[True,False,False,True]));ls2d=gdstk.Label("SET2_D",(gorx,yct),layer=LAYER_LABEL);cell.add(ls2d);connection_points["SET2_D"]=ls2d.origin;lsre=lscpgx+shw;smxm=lsre+sd_gap_to_gate;rsle=rscpgx-shw;smxa=rsle-sd_gap_to_gate;msdp=cr(smxm,yct-sd_height/2,smxa,yct+sd_height/2,[True,True,True,True]);msdcx=(smxm+smxa)/2;mlye=yct+lead_length_top;msdl=gdstk.rectangle((msdcx-sd_lead_width/2,yct),(msdcx+sd_lead_width/2,mlye),layer=LAYER_SD);sd_shps.extend(gdstk.boolean([msdp,msdl],[],'or',layer=LAYER_SD));lsets=gdstk.Label("SET_S",(msdcx,mlye),layer=LAYER_LABEL);cell.add(lsets);connection_points["SET_S"]=lsets.origin;all_shapes.extend(sd_shps);cell.add(*all_shapes)
    return lib, cell, connection_points

def create_rect_wire_layout(
    N=6, layout_width=1400, layout_height=1400, pad_width=100, pad_height=100, pad_spacing=30, edge_margin=50, active_size=180, trace_width=10, trace_spacing=10, taper_length=40, active_entry_len=40
):
    lib=gdstk.Library();cell=lib.new_cell('PAD_FRAME_WITH_WIRES');connection_points={};all_pads=[];
    def nsk(s):return[int(t)if t.isdigit()else t.lower()for t in re.split('([0-9]+)',s)]
    pgl=sorted([f"QD_PG{i+1}"for i in range(N)],key=nsk);bgl=sorted([f"QD_B{i+1}"for i in range(N+1)],key=nsk);gs=[];[gs.extend([bgl[i],pgl[i]])for i in range(N)];gs.append(bgl[N]);tp=[];ns=math.ceil(N/3);
    if ns>=1:tp.extend(["SET1_B1","SET1_B2","SET1_G"]);tp.append("SET_S");
    if ns>=2:tp.extend(["SET2_G","SET2_B2","SET2_B1"])
    lf=["SG1"];
    if ns>=1:lf.append("SET1_D");lf.extend(["SG2","QD_D"]);rf=[]
    if ns>=2:rf.append("SET2_D");rf.append("QD_S");rf.extend(["GND","SG3","GND"])
    ls=len(gs);ll=len(lf);lr=len(rf);lt=len(tp);aw=layout_width-2*edge_margin;bc=(0,ls,0);md=float('inf');ss=max(ll,lr)
    for ts in range(ss,ss+ls):
        al=ts-ll;ar=ts-lr;
        if al+ar>ls:break
        nb=ls-(al+ar)
        if nb>0 and(nb*pad_width+(nb-1)*pad_spacing>aw):continue
        d=abs(nb-lt);
        if d<md:md=d;bc=(al,nb,ar)
    nl,nb,nr=bc;ple=gs[:nl];pb=gs[nl:nl+nb];pre=gs[nl+nb:];lp=lf+ple;rp=rf+list(reversed(pre))
    def ac(pads,edge):
        c=len(pads)
        if c==0:return
        if edge in['top','bottom']:
            ts=c*pad_width+(c-1)*pad_spacing;sx=(layout_width-ts)/2;y=layout_height-edge_margin-pad_height if edge=='top'else edge_margin
            for i,lbl in enumerate(pads):x=sx+i*(pad_width+pad_spacing);cx,cy=x+pad_width/2,y+pad_height/2;all_pads.append({'label':lbl,'rect':(x,y,pad_width,pad_height),'center':(cx,cy),'edge':edge})
        else:
            ts=c*pad_height+(c-1)*pad_spacing;ty=(layout_height+ts)/2;x=edge_margin if edge=='left'else layout_width-edge_margin-pad_width
            for i,lbl in enumerate(pads):y=ty-(i+1)*pad_height-i*pad_spacing;cx,cy=x+pad_height/2,y+pad_height/2;all_pads.append({'label':lbl,'rect':(x,y,pad_height,pad_width),'center':(cx,cy),'edge':edge})
    ac(tp,'top');ac(pb,'bottom');ac(lp,'left');ac(rp,'right')
    cx=layout_width/2;cy=layout_height/2;acenter=(cx,cy);ax1=cx-active_size/2;ax2=cx+active_size/2;ay1=cy-active_size/2;ay2=cy+active_size/2;pbe={'top':[],'bottom':[],'left':[],'right':[]};[pbe[p['edge']].append(p)for p in all_pads]
    pbe['top'].sort(key=lambda p:p['center'][0]);pbe['bottom'].sort(key=lambda p:p['center'][0]);pbe['left'].sort(key=lambda p:p['center'][1],reverse=True);pbe['right'].sort(key=lambda p:p['center'][1],reverse=True);pitch=trace_width+trace_spacing
    def at(pads,edge):
        n=len(pads)
        if n==0:return
        ind=np.arange(n)-(n-1)/2;offs=ind*pitch
        if edge=='top':[p.update({'target':(cx+offs[i],ay2),'target_edge':'top'})for i,p in enumerate(pads)]
        elif edge=='bottom':[p.update({'target':(cx+offs[i],ay1),'target_edge':'bottom'})for i,p in enumerate(pads)]
        elif edge=='left':[p.update({'target':(ax1,cy+offs[-(i+1)]),'target_edge':'left'})for i,p in enumerate(pads)]
        elif edge=='right':[p.update({'target':(ax2,cy+offs[-(i+1)]),'target_edge':'right'})for i,p in enumerate(pads)]
    [at(pbe[edge],edge)for edge in['top','bottom','left','right']]
    LS=0;LP=1;LT=1;LTXT=2;LA=10;LSIO2=11;sm=10
    cell.add(gdstk.rectangle((0,0),(layout_width,layout_height),layer=LS));cell.add(gdstk.rectangle((ax1,ay1),(ax2,ay2),layer=LA))
    def ct(pr,edge,tl,tw):
        x,y,w,h=pr;cx,cy=x+w/2,y+h/2
        if edge=='top':pts=[(x,y),(x+w,y),(cx+tw/2,y-tl),(cx-tw/2,y-tl)];start=(cx,y-tl)
        elif edge=='bottom':pts=[(x,y+h),(x+w,y+h),(cx+tw/2,y+h+tl),(cx-tw/2,y+h+tl)];start=(cx,y+h+tl)
        elif edge=='left':pts=[(x+w,y),(x+w,y+h),(x+w+tl,cy+tw/2),(x+w+tl,cy-tw/2)];start=(x+w+tl,cy)
        else:pts=[(x,y),(x,y+h),(x-tl,cy+tw/2),(x-tl,cy-tw/2)];start=(x-tl,cy)
        return gdstk.Polygon(pts,layer=LT),start
    def rsz(sp,tp,edge,width):
        x1,y1=sp;x2,y2=tp;sm=20;pts=[(x1,y1)]
        if edge in['top','bottom']:
            sy=-1 if edge=='top'else 1;ys=y1+sy*sm;pts.append((x1,ys));ye=y2+sy*active_entry_len;ah=abs(ye-ys);nw=abs(x2-x1);sx=1 if x2>x1 else -1
            if nw>ah:hl=ah/2;ym=ys+sy*hl;xm1=x1+sx*hl;xm2=x2-sx*hl;pts.extend([(xm1,ym),(xm2,ym)])
            else:dyd=nw;yte=ys+sy*dyd;pts.append((x2,yte))
            pts.append((x2,y2))
        else:
            sx=1 if edge=='left'else -1;xs=x1+sx*sm;pts.append((xs,y1));xe=x2+sx*active_entry_len;aw=abs(xe-xs);nh=abs(y2-y1);sy=1 if y2>y1 else -1
            if nh>aw:wl=aw/2;xm=xs+sx*wl;ym1=y1+sy*wl;ym2=y2-sy*wl;pts.extend([(xm,ym1),(xm,ym2)])
            else:dxd=nh;xte=xs+sx*dxd;pts.append((xte,y2))
            pts.append((x2,y2))
        return gdstk.FlexPath(pts,width,ends="flush",joins="round",layer=LT),pts
    for p in all_pads:
        x,y,w,h=p['rect'];lbl=p['label'];edge=p['edge']
        cell.add(gdstk.rectangle((x-sm,y-sm),(x+w+sm,y+h+sm),layer=LSIO2));cell.add(gdstk.rectangle((x,y),(x+w,y+h),layer=LP));cx,cy=x+w/2,y+h/2;rot=0 if w>h else math.pi/2
        cell.add(gdstk.Label(lbl,(cx,cy),layer=LTXT,magnification=10,rotation=rot));poly_t,t_start=ct(p['rect'],edge,taper_length,trace_width);cell.add(poly_t);path_obj,path_pts=rsz(t_start,p['target'],edge,trace_width);cell.add(path_obj);p['trace_points']=path_pts;connection_points[p['label']]=p['target']
    return lib, cell, connection_points, acenter, all_pads

# ===================================================================
#                      分步执行与可视化逻辑
# ===================================================================
def visualize_and_save_pad_frame(lib, cell_to_show, title, gds_filename):
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
#           [核心修改] 全新 Fan-out 结构生成与布线算法
# ===================================================================
class RoutingGenerator:
    """
    一个全新的布线管理器，负责全局规划和DRC检查。
    """
    def __init__(self, chip_center, active_region_bbox, quantum_dot_polygons, layer, pad_connection_extension_len=5.0):
        self.chip_center = chip_center
        self.active_region_bbox = active_region_bbox
        self.quantum_dot_polygons = quantum_dot_polygons
        self.layer = layer
        self.route_requests = []
        self.routes = [] # Stores {'geo': gdstk.Polygon, 'label': str, 'path_points': list, 'path_widths': list, 'max_width': float}

        # 线宽定义 (nm -> um)
        self.W = {
            'tier1': 0.040,  # 40nm
            'tier2': 0.100,  # 100nm
            'tier3': 0.500   # 500nm
        }
        # 间距定义 (间距必须大于线宽)
        self.S = {
            'tier1': 0.200,  # Increased spacing for 40nm lines
            'tier2': 0.400,  # Increased spacing for 100nm lines
            'tier3': 1.500   # Increased spacing for 500nm lines
        }
        # 初始引出线长度和pad连接线长度
        self.initial_extension_len = 0.500  # 500nm长度的40nm宽度线
        self.pad_connection_extension_len = pad_connection_extension_len  # 500nm线连接pad引线的部分，垂直引出10um (default 5um)

        # 用于布线规划的参数 (可以根据需要调整)
        # These radial distances define the "bands" where each tier of wiring will primarily run.
        # They are relative to the chip_center.
        self.tier_radial_distances = {
            'tier1': 150.0, # Increased for better spacing
            'tier2': 200.0, # Increased for better spacing
            'tier3': 300.0  # Increased for better spacing
        }
        self.min_segment_length = 0.1 # 最小线段长度，避免生成过短的线段

    def add_route_request(self, label, start_pt, end_pt, side_info, electrode_direction=None):
        self.route_requests.append({
            'label': label,
            'start_pt': start_pt,
            'end_pt': end_pt,
            'side_info': side_info,
            'electrode_direction': electrode_direction
        })

    def plan_and_generate_routes(self):
        print("\n--- 开始全局布线规划 ---")
        all_geometries = []
        
        requests_by_side = {'top': [], 'bottom': [], 'left': [], 'right': []}
        for req in self.route_requests:
            requests_by_side[req['side_info'][0]].append(req)

        for side, requests in requests_by_side.items():
            if not requests:
                continue
            if side in ['top', 'bottom']:
                requests.sort(key=lambda r: r['start_pt'][0])
            else:
                requests.sort(key=lambda r: r['start_pt'][1], reverse=True)

            for i, req in enumerate(requests):
                label = req['label']
                start_pt = np.array(req['start_pt'])
                end_pt = np.array(req['end_pt'])
                electrode_direction = req['electrode_direction']

                target_waypoints = self._calculate_target_waypoints(start_pt, end_pt, side, i, len(requests), electrode_direction)

                full_path_points = []
                full_path_widths = []

                # Segment 1: Initial extension (40nm)
                segment1_points = self._create_orthogonal_path(target_waypoints['p0'], target_waypoints['p0_ext'])
                full_path_points.extend(segment1_points)
                full_path_widths.extend([self.W['tier1']] * len(segment1_points))

                # Verify initial lead width
                if len(segment1_points) > 1:
                    first_segment_length = np.linalg.norm(np.array(segment1_points[1]) - np.array(segment1_points[0]))
                    if not np.isclose(first_segment_length, self.initial_extension_len, atol=1e-3):
                        print(f"Warning: Route {label} initial segment length is {first_segment_length:.3f}um, expected {self.initial_extension_len}um.")
                    if not np.isclose(full_path_widths[0], self.W['tier1'], atol=1e-3):
                        print(f"Warning: Route {label} initial segment width is {full_path_widths[0]:.3f}um, expected {self.W['tier1']}um.")

                # Segment 2: 40nm to 100nm transition
                full_path_points.append(target_waypoints['p0_ext_to_p1'])
                full_path_widths.append(self.W['tier1'])
                full_path_points.append(target_waypoints['p0_ext_to_p1'])
                full_path_widths.append(self.W['tier2'])

                segment2_points = self._create_orthogonal_path(target_waypoints['p0_ext_to_p1'], target_waypoints['p1_to_p2'])
                full_path_points.extend(segment2_points[1:])
                full_path_widths.extend([self.W['tier2']] * (len(segment2_points) - 1))

                # Segment 3: 100nm to 500nm transition
                full_path_points.append(target_waypoints['p1_to_p2'])
                full_path_widths.append(self.W['tier2'])
                full_path_points.append(target_waypoints['p1_to_p2'])
                full_path_widths.append(self.W['tier3'])

                segment3_points = self._create_orthogonal_path(target_waypoints['p1_to_p2'], target_waypoints['p2_to_p3'])
                full_path_points.extend(segment3_points[1:])
                full_path_widths.extend([self.W['tier3']] * (len(segment3_points) - 1))

                # Segment 4: 500nm to pad connection
                segment4_points = self._create_orthogonal_path(target_waypoints['p2_to_p3'], target_waypoints['p3_ext'])
                full_path_points.extend(segment4_points[1:])
                full_path_widths.extend([self.W['tier3']] * (len(segment4_points) - 1))

                # Final connection to pad
                segment_final_points = self._create_orthogonal_path(target_waypoints['p3_ext'], target_waypoints['p_final'])
                full_path_points.extend(segment_final_points[1:])
                full_path_widths.extend([self.W['tier3']] * (len(segment_final_points) - 1))

                if len(full_path_points) > 1 and len(full_path_points) == len(full_path_widths):
                    full_route_path = gdstk.FlexPath(full_path_points, full_path_widths, joins="miter", ends="flush", layer=self.layer)
                    merged_route_geo = gdstk.boolean(full_route_path.to_polygons(), [], 'or', layer=self.layer)
                    if merged_route_geo:
                        max_width_for_route = max(full_path_widths) if full_path_widths else 0.0
                        self.routes.append({
                            'geo': merged_route_geo[0], 
                            'label': label, 
                            'path_points': full_path_points, 
                            'path_widths': full_path_widths,
                            'max_width': max_width_for_route # Store max width here
                        })
                        all_geometries.extend(merged_route_geo)
                else:
                    print(f"Warning: Skipping route {label} due to insufficient points or width mismatch.")

        print("全局规划和几何生成完毕。")
        return all_geometries

    def _calculate_target_waypoints(self, p0, p_final, side, index, total_on_side, electrode_direction):
        waypoints = {'p0': p0, 'p_final': p_final}
        center_x, center_y = self.chip_center
        active_min_x, active_min_y = self.active_region_bbox[0]
        active_max_x, active_max_y = self.active_region_bbox[1]

        # 1. Initial orthogonal extension (40nm width, 500nm length)
        p0_ext = np.array(p0)
        initial_extension_len_actual = self.initial_extension_len

        # Ensure initial extension is orthogonal and respects electrode direction
        if electrode_direction == 'down':
            p0_ext[1] -= initial_extension_len_actual
        elif electrode_direction == 'up':
            p0_ext[1] += initial_extension_len_actual
        elif electrode_direction == 'left':
            p0_ext[0] -= initial_extension_len_actual
        elif electrode_direction == 'right':
            p0_ext[0] += initial_extension_len_actual
        else: # Default behavior if no specific electrode direction
            if side == 'top':
                p0_ext[1] += initial_extension_len_actual
            elif side == 'bottom':
                p0_ext[1] -= initial_extension_len_actual
            elif side == 'left':
                p0_ext[0] -= initial_extension_len_actual
            elif side == 'right':
                p0_ext[0] += initial_extension_len_actual
        waypoints['p0_ext'] = p0_ext

        # Calculate lateral offsets for each tier based on its width and spacing
        # This ensures lines within a tier are properly spaced and centered.
        
        # Calculate the total width/height needed for all routes on this side for a given tier
        # This is used to center the group of routes.
        total_width_needed_tier1 = total_on_side * self.W['tier1'] + (total_on_side - 1) * self.S['tier1']
        total_width_needed_tier2 = total_on_side * self.W['tier2'] + (total_on_side - 1) * self.S['tier2']
        total_width_needed_tier3 = total_on_side * self.W['tier3'] + (total_on_side - 1) * self.S['tier3']

        # Calculate the starting coordinate for the first route on this side, centered
        if side in ['top', 'bottom']: # Lateral spreading along X-axis
            start_lateral_coord_tier1 = center_x - total_width_needed_tier1 / 2.0
            start_lateral_coord_tier2 = center_x - total_width_needed_tier2 / 2.0
            start_lateral_coord_tier3 = center_x - total_width_needed_tier3 / 2.0
            
            # Calculate the current route's lateral position (center of the wire)
            current_lateral_pos_tier1 = start_lateral_coord_tier1 + index * (self.W['tier1'] + self.S['tier1']) + self.W['tier1'] / 2.0
            current_lateral_pos_tier2 = start_lateral_coord_tier2 + index * (self.W['tier2'] + self.S['tier2']) + self.W['tier2'] / 2.0
            current_lateral_pos_tier3 = start_lateral_coord_tier3 + index * (self.W['tier3'] + self.S['tier3']) + self.W['tier3'] / 2.0

            # Tier 1 (40nm) waypoint: Move radially out, then laterally
            p0_ext_to_p1 = np.array([current_lateral_pos_tier1, (active_min_y if side == 'top' else active_max_y) + (self.tier_radial_distances['tier1'] * (1 if side == 'top' else -1))])
            
            # Tier 2 (100nm) waypoint: Move further radially out, then adjust lateral position
            p1_to_p2 = np.array([current_lateral_pos_tier2, (active_min_y if side == 'top' else active_max_y) + (self.tier_radial_distances['tier2'] * (1 if side == 'top' else -1))])

            # Tier 3 (500nm) waypoint: Move even further radially out, then adjust lateral position
            p2_to_p3 = np.array([current_lateral_pos_tier3, (active_min_y if side == 'top' else active_max_y) + (self.tier_radial_distances['tier3'] * (1 if side == 'top' else -1))])

        else: # 'left' or 'right' - Lateral spreading along Y-axis
            start_lateral_coord_tier1 = center_y - total_width_needed_tier1 / 2.0
            start_lateral_coord_tier2 = center_y - total_width_needed_tier2 / 2.0
            start_lateral_coord_tier3 = center_y - total_width_needed_tier3 / 2.0

            current_lateral_pos_tier1 = start_lateral_coord_tier1 + index * (self.W['tier1'] + self.S['tier1']) + self.W['tier1'] / 2.0
            current_lateral_pos_tier2 = start_lateral_coord_tier2 + index * (self.W['tier2'] + self.S['tier2']) + self.W['tier2'] / 2.0
            current_lateral_pos_tier3 = start_lateral_coord_tier3 + index * (self.W['tier3'] + self.S['tier3']) + self.W['tier3'] / 2.0

            # Tier 1 (40nm) waypoint: Move radially out, then laterally
            p0_ext_to_p1 = np.array([(active_min_x if side == 'left' else active_max_x) + (self.tier_radial_distances['tier1'] * (1 if side == 'left' else -1)), current_lateral_pos_tier1])
            
            # Tier 2 (100nm) waypoint: Move further radially out, then adjust lateral position
            p1_to_p2 = np.array([(active_min_x if side == 'left' else active_max_x) + (self.tier_radial_distances['tier2'] * (1 if side == 'left' else -1)), current_lateral_pos_tier2])

            # Tier 3 (500nm) waypoint: Move even further radially out, then adjust lateral position
            p2_to_p3 = np.array([(active_min_x if side == 'left' else active_max_x) + (self.tier_radial_distances['tier3'] * (1 if side == 'left' else -1)), current_lateral_pos_tier3])

        waypoints['p0_ext_to_p1'] = p0_ext_to_p1
        waypoints['p1_to_p2'] = p1_to_p2
        waypoints['p2_to_p3'] = p2_to_p3

        # 5. Orthogonal extension from p2_to_p3 towards the pad (p3_ext)
        p3_ext = np.array(p_final)
        if side == 'top':
            p3_ext[1] -= self.pad_connection_extension_len
        elif side == 'bottom':
            p3_ext[1] += self.pad_connection_extension_len
        elif side == 'left':
            p3_ext[0] += self.pad_connection_extension_len
        elif side == 'right':
            p3_ext[0] -= self.pad_connection_extension_len
        waypoints['p3_ext'] = p3_ext

        return waypoints

    def _create_orthogonal_path(self, p_start, p_end):
        # This function generates a path segment between p_start and p_end
        # ensuring only orthogonal turns (L-bends).
        # It returns a list of points including p_start and p_end, and one intermediate point if needed.

        path_segment = [p_start]
        
        dx = p_end[0] - p_start[0]
        dy = p_end[1] - p_start[1]

        # If points are already aligned horizontally or vertically
        if np.isclose(dx, 0) and np.isclose(dy, 0): # Points are identical
            return path_segment
        elif np.isclose(dx, 0): # Already vertical
            path_segment.append(p_end)
            return path_segment
        elif np.isclose(dy, 0): # Already horizontal
            path_segment.append(p_end)
            return path_segment
        
        # Create an L-bend. Prioritize horizontal then vertical.
        mid_pt = (p_end[0], p_start[1])
        path_segment.append(mid_pt)
        path_segment.append(p_end)
        
        return path_segment

    def check_drc(self):
        print("\n--- 开始执行DRC检查 ---")
        violations_found = False
        num_routes = len(self.routes)

        if num_routes < 1:
            print("没有路由，无需检查。")
            return False

        # 1. Angle Check (Orthogonal)
        print("进行角度检查...")
        for route_info in self.routes:
            path_points = route_info['path_points']
            label = route_info['label']
            
            if len(path_points) < 2:
                continue

            for k in range(1, len(path_points) - 1):
                p1 = np.array(path_points[k-1])
                p2 = np.array(path_points[k])
                p3 = np.array(path_points[k+1])
                
                v1 = p2 - p1
                v2 = p3 - p2
                
                # Avoid division by zero for zero-length segments
                if np.linalg.norm(v1) < 1e-9 or np.linalg.norm(v2) < 1e-9:
                    continue

                # Calculate angle between segments
                dot_product = np.dot(v1, v2)
                magnitudes = np.linalg.norm(v1) * np.linalg.norm(v2)
                
                if magnitudes < 1e-9:
                    continue

                angle_rad = np.arccos(np.clip(dot_product / magnitudes, -1.0, 1.0))
                angle_deg = np.degrees(angle_rad)

                is_valid_angle = False
                for target_angle in [0, 90, 180]:
                    if np.isclose(angle_deg, target_angle, atol=5):
                        is_valid_angle = True
                        break
                
                if not is_valid_angle:
                    print(f"[DRC角度违规] 路由 '{label}' 在点 {p2} 发现一个角度为 {angle_deg:.2f}° 的非正交转角。")
                    violations_found = True

        # 2. Spacing Check (between different lines)
        print("进行间距检查...")
        for i in range(num_routes):
            for j in range(i + 1, num_routes):
                route_i_info = self.routes[i]
                route_j_info = self.routes[j]
                route_i_geo = route_i_info['geo']
                route_j_geo = route_j_info['geo']

                max_width_i = route_i_info['max_width']
                max_width_j = route_j_info['max_width']
                
                required_clearance = max(max_width_i, max_width_j) + 0.001
                
                buffered_route_i = gdstk.offset([route_i_geo], required_clearance / 2.0, join="miter", layer=self.layer + 1)
                
                if buffered_route_i:
                    overlap_check = gdstk.boolean(buffered_route_i, [route_j_geo], 'and', layer=self.layer + 2)
                    if overlap_check:
                        violations_found = True
                        print(f"[DRC间距违规] 路由 '{self.routes[i]['label']}' 和 '{self.routes[j]['label']}' 之间间距不足或重叠。")

        # 3. Active Region and Quantum Dot Overlap Check
        print("进行区域和量子点重叠检查...")
        active_region_poly = gdstk.rectangle(self.active_region_bbox[0], self.active_region_bbox[1], layer=self.layer + 3)
        
        for route_info in self.routes:
            route_geo = route_info['geo']
            
            intersection_with_active = gdstk.boolean([route_geo], [active_region_poly], 'and', layer=self.layer + 4)
            if not intersection_with_active or abs(route_geo.area() - intersection_with_active[0].area()) > 1e-9:
                violations_found = True
                print(f"[DRC区域违规] 路由 '{route_info['label']}' 未完全在Active Region内。")
            
            for qd_poly in self.quantum_dot_polygons:
                overlap_qd = gdstk.boolean([route_geo], [qd_poly], 'and', layer=self.layer + 5)
                if overlap_qd:
                    violations_found = True
                    print(f"[DRC重叠违规] 路由 '{route_info['label']}' 与量子点重叠。")

        if not violations_found:
            print("[DRC通过] 未发现违规。")
        return violations_found

    def _get_route_max_width(self, route_geo):
        for route_entry in self.routes:
            if route_entry['geo'] == route_geo:
                if 'max_width' in route_entry:
                    return route_entry['max_width']
        
        if isinstance(route_geo, gdstk.Polygon):
            bbox = route_geo.bounding_box()
            if bbox:
                (min_x, min_y), (max_x, max_y) = bbox
                return min(max_x - min_x, max_y - min_y)
        return 0.0

# ===================================================================
#                           主程序入口
# ===================================================================
if __name__ == "__main__":

    ACTIVE_SIZE = 180 # Reverted to original size
    ROUTING_LAYER = 5 # Using layer 5 for routing as per existing code

    lib_device, cell_device, device_points = create_6qd_layout_with_labels()
    # plot_gds(cell_device, title="步骤 1: 仅带长引线的量子点器件")
    lib_device.write_gds("step1_device_with_leads.gds")

    # Set pad_connection_extension_len to 10um as requested (default is 5um in class)
    PAD_CONNECTION_EXTENSION_LEN_UM = 10.0 

    lib_pads, cell_pads, pad_points, active_center, all_pads_info = create_rect_wire_layout(active_size=ACTIVE_SIZE, edge_margin=50, active_entry_len=40)
    # visualize_and_save_pad_frame(lib_pads, cell_pads, "步骤 2: 仅带引线的焊盘框架", "step2_pad_frame.gds")

    label_to_side_info = {}
    pads_by_edge = {'top': [], 'bottom': [], 'left': [], 'right': []}

    for pad_info in all_pads_info:
        pads_by_edge[pad_info['edge']].append(pad_info)

    for side, pads in pads_by_edge.items():
        if side in ['top', 'bottom']:
            pads.sort(key=lambda p: p['center'][0])
        else: # 'left', 'right'
            pads.sort(key=lambda p: p['center'][1], reverse=True)
        
        for index, pad_info in enumerate(pads):
            label_to_side_info[pad_info['label']] = (side, index, len(pads))

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
    
    print("--- 步骤 4 & 5: 正在使用全新 RoutingGenerator 连接所有线... ---")

    # Extract active region bounding box
    active_region_bbox = (
        (active_center[0] - ACTIVE_SIZE / 2, active_center[1] - ACTIVE_SIZE / 2),
        (active_center[0] + ACTIVE_SIZE / 2, active_center[1] + ACTIVE_SIZE / 2)
    )

    # Extract quantum dot polygons from cell_device
    quantum_dot_polygons = []
    for layer in [0, 2, 3]: # S/D, BG, PG layers
        quantum_dot_polygons.extend(cell_device.get_polygons(layer=layer, datatype=0))
    
    # Translate quantum dot polygons to their global position
    translated_qd_polygons = []
    for poly in quantum_dot_polygons:
        translated_poly = gdstk.Polygon(poly.points) # Create a copy
        translated_poly.translate(active_center[0], active_center[1])
        translated_qd_polygons.append(translated_poly)

    manager = RoutingGenerator(
        chip_center=active_center,
        active_region_bbox=active_region_bbox,
        quantum_dot_polygons=translated_qd_polygons,
        layer=ROUTING_LAYER,
        pad_connection_extension_len=PAD_CONNECTION_EXTENSION_LEN_UM # Use 10um as requested
    )

    # Define electrode directions for specific QDs
    electrode_directions = {
        "QD_B1": "down", "QD_PG1": "down", "QD_B2": "down",
        "QD_PG2": "down", "QD_B3": "down", "QD_PG3": "down",
        "QD_B4": "down", "QD_PG4": "down", "QD_B5": "down",
        "QD_PG5": "down", "QD_B6": "down", "QD_PG6": "down",
        "QD_B7": "down"
    }
    
    for label, local_dev_coord in device_points.items():
        if label in pad_points and label in label_to_side_info:
            pad_coord = pad_points[label]
            global_dev_coord = (local_dev_coord[0] + active_center[0], local_dev_coord[1] + active_center[1])
            side_info = label_to_side_info[label]
            
            # Get electrode direction, default to None if not specified
            direction = electrode_directions.get(label, None)
            
            manager.add_route_request(
                label=label,
                start_pt=global_dev_coord,
                end_pt=pad_coord,
                side_info=side_info,
                electrode_direction=direction
            )

    all_final_geometries = manager.plan_and_generate_routes()
    assembly_cell.add(*all_final_geometries)
            
    manager.check_drc()
    visualize_and_save_assembly(assembly_lib, assembly_cell, "步骤 5: 最终完整布线", "step5_fully_routed.gds")

    print("全部分步流程执行完毕！")

</final_file_content>

IMPORTANT: For any future changes to this file, use the final_file_content shown above as your reference. This content reflects the current state of the file, including any auto-formatting (e.g., if you used single quotes but the formatter converted them to double quotes). Always base your SEARCH/REPLACE operations on this final version to ensure accuracy.



# TODO LIST RECOMMENDED
When starting a new task, it is recommended to create a todo list.




1. Include the task_progress parameter in your next tool call

2. Create a comprehensive checklist of all steps needed

3. Use markdown format: - [ ] for incomplete, - [x] for complete



**Benefits of creating a todo list now:**

	- Clear roadmap for implementation

	- Progress tracking throughout the task

	- Nothing gets forgotten or missed

	- Users can see, monitor, and edit the plan



**Example structure:**
```

- [ ] Analyze requirements

- [ ] Set up necessary files

- [ ] Implement main functionality

- [ ] Handle edge cases

- [ ] Test the implementation

- [ ] Verify results
```



Keeping the todo list updated helps track progress and ensures nothing is missed.
<environment_details>
# VSCode Visible Files
GDS_Drawer/Prefab-6-dot/20251228/test.ipynb
GDS_Drawer/Prefab-6-dot/20251229/new_routing_code.py

# VSCode Open Tabs
GDS_Drawer/Prefab-6-dot/20251227/test.ipynb
GDS_Drawer/Prefab-6-dot/20251228/test.ipynb
GDS_Drawer/Prefab-6-dot/20251229/new_routing_code.py
GDS_Drawer/Prefab-6-dot/20251229/code.py

# Recently Modified Files
These files have been modified since you last accessed them (file was just edited so you may need to re-read it before editing):
GDS_Drawer/Prefab-6-dot/20251229/new_routing_code.py

# Current Time
12/29/2025, 2:16:03 PM (Asia/Shanghai, UTC+8:00)

# Context Window Usage
914,781 / 1,048.576K tokens used (87%)

# Current Mode
ACT MODE
</environment_details>
