import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import math
import re
import os

# ===================================================================
#                      解决中文显示问题的代码
# ===================================================================
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
# ===================================================================

class QuantumDotPadGenerator18Dot:
    def __init__(self):
        # Default parameters for quantum dot generation
        self.qd_params = {
            'pg_max_width': 0.120, 'pg_vert_side_len': 0.040, 'pg_chamfer_h': 0.040,
            'pg_bot_flat_w': 0.042, 'pg_top_flat_w': 0.040, 'bg_max_width': 0.060,
            'bg_vert_side_len': 0.040, 'bg_top_flat_w': 0.020, 'bg_bot_flat_w': 0.042,
            'gap_pg_bg': 0.002, 'd1_gap': 0.020, 'd2_gap': 0.020,
            'sd_gap_to_gate': 0.01, 'gap_gate_outer_sg': 0.030, 'sg_mid_thick': 0.100,
            'sg_top_thick': 0.300, 'sg_bot_thick': 0.300, 'sg_extension': 0.200,
            'lead_width': 0.042, 'lead_length_bot': 0.4, 'lead_length_top': 0.5,
            'lead_overlap': 0.025, 'sd_height': 0.10, 'sd_width': 0.35,
            'sd_lead_width': 0.060, 'sd_outer_length': 0.300
        }

        # Default parameters for pad frame generation
        self.pad_params = {
            'N': 18, 'layout_width': 2500, 'layout_height': 2500,
            'pad_width': 100, 'pad_height': 100, 'pad_spacing': 30,
            'edge_margin': 100, 'active_size': 500, 'trace_width': 10,
            'trace_spacing': 10, 'taper_length': 40, 'active_entry_len': 40
        }

        self.device_lib = None
        self.device_cell = None
        self.device_connection_points = None
        self.pad_lib = None
        self.pad_cell = None
        self.pad_connection_points = None
        self.active_center = None
        self.all_pads_info = None

    def _plot_gds(self, cell, title="18-dot 量子点器件布局", show_plot=True):
        if not show_plot:
            return

        fig, ax = plt.subplots(figsize=(16, 12)); ax.set_title(title, fontsize=18)
        layer_config = {0:{'color':'#87CEEB','alpha':0.6,'label':'S/D (L0)','zorder':2}, 1:{'color':'#D3D3D3','alpha':0.5,'label':'SG (L1)','hatch':'///','zorder':1}, 2:{'color':'#FF1493','alpha':0.9,'label':'BG (L2)','zorder':3}, 3:{'color':'#8A2BE2','alpha':0.9,'label':'PG (L3)','zorder':4}}
        drawn_labels = set()
        for layer in [1, 0, 2, 3]: 
            polygons = cell.get_polygons(layer=layer, datatype=0)
            if not polygons: continue
            merged_polys = gdstk.boolean(polygons, [], 'or', layer=layer)
            cfg = layer_config.get(layer); label = cfg['label']
            for gds_poly in merged_polys:
                pts = gds_poly.points
                if pts.ndim != 2 or len(pts) < 3: continue
                current_label = label if label not in drawn_labels else None
                poly = MplPolygon(pts, closed=True, facecolor=cfg['color'], edgecolor='none', linewidth=0, alpha=cfg['alpha'], label=current_label, zorder=cfg['zorder'])
                ax.add_patch(poly); drawn_labels.add(label)
        handles, labels = ax.get_legend_handles_labels()
        if handles: by_label = dict(zip(labels, handles)); ax.legend(by_label.values(), by_label.keys(), loc='upper right')
        ax.set_aspect('equal')
        bbox = cell.bounding_box()
        if bbox is not None:
            min_x, min_y = bbox[0]; max_x, max_y = bbox[1]; width = max_x - min_x; height = max_y - min_y; margin_x = width * 0.1 if width > 0 else 1; margin_y = height * 0.1 if height > 0 else 1
            ax.set_xlim(min_x - margin_x, max_x + margin_x); ax.set_ylim(min_y - margin_y, max_y + margin_y)
        else: ax.autoscale_view()
        ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); plt.grid(True, which='both', linestyle='--', alpha=0.3); plt.tight_layout(); plt.show()

    def _create_18qd_layout_with_labels(self, **kwargs):
        params = {**self.qd_params, **kwargs}
        pg_max_width = params['pg_max_width']; pg_vert_side_len = params['pg_vert_side_len']; pg_chamfer_h = params['pg_chamfer_h']
        pg_bot_flat_w = params['pg_bot_flat_w']; pg_top_flat_w = params['pg_top_flat_w']; bg_max_width = params['bg_max_width']
        bg_vert_side_len = params['bg_vert_side_len']; bg_top_flat_w = params['bg_top_flat_w']; bg_bot_flat_w = params['bg_bot_flat_w']
        gap_pg_bg = params['gap_pg_bg']; d1_gap = params['d1_gap']; d2_gap = params['d2_gap']
        sd_gap_to_gate = params['sd_gap_to_gate']; gap_gate_outer_sg = params['gap_gate_outer_sg']; sg_mid_thick = params['sg_mid_thick']
        sg_top_thick = params['sg_top_thick']; sg_bottom_thick = params['sg_bot_thick']; sg_extension = params['sg_extension']
        lead_width = params['lead_width']; lead_length_bot = params['lead_length_bot']; lead_length_top = params['lead_length_top']
        lead_overlap = params['lead_overlap']; sd_height = params['sd_height']; sd_width = params['sd_width']
        sd_lead_width = params['sd_lead_width']; sd_outer_length = params['sd_outer_length']

        LAYER_SD=0;LAYER_SG=1;LAYER_BG=2;LAYER_PG=3;LAYER_LABEL=100;gate_half_height=pg_vert_side_len/2+pg_chamfer_h;row_pitch=gate_half_height+d1_gap+sg_mid_thick+d2_gap+gate_half_height;lib=gdstk.Library();cell=lib.new_cell('18QD_DEVICE_WITH_LEADS');connection_points={};
        
        def get_pg_points(c): cx,cy=c;yvt=pg_vert_side_len/2;yvb=-pg_vert_side_len/2;yt=yvt+pg_chamfer_h;yb=yvb-pg_chamfer_h;xm=pg_max_width/2;xtf=pg_top_flat_w/2;xbf=pg_bot_flat_w/2;return [(cx+xtf,cy+yt),(cx+xm,cy+yvt),(cx+xm,cy+yvb),(cx+xbf,cy+yb),(cx-xbf,cy+yb),(cx-xm,cy+yvb),(cx-xm,cy+yvt),(cx-xtf,cy+yt)],yt,yb
        def get_bg_points(c): cx,cy=c;yvt=bg_vert_side_len/2;yvb=-bg_vert_side_len/2;bcht=(bg_max_width-bg_top_flat_w)/2;bchb=(bg_max_width-bg_bot_flat_w)/2;yt=yvt+bcht;yb=yvb-bchb;xm=bg_max_width/2;xtf=bg_top_flat_w/2;xbf=bg_bot_flat_w/2;return [(cx+xtf,cy+yt),(cx+xm,cy+yvt),(cx+xm,cy+yvb),(cx+xbf,cy+yb),(cx-xbf,cy+yb),(cx-xm,cy+yvb),(cx-xm,cy+yvt),(cx-xtf,cy+yt)],yt,yb
        
        pglb=-(pg_vert_side_len/2+pg_chamfer_h);bglb=-(bg_vert_side_len/2+(bg_max_width-bg_bot_flat_w)/2);fbe=min(pglb,bglb)-lead_length_bot;yct=row_pitch;pglt=(pg_vert_side_len/2+pg_chamfer_h);bglt=(bg_vert_side_len/2+(bg_max_width-bg_top_flat_w)/2);fte=yct+max(pglt,bglt)+lead_length_top
        
        all_shapes=[];ycb=0;curr_x=0;bgn=[];[bgn.extend([f'QD_B{i}',f'QD_PG{i}'])for i in range(1,19)];bgn.append('QD_B19');bot_tips_y=[]
        pg_centers_x = []
        
        for i,n in enumerate(bgn):
            gt='PG' if 'PG' in n else 'BG';w=pg_max_width if gt=='PG' else bg_max_width;lyr=LAYER_PG if gt=='PG' else LAYER_BG;
            if i>0:pw=pg_max_width if 'PG' in bgn[i-1] else bg_max_width;curr_x+=(pw/2)+gap_pg_bg+(w/2)
            cen=(curr_x,ycb);pts,tip_y,by=get_pg_points(cen)if gt=='PG'else get_bg_points(cen);bot_tips_y.append(tip_y);poly=gdstk.Polygon(pts,layer=lyr);lys=by+lead_overlap;lye=fbe;lead=gdstk.rectangle((curr_x-lead_width/2,lys),(curr_x+lead_width/2,lye),layer=lyr)
            all_shapes.extend(gdstk.boolean([poly,lead],[],'or',layer=lyr));lbl=gdstk.Label(n,(curr_x,lye),layer=LAYER_LABEL);cell.add(lbl);connection_points[n]=lbl.origin
            if gt == 'PG': pg_centers_x.append(curr_x)

        bmx=-bg_max_width/2;bax=curr_x+bg_max_width/2;cxg=(bmx+bax)/2;
        set_centers_x = [pg_centers_x[1], pg_centers_x[4], pg_centers_x[7], pg_centers_x[10], pg_centers_x[13], pg_centers_x[16]]
        
        def create_set_at(cpgx, set_idx):
            shps=[];offs=[-(pg_max_width/2+gap_pg_bg+bg_max_width/2),0,(pg_max_width/2+gap_pg_bg+bg_max_width/2)];
            sgn=[f'SET{set_idx}_B1', f'SET{set_idx}_G', f'SET{set_idx}_B2']
            for off,n in zip(offs,sgn):
                cx=cpgx+off;gt='PG' if '_G' in n else 'BG';lyr=LAYER_PG if gt=='PG' else LAYER_BG;pts,_,_=get_pg_points((cx,yct))if gt=='PG'else get_bg_points((cx,yct));poly=gdstk.Polygon(pts,layer=lyr);Wty=np.max([pt[1]for pt in pts]);lys=Wty-lead_overlap;lye=fte;lead=gdstk.rectangle((cx-lead_width/2,lys),(cx+lead_width/2,lye),layer=lyr)
                shps.extend(gdstk.boolean([poly,lead],[],'or',layer=lyr));lbl=gdstk.Label(n,(cx,lye),layer=LAYER_LABEL);cell.add(lbl);connection_points[n]=lbl.origin
            return shps

        for i in range(1, 7):
            all_shapes.extend(create_set_at(set_centers_x[i-1], i))

        sgxs=bmx-sg_extension;sgxe=bax+sg_extension;ysgmb=max(bot_tips_y)+d1_gap;ysgmt=ysgmb+sg_mid_thick;
        sg_mid=gdstk.rectangle((sgxs,ysgmb),(sgxe,ysgmt),layer=LAYER_SG);ysgtb=yct+gate_half_height+gap_gate_outer_sg;ysgtt=ysgtb+sg_top_thick;
        sg_top=gdstk.rectangle((sgxs,ysgtb),(sgxe,ysgtt),layer=LAYER_SG);ysgbt=ycb-gate_half_height-gap_gate_outer_sg;ysgbb=ysgbt-sg_bottom_thick;
        sg_bot=gdstk.rectangle((sgxs,ysgbb),(sgxe,ysgbt),layer=LAYER_SG);all_shapes.extend([sg_top,sg_mid,sg_bot])
        
        lsg1=gdstk.Label("SG1",(sgxs,(ysgtb+ysgtt)/2),layer=LAYER_LABEL);cell.add(lsg1);connection_points["SG1"]=lsg1.origin;
        lsg2=gdstk.Label("SG2",(sgxs,(ysgmb+ysgmt)/2),layer=LAYER_LABEL);cell.add(lsg2);connection_points["SG2"]=lsg2.origin;
        lsg3=gdstk.Label("SG3",(sgxe,(ysgbb+ysgbt)/2),layer=LAYER_LABEL);cell.add(lsg3);connection_points["SG3"]=lsg3.origin

        sd_shps=[];sdcr=0.020
        def cr(xmin,ymin,xmax,ymax,crnrs):poly=gdstk.Polygon([(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)],layer=LAYER_SD);poly.fillet([sdcr if c else 0 for c in crnrs],tolerance=0.001);return poly
        
        blsix=bmx-sd_gap_to_gate;brsix=bax+sd_gap_to_gate;golx=blsix-sd_outer_length;gorx=brsix+sd_outer_length;
        sd_shps.append(cr(golx,ycb-sd_height/2,blsix,ycb+sd_height/2,[False,True,True,False]));lqdd=gdstk.Label("QD_D",(golx,ycb),layer=LAYER_LABEL);cell.add(lqdd);connection_points["QD_D"]=lqdd.origin;
        sd_shps.append(cr(brsix,ycb-sd_height/2,gorx,ycb+sd_height/2,[True,False,False,True]));lqds=gdstk.Label("QD_S",(gorx,ycb),layer=LAYER_LABEL);cell.add(lqds);connection_points["QD_S"]=lqds.origin;
        
        shw=bg_max_width+gap_pg_bg+pg_max_width/2;
        lsle=set_centers_x[0]-shw;tlsix=lsle-sd_gap_to_gate;
        sd_shps.append(cr(golx,yct-sd_height/2,tlsix,yct+sd_height/2,[False,True,True,False]));
        ls1d=gdstk.Label("SET1_D",(golx,yct),layer=LAYER_LABEL);cell.add(ls1d);connection_points["SET1_D"]=ls1d.origin;
        
        rsre=set_centers_x[5]+shw;trsix=rsre+sd_gap_to_gate;
        sd_shps.append(cr(trsix,yct-sd_height/2,gorx,yct+sd_height/2,[True,False,False,True]));
        ls6d=gdstk.Label("SET6_D",(gorx,yct),layer=LAYER_LABEL);cell.add(ls6d);connection_points["SET6_D"]=ls6d.origin;
        
        for i in range(1, 6):
            smxm = set_centers_x[i-1] + shw + sd_gap_to_gate
            smxa = set_centers_x[i] - shw - sd_gap_to_gate
            msdp = cr(smxm, yct-sd_height/2, smxa, yct+sd_height/2, [True,True,True,True])
            msdcx = (smxm+smxa)/2
            mlye = yct + lead_length_top
            msdl = gdstk.rectangle((msdcx-sd_lead_width/2, yct), (msdcx+sd_lead_width/2, mlye), layer=LAYER_SD)
            sd_shps.extend(gdstk.boolean([msdp, msdl], [], 'or', layer=LAYER_SD))
            label_name = f"SET{i}{i+1}_S" if i % 2 != 0 else f"SET{i}{i+1}_D"
            lbl = gdstk.Label(label_name, (msdcx, mlye), layer=LAYER_LABEL)
            cell.add(lbl); connection_points[label_name] = lbl.origin

        all_shapes.extend(sd_shps);cell.add(*all_shapes)
        return lib, cell, connection_points

    def _create_rect_wire_layout(self, **kwargs):
        params = {**self.pad_params, **kwargs}
        N = params['N']; layout_width = params['layout_width']; layout_height = params['layout_height']
        pad_width = params['pad_width']; pad_height = params['pad_height']; pad_spacing = params['pad_spacing']
        edge_margin = params['edge_margin']; active_size = params['active_size']; trace_width = params['trace_width']
        trace_spacing = params['trace_spacing']; taper_length = params['taper_length']; active_entry_len = params['active_entry_len']

        lib=gdstk.Library();cell=lib.new_cell('PAD_FRAME_WITH_WIRES');connection_points={};all_pads=[];
        def nsk(s):return[int(t)if t.isdigit()else t.lower()for t in re.split('([0-9]+)',s)]
        
        pgl=sorted([f"QD_PG{i+1}"for i in range(N)],key=nsk);bgl=sorted([f"QD_B{i+1}"for i in range(N+1)],key=nsk);gs=[];[gs.extend([bgl[i],pgl[i]])for i in range(N)];gs.append(bgl[N]);
        
        set_electrodes = []
        for i in range(1, 7):
            if i == 1: set_electrodes.append("SET1_D")
            set_electrodes.extend([f"SET{i}_B1", f"SET{i}_G", f"SET{i}_B2"])
            if i < 6:
                label_name = f"SET{i}{i+1}_S" if i % 2 != 0 else f"SET{i}{i+1}_D"
                set_electrodes.append(label_name)
            if i == 6: set_electrodes.append("SET6_D")
        
        pb = gs[10:27]
        lp = set_electrodes[1:4][::-1] + ["SG1"] + set_electrodes[:1] + ["SG2", "QD_D"] + gs[:10]
        rp = set_electrodes[22:25] + ["QD_S", "GND", "SG3", "GND"] + list(reversed(gs[27:]))
        tp = set_electrodes[4:22]
        
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
            else: [p.update({'target':(ax2,cy+offs[-(i+1)]),'target_edge':'right'})for i,p in enumerate(pads)]
        
        [at(pbe[edge],edge)for edge in['top','bottom','left','right']]
        # Assign different layer numbers to avoid overlap with QD layers (0, 1, 2, 3, 100)
        # Spread them out to help KLayout assign distinct default colors
        LS=30;LP=31;LT=31;LTXT=40;LA=50;LSIO2=60;sm=10
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
                if nw>ah:hl=ah/2;ym = ys+sy*hl;xm1=x1+sx*hl;xm2=x2-sx*hl;pts.extend([(xm1,ym),(xm2,ym)])
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

    def generate_quantum_dot_layout(self, **kwargs):
        self.device_lib, self.device_cell, self.device_connection_points = self._create_18qd_layout_with_labels(**kwargs)
        return self.device_lib, self.device_cell, self.device_connection_points

    def generate_pad_frame_layout(self, **kwargs):
        self.pad_lib, self.pad_cell, self.pad_connection_points, self.active_center, self.all_pads_info = self._create_rect_wire_layout(**kwargs)
        return self.pad_lib, self.pad_cell, self.pad_connection_points, self.active_center, self.all_pads_info

    def assemble_device_and_pad_frame(self):
        if self.device_cell is None or self.pad_cell is None:
            print("Error: Device cell or Pad cell not generated.")
            return None, None

        # Create a new library for the assembly with standard units
        assembly_lib = gdstk.Library(name="ASSEMBLY", unit=1e-6, precision=1e-9)
        # Use a very simple name for the top cell to ensure KLayout picks it up
        assembly_cell = assembly_lib.new_cell('18DOT_FULL_ASSEMBLY')
        
        # Brute-force flattening: Copy all geometry directly into the assembly cell
        # Use depth=None to get all geometry recursively
        assembly_cell.add(*self.pad_cell.get_polygons(depth=None))
        assembly_cell.add(*self.pad_cell.get_paths(depth=None))
        assembly_cell.add(*self.pad_cell.get_labels(depth=None))
            
        device_polys = self.device_cell.get_polygons(depth=None)
        device_paths = self.device_cell.get_paths(depth=None)
        device_labels = self.device_cell.get_labels(depth=None)
        
        bbox = self.device_cell.bounding_box()
        if bbox is not None:
            (min_x, min_y), (max_x, max_y) = bbox
            cx, cy = (min_x + max_x) / 2, (min_y + max_y) / 2
            tx, ty = self.active_center[0] - cx, self.active_center[1] - cy
            
            for p in device_polys: p.translate(tx, ty)
            for p in device_paths: p.translate(tx, ty)
            for l in device_labels: l.origin = (l.origin[0] + tx, l.origin[1] + ty)
            
        assembly_cell.add(*device_polys)
        assembly_cell.add(*device_paths)
        assembly_cell.add(*device_labels)
        
        return assembly_lib, assembly_cell

    def visualize_and_save_pad_frame(self, lib, cell_to_show, title, gds_filename, show_plot=True):
        print(f"--- {title} ---"); print(f"正在可视化: {title}...")
        if show_plot:
            fig, ax = plt.subplots(figsize=(12, 12), dpi=100)
            # Use depth=None to get all geometry
            polys = cell_to_show.get_polygons(depth=None)
            for p in polys:
                if p.layer == 31: color, alpha, zorder = '#f1c40f', 0.8, 5 # Gold
                elif p.layer == 50: color, alpha, zorder = '#1abc9c', 0.3, 2 # Turquoise
                elif p.layer == 60: color, alpha, zorder = '#e67e22', 0.2, 1 # Orange
                elif p.layer == 30: color, alpha, zorder = '#ecf0f1', 0.2, 0 # Light Gray
                else: continue
                ax.add_patch(MplPolygon(p.points, closed=True, facecolor=color, alpha=alpha, edgecolor='none', zorder=zorder))
            
            # Add labels
            labels = cell_to_show.get_labels(depth=None)
            for lbl in labels:
                ax.text(lbl.origin[0], lbl.origin[1], lbl.text, fontsize=9, ha='center', va='center', color='blue', fontweight='bold', alpha=0.8)
                
            ax.set_aspect('equal')
            ax.set_title(title, fontsize=14)
            bbox = cell_to_show.bounding_box()
            if bbox:
                (min_x, min_y), (max_x, max_y) = bbox
                margin = (max_x - min_x) * 0.05
                ax.set_xlim(min_x - margin, max_x + margin)
                ax.set_ylim(min_y - margin, max_y + margin)
            plt.tight_layout()
            plt.show()
            
        lib.write_gds(gds_filename); print(f"GDS 文件已保存到: '{gds_filename}'\n")

    def visualize_and_save_assembly(self, lib, cell_to_show, title, gds_filename, show_plot=True):
        print(f"--- {title} ---"); print(f"正在可视化: {title}...")
        
        if os.path.exists(gds_filename): os.remove(gds_filename)
        lib.write_gds(gds_filename)
        print(f"GDS 文件已保存到: '{gds_filename}'")

        if not show_plot:
            return

        # Revised layer config with more distinct colors and updated layer numbers
        layer_config = {
            0:  {'color': '#3498db', 'alpha': 0.8, 'label': 'S/D (L0)', 'zorder': 10},
            1:  {'color': '#bdc3c7', 'alpha': 0.3, 'label': 'SG (L1)', 'zorder': 1},
            2:  {'color': '#e74c3c', 'alpha': 0.8, 'label': 'BG (L2)', 'zorder': 11},
            3:  {'color': '#9b59b6', 'alpha': 0.8, 'label': 'PG (L3)', 'zorder': 12},
            5:  {'color': '#2ecc71', 'alpha': 1.0, 'label': 'Routing (L5)', 'zorder': 13},
            30: {'color': '#ecf0f1', 'alpha': 0.2, 'label': 'Substrate (L30)', 'zorder': 0},
            31: {'color': '#f1c40f', 'alpha': 0.7, 'label': 'Pads (L31)', 'zorder': 5},
            40: {'color': '#e67e22', 'alpha': 0.8, 'label': 'Labels (L40)', 'zorder': 15},
            50: {'color': '#1abc9c', 'alpha': 0.4, 'label': 'Active (L50)', 'zorder': 2},
            60: {'color': '#e67e22', 'alpha': 0.2, 'label': 'SiO2 (L60)', 'zorder': 1},
        }

        def draw_on_ax(ax_obj, is_zoomed=False):
            drawn_labels = set()
            for layer in sorted(layer_config.keys(), key=lambda x: layer_config[x]['zorder']):
                cfg = layer_config[layer]
                # Fix: Filtering requires both layer and datatype/texttype to be set in some gdstk versions
                polys = cell_to_show.get_polygons(depth=None, layer=layer, datatype=0)
                if not polys: continue
                merged_polys = gdstk.boolean(polys, [], 'or', layer=layer)
                for gds_poly in merged_polys:
                    label = cfg['label'] if cfg['label'] not in drawn_labels else None
                    is_fill = True # All solid fills for clarity
                    ax_obj.add_patch(MplPolygon(gds_poly.points, closed=True, fill=is_fill, 
                                              facecolor=cfg['color'], edgecolor='none', 
                                              alpha=cfg['alpha'], label=label, zorder=cfg['zorder']))
                    if label: drawn_labels.add(label)
            
            labels = cell_to_show.get_labels(depth=None)
            xmin, xmax = ax_obj.get_xlim(); ymin, ymax = ax_obj.get_ylim()
            for lbl in labels:
                if xmin <= lbl.origin[0] <= xmax and ymin <= lbl.origin[1] <= ymax:
                    if is_zoomed or len(lbl.text) > 4:
                        fs = 10 if is_zoomed else 8
                        ax_obj.text(lbl.origin[0], lbl.origin[1], lbl.text, fontsize=fs, 
                                  ha='center', va='center', color='blue', fontweight='bold', alpha=0.8)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), dpi=150)
        bbox = cell_to_show.bounding_box()
        if bbox:
            (min_x, min_y), (max_x, max_y) = bbox; margin = (max_x - min_x) * 0.05
            ax1.set_xlim(min_x - margin, max_x + margin); ax1.set_ylim(min_y - margin, max_y + margin)
        
        if self.active_center:
            acx, acy = self.active_center; zs = self.pad_params.get('active_size', 500) * 0.6
            ax2.set_xlim(acx - zs, acx + zs); ax2.set_ylim(acy - zs, acy + zs)

        draw_on_ax(ax1, is_zoomed=False); ax1.set_aspect('equal'); ax1.set_title(f"{title} (Full View)")
        draw_on_ax(ax2, is_zoomed=True); ax2.set_aspect('equal'); ax2.set_title(f"{title} (Zoomed Center)")
        plt.tight_layout(); plt.show()
