import numpy as np
import math
import re
from typing import List, Tuple, Dict

import gdsfactory as gf
from gdsfactory.component import Component
from gdsfactory.port import Port

class PadLeadGenerator:
    def __init__(self):
        self.LAYER_LS = 0 # Layout Structure
        self.LAYER_LP = 1 # Pads
        self.LAYER_LT = 1 # Traces (same as pads for now, will be routed on ROUTING_LAYER later)
        self.LAYER_LTXT = 2 # Labels
        self.LAYER_LA = 10 # Active Area
        self.LAYER_LSIO2 = 11 # SiO2

    def create_rect_wire_layout(
        self,
        N: int = 6, layout_width: float = 1400, layout_height: float = 1400, pad_width: float = 100, pad_height: float = 100, pad_spacing: float = 30,
        edge_margin: float = 50, active_size: float = 180, trace_width: float = 10, trace_spacing: float = 10, taper_length: float = 40,
        active_entry_len: float = 40
    ) -> Tuple[Component, Dict[str, Port], Tuple[float, float], List[Dict]]:
        c = gf.Component('PAD_FRAME_WITH_WIRES')
        connection_points = {}
        all_pads = []

        def nsk(s):
            return [int(t) if t.isdigit() else t.lower() for t in re.split('([0-9]+)', s)]

        pgl = sorted([f"QD_PG{i + 1}" for i in range(N)], key=nsk)
        bgl = sorted([f"QD_B{i + 1}" for i in range(N + 1)], key=nsk)
        gs = []
        [gs.extend([bgl[i], pgl[i]]) for i in range(N)]
        gs.append(bgl[N])
        tp = []
        ns = math.ceil(N / 3)

        if ns >= 1:
            tp.extend(["SET1_B1", "SET1_B2", "SET1_G"])
            tp.append("SET_S")
        if ns >= 2:
            tp.extend(["SET2_G", "SET2_B2", "SET2_B1"])

        lf = ["SG1"]
        if ns >= 1:
            lf.append("SET1_D")
            lf.extend(["SG2", "QD_D"])
        rf = []
        if ns >= 2:
            rf.append("SET2_D")
            rf.append("QD_S")
            rf.extend(["GND", "SG3", "GND"])

        ls = len(gs)
        ll = len(lf)
        lr = len(rf)
        lt = len(tp)
        aw = layout_width - 2 * edge_margin
        bc = (0, ls, 0)
        md = float('inf')
        ss = max(ll, lr)

        for ts in range(ss, ss + ls):
            al = ts - ll
            ar = ts - lr
            if al + ar > ls:
                break
            nb = ls - (al + ar)
            if nb > 0 and (nb * pad_width + (nb - 1) * pad_spacing > aw):
                continue
            d = abs(nb - lt)
            if d < md:
                md = d
                bc = (al, nb, ar)

        nl, nb, nr = bc
        ple = gs[:nl]
        pb = gs[nl:nl + nb]
        pre = gs[nl + nb:]
        lp = lf + ple
        rp = rf + list(reversed(pre))

        def ac(pads, edge):
            num_pads = len(pads)
            if num_pads == 0:
                return
            if edge in ['top', 'bottom']:
                total_span = num_pads * pad_width + (num_pads - 1) * pad_spacing
                start_x = (layout_width - total_span) / 2
                y_coord = layout_height - edge_margin - pad_height if edge == 'top' else edge_margin
                for i, lbl in enumerate(pads):
                    x_coord = start_x + i * (pad_width + pad_spacing)
                    center_x, center_y = x_coord + pad_width / 2, y_coord + pad_height / 2
                    all_pads.append({'label': lbl, 'rect': (x_coord, y_coord, pad_width, pad_height), 'center': (center_x, center_y), 'edge': edge})
            else: # 'left', 'right'
                total_span = num_pads * pad_height + (num_pads - 1) * pad_spacing
                start_y = (layout_height + total_span) / 2
                x_coord = edge_margin if edge == 'left' else layout_width - edge_margin - pad_width
                for i, lbl in enumerate(pads):
                    y_coord = start_y - (i + 1) * pad_height - i * pad_spacing
                    center_x, center_y = x_coord + pad_height / 2, y_coord + pad_height / 2
                    all_pads.append({'label': lbl, 'rect': (x_coord, y_coord, pad_height, pad_width), 'center': (center_x, center_y), 'edge': edge})

        ac(tp, 'top')
        ac(pb, 'bottom')
        ac(lp, 'left')
        ac(rp, 'right')

        center_x = layout_width / 2
        center_y = layout_height / 2
        acenter = (center_x, center_y)
        ax1 = center_x - active_size / 2
        ax2 = center_x + active_size / 2
        ay1 = center_y - active_size / 2
        ay2 = center_y + active_size / 2
        pbe = {'top': [], 'bottom': [], 'left': [], 'right': []}
        [pbe[p['edge']].append(p) for p in all_pads]

        pbe['top'].sort(key=lambda p: p['center'][0])
        pbe['bottom'].sort(key=lambda p: p['center'][0])
        pbe['left'].sort(key=lambda p: p['center'][1], reverse=True)
        pbe['right'].sort(key=lambda p: p['center'][1], reverse=True)
        pitch = trace_width + trace_spacing

        def at(pads, edge):
            num_pads = len(pads)
            if num_pads == 0:
                return
            indices = np.arange(num_pads) - (num_pads - 1) / 2
            offsets = indices * pitch
            if edge == 'top':
                [p.update({'target': (center_x + offsets[i], ay2), 'target_edge': 'top'}) for i, p in enumerate(pads)]
            elif edge == 'bottom':
                [p.update({'target': (center_x + offsets[i], ay1), 'target_edge': 'bottom'}) for i, p in enumerate(pads)]
            elif edge == 'left':
                [p.update({'target': (ax1, center_y + offsets[-(i + 1)]), 'target_edge': 'left'}) for i, p in enumerate(pads)]
            elif edge == 'right':
                [p.update({'target': (ax2, center_y + offsets[-(i + 1)]), 'target_edge': 'right'}) for i, p in enumerate(pads)]

        [at(pbe[edge], edge) for edge in ['top', 'bottom', 'left', 'right']]

        sm = 10
        layout_rect = gf.components.rectangle(size=(layout_width, layout_height), layer=self.LAYER_LS)
        layout_ref = c.add_ref(layout_rect)
        layout_ref.center = (layout_width/2, layout_height/2)

        active_rect = gf.components.rectangle(size=(active_size, active_size), layer=self.LAYER_LA)
        active_ref = c.add_ref(active_rect)
        active_ref.center = (center_x, center_y)

        def ct_gf(pr, edge, tl, tw, layer):
            x, y, w, h = pr
            cx, cy = x + w / 2, y + h / 2
            if edge == 'top':
                pts = [(x, y), (x + w, y), (cx + tw / 2, y - tl), (cx - tw / 2, y - tl)]
                start = (cx, y - tl)
                orientation = 270 # South
            elif edge == 'bottom':
                pts = [(x, y + h), (x + w, y + h), (cx + tw / 2, y + h + tl), (cx - tw / 2, y + h + tl)]
                start = (cx, y + h + tl)
                orientation = 90 # North
            elif edge == 'left':
                pts = [(x + w, y), (x + w, y + h), (x + w + tl, cy + tw / 2), (x + w + tl, cy - tw / 2)]
                start = (x + w + tl, cy)
                orientation = 0 # East
            else: # right
                pts = [(x, y), (x, y + h), (x - tl, cy + tw / 2), (x - tl, cy - tw / 2)]
                start = (x - tl, cy)
                orientation = 180 # West
            poly_comp = gf.Component()
            poly_comp.add_polygon(pts, layer=layer)
            return poly_comp, start, orientation

        def rsz_gf(sp, tp, edge, width, layer):
            x1, y1 = sp
            x2, y2 = tp
            sm = 20
            path_pts = [(x1, y1)]
            if edge in ['top', 'bottom']:
                sy = -1 if edge == 'top' else 1
                ys = y1 + sy * sm
                path_pts.append((x1, ys))
                ye = y2 + sy * active_entry_len
                ah = abs(ye - ys)
                nw = abs(x2 - x1)
                sx = 1 if x2 > x1 else -1
                if nw > ah:
                    hl = ah / 2
                    ym = ys + sy * hl
                    xm1 = x1 + sx * hl
                    xm2 = x2 - sx * hl
                    path_pts.extend([(xm1, ym), (xm2, ym)])
                else:
                    dyd = nw
                    yte = ys + sy * dyd
                    path_pts.append((x2, yte))
                path_pts.append((x2, y2))
            else: # 'left', 'right'
                sx = 1 if edge == 'left' else -1
                xs = x1 + sx * sm
                path_pts.append((xs, y1))
                xe = x2 + sx * active_entry_len
                aw = abs(xe - xs)
                nh = abs(y2 - y1)
                sy = 1 if y2 > y1 else -1
                if nh > aw:
                    wl = aw / 2
                    xm = xs + sx * wl
                    ym1 = y1 + sy * wl
                    ym2 = y2 - sy * wl
                    path_pts.extend([(xm, ym1), (xm, ym2)])
                else:
                    dxd = nh
                    xte = xs + sx * dxd
                    path_pts.append((xte, y2))
                path_pts.append((x2, y2))
            
            # Create a Path object and extrude it
            path = gf.path.Path(path_pts)
            return path.extrude(width=width, layer=layer), path_pts

        for p in all_pads:
            x, y, w, h = p['rect']
            lbl = p['label']
            edge = p['edge']
            
            # SiO2 rectangle
            sio2_rect_comp = gf.components.rectangle(size=(w + 2 * sm, h + 2 * sm), layer=self.LAYER_LSIO2)
            sio2_rect_ref = c.add_ref(sio2_rect_comp)
            sio2_rect_ref.center = (x + w / 2, y + h / 2)

            # Pad rectangle
            pad_rect_comp = gf.components.rectangle(size=(w, h), layer=self.LAYER_LP)
            pad_rect_ref = c.add_ref(pad_rect_comp)
            pad_rect_ref.center = (x + w / 2, y + h / 2)

            # Label
            cx, cy = x + w / 2, y + h / 2
            # gdsfactory's add_label does not directly support magnification or rotation as keyword arguments.
            # Rotation can be applied to the component reference if needed, but for a simple label, it's usually not.
            c.add_label(text=lbl, position=(cx, cy), layer=(self.LAYER_LTXT, 0))

            # Taper
            poly_t_comp, t_start_coords, t_start_orientation = ct_gf(p['rect'], edge, taper_length, trace_width, self.LAYER_LT)
            c.add_ref(poly_t_comp)
            
            # Routing path
            path_comp, path_pts = rsz_gf(t_start_coords, p['target'], edge, trace_width, self.LAYER_LT)
            c.add_ref(path_comp)
            p['trace_points'] = path_pts
            
            # Create a Port at the target for routing
            target_x, target_y = p['target']
            target_orientation = 0 # Default to East, will refine based on target_edge
            if p['target_edge'] == 'top':
                target_orientation = 90
            elif p['target_edge'] == 'bottom':
                target_orientation = 270
            elif p['target_edge'] == 'left':
                target_orientation = 180
            elif p['target_edge'] == 'right':
                target_orientation = 0

            connection_points[p['label']] = gf.Port(
                name=p['label'],
                center=(target_x, target_y),
                width=trace_width,
                orientation=target_orientation,
                layer=self.LAYER_LT
            )
        return c, connection_points, acenter, all_pads
