import gdsfactory as gf
from gdsfactory.component import Component
from gdsfactory.port import Port
import numpy as np
from typing import List, Tuple

import gdsfactory as gf
from gdsfactory.component import Component
from gdsfactory.port import Port


class QuantumDotGenerator:
    def __init__(self):
        self.LAYER_SD = 0
        self.LAYER_SG = 1
        self.LAYER_BG = 2
        self.LAYER_PG = 3
        self.LAYER_LABEL = 100

    def create_6qd_layout_with_labels(
        self,
        pg_max_width=0.120, pg_vert_side_len=0.040, pg_chamfer_h=0.040, pg_bot_flat_w=0.042, pg_top_flat_w=0.040,
        bg_max_width=0.060, bg_vert_side_len=0.040, bg_top_flat_w=0.020, bg_bot_flat_w=0.042,
        gap_pg_bg=0.002, d1_gap=0.020, d2_gap=0.020, sd_gap_to_gate=0.01, gap_gate_outer_sg=0.030,
        sg_mid_thick=0.100, sg_top_thick=0.300, sg_bot_thick=0.300, sg_extension=0.200,
        lead_width=0.042, lead_length_bot=0.4, lead_length_top=0.5, lead_overlap=0.025,
        sd_height=0.10, sd_width=0.35, sd_lead_width=0.060, sd_outer_length=0.300
    ) -> Tuple[Component, dict]:
        gate_half_height = pg_vert_side_len / 2 + pg_chamfer_h
        row_pitch = gate_half_height + d1_gap + sg_mid_thick + d2_gap + gate_half_height
        c = gf.Component('6QD_DEVICE_WITH_LEADS')
        connection_points = {}

        def get_pg_polygon(c_center, layer):
            cx, cy = c_center
            yvt = pg_vert_side_len / 2
            yvb = -pg_vert_side_len / 2
            yt = yvt + pg_chamfer_h
            yb = yvb - pg_chamfer_h
            xm = pg_max_width / 2
            xtf = pg_top_flat_w / 2
            xbf = pg_bot_flat_w / 2
            pts = [(cx + xtf, cy + yt), (cx + xm, cy + yvt), (cx + xm, cy + yvb), (cx + xbf, cy + yb),
                   (cx - xbf, cy + yb), (cx - xm, cy + yvb), (cx - xm, cy + yvt), (cx - xtf, cy + yt)]
            # Create a component from the polygon points
            poly_comp = gf.Component()
            poly_comp.add_polygon(pts, layer=layer)
            return poly_comp, yt, yb

        def get_bg_polygon(c_center, layer):
            cx, cy = c_center
            yvt = bg_vert_side_len / 2
            yvb = -bg_vert_side_len / 2
            bcht = (bg_max_width - bg_top_flat_w) / 2
            bchb = (bg_max_width - bg_bot_flat_w) / 2
            yt = yvt + bcht
            yb = yvb - bchb
            xm = bg_max_width / 2
            xtf = bg_top_flat_w / 2
            xbf = bg_bot_flat_w / 2
            pts = [(cx + xtf, cy + yt), (cx + xm, cy + yvt), (cx + xm, cy + yvb), (cx + xbf, cy + yb),
                   (cx - xbf, cy + yb), (cx - xm, cy + yvb), (cx - xm, cy + yvt), (cx - xtf, cy + yt)]
            # Create a component from the polygon points
            poly_comp = gf.Component()
            poly_comp.add_polygon(pts, layer=layer)
            return poly_comp, yt, yb

        pglb = -(pg_vert_side_len / 2 + pg_chamfer_h)
        bglb = -(bg_vert_side_len / 2 + (bg_max_width - bg_bot_flat_w) / 2)
        fbe = min(pglb, bglb) - lead_length_bot
        yct = row_pitch
        pglt = (pg_vert_side_len / 2 + pg_chamfer_h)
        bglt = (bg_vert_side_len / 2 + (bg_max_width - bg_top_flat_w) / 2)
        fte = yct + max(pglt, bglt) + lead_length_top

        ycb = 0
        curr_x = 0
        bgn = []
        [bgn.extend([f'QD_B{i}', f'QD_PG{i}']) for i in range(1, 7)]
        bgn.append('QD_B7')
        bot_tips_y = []

        for i, n in enumerate(bgn):
            gt = 'PG' if 'PG' in n else 'BG'
            w = pg_max_width if gt == 'PG' else bg_max_width
            lyr = self.LAYER_PG if gt == 'PG' else self.LAYER_BG
            if i > 0:
                pw = pg_max_width if 'PG' in bgn[i - 1] else bg_max_width
                curr_x += (pw / 2) + gap_pg_bg + (w / 2)
            cen = (curr_x, ycb)
            poly_comp, tip_y, by = get_pg_polygon(cen, lyr) if gt == 'PG' else get_bg_polygon(cen, lyr)
            bot_tips_y.append(tip_y)
            
            lys = by + lead_overlap
            lye = fbe
            lead_comp = gf.components.rectangle(size=(lead_width, abs(lye - lys)), layer=lyr)
            lead_ref = c.add_ref(lead_comp)
            lead_ref.center = (curr_x, (lys + lye) / 2)
            
            # Boolean operation to combine polygon and lead
            # gf.boolean operates on components, so we need to ensure poly_comp is a component.
            # Since get_pg_polygon/get_bg_polygon now return components, this should be fine.
            # We need to add poly_comp as a reference to c before boolean operation if it's not already.
            poly_ref = c.add_ref(poly_comp)
            combined_shape = gf.boolean(poly_ref, lead_ref, 'or', layer=lyr)
            c.add_ref(combined_shape)

            # Add label and port
            c.add_label(text=n, position=(curr_x, lye), layer=(self.LAYER_LABEL, 0)) # gdsfactory labels use (layer, datatype)
            connection_points[n] = gf.Port(name=n, center=(curr_x, lye), width=lead_width, orientation=270, layer=lyr) # Assuming South orientation for bottom leads

        bmx = -bg_max_width / 2
        bax = curr_x + bg_max_width / 2
        cxg = (bmx + bax) / 2
        igr = 2 * sd_gap_to_gate + sd_width
        off_cen = igr / 2 + bg_max_width / 2 + gap_pg_bg + pg_max_width / 2
        lscpgx = cxg - off_cen
        rscpgx = cxg + off_cen

        def create_set_at(cpgx, side):
            offs = [-(pg_max_width / 2 + gap_pg_bg + bg_max_width / 2), 0, (pg_max_width / 2 + gap_pg_bg + bg_max_width / 2)]
            sgn = ['SET1_B1', 'SET1_G', 'SET1_B2'] if side == 'L' else ['SET2_B1', 'SET2_G', 'SET2_B2']
            for off, n in zip(offs, sgn):
                cx = cpgx + off
                gt = 'PG' if '_G' in n else 'BG'
                lyr = self.LAYER_PG if gt == 'PG' else self.LAYER_BG
                poly_comp, _, _ = get_pg_polygon((cx, yct), lyr) if gt == 'PG' else get_bg_polygon((cx, yct), lyr)
                poly_ref = c.add_ref(poly_comp)
                Wty = poly_ref.ymax # Access ymax from the component reference
                lys = Wty - lead_overlap
                lye = fte
                lead_comp = gf.components.rectangle(size=(lead_width, abs(lye - lys)), layer=lyr)
                lead_ref = c.add_ref(lead_comp)
                lead_ref.center = (cx, (lys + lye) / 2)
                
                combined_shape = gf.boolean(poly_ref, lead_ref, 'or', layer=lyr)
                c.add_ref(combined_shape)

                c.add_label(text=n, position=(cx, lye), layer=(self.LAYER_LABEL, 0))
                connection_points[n] = gf.Port(name=n, center=(cx, lye), width=lead_width, orientation=90, layer=lyr) # Assuming North orientation for top leads

        create_set_at(lscpgx, side='L')
        create_set_at(rscpgx, side='R')

        sgxs = bmx - sg_extension
        sgxe = bax + sg_extension
        ysgmb = max(bot_tips_y) + d1_gap
        ysgmt = ysgmb + sg_mid_thick
        sg_mid_comp = gf.components.rectangle(size=(sgxe - sgxs, ysgmt - ysgmb), layer=self.LAYER_SG)
        sg_mid_ref = c.add_ref(sg_mid_comp)
        sg_mid_ref.center = ((sgxs + sgxe) / 2, (ysgmb + ysgmt) / 2)
        ysgtb = yct + gate_half_height + gap_gate_outer_sg
        ysgtt = ysgtb + sg_top_thick
        sg_top_comp = gf.components.rectangle(size=(sgxe - sgxs, ysgtt - ysgtb), layer=self.LAYER_SG)
        sg_top_ref = c.add_ref(sg_top_comp)
        sg_top_ref.center = ((sgxs + sgxe) / 2, (ysgtb + ysgtt) / 2)
        ysgbt = ycb - gate_half_height - gap_gate_outer_sg
        ysgbb = ysgbt - sg_bot_thick
        sg_bot_comp = gf.components.rectangle(size=(sgxe - sgxs, ysgbt - ysgbb), layer=self.LAYER_SG)
        sg_bot_ref = c.add_ref(sg_bot_comp)
        sg_bot_ref.center = ((sgxs + sgxe) / 2, (ysgbb + ysgbt) / 2)

        c.add_label(text="SG1", position=(sgxs, (ysgtb + ysgtt) / 2), layer=(self.LAYER_LABEL, 0))
        connection_points["SG1"] = gf.Port(name="SG1", center=(sgxs, (ysgtb + ysgtt) / 2), width=sg_top_thick, orientation=180, layer=self.LAYER_SG)
        c.add_label(text="SG2", position=(sgxs, (ysgmb + ysgmt) / 2), layer=(self.LAYER_LABEL, 0))
        connection_points["SG2"] = gf.Port(name="SG2", center=(sgxs, (ysgmb + ysgmt) / 2), width=sg_mid_thick, orientation=180, layer=self.LAYER_SG)
        c.add_label(text="SG3", position=(sgxs, (ysgbb + ysgbt) / 2), layer=(self.LAYER_LABEL, 0))
        connection_points["SG3"] = gf.Port(name="SG3", center=(sgxs, (ysgbb + ysgbt) / 2), width=sg_bot_thick, orientation=180, layer=self.LAYER_SG)

        sdcr = 0.020

        def cr_gf(xmin, ymin, xmax, ymax, crnrs, layer):
            # gdsfactory does not have a direct fillet function for polygons like gdstk
            # For rounded rectangles, gf.components.rounded() can be used, but it's for a full component.
            # For individual polygons, we might need to construct it from paths or use boolean operations.
            # For now, let's create a simple rectangle and add it.
            # If specific rounding is critical, a custom component would be needed.
            rect = gf.components.rectangle(size=(xmax - xmin, ymax - ymin), layer=layer)
            rect_ref = gf.Component()
            rect_ref.add_ref(rect).center = ((xmin + xmax) / 2, (ymin + ymax) / 2)
            # Applying fillets to individual polygon points is not straightforward with gdsfactory's high-level API.
            # For this conversion, we will approximate with a simple rectangle.
            return rect_ref

        blbe = bmx
        blsix = blbe - sd_gap_to_gate
        brbe = bax
        brsix = brbe + sd_gap_to_gate
        golx = blsix - sd_outer_length
        gorx = brsix + sd_outer_length
        
        # QD_D
        qd_d_comp = cr_gf(golx, ycb - sd_height / 2, blsix, ycb + sd_height / 2, [False, True, True, False], self.LAYER_SD)
        qd_d_ref = c.add_ref(qd_d_comp)
        c.add_label(text="QD_D", position=(golx, ycb), layer=(self.LAYER_LABEL, 0))
        connection_points["QD_D"] = gf.Port(name="QD_D", center=(golx, ycb), width=sd_height, orientation=180, layer=self.LAYER_SD) # West

        # QD_S
        qd_s_comp = cr_gf(brsix, ycb - sd_height / 2, gorx, ycb + sd_height / 2, [True, False, False, True], self.LAYER_SD)
        qd_s_ref = c.add_ref(qd_s_comp)
        c.add_label(text="QD_S", position=(gorx, ycb), layer=(self.LAYER_LABEL, 0))
        connection_points["QD_S"] = gf.Port(name="QD_S", center=(gorx, ycb), width=sd_height, orientation=0, layer=self.LAYER_SD) # East

        shw = bg_max_width + gap_pg_bg + pg_max_width / 2
        lsle = lscpgx - shw
        tlsix = lsle - sd_gap_to_gate
        rsre = rscpgx + shw
        trsix = rsre - sd_gap_to_gate
        
        # SET1_D
        set1_d_comp = cr_gf(golx, yct - sd_height / 2, tlsix, yct + sd_height / 2, [False, True, True, False], self.LAYER_SD)
        set1_d_ref = c.add_ref(set1_d_comp)
        c.add_label(text="SET1_D", position=(golx, yct), layer=(self.LAYER_LABEL, 0))
        connection_points["SET1_D"] = gf.Port(name="SET1_D", center=(golx, yct), width=sd_height, orientation=180, layer=self.LAYER_SD) # West

        # SET2_D
        set2_d_comp = cr_gf(trsix, yct - sd_height / 2, gorx, yct + sd_height / 2, [True, False, False, True], self.LAYER_SD)
        set2_d_ref = c.add_ref(set2_d_comp)
        c.add_label(text="SET2_D", position=(gorx, yct), layer=(self.LAYER_LABEL, 0))
        connection_points["SET2_D"] = gf.Port(name="SET2_D", center=(gorx, yct), width=sd_height, orientation=0, layer=self.LAYER_SD) # East

        lsre = lscpgx + shw
        smxm = lsre + sd_gap_to_gate
        rsle = rscpgx - shw
        smxa = rsle - sd_gap_to_gate
        
        # SET_S
        msdp_comp = cr_gf(smxm, yct - sd_height / 2, smxa, yct + sd_height / 2, [True, True, True, True], self.LAYER_SD)
        msdcx = (smxm + smxa) / 2
        mlye = yct + lead_length_top
        msdl_comp = gf.components.rectangle(size=(sd_lead_width, abs(mlye - yct)), layer=self.LAYER_SD)
        msdl_ref = c.add_ref(msdl_comp)
        msdl_ref.center = (msdcx, (yct + mlye) / 2)
        
        msdp_ref = c.add_ref(msdp_comp)
        combined_sd_shape = gf.boolean(msdp_ref, msdl_ref, 'or', layer=self.LAYER_SD)
        c.add_ref(combined_sd_shape)

        c.add_label(text="SET_S", position=(msdcx, mlye), layer=(self.LAYER_LABEL, 0))
        connection_points["SET_S"] = gf.Port(name="SET_S", center=(msdcx, mlye), width=sd_lead_width, orientation=90, layer=self.LAYER_SD) # North
        
        return c, connection_points
