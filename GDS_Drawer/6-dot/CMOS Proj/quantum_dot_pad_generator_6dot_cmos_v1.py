import os
import re
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

# Reuse the 6-dot CMOS device generator (device core + labels/connection points)
from six_dot_generator_SG1_CMOS_20260312 import create_6qd_device_core, plot_gds


# ===================================================================
# 解决中文显示问题
# ===================================================================
plt.rcParams["font.sans-serif"] = ["SimHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False


Point = Tuple[float, float]


class QuantumDotPadGenerator6DotCMOS:
    """
    6-dot CMOS: device(core with leads) + pad frame generator.

    Design goals (aligned with Prefab-18-dot 20260109 v3 workflow):
    - Generate device GDS (6-dot CMOS core) with labels/connection_points
    - Generate pad_frame GDS using the same style as 18dot_v3: pads on 4 edges + 4 active regions + multi-stage inner leads
    - Assemble (device + pad_frame) into a single GDS by translating device to active_center
    - Do NOT auto-route between device and pads in this script (routing is left for the notebook / manual wiring)
    """

    def __init__(self):
        # ---------------------------
        # Pad-frame default parameters (copied from 18dot_v3)
        # ---------------------------
        self.pad_params = {
            "layout_width": 2500,
            "layout_height": 2500,
            "pad_width": 100,
            "pad_height": 100,
            "pad_spacing": 30,
            "edge_margin": 100,
            
            "active_size": 400,
            "active_size_2": 150,
            "active_size_3": 60,
            "active_width_4": 30,
            "active_height_4": 15,
            
            "trace_width": 10,
            "trace_spacing": 10,
            "taper_length": 40,

            "active_entry_len": 40,
            "active_entry_len_2": 10,
            "active_entry_len_3": 5,
            "active_entry_len_4": 2,

            "inner_trace_width": 1.0,
            "inner_trace_spacing": 2.0,
            "inner_trace_width_3": 0.5,
            "inner_trace_spacing_3": 1.0,
            "inner_trace_width_4": 0.1,
            "inner_trace_spacing_4": 0.5,
        }

        # Runtime generated artifacts
        self.device_lib = None
        self.device_cell = None
        self.device_connection_points: Optional[Dict[str, Point]] = None

        self.pad_lib = None
        self.pad_cell = None
        self.pad_connection_points: Optional[Dict[str, Point]] = None

        self.active_center: Optional[Point] = None
        self.all_pads_info = None

    # ============================================================
    # Visualization helpers (optional)
    # ============================================================
    def plot_pad_frame(self, cell_to_show, title="Pad Frame", show_plot=True):
        if not show_plot:
            return

        fig, ax = plt.subplots(figsize=(12, 12), dpi=120)
        ax.set_title(title)

        # Match 18dot_v3 visualization colors for pad frame
        polys = cell_to_show.get_polygons(depth=None)
        for p in polys:
            if p.layer == 31:
                color, alpha, zorder = "#f1c40f", 0.8, 5  # Pads
            elif p.layer == 50:
                color, alpha, zorder = "#1abc9c", 0.3, 2  # Active 1
            elif p.layer == 51:
                color, alpha, zorder = "#3498db", 0.3, 3  # Active 2
            elif p.layer == 52:
                color, alpha, zorder = "#9b59b6", 0.3, 4  # Active 3
            elif p.layer == 53:
                color, alpha, zorder = "#e74c3c", 0.3, 5  # Active 4
            elif p.layer == 60:
                color, alpha, zorder = "#e67e22", 0.2, 1  # SiO2
            elif p.layer == 30:
                color, alpha, zorder = "#ecf0f1", 0.2, 0  # Substrate
            else:
                continue
            ax.add_patch(MplPolygon(p.points, closed=True, facecolor=color, alpha=alpha, edgecolor="none", zorder=zorder))

        labels = cell_to_show.get_labels(depth=None)
        for lbl in labels:
            ax.text(lbl.origin[0], lbl.origin[1], lbl.text, fontsize=8, ha="center", va="center", color="blue", alpha=0.75)

        ax.set_aspect("equal")
        bbox = cell_to_show.bounding_box()
        if bbox:
            (min_x, min_y), (max_x, max_y) = bbox
            margin = (max_x - min_x) * 0.05
            ax.set_xlim(min_x - margin, max_x + margin)
            ax.set_ylim(min_y - margin, max_y + margin)
        plt.tight_layout()
        plt.show()

    # ============================================================
    # Device generation
    # ============================================================
    def generate_quantum_dot_layout(self, **kwargs):
        """
        Generate 6-dot CMOS device core. Returns (lib, cell, connection_points).

        If kwargs is empty, use the same default parameters as
        `six_dot_generator_SG1_CMOS_20260312.py` main() for consistency.
        """
        if not kwargs:
            kwargs = dict(
                num_dots=6,
                pg_max_width=0.120,
                pg_vert_side_len=0.040,
                pg_chamfer_h=0.040,
                bg_max_width=0.060,
                bg_vert_side_len=0.040,
                gap_pg_bg=0.002,
                d1_gap=0.020,
                d2_gap=0.020,
                sg_mid_thick=0.100,
                sg_top_thick=0.100,
                sg_bot_thick=0.100,
                anneal_sd_enable=True,
                anneal_sd_gap_to_accum=0.50,
                anneal_sd_rect_w=0.10,
                anneal_sd_rect_h=0.30,
                anneal_sd_dx=0.30,
                anneal_sd_dy=0.20,
                anneal_central_set_dx=0.0,
                anneal_central_set_dy=0.5,
                anneal_sd_gap_override=False,
                accum_sd_lead_extend_to_anneal=0.20,
                anneal_set_s_enable=True,
                lead_width=0.042,
                lead_length_bot=0.4,
                lead_length_top=1,
            )

        self.device_lib, self.device_cell, self.device_connection_points = create_6qd_device_core(
            num_dots=6,
            pg_max_width=0.120,
            pg_vert_side_len=0.040,
            pg_chamfer_h=0.040,
            bg_max_width=0.060,
            bg_vert_side_len=0.040,
            gap_pg_bg=0.002,
            d1_gap=0.020,
            d2_gap=0.020,
            sg_mid_thick=0.100,
            sg_top_thick=0.100,
            sg_bot_thick=0.100,
            
            anneal_sd_enable=True,
            anneal_sd_gap_to_accum=0.50,   # 与对应 accumulation S/D 的边缘间距（500nm）- 用于计算 lead 延长/重叠
            anneal_sd_rect_w=0.10,         # 竖直长方形 x 宽度（100nm）
            anneal_sd_rect_h=0.30,         # 竖直长方形 y 高度（400nm）
            anneal_sd_dx=0.30,             # 相对 accumulation 的 x 偏移量（±500nm）
            anneal_sd_dy=0.20,  
            anneal_central_set_dx = 0.0,           # 相对 accumulation 的 y 偏移量（±500nm）
            anneal_central_set_dy= 0.5,

            anneal_sd_gap_override=False,  # True 时强制用“边缘间距”而不是 (dx,dy) 放置
            accum_sd_lead_extend_to_anneal=0.20,  # 蓝色 accumulation 的 lead 额外延长（朝向 anneal），用于与 anneal 重叠（um）
            anneal_set_s_enable=True,
            
            lead_width=0.042,
            lead_length_bot=0.4,
            lead_length_top=1
        )

        return self.device_lib, self.device_cell, self.device_connection_points

    # ============================================================
    # Pad frame generation (ported from 18dot_v3 with customization)
    # ============================================================
    def generate_pad_frame_layout(
        self,
        tp: List[str],
        bp: List[str],
        lp: List[str],
        rp: List[str],
        **kwargs,
    ):
        """
        Create pad frame with pads distributed on 4 edges.
        Order is exactly the order of the provided lists (like 18dot notebook manual lists).

        Returns:
            pad_lib, pad_cell, pad_connection_points, active_center, all_pads_info
        """
        self.pad_params.update(kwargs)
        (
            self.pad_lib,
            self.pad_cell,
            self.pad_connection_points,
            self.active_center,
            self.all_pads_info,
        ) = self._create_rect_wire_layout(tp=tp, bp=bp, lp=lp, rp=rp, **self.pad_params)
        return self.pad_lib, self.pad_cell, self.pad_connection_points, self.active_center, self.all_pads_info

    def _create_rect_wire_layout(
        self,
        tp: List[str],
        bp: List[str],
        lp: List[str],
        rp: List[str],
        **params,
    ):
        layout_width = params["layout_width"]
        layout_height = params["layout_height"]
        pad_width = params["pad_width"]
        pad_height = params["pad_height"]
        pad_spacing = params["pad_spacing"]
        edge_margin = params["edge_margin"]

        # Active Region No.1 (support independent width/height)
        active_size = params.get("active_size", 400)
        active_width = params.get("active_width", active_size)
        active_height = params.get("active_height", active_size)

        trace_width = params["trace_width"]
        trace_spacing = params["trace_spacing"]
        taper_length = params["taper_length"]
        active_entry_len = params["active_entry_len"]

        active_size_2 = params.get("active_size_2", 150)
        inner_trace_width = params.get("inner_trace_width", 1.0)
        inner_trace_spacing = params.get("inner_trace_spacing", 2.0)
        active_entry_len_2 = params.get("active_entry_len_2", 10)

        active_size_3 = params.get("active_size_3", 60)
        inner_trace_width_3 = params.get("inner_trace_width_3", 0.5)
        inner_trace_spacing_3 = params.get("inner_trace_spacing_3", 1.0)
        active_entry_len_3 = params.get("active_entry_len_3", 5)

        active_width_4 = params.get("active_width_4", 30)
        active_height_4 = params.get("active_height_4", 15)
        inner_trace_width_4 = params.get("inner_trace_width_4", 0.1)
        inner_trace_spacing_4 = params.get("inner_trace_spacing_4", 0.5)
        active_entry_len_4 = params.get("active_entry_len_4", 2)

        lib = gdstk.Library()
        cell = lib.new_cell("PAD_FRAME_WITH_WIRES")
        connection_points: Dict[str, Point] = {}
        all_pads = []
        global_pad_counter = 0

        # Collect pads in the order user specified per edge
        pbe = {"top": [], "bottom": [], "left": [], "right": []}

        def add_edge_pads(pads: List[str], edge: str):
            nonlocal global_pad_counter
            c = len(pads)
            if c == 0:
                return

            if edge in ["top", "bottom"]:
                total_span = c * pad_width + (c - 1) * pad_spacing
                sx = (layout_width - total_span) / 2
                y = layout_height - edge_margin - pad_height if edge == "top" else edge_margin
                for i, lbl in enumerate(pads):
                    x = sx + i * (pad_width + pad_spacing)
                    cx, cy = x + pad_width / 2, y + pad_height / 2
                    rec = {"label": lbl, "rect": (x, y, pad_width, pad_height), "center": (cx, cy), "edge": edge, "global_index": global_pad_counter}
                    global_pad_counter += 1
                    all_pads.append(rec)
                    pbe[edge].append(rec)
            else:
                # left/right edges: pads are rotated by 90deg (w = pad_height, h = pad_width)
                # IMPORTANT: The pad-to-pad edge gap on ALL sides is controlled by `pad_spacing`.
                # Since left/right pads have y-size = pad_width, the vertical pitch must be (pad_width + pad_spacing).
                total_span = c * pad_width + (c - 1) * pad_spacing
                start_y = (layout_height - total_span) / 2
                x = edge_margin if edge == "left" else layout_width - edge_margin - pad_width
                for i, lbl in enumerate(pads):
                    # Keep list order visually from top -> bottom, like the old implementation
                    y = start_y + (c - 1 - i) * (pad_width + pad_spacing)

                    # rect is (x, y, w, h) where w = pad_height and h = pad_width for left/right
                    cx, cy = x + pad_height / 2, y + pad_width / 2

                    rec = {
                        "label": lbl,
                        "rect": (x, y, pad_height, pad_width),
                        "center": (cx, cy),
                        "edge": edge,
                        "global_index": global_pad_counter,
                    }
                    global_pad_counter += 1
                    all_pads.append(rec)
                    pbe[edge].append(rec)

        add_edge_pads(tp, "top")
        add_edge_pads(bp, "bottom")
        add_edge_pads(lp, "left")
        add_edge_pads(rp, "right")

        # Active regions
        cx = layout_width / 2
        cy = layout_height / 2
        acenter = (cx, cy)

        ax1, ax2 = cx - active_width / 2, cx + active_width / 2
        ay1, ay2 = cy - active_height / 2, cy + active_height / 2

        ax1_2, ax2_2 = cx - active_size_2 / 2, cx + active_size_2 / 2
        ay1_2, ay2_2 = cy - active_size_2 / 2, cy + active_size_2 / 2

        ax1_3, ax2_3 = cx - active_size_3 / 2, cx + active_size_3 / 2
        ay1_3, ay2_3 = cy - active_size_3 / 2, cy + active_size_3 / 2

        ax1_4, ax2_4 = cx - active_width_4 / 2, cx + active_width_4 / 2
        ay1_4, ay2_4 = cy - active_height_4 / 2, cy + active_height_4 / 2

        # Assign targets on boundary of Active Region 1
        pitch = trace_width + trace_spacing

        # For stable layout:
        # - top/bottom sort by x
        # - left/right sort by y (descending)
        pbe["top"].sort(key=lambda p: p["center"][0])
        pbe["bottom"].sort(key=lambda p: p["center"][0])
        pbe["left"].sort(key=lambda p: p["center"][1], reverse=True)
        pbe["right"].sort(key=lambda p: p["center"][1], reverse=True)

        def assign_targets(pads, edge):
            n = len(pads)
            if n == 0:
                return
            ind = np.arange(n) - (n - 1) / 2
            offs = ind * pitch
            if edge == "top":
                for i, p in enumerate(pads):
                    p.update({"target": (cx + offs[i], ay2), "target_edge": "top"})
            elif edge == "bottom":
                for i, p in enumerate(pads):
                    p.update({"target": (cx + offs[i], ay1), "target_edge": "bottom"})
            elif edge == "left":
                for i, p in enumerate(pads):
                    p.update({"target": (ax1, cy + offs[-(i + 1)]), "target_edge": "left"})
            else:
                for i, p in enumerate(pads):
                    p.update({"target": (ax2, cy + offs[-(i + 1)]), "target_edge": "right"})

        for edge in ["top", "bottom", "left", "right"]:
            assign_targets(pbe[edge], edge)

        # Inner targets for regions 2/3/4
        def assign_targets_inner(pads, edge, pitch_inner, tx_top, tx_bot, tx_left, tx_right, suffix):
            n = len(pads)
            if n == 0:
                return
            ind = np.arange(n) - (n - 1) / 2
            offs = ind * pitch_inner
            if edge == "top":
                for i, p in enumerate(pads):
                    p.update({f"target_{suffix}": (cx + offs[i], tx_top)})
            elif edge == "bottom":
                for i, p in enumerate(pads):
                    p.update({f"target_{suffix}": (cx + offs[i], tx_bot)})
            elif edge == "left":
                for i, p in enumerate(pads):
                    p.update({f"target_{suffix}": (tx_left, cy + offs[-(i + 1)])})
            else:
                for i, p in enumerate(pads):
                    p.update({f"target_{suffix}": (tx_right, cy + offs[-(i + 1)])})

        inner_pitch_2 = inner_trace_width + inner_trace_spacing
        inner_pitch_3 = inner_trace_width_3 + inner_trace_spacing_3
        inner_pitch_4 = inner_trace_width_4 + inner_trace_spacing_4

        for edge in ["top", "bottom", "left", "right"]:
            assign_targets_inner(pbe[edge], edge, inner_pitch_2, ay2_2, ay1_2, ax1_2, ax2_2, "2")
            assign_targets_inner(pbe[edge], edge, inner_pitch_3, ay2_3, ay1_3, ax1_3, ax2_3, "3")
            assign_targets_inner(pbe[edge], edge, inner_pitch_4, ay2_4, ay1_4, ax1_4, ax2_4, "4")

        # Layers (same as 18dot_v3 pad frame)
        LS = 30
        LP = 31  # pad metal
        LT = 31  # traces are also on pad metal in v3
        LTXT = 40
        LA1, LA2, LA3, LA4 = 50, 51, 52, 53
        LSIO2 = 60
        sm = 10

        # Substrate outline
        cell.add(gdstk.rectangle((0, 0), (layout_width, layout_height), layer=LS))

        # Active regions
        cell.add(gdstk.rectangle((ax1, ay1), (ax2, ay2), layer=LA1))
        cell.add(gdstk.Label("Active Region No.1", (ax1, ay2), layer=LTXT, magnification=10))

        cell.add(gdstk.rectangle((ax1_2, ay1_2), (ax2_2, ay2_2), layer=LA2))
        cell.add(gdstk.Label("Active Region No.2", (ax1_2, ay2_2), layer=LTXT, magnification=10))

        cell.add(gdstk.rectangle((ax1_3, ay1_3), (ax2_3, ay2_3), layer=LA3))
        cell.add(gdstk.Label("Active Region No.3", (ax1_3, ay2_3), layer=LTXT, magnification=5))

        cell.add(gdstk.rectangle((ax1_4, ay1_4), (ax2_4, ay2_4), layer=LA4))
        cell.add(gdstk.Label("Active Region No.4", (ax1_4, ay2_4), layer=LTXT, magnification=2))

        # Taper polygon helper
        def create_taper(pad_rect, edge, taper_len, taper_w):
            x, y, w, h = pad_rect
            cxp, cyp = x + w / 2, y + h / 2
            if edge == "top":
                pts = [(x, y), (x + w, y), (cxp + taper_w / 2, y - taper_len), (cxp - taper_w / 2, y - taper_len)]
                start = (cxp, y - taper_len)
            elif edge == "bottom":
                pts = [(x, y + h), (x + w, y + h), (cxp + taper_w / 2, y + h + taper_len), (cxp - taper_w / 2, y + h + taper_len)]
                start = (cxp, y + h + taper_len)
            elif edge == "left":
                pts = [(x + w, y), (x + w, y + h), (x + w + taper_len, cyp + taper_w / 2), (x + w + taper_len, cyp - taper_w / 2)]
                start = (x + w + taper_len, cyp)
            else:
                pts = [(x, y), (x, y + h), (x - taper_len, cyp + taper_w / 2), (x - taper_len, cyp - taper_w / 2)]
                start = (x - taper_len, cyp)
            return gdstk.Polygon(pts, layer=LT), start

        # Routing helper (same logic as 18dot_v3 rsz)
        def route_straight_zigzag(sp, tp, edge, width, custom_sm=None, custom_entry=None):
            x1, y1 = sp
            x2, y2 = tp
            smooth = custom_sm if custom_sm is not None else 20
            entry = custom_entry if custom_entry is not None else active_entry_len
            pts = [(x1, y1)]

            if edge in ["top", "bottom"]:
                sy = -1 if edge == "top" else 1
                ys = y1 + sy * smooth
                pts.append((x1, ys))
                ye = y2 + sy * entry
                ah = abs(ye - ys)
                nw = abs(x2 - x1)
                sx = 1 if x2 > x1 else -1

                if nw > ah:
                    hl = ah / 2
                    ym = ys + sy * hl
                    xm1 = x1 + sx * hl
                    xm2 = x2 - sx * hl
                    pts.extend([(xm1, ym), (xm2, ym)])
                else:
                    dyd = nw
                    yte = ys + sy * dyd
                    pts.append((x2, yte))

                pts.append((x2, y2))
            else:
                sx = 1 if edge == "left" else -1
                xs = x1 + sx * smooth
                pts.append((xs, y1))
                xe = x2 + sx * entry
                aw = abs(xe - xs)
                nh = abs(y2 - y1)
                sy = 1 if y2 > y1 else -1

                if nh > aw:
                    wl = aw / 2
                    xm = xs + sx * wl
                    ym1 = y1 + sy * wl
                    ym2 = y2 - sy * wl
                    pts.extend([(xm, ym1), (xm, ym2)])
                else:
                    dxd = nh
                    xte = xs + sx * dxd
                    pts.append((xte, y2))

                pts.append((x2, y2))

            return gdstk.FlexPath(pts, width, ends="flush", joins="round", layer=LT), pts

        # Create each pad + multi-stage leads
        for p in all_pads:
            x, y, w, h = p["rect"]
            lbl = p["label"]
            edge = p["edge"]

            # SiO2 around pad
            cell.add(gdstk.rectangle((x - sm, y - sm), (x + w + sm, y + h + sm), layer=LSIO2))
            # Pad metal
            cell.add(gdstk.rectangle((x, y), (x + w, y + h), layer=LP))

            # Label text on pads
            cxp, cyp = x + w / 2, y + h / 2
            rot = 0 if w > h else math.pi / 2
            cell.add(gdstk.Label(str(lbl), (cxp, cyp), layer=LTXT, magnification=10, rotation=rot))

            # Taper + trace to Active region 1 boundary
            poly_t, t_start = create_taper(p["rect"], edge, taper_length, trace_width)
            cell.add(poly_t)
            path_obj, path_pts = route_straight_zigzag(t_start, p["target"], edge, trace_width)
            cell.add(path_obj)
            p["trace_points"] = path_pts

            # Region 1 -> Region 2
            inner_path_obj_2, _ = route_straight_zigzag(p["target"], p["target_2"], edge, inner_trace_width, custom_sm=10, custom_entry=active_entry_len_2)
            cell.add(inner_path_obj_2)

            # Region 2 -> Region 3
            inner_path_obj_3, _ = route_straight_zigzag(p["target_2"], p["target_3"], edge, inner_trace_width_3, custom_sm=5, custom_entry=active_entry_len_3)
            cell.add(inner_path_obj_3)

            # Region 3 -> Region 4
            inner_path_obj_4, _ = route_straight_zigzag(p["target_3"], p["target_4"], edge, inner_trace_width_4, custom_sm=2, custom_entry=active_entry_len_4)
            cell.add(inner_path_obj_4)

            connection_points[str(lbl)] = p["target_4"]

        return lib, cell, connection_points, acenter, all_pads

    # ============================================================
    # Assembly (device + pad frame)
    # ============================================================
    def assemble_device_and_pad_frame(self):
        if self.device_cell is None or self.pad_cell is None:
            print("Error: Device cell or Pad cell not generated.")
            return None, None

        if self.active_center is None:
            print("Error: active_center is None. Generate pad frame first.")
            return None, None

        assembly_lib = gdstk.Library(name="ASSEMBLY", unit=1e-6, precision=1e-9)
        assembly_cell = assembly_lib.new_cell("6DOT_FULL_ASSEMBLY")

        # Flatten pad frame into assembly
        assembly_cell.add(*self.pad_cell.get_polygons(depth=None))
        assembly_cell.add(*self.pad_cell.get_paths(depth=None))
        assembly_cell.add(*self.pad_cell.get_labels(depth=None))

        # Copy device geometry and translate to active_center
        device_polys = self.device_cell.get_polygons(depth=None)
        device_paths = self.device_cell.get_paths(depth=None)
        device_labels = self.device_cell.get_labels(depth=None)

        bbox = self.device_cell.bounding_box()
        if bbox is not None:
            (min_x, min_y), (max_x, max_y) = bbox
            dev_cx, dev_cy = (min_x + max_x) / 2, (min_y + max_y) / 2
            tx, ty = self.active_center[0] - dev_cx, self.active_center[1] - dev_cy

            for p in device_polys:
                p.translate(tx, ty)
            for p in device_paths:
                p.translate(tx, ty)
            for l in device_labels:
                l.origin = (l.origin[0] + tx, l.origin[1] + ty)

        assembly_cell.add(*device_polys)
        assembly_cell.add(*device_paths)
        assembly_cell.add(*device_labels)

        return assembly_lib, assembly_cell

    def plot_assembly(self, cell_to_show, title="6-dot CMOS Full Assembly", show_plot=True, zoom_factor=0.06):
        """
        Show 2 subplots:
        - Full view
        - Zoomed view around active_center (device + inner pad endpoints)

        zoom_factor is relative to Active Region No.1 size (active_width/height).
        """
        if not show_plot:
            return
        if self.active_center is None:
            raise ValueError("active_center is None. Generate pad frame first.")

        layer_config = {
            # device (CMOS)
            1: {"color": "#FF8C00", "alpha": 0.85, "label": "Anneal S/D (L1)", "zorder": 10},
            2: {"color": "#bdc3c7", "alpha": 0.35, "label": "SG (L2)", "zorder": 2},
            3: {"color": "#7B68EE", "alpha": 0.65, "label": "PG/Accum (L3)", "zorder": 11},
            4: {"color": "#FF69B4", "alpha": 0.80, "label": "BG (L4)", "zorder": 12},
            5: {"color": "#2ecc71", "alpha": 0.95, "label": "Wires (L5)", "zorder": 30},
            100: {"color": "#2980b9", "alpha": 1.0, "label": "Device Labels (L100)", "zorder": 20},
            # pad frame (18dot_v3 style layers)
            30: {"color": "#ecf0f1", "alpha": 0.2, "label": "Substrate (L30)", "zorder": 0},
            31: {"color": "#f1c40f", "alpha": 0.7, "label": "Pads/Traces (L31)", "zorder": 5},
            40: {"color": "#e67e22", "alpha": 0.9, "label": "Pad Labels (L40)", "zorder": 21},
            50: {"color": "#1abc9c", "alpha": 0.25, "label": "Active1 (L50)", "zorder": 1},
            51: {"color": "#3498db", "alpha": 0.25, "label": "Active2 (L51)", "zorder": 2},
            52: {"color": "#9b59b6", "alpha": 0.25, "label": "Active3 (L52)", "zorder": 3},
            53: {"color": "#e74c3c", "alpha": 0.25, "label": "Active4 (L53)", "zorder": 4},
            60: {"color": "#e67e22", "alpha": 0.18, "label": "SiO2 (L60)", "zorder": 0},
        }

        def draw_on_ax(ax_obj, is_zoomed=False):
            drawn_labels = set()

            # polygons
            for layer in sorted(layer_config.keys(), key=lambda x: layer_config[x]["zorder"]):
                cfg = layer_config[layer]
                polys = cell_to_show.get_polygons(depth=None, layer=layer, datatype=0)
                if not polys:
                    continue
                merged = gdstk.boolean(polys, [], "or", layer=layer)
                for gds_poly in merged:
                    label = cfg["label"] if cfg["label"] not in drawn_labels else None
                    ax_obj.add_patch(
                        MplPolygon(
                            gds_poly.points,
                            closed=True,
                            facecolor=cfg["color"],
                            edgecolor="none",
                            alpha=cfg["alpha"],
                            label=label,
                            zorder=cfg["zorder"],
                        )
                    )
                    if label:
                        drawn_labels.add(cfg["label"])

            # paths (e.g., wires created by create_gds_wire -> gdstk.FlexPath)
            for layer in sorted(layer_config.keys(), key=lambda x: layer_config[x]["zorder"]):
                cfg = layer_config[layer]
                paths = cell_to_show.get_paths(depth=None, layer=layer, datatype=0)
                if not paths:
                    continue

                for path in paths:
                    # Prefer rendering paths as filled polygons for consistent visualization
                    poly_pts_list = []
                    if hasattr(path, "to_polygons"):
                        try:
                            poly_pts_list = path.to_polygons()
                        except Exception:
                            poly_pts_list = []

                    # Filter out empty / malformed polygon outputs (prevents Matplotlib ValueError)
                    clean_poly_pts_list = []
                    for poly_pts in poly_pts_list or []:
                        arr = np.asarray(poly_pts)
                        if arr.ndim != 2 or arr.shape[0] < 3 or arr.shape[1] != 2:
                            continue
                        clean_poly_pts_list.append(arr)

                    if not clean_poly_pts_list:
                        # Fallback: if no polygon data is available, try plotting centerline
                        if hasattr(path, "points") and path.points is not None:
                            pts = np.asarray(path.points)
                            if pts.ndim == 2 and pts.shape[0] >= 2 and pts.shape[1] == 2:
                                ax_obj.plot(
                                    pts[:, 0],
                                    pts[:, 1],
                                    color=cfg["color"],
                                    alpha=cfg["alpha"],
                                    linewidth=1.2,
                                    zorder=cfg["zorder"],
                                )
                        continue

                    for poly_pts in clean_poly_pts_list:
                        label = cfg["label"] if cfg["label"] not in drawn_labels else None
                        ax_obj.add_patch(
                            MplPolygon(
                                poly_pts,
                                closed=True,
                                facecolor=cfg["color"],
                                edgecolor="none",
                                alpha=cfg["alpha"],
                                label=label,
                                zorder=cfg["zorder"],
                            )
                        )
                        if label:
                            drawn_labels.add(cfg["label"])

            # labels
            labels = cell_to_show.get_labels(depth=None)
            xmin, xmax = ax_obj.get_xlim()
            ymin, ymax = ax_obj.get_ylim()
            for lbl in labels:
                if xmin <= lbl.origin[0] <= xmax and ymin <= lbl.origin[1] <= ymax:
                    # zoomed: show more labels; full view: show shorter labels only
                    if is_zoomed or len(lbl.text) <= 6:
                        fs = 9 if is_zoomed else 7
                        ax_obj.text(
                            lbl.origin[0],
                            lbl.origin[1],
                            lbl.text,
                            fontsize=fs,
                            ha="center",
                            va="center",
                            color="blue",
                            alpha=0.8,
                        )

            ax_obj.set_aspect("equal")

        # figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), dpi=150)
        bbox = cell_to_show.bounding_box()
        if bbox:
            (min_x, min_y), (max_x, max_y) = bbox
            margin = (max_x - min_x) * 0.05
            ax1.set_xlim(min_x - margin, max_x + margin)
            ax1.set_ylim(min_y - margin, max_y + margin)

        # zoom around active center
        acx, acy = self.active_center
        aw = self.pad_params.get("active_width", self.pad_params.get("active_size", 400))
        ah = self.pad_params.get("active_height", self.pad_params.get("active_size", 400))
        zw = aw * zoom_factor
        zh = ah * zoom_factor
        ax2.set_xlim(acx - zw, acx + zw)
        ax2.set_ylim(acy - zh, acy + zh)

        draw_on_ax(ax1, is_zoomed=False)
        ax1.set_title(f"{title} (Full View)")
        draw_on_ax(ax2, is_zoomed=True)
        ax2.set_title(f"{title} (Zoomed Center)")

        # legend only on full view
        handles, labels = ax1.get_legend_handles_labels()
        if handles:
            by_label = dict(zip(labels, handles))
            ax1.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=8)

        plt.tight_layout()
        plt.show()

    # ============================================================
    # Manual wiring helpers (for notebook)
    # ============================================================
    def create_gds_wire(self, pts: List[Point], width=0.040, layer=5):
        return gdstk.FlexPath(pts, width, layer=layer, ends="flush", joins="round")

    def compute_device_translation_to_active_center(self) -> Optional[Tuple[float, float]]:
        if self.device_cell is None or self.active_center is None:
            return None
        bbox = self.device_cell.bounding_box()
        if bbox is None:
            return None
        dev_cx, dev_cy = (bbox[0][0] + bbox[1][0]) / 2, (bbox[0][1] + bbox[1][1]) / 2
        tx, ty = self.active_center[0] - dev_cx, self.active_center[1] - dev_cy
        return tx, ty

    def add_connection_wire(self, cell, label: str, width=0.040, layer=5):
        """
        Straight line wiring between pad_connection_points[label] and device_connection_points[label] (after translation).
        This is only a helper; for real manual wiring, you can build your own polyline.
        """
        if self.device_connection_points is None or self.pad_connection_points is None:
            print("Error: Connection points not generated.")
            return None

        if label not in self.device_connection_points or label not in self.pad_connection_points:
            return None

        t = self.compute_device_translation_to_active_center()
        if t is None:
            return None
        tx, ty = t

        p_pad = self.pad_connection_points[label]
        p_dev = (self.device_connection_points[label][0] + tx, self.device_connection_points[label][1] + ty)

        wire = self.create_gds_wire([p_pad, p_dev], width=width, layer=layer)
        cell.add(wire)
        return wire

    # ============================================================
    # Save helpers
    # ============================================================
    @staticmethod
    def save_gds(lib: gdstk.Library, filepath: str):
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        if os.path.exists(filepath):
            os.remove(filepath)
        lib.write_gds(filepath)
        print(f"GDS saved: {filepath}")
