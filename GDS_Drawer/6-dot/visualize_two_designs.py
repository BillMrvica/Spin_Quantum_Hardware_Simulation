"""
Visualization: compare two 6-dot device designs side-by-side

Design A:
- create_6qd_device_core (layers L0-L3)

Design B:
- create_6qd_ohmic_design (layers L4-L6)

Outputs (default):
- GDS_Drawer/6-dot/output_6dot/compare_two_designs.png

Optional:
- Save GDS for each design as well (via --save-gds)

Run:
  python GDS_Drawer/6-dot/visualize_two_designs.py
  python GDS_Drawer/6-dot/visualize_two_designs.py --save-gds

Notes:
- Uses the generator functions from six_dot_generator.py (same directory).
- Uses a single matplotlib figure with 2 subplots.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

# Import generator functions
from six_dot_generator import (
    create_6qd_device_core,
    create_6qd_ohmic_design,
)

# Layer constants (copied for clarity; these match six_dot_generator.py)
LAYER_SD = 0
LAYER_SG = 1
LAYER_BG = 2
LAYER_PG = 3

LAYER_OHMIC = 4
LAYER_ACCUM = 5
LAYER_SCREEN = 6


def _set_axis_to_cell_bbox(ax, cell, margin_frac: float = 0.10) -> None:
    bbox = cell.bounding_box()
    if bbox is None:
        ax.autoscale_view()
        return
    (min_x, min_y), (max_x, max_y) = bbox
    width = max_x - min_x
    height = max_y - min_y
    margin_x = width * margin_frac if width > 0 else 1
    margin_y = height * margin_frac if height > 0 else 1
    ax.set_xlim(min_x - margin_x, max_x + margin_x)
    ax.set_ylim(min_y - margin_y, max_y + margin_y)


def plot_gds_on_axis(ax, cell, title: str = "") -> None:
    """
    Plot device core design (layers 0-3) onto an existing matplotlib axis.
    """
    ax.set_title(title, fontsize=14)

    layer_config = {
        0: {"color": "#87CEEB", "alpha": 0.6, "label": "S/D (L0)", "zorder": 2},
        1: {"color": "#D3D3D3", "alpha": 0.5, "label": "SG (L1)", "hatch": "///", "zorder": 1},
        2: {"color": "#FF69B4", "alpha": 0.9, "label": "BG (L2)", "zorder": 3},
        3: {"color": "#9370DB", "alpha": 0.9, "label": "PG (L3)", "zorder": 4},
    }

    drawn_labels = set()
    sorted_layers = [1, 0, 2, 3]  # match generator z-order

    for layer in sorted_layers:
        polygons = cell.get_polygons(layer=layer, datatype=0)
        if not polygons:
            continue

        cfg = layer_config[layer]
        for gds_poly in polygons:
            pts = gds_poly.points
            if pts.ndim != 2 or len(pts) < 3:
                continue

            label = cfg["label"] if cfg["label"] not in drawn_labels else None
            poly = MplPolygon(
                pts,
                closed=True,
                facecolor=cfg["color"],
                edgecolor="black",
                linewidth=0.5,
                alpha=cfg["alpha"],
                hatch=cfg.get("hatch", ""),
                label=label,
                zorder=cfg.get("zorder", 1),
            )
            ax.add_patch(poly)
            if label:
                drawn_labels.add(cfg["label"])

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=10)

    ax.set_aspect("equal")
    _set_axis_to_cell_bbox(ax, cell)

    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, which="both", linestyle="--", alpha=0.25)


def plot_ohmic_on_axis(ax, cell, title: str = "") -> None:
    """
    Plot ohmic design (layers 4-6) onto an existing matplotlib axis.
    """
    ax.set_title(title, fontsize=14)

    layer_config = {
        4: {"color": "#32CD32", "alpha": 0.9, "label": "Ohmic (L4)"},
        5: {"color": "#1E90FF", "alpha": 0.9, "label": "Accum (L5)"},
        6: {"color": "#FF8C00", "alpha": 0.9, "label": "Screen (L6)"},
    }

    drawn_labels = set()
    for layer in [4, 5, 6]:
        polygons = cell.get_polygons(layer=layer, datatype=0)
        if not polygons:
            continue

        cfg = layer_config[layer]
        for gds_poly in polygons:
            pts = gds_poly.points
            if pts.ndim != 2 or len(pts) < 3:
                continue

            label = cfg["label"] if cfg["label"] not in drawn_labels else None
            poly = MplPolygon(
                pts,
                closed=True,
                facecolor=cfg["color"],
                edgecolor="black",
                linewidth=0.5,
                alpha=cfg["alpha"],
                label=label,
            )
            ax.add_patch(poly)
            if label:
                drawn_labels.add(cfg["label"])

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=10)

    ax.set_aspect("equal")
    _set_axis_to_cell_bbox(ax, cell)

    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, which="both", linestyle="--", alpha=0.25)


def build_designs() -> Tuple[Tuple[object, object, Dict], Tuple[object, object, Dict]]:
    """
    Build the two designs using the same parameters as the __main__ section
    in six_dot_generator.py (as of 2026-03-11).
    """
    # Design A: 6qd device core
    lib_a, cell_a, points_a = create_6qd_device_core(
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
        sg_top_thick=0.300,
        sg_bot_thick=0.300,
        lead_width=0.042,
        lead_length_bot=0.4,
        lead_length_top=0.5,
    )

    # Design B: ohmic design (non-reshaped screening) - avoids covering active regions
    lib_b, cell_b, points_b = create_6qd_ohmic_design(
        num_dots=4,
        ohmic_width=0.080,
        ohmic_height=0.150,
        ohmic_extension=0.300,
        accum_width=0.120,
        accum_height=0.180,
        accum_extension=0.250,
        screen_width=0.060,
        screen_height=0.140,
    )

    return (lib_a, cell_a, points_a), (lib_b, cell_b, points_b)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare two 6-dot device designs side-by-side.")
    parser.add_argument(
        "--outdir",
        default=os.path.join(os.path.dirname(__file__), "output_6dot"),
        help="Output directory for saved PNG/GDS (default: ./output_6dot next to this script).",
    )
    parser.add_argument(
        "--save-gds",
        action="store_true",
        help="Also save the GDS files for each design into outdir.",
    )
    parser.add_argument(
        "--png-name",
        default="compare_two_designs.png",
        help="Output PNG filename (default: compare_two_designs.png).",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    (lib_a, cell_a, _), (lib_b, cell_b, _) = build_designs()

    fig, axes = plt.subplots(1, 2, figsize=(20, 9), constrained_layout=True)

    plot_gds_on_axis(axes[0], cell_a, title="Design A: 6QD device core (L0-L3)")
    plot_ohmic_on_axis(axes[1], cell_b, title="Design B: Ohmic design (L4-L6)")

    png_path = os.path.join(args.outdir, args.png_name)
    fig.savefig(png_path, dpi=200)
    print(f"[OK] Saved comparison PNG: {png_path}")

    if args.save_gds:
        gds_a = os.path.join(args.outdir, "6qd_device_core.gds")
        gds_b = os.path.join(args.outdir, "6qd_ohmic_design.gds")
        lib_a.write_gds(gds_a)
        lib_b.write_gds(gds_b)
        print(f"[OK] Saved GDS A: {gds_a}")
        print(f"[OK] Saved GDS B: {gds_b}")

    # Show interactive window at the end (useful in notebooks / local runs)
    plt.show()


if __name__ == "__main__":
    main()
