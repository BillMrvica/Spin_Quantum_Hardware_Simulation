import math
from dataclasses import dataclass
from pathlib import Path

import gdstk
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class LayerMap:
    # Mask aligner Al electrodes
    MA_GREY: tuple[int, int] = (10, 0)
    MA_ORANGE: tuple[int, int] = (11, 0)

    # EBL Al electrodes
    EBL_GREY: tuple[int, int] = (12, 0)
    EBL_ORANGE: tuple[int, int] = (13, 0)

    # Al2O3 formed by O2 plasma ashing on top of Al (oxide area)
    ALOX: tuple[int, int] = (20, 0)

    # Annotation / labels
    TEXT: tuple[int, int] = (30, 0)


def add_text(cell: gdstk.Cell, text: str, origin: tuple[float, float], size_um: float, layer: tuple[int, int]):
    # gdstk.text returns a list of polygons
    polys = gdstk.text(text, size_um, origin, layer=layer[0], datatype=layer[1])
    cell.add(*polys)


def rect(center: tuple[float, float], size: tuple[float, float], layer: tuple[int, int]) -> gdstk.Polygon:
    cx, cy = center
    w, h = size
    return gdstk.rectangle(
        (cx - w / 2, cy - h / 2),
        (cx + w / 2, cy + h / 2),
        layer=layer[0],
        datatype=layer[1],
    )


def make_overlap_capacitor(
    cell: gdstk.Cell,
    origin: tuple[float, float],
    *,
    w: float,
    l: float,
    overlap: float,
    gap: float,
    alox_enclosure: float,
    layers: LayerMap,
    label: str,
):
    """
    Simple MIM-like overlap test:
      - bottom Al rectangle
      - top Al rectangle offset in y, producing an overlap region (controls effective AlOx area)
      - AlOx polygon covering the overlap area (+ optional enclosure)
    This is useful to probe oxide quality via leakage/capacitance vs area.
    """
    ox, oy = origin

    # Bottom electrode centered at origin
    bottom_center = (ox, oy)
    bottom = rect(bottom_center, (w, l), layers.AL)
    cell.add(bottom)

    # Top electrode shifted upward: create overlap region of length `overlap`
    # distance between centers along y:
    # bottom half-length + top half-length - overlap + gap
    dy = (l / 2) + (l / 2) - overlap + gap
    top_center = (ox, oy + dy)
    top = rect(top_center, (w, l), layers.AL)
    cell.add(top)

    # Overlap window bounds in y:
    bottom_top_y = oy + l / 2
    top_bottom_y = (oy + dy) - l / 2
    overlap_len = max(0.0, bottom_top_y - top_bottom_y)

    # Oxide area: cover the overlap region with optional enclosure margin.
    # If overlap_len is 0, still place a tiny marker rectangle for visibility.
    ox_w = w + 2 * alox_enclosure
    ox_l = overlap_len + 2 * alox_enclosure
    if ox_l <= 0:
        ox_l = 0.2  # 200 nm fallback

    overlap_center_y = (bottom_top_y + top_bottom_y) / 2 if overlap_len > 0 else (oy + dy / 2)
    alox = rect((ox, overlap_center_y), (ox_w, ox_l), layers.ALOX)
    cell.add(alox)

    add_text(cell, label, (ox - w / 2, oy - l / 2 - 40), 12, layers.TEXT)


def make_leakage_bridge(
    cell: gdstk.Cell,
    origin: tuple[float, float],
    *,
    pad_size: float,
    line_w: float,
    gap: float,
    bridge_len: float,
    alox_strip_w: float,
    alox_strip_l: float,
    layers: LayerMap,
    label: str,
):
    """
    Two facing Al pads with a narrow gap. Add an AlOx strip spanning the gap region.
    Useful for stressing oxide in a defined field region and probing leakage/BD.
    """
    ox, oy = origin

    # Left pad and lead
    left_center = (ox - (gap / 2 + pad_size / 2), oy)
    right_center = (ox + (gap / 2 + pad_size / 2), oy)

    cell.add(rect(left_center, (pad_size, pad_size), layers.AL))
    cell.add(rect(right_center, (pad_size, pad_size), layers.AL))

    # Narrow "bridge leads" pointing into the gap to define field concentration
    # Left lead goes right; right lead goes left
    lead_len = bridge_len
    cell.add(rect((left_center[0] + pad_size / 2 + lead_len / 2, oy), (lead_len, line_w), layers.AL))
    cell.add(rect((right_center[0] - pad_size / 2 - lead_len / 2, oy), (lead_len, line_w), layers.AL))

    # Oxide strip centered on the gap (spans slightly into both leads)
    cell.add(rect((ox, oy), (alox_strip_l, alox_strip_w), layers.ALOX))

    add_text(cell, label, (ox - pad_size, oy - pad_size / 2 - 40), 12, layers.TEXT)


def make_via_chain_oxide_windows(
    cell: gdstk.Cell,
    origin: tuple[float, float],
    *,
    n: int,
    pitch: float,
    al_size: tuple[float, float],
    alox_size: tuple[float, float],
    layers: LayerMap,
    label: str,
):
    """
    Array of Al rectangles each covered by a slightly smaller (or larger) AlOx window.
    Can be used to evaluate oxide uniformity/density across many sites.
    """
    ox, oy = origin
    w_al, h_al = al_size
    w_ox, h_ox = alox_size

    for i in range(n):
        x = ox + i * pitch
        cell.add(rect((x, oy), (w_al, h_al), layers.AL))
        cell.add(rect((x, oy), (w_ox, h_ox), layers.ALOX))

    add_text(cell, label, (ox, oy - max(h_al, h_ox) / 2 - 40), 12, layers.TEXT)


def make_alignment_frame(cell: gdstk.Cell, *, size: float, border: float, layers: LayerMap):
    # Simple frame + crosshair at origin
    half = size / 2
    outer = gdstk.rectangle((-half, -half), (half, half), layer=layers.MARK[0], datatype=layers.MARK[1])
    inner = gdstk.rectangle(
        (-(half - border), -(half - border)),
        (half - border, half - border),
        layer=layers.MARK[0],
        datatype=layers.MARK[1],
    )
    # subtract inner from outer by boolean to create border ring
    ring = gdstk.boolean(outer, inner, "not", layer=layers.MARK[0], datatype=layers.MARK[1])
    if ring:
        cell.add(*ring)

    cell.add(rect((0, 0), (size * 0.12, border), layers.MARK))
    cell.add(rect((0, 0), (border, size * 0.12), layers.MARK))


def build_library() -> gdstk.Library:
    """
    Rebuild layout from provided draw.io mxGraphModel.

    Interpreting the draw.io geometry:
      - The mxGeometry x/y/width/height are in "drawing units" (px-like).
      - The dimension arrows label 150 um, 200 um, 100 um, 10 um in the figure.
      - We map 1 drawing unit -> 1 um, because the big blocks have sizes 160/320
        and the labels indicate 150/200 um scale (close enough for a direct mapping).
        If you want a strict scaling factor, tell me one anchor pair and I will update.

    Shapes extracted (rectangles only):
      Mask Aligner:
        - Grey (#647687)  -> MA_GREY
        - Orange (#f0a30a)-> MA_ORANGE
      EBL:
        - Grey (#76608a)  -> EBL_GREY
        - Orange (#e51400)-> EBL_ORANGE
      Oxide (AlOx): generated as overlap(MA_GREY, MA_ORANGE) + overlap(EBL_GREY, EBL_ORANGE),
        with a small enclosure, on ALOX layer.
    """
    layers = LayerMap()
    lib = gdstk.Library(unit=1e-6, precision=1e-9)  # um units, nm precision
    top = lib.new_cell("ALOX_TEST_STRUCTURES")

    # --- Parameters ---
    scale = 1.0  # um per draw.io unit (see docstring)
    alox_enclosure = 1.0  # um oxide enclosure around overlap areas

    # EBL line width requirement: 50 nm = 0.05 um
    ebl_linewidth_um = 0.05

    def add_drawio_rect(
        x: float,
        y: float,
        w: float,
        h: float,
        layer: tuple[int, int],
        *,
        rotation_deg: float = 0.0,
    ):
        """
        draw.io rectangle:
          - geometry (x,y,w,h) is the axis-aligned bounding box in draw.io canvas units
          - style may include rotation (±90 deg in this design)

        We build a rectangle centered at the bbox center, then rotate around the center.
        """
        # draw.io origin is top-left, y increases downward.
        # For GDS, use y up by flipping y.
        cx = (x + w / 2) * scale
        cy = -(y + h / 2) * scale

        ww = w * scale
        hh = h * scale

        poly = gdstk.rectangle(
            (cx - ww / 2, cy - hh / 2),
            (cx + ww / 2, cy + hh / 2),
            layer=layer[0],
            datatype=layer[1],
        )

        if rotation_deg:
            poly.rotate(math.radians(rotation_deg), (cx, cy))

        return poly

    # --- Rectangles from mxGraphModel (x, y, w, h) ---
    # Per your confirmation:
    #   - Purple (#76608a) = EBL grey
    #   - Red (#e51400)    = EBL orange
    #   - Other grey/orange are Mask Aligner grey/orange
    ma_grey_rects = [
        (240, 190, 320, 20, 0),    # kfo...-1
        (80, 120, 160, 160, 0),    # kfo...-2
        
        (980, 197.5, 16, 5, 0),
        (800, 190, 180, 20, 0),    # alj...-6
        (640, 120, 160, 160, 0),   # alj...-7
        # (620, 460, 20, 40, 0),     # alj...-30 (legend grey for MA)
    ]
    ma_orange_rects = [
        (280, 230, 320, 20, -90),  # kfo...-3
        (360, 400, 160, 160, 0),   # kfo...-4
        
        (992, 209.5, 16, 5, +90),  # alj...-8
        (910, 300, 180, 20, +90),  # alj...-8
        (920, 400, 160, 160, 0),   # alj...-9
        # (590, 470, 40, 20, -90),   # alj...-29 (legend orange for MA)
    ]
    # For EBL bars, enforce linewidth = 50 nm by overriding the "thin" dimension to 0.05 um.
    # draw.io provides nominal 20 units thickness; we replace that with ebl_linewidth_um.
    ebl_grey_rects = [
        # alj...-14 (purple bar): nominal 120 x 20 -> 120 x 0.05
        (990, 200, 20, ebl_linewidth_um, 0),
        (990, 200, 20, 0.5, 0),
        # alj...-27 (legend purple for EBL): keep as-is (legend, not device)
        # (620, 410, 20, 40, 0),
    ]
    ebl_orange_rects = [
        # alj...-15 (red bar, rot=-90): nominal 120 x 20 -> 120 x 0.05 then rotate
        (990, 200, 20, ebl_linewidth_um, -90),
        # alj...-26 (legend red for EBL): keep as-is (legend, not device)
        # (590, 420, 40, 20, -90),
    ]

    # Place metals on separate layers
    for x, y, w, h, rot in ma_grey_rects:
        top.add(add_drawio_rect(x, y, w, h, layers.MA_GREY, rotation_deg=rot))
    for x, y, w, h, rot in ma_orange_rects:
        top.add(add_drawio_rect(x, y, w, h, layers.MA_ORANGE, rotation_deg=rot))
    for x, y, w, h, rot in ebl_grey_rects:
        top.add(add_drawio_rect(x, y, w, h, layers.EBL_GREY, rotation_deg=rot))
    for x, y, w, h, rot in ebl_orange_rects:
        top.add(add_drawio_rect(x, y, w, h, layers.EBL_ORANGE, rotation_deg=rot))

    # --- Oxide: overlaps within each process pair ---
    ma_grey_polys = [add_drawio_rect(x, y, w, h, layers.MA_GREY, rotation_deg=rot) for x, y, w, h, rot in ma_grey_rects]
    ma_orange_polys = [add_drawio_rect(x, y, w, h, layers.MA_ORANGE, rotation_deg=rot) for x, y, w, h, rot in ma_orange_rects]
    ebl_grey_polys = [add_drawio_rect(x, y, w, h, layers.EBL_GREY, rotation_deg=rot) for x, y, w, h, rot in ebl_grey_rects]
    ebl_orange_polys = [add_drawio_rect(x, y, w, h, layers.EBL_ORANGE, rotation_deg=rot) for x, y, w, h, rot in ebl_orange_rects]

    overlaps_ma = gdstk.boolean(ma_grey_polys, ma_orange_polys, "and", layer=layers.ALOX[0], datatype=layers.ALOX[1])
    overlaps_ebl = gdstk.boolean(ebl_grey_polys, ebl_orange_polys, "and", layer=layers.ALOX[0], datatype=layers.ALOX[1])
    overlaps = []
    if overlaps_ma:
        overlaps.extend(overlaps_ma)
    if overlaps_ebl:
        overlaps.extend(overlaps_ebl)

    if overlaps:
        grown = gdstk.offset(overlaps, alox_enclosure, join="round", layer=layers.ALOX[0], datatype=layers.ALOX[1])
        top.add(*grown)

    # --- Labels ---
    # add_text(top, "MA_GREY  L10/0 (grey #647687)", (60, -40), 12, layers.TEXT)
    # add_text(top, "MA_ORANGE L11/0 (orange #f0a30a)", (60, -60), 12, layers.TEXT)
    # add_text(top, "EBL_GREY L12/0 (purple #76608a)", (60, -80), 12, layers.TEXT)
    # add_text(top, "EBL_ORANGE L13/0 (red #e51400)", (60, -100), 12, layers.TEXT)
    # add_text(top, "AlOx (overlap within MA/EBL pairs) L20/0", (60, -120), 12, layers.TEXT)
    # add_text(top, "Imported from draw.io mxGraphModel (scale 1 unit = 1 um)", (60, -140), 12, layers.TEXT)

    return lib


def plot_cell(cell: gdstk.Cell, out_png: Path):
    # Quick visualization using matplotlib (similar style to other scripts in this repo)
    fig, ax = plt.subplots(figsize=(14, 10))

    layer_colors: dict[int, tuple[str, float]] = {
        10: ("#647687", 0.85),  # MA grey
        11: ("#f0a30a", 0.85),  # MA orange
        12: ("#76608a", 0.85),  # EBL grey (purple)
        13: ("#e51400", 0.85),  # EBL orange (red)
        20: ("#808080", 0.50),  # AlOx
        30: ("#000000", 1.00),  # text
    }

    polys = cell.get_polygons()
    for poly in polys:
        layer = poly.layer
        color, alpha = layer_colors.get(layer, ("#cccccc", 0.4))
        patch = plt.Polygon(poly.points, facecolor=color, alpha=alpha, edgecolor="black", linewidth=0.3)
        ax.add_patch(patch)

    ax.set_aspect("equal", "box")
    ax.autoscale_view()
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, linestyle="--", linewidth=0.3)
    ax.set_title(cell.name)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main():
    out_dir = Path(__file__).resolve().parent
    out_gds = out_dir / "AlOxide_test_stru.gds"
    out_png = out_dir / "AlOxide_test_stru_preview.png"

    lib = build_library()
    lib.write_gds(str(out_gds))
    print(f"Wrote: {out_gds}")

    top = lib.cells[0] if len(lib.cells) == 1 else lib.cells[-1]
    # We created the top cell as "ALOX_TEST_STRUCTURES" in build_library(); select by name if present.
    for c in lib.cells:
        if c.name == "ALOX_TEST_STRUCTURES":
            top = c
            break
    plot_cell(top, out_png)
    print(f"Wrote: {out_png}")


if __name__ == "__main__":
    main()
