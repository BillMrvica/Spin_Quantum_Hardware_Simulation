import gdstk
import math
import re

def create_quantum_dot_cell(N=28):
    """
    Creates a GDSII Cell containing the quantum dot layout. This function is correct.
    """
    # 1. Constants (units are now in micrometers, um)
    layout_width, layout_height = 2500, 2500
    pad_width, pad_height = 100, 150
    pad_spacing, row_spacing = 30, 40
    horizontal_margin, vertical_margin = 300, 350
    set_group_spacing, set_horizontal_margin = -30, 200

    qd_cell = gdstk.Cell("QD_Chip")
    
    # Layer 100: Cell boundary
    qd_cell.add(gdstk.rectangle((0, 0), (layout_width, layout_height), layer=100))
    
    # --- This entire placement logic is confirmed to be working correctly ---
    def natural_sort_key(s): return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]
    
    num_sets = math.ceil(N / 3)
    pg_labels = sorted([f"QD_PG{i+1}" for i in range(N)], key=natural_sort_key)
    bg_labels = sorted([f"QD_B{i+1}" for i in range(N + 1)], key=natural_sort_key)
    set_d_labels = [f"SET{i+1}_D" for i in range(num_sets)]
    set_s_label = ["SET_S"]
    set_groups = [{'b1': f'SET{i}_B1', 'b2': f'SET{i}_B2', 'g': f'SET{i}_G'} for i in range(1, num_sets + 1)]
    
    pad_placements = []
    pitch = pad_width + pad_spacing
    available_horizontal_width = layout_width - 2 * horizontal_margin

    def place_row(labels, y, x_start):
        num_fit = math.floor((layout_width - horizontal_margin - x_start + pad_spacing) / pitch) if labels else 0
        on_row, overflow = labels[:num_fit], labels[num_fit:]
        for i, lbl in enumerate(on_row):
            pad_placements.append({'label': lbl, 'position': (x_start + i * pitch, y), 'dimensions': (pad_width, pad_height)})
        return overflow

    num_bg_bot = math.floor((available_horizontal_width + pad_spacing) / pitch)
    start_x_bg = horizontal_margin + (available_horizontal_width - (num_bg_bot * pad_width + (num_bg_bot - 1) * pad_spacing)) / 2
    bg_overflow = place_row(bg_labels, 0, start_x_bg)
    pg_overflow = place_row(pg_labels, pad_height + row_spacing, start_x_bg + pitch / 2)
    
    t1_ov, t2_ov, t3_ov = [], [], []
    group_w, eff_space = 2 * pad_width + pad_spacing, set_group_spacing
    unit_w = group_w + eff_space
    num_groups_fit = math.floor((layout_width - 2 * set_horizontal_margin + eff_space) / unit_w)
    groups_top, overflow_groups = set_groups[:num_groups_fit], set_groups[num_groups_fit:]
    set_d_top, set_d_ov = set_d_labels[:num_groups_fit], set_d_labels[num_groups_fit:]
    total_w = num_groups_fit * group_w + (num_groups_fit - 1) * eff_space
    start_x_groups = horizontal_margin + (available_horizontal_width - total_w) / 2
    y1, y2, y3 = layout_height - pad_height, layout_height - 2 * pad_height - row_spacing, layout_height - 3 * pad_height - 2 * row_spacing
    cx = start_x_groups
    
    for i, grp in enumerate(groups_top):
        gx = cx + (group_w - pad_width) / 2
        if i % 2 == 0:
            pad_placements.extend([
                {'label': grp['b1'], 'position': (cx, y1), 'dimensions': (pad_width, pad_height)},
                {'label': grp['b2'], 'position': (cx + pitch, y1), 'dimensions': (pad_width, pad_height)},
                {'label': grp['g'], 'position': (gx, y2), 'dimensions': (pad_width, pad_height)}])
        else:
            pad_placements.extend([
                {'label': grp['g'], 'position': (gx, y1), 'dimensions': (pad_width, pad_height)},
                {'label': grp['b1'], 'position': (cx, y2), 'dimensions': (pad_width, pad_height)},
                {'label': grp['b2'], 'position': (cx + pitch, y2), 'dimensions': (pad_width, pad_height)}])
        dx = cx + group_w + (eff_space - pad_width) / 2
        pad_placements.append({'label': set_d_top[i], 'position': (dx, y3), 'dimensions': (pad_width, pad_height)})
        if i == 0 and set_s_label:
            pad_placements.append({'label': set_s_label[0], 'position': (dx - unit_w, y3), 'dimensions': (pad_width, pad_height)})
        cx += unit_w
        
    for grp in overflow_groups: t1_ov.extend([grp['b1'], grp['b2']]); t2_ov.append(grp['g'])
    t3_ov.extend(set_d_ov)
    if num_groups_fit < 1 and set_s_label: t3_ov.extend(set_s_label)
    
    w_s, h_s, p_s, base_h, cross_h = pad_height, pad_width, pad_width + pad_spacing, vertical_margin, pad_width * 0.5 + pad_spacing
    o_l_man, i_l_man, qds_g, sg2_g = ['non', 'non', 'QD_D', 'SG1'], ['GND', 'SG3', 'GND'], ['QD_S'], ['SG2']
    man_l, man_r = len(o_l_man) + len(i_l_man), len(qds_g) + len(sg2_g)
    bot_ov, top_ov = pg_overflow + bg_overflow, t1_ov + t2_ov + t3_ov
    target_l_bot = math.floor((len(bot_ov) - (man_l - man_r)) / 2)
    num_pg_l = min(len(pg_overflow), math.ceil(target_l_bot / 2))
    num_bg_l = target_l_bot - num_pg_l
    pg_ov_l, pg_ov_r = pg_overflow[:num_pg_l], pg_overflow[num_pg_l:]
    bg_ov_l, bg_ov_r = bg_overflow[:num_bg_l], bg_overflow[num_bg_l:]
    num_t_l = math.ceil(len(top_ov) / 2)
    t_ov_l, t_ov_r = top_ov[:num_t_l], top_ov[num_t_l:]

    def place_side(ovfl, side, outer, top):
        y_b = base_h if not top else layout_height - base_h
        x_p = 0 if outer else w_s + row_spacing
        if side == 'right': x_p = layout_width - w_s if outer else layout_width - 2 * w_s - row_spacing
        for i, lbl in enumerate(ovfl):
            y_p = (layout_height - h_s - i * p_s - base_h - (0 if outer else cross_h)) if top else (base_h + i * p_s + (0 if not outer else cross_h))
            y_b = y_p if top else y_p + h_s
            pad_placements.append({'label': lbl, 'position': (x_p, y_p), 'dimensions': (w_s, h_s)})
        return y_b

    y_end_l_b_o, y_end_l_b_i = place_side(bg_ov_l, 'left', True, False), place_side(pg_ov_l, 'left', False, False)
    y_end_r_b_o, y_end_r_b_i = place_side(bg_ov_r, 'right', True, False), place_side(pg_ov_r, 'right', False, False)
    place_side(sorted([p for p in t_ov_l if '_B' in p]), 'left', True, True)
    place_side(sorted([p for p in t_ov_l if '_G' in p or '_D' in p or '_S' in p]), 'left', False, True)
    place_side(sorted([p for p in t_ov_r if '_B' in p]), 'right', True, True)
    place_side(sorted([p for p in t_ov_r if '_G' in p or '_D' in p or '_S' in p]), 'right', False, True)
    
    y_start_ol = y_end_l_b_o + pad_spacing if bg_ov_l else y_end_l_b_o
    for i, lbl in enumerate(o_l_man): pad_placements.append({'label': lbl, 'position': (0, y_start_ol + i * p_s), 'dimensions': (w_s, h_s)})
    y_start_il = y_end_l_b_i + pad_spacing if pg_ov_l else y_end_l_b_i
    for i, lbl in enumerate(i_l_man): pad_placements.append({'label': lbl, 'position': (w_s + row_spacing, y_start_il + i * p_s), 'dimensions': (w_s, h_s)})
    y_start_or = y_end_r_b_o + pad_spacing if bg_ov_r else y_end_r_b_o
    for i, lbl in enumerate(sg2_g): pad_placements.append({'label': lbl, 'position': (layout_width - w_s, y_start_or + i * p_s), 'dimensions': (w_s, h_s)})
    y_start_ir = y_end_r_b_i + pad_spacing if pg_ov_r else y_end_r_b_i
    for i, lbl in enumerate(qds_g): pad_placements.append({'label': lbl, 'position': (layout_width - 2 * w_s - row_spacing, y_start_ir + i * p_s), 'dimensions': (w_s, h_s)})

    for pad in pad_placements:
        if "non" in pad['label']: continue
        pos, dims = pad['position'], pad['dimensions']
        rect = gdstk.rectangle(pos, (pos[0] + dims[0], pos[1] + dims[1]), layer=1)
        qd_cell.add(rect)
        label_pos = (pos[0] + dims[0] / 2, pos[1] + dims[1] / 2)
        text = gdstk.Label(pad['label'], label_pos, layer=10)
        qd_cell.add(text)
    return qd_cell

def generate_wafer_layout_with_flat(cell, cell_distance_mm=0.5, wafer_flat_length_mm=32.5):
    """
    Creates a GDSII cell for a 4-inch wafer with a bottom flat, and places the
    cell only where it fits completely inside the wafer boundary.
    """
    wafer_diameter_um = 101600
    wafer_radius_um = wafer_diameter_um / 2
    flat_length_um = wafer_flat_length_mm * 1000
    
    cell_width, cell_height = cell.bounding_box()[1]
    cell_dist_um = cell_distance_mm * 1000
    pitch_x, pitch_y = cell_width + cell_dist_um, cell_height + cell_dist_um
    nx, ny = math.floor(wafer_radius_um / pitch_x), math.floor(wafer_radius_um / pitch_y)
    
    print(f"Placing cell '{cell.name}' with size {cell_width}x{cell_height} um on a {wafer_diameter_um} um wafer.")
    print(f"Cell-to-cell distance: {cell_dist_um} um. Grid pitch: {pitch_x}x{pitch_y} um.")
    print(f"Wafer flat length: {flat_length_um} um.")
    
    wafer_cell = gdstk.Cell("Wafer_Layout")
    
    # --- Create the D-shaped wafer boundary ---
    # Calculate the y-position of the horizontal flat line using pythagorean theorem
    # r^2 = (L/2)^2 + d^2  => d = sqrt(r^2 - (L/2)^2)
    flat_y_pos = -math.sqrt(wafer_radius_um**2 - (flat_length_um / 2)**2)
    
    # Create the full circle
    wafer_circle = gdstk.ellipse((0, 0), wafer_radius_um)
    # Create a large rectangle to cut off the bottom
    cutter = gdstk.rectangle((-wafer_radius_um, -wafer_radius_um), (wafer_radius_um, flat_y_pos))
    # Subtract the cutter from the circle
    wafer_shape = gdstk.boolean(wafer_circle, cutter, "not")
    wafer_cell.add(*wafer_shape) # Add the resulting polygon(s) to the cell
    
    # --- Place chips based on the new boundary ---
    for i in range(-nx, nx + 1):
        for j in range(-ny, ny + 1):
            center_x = i * pitch_x
            center_y = j * pitch_y
            
            # Define the four corners of the chip's bounding box
            c1 = (center_x - cell_width/2, center_y - cell_height/2) # Bottom-left
            c2 = (center_x + cell_width/2, center_y - cell_height/2) # Bottom-right
            c3 = (center_x + cell_width/2, center_y + cell_height/2) # Top-right
            c4 = (center_x - cell_width/2, center_y + cell_height/2) # Top-left
            
            # Condition 1: All four corners must be inside the circular part of the wafer
            in_circle = all(corner[0]**2 + corner[1]**2 < wafer_radius_um**2 for corner in [c1, c2, c3, c4])
            # Condition 2: The chip's bottom edge must be above the flat line
            above_flat = c1[1] > flat_y_pos
            
            if in_circle and above_flat:
                ref_origin_x = center_x - cell_width / 2
                ref_origin_y = center_y - cell_height / 2
                ref = gdstk.Reference(cell, (ref_origin_x, ref_origin_y))
                wafer_cell.add(ref)

    return wafer_cell

# --- Main script execution ---
if __name__ == "__main__":
    print("Starting GDSII generation process...")
    
    # --- STAGE 1: Create and save the single quantum dot chip ---
    print("\n[Stage 1] Creating the single chip cell...")
    quantum_dot_chip = create_quantum_dot_cell(N=28)
    chip_lib = gdstk.Library()
    chip_lib.add(quantum_dot_chip)
    chip_output_file = "quantum_dot_chip.gds"
    chip_lib.write_gds(chip_output_file)
    print(f"Successfully generated single chip GDS file: {chip_output_file}")

    # --- STAGE 2: Create and save the full wafer layout ---
    print("\n[Stage 2] Creating the full 4-inch wafer layout with flat...")
    full_wafer = generate_wafer_layout_with_flat(quantum_dot_chip, cell_distance_mm=2, wafer_flat_length_mm=32.5)
    wafer_lib = gdstk.Library()
    wafer_lib.add(full_wafer, quantum_dot_chip)
    wafer_output_file = "4in_30QD_2mm_20251010.gds"
    wafer_lib.write_gds(wafer_output_file)
    
    print(f"Successfully generated wafer GDS file: {wafer_output_file}")
    print("\nProcess complete.")