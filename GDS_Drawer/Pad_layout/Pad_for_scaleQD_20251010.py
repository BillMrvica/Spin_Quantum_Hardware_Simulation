import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import re

def generate_quantum_dot_layout_v22_modified(N, set_group_spacing=50, set_horizontal_margin=300):
    """
    Generates a layout with locally grouped triangular SETs.
    - Places SET_D pads in the space between the top-row SET groups.
    - Places SET_S to the left of the SET1_D pad.
    """
    # 1. Constants
    layout_width = 2500
    layout_height = 2500
    pad_width = 100
    pad_height = 150
    pad_spacing = 30
    row_spacing = 40
    horizontal_margin = 300
    vertical_margin = 350

    def natural_sort_key(s):
        return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

    # 2. Label Generation
    num_sets = math.ceil(N / 3)
    sg_labels = [f"SG{i+1}" for i in range(3)]
    pg_labels = sorted([f"QD_PG{i+1}" for i in range(N)], key=natural_sort_key)
    bg_labels = sorted([f"QD_B{i+1}" for i in range(N + 1)], key=natural_sort_key)
    
    set_d_labels = [f"SET{i+1}_D" for i in range(num_sets)]
    set_s_label = ["SET_S"]

    set_groups = []
    for i in range(1, num_sets + 1):
        set_groups.append({
            'b1': f'SET{i}_B1',
            'b2': f'SET{i}_B2',
            'g': f'SET{i}_G'
        })
        
    pad_placements = []
    pitch = pad_width + pad_spacing
    available_horizontal_width = layout_width - 2 * horizontal_margin

    # 3. Placement Logic
    def place_row_and_get_overflow(labels, y_pos, start_x_abs):
        num_can_fit = math.floor((layout_width - horizontal_margin - start_x_abs + pad_spacing) / pitch) if labels else 0
        pads_on_row = labels[:num_can_fit]
        overflow_pads = labels[num_can_fit:]
        for i, label in enumerate(pads_on_row):
            x = start_x_abs + i * pitch
            pad_placements.append({'label': label, 'position': (x, y_pos), 'dimensions': (pad_width, pad_height), 'rotation': 90})
        return overflow_pads

    # --- Bottom Rows ---
    num_bg_on_bottom = math.floor((available_horizontal_width + pad_spacing) / pitch) if bg_labels else 0
    fit_width_bg = num_bg_on_bottom * pad_width + (num_bg_on_bottom - 1) * pad_spacing if num_bg_on_bottom > 0 else 0
    start_x_bg_abs = horizontal_margin + (available_horizontal_width - fit_width_bg) / 2
    bg_overflow = place_row_and_get_overflow(bg_labels, 0, start_x_bg_abs)

    start_x_pg_abs = start_x_bg_abs + pitch / 2
    pg_overflow = place_row_and_get_overflow(pg_labels, pad_height + row_spacing, start_x_pg_abs)
    
    # --- Top Row Logic ---
    t1_overflow, t2_overflow, t3_overflow = [], [], []
    group_width = 2 * pad_width + pad_spacing
    
    min_spacing_for_drain = pad_width + pad_spacing
    # effective_spacing = max(set_group_spacing, min_spacing_for_drain)
    effective_spacing = set_group_spacing
    
    unit_width = group_width + effective_spacing
    num_groups_can_fit = math.floor((layout_width - 2 * set_horizontal_margin + effective_spacing) / unit_width) if set_groups else 0
    groups_on_top = set_groups[:num_groups_can_fit]
    overflow_groups = set_groups[num_groups_can_fit:]
    
    set_d_on_top = set_d_labels[:num_groups_can_fit]
    set_d_overflow = set_d_labels[num_groups_can_fit:]

    total_fit_width = num_groups_can_fit * group_width + (num_groups_can_fit - 1) * effective_spacing if num_groups_can_fit > 0 else 0
    start_x_groups = horizontal_margin + (available_horizontal_width - total_fit_width) / 2
    
    y1_top, y2_top, y3_top = layout_height - pad_height, layout_height - 2*pad_height - row_spacing, layout_height - 3*pad_height - 2*row_spacing
    
    current_x = start_x_groups
    for i, group in enumerate(groups_on_top):
        g_x = current_x + (group_width - pad_width) / 2
        if (i % 2) == 0: # Up-triangle
            pad_placements.extend([
                {'label': group['b1'], 'position': (current_x, y1_top), 'dimensions': (pad_width, pad_height), 'rotation': 90},
                {'label': group['b2'], 'position': (current_x + pad_width + pad_spacing, y1_top), 'dimensions': (pad_width, pad_height), 'rotation': 90},
                {'label': group['g'], 'position': (g_x, y2_top), 'dimensions': (pad_width, pad_height), 'rotation': 90}
            ])
        else: # Down-triangle
            pad_placements.extend([
                {'label': group['g'], 'position': (g_x, y1_top), 'dimensions': (pad_width, pad_height), 'rotation': 90},
                {'label': group['b1'], 'position': (current_x, y2_top), 'dimensions': (pad_width, pad_height), 'rotation': 90},
                {'label': group['b2'], 'position': (current_x + pad_width + pad_spacing, y2_top), 'dimensions': (pad_width, pad_height), 'rotation': 90}
            ])
        
        # Place SET_D in the gap *after* the current group
        d_pad_x = current_x + group_width + (effective_spacing - pad_width) / 2
        pad_placements.append({'label': set_d_on_top[i], 'position': (d_pad_x, y3_top), 'dimensions': (pad_width, pad_height), 'rotation': 90})
        
        # MODIFICATION: If this is the first group (SET1), place SET_S to its left
        if i == 0 and set_s_label:
            s_pad_x = d_pad_x - unit_width
            pad_placements.append({'label': set_s_label[0], 'position': (s_pad_x, y3_top), 'dimensions': (pad_width, pad_height), 'rotation': 90})

        current_x += unit_width

    # --- Collate Overflows ---
    for group in overflow_groups:
        t1_overflow.extend([group['b1'], group['b2']])
        t2_overflow.append(group['g'])
    t3_overflow.extend(set_d_overflow)

    # If SET1 didn't fit on top, SET_S must also go to overflow
    if num_groups_can_fit < 1 and set_s_label:
        t3_overflow.extend(set_s_label)

    # --- Side Placement (copied from previous version) ---
    w_side, h_side = pad_height, pad_width
    pitch_side = h_side + pad_spacing
    base_height = vertical_margin
    cross_height = pad_width*0.5+pad_spacing
    
    bg_overflow_left, bg_overflow_right = bg_overflow[:len(bg_overflow)//2], bg_overflow[len(bg_overflow)//2:]
    pg_overflow_left, pg_overflow_right = pg_overflow[:len(pg_overflow)//2], pg_overflow[len(pg_overflow)//2:]
    t1_overflow_left, t1_overflow_right = sorted(t1_overflow, key=natural_sort_key)[:len(t1_overflow)//2], sorted(t1_overflow, key=natural_sort_key)[len(t1_overflow)//2:]
    t2_overflow_left, t2_overflow_right = sorted(t2_overflow, key=natural_sort_key)[:len(t2_overflow)//2], sorted(t2_overflow, key=natural_sort_key)[len(t2_overflow)//2:]
    t3_overflow_left, t3_overflow_right = sorted(t3_overflow, key=natural_sort_key)[:len(t3_overflow)//2], sorted(t3_overflow, key=natural_sort_key)[len(t3_overflow)//2:]

    def place_side_overflows_and_get_boundary(overflows, side, is_outer, is_top):
        y_boundary = base_height if not is_top else layout_height - base_height
        x_pos = 0 if is_outer else w_side + row_spacing
        if side == 'right': x_pos = layout_width - w_side if is_outer else layout_width - 2*w_side - row_spacing
        for i, label in enumerate(overflows):
            if is_top:
                y_pos = layout_height - h_side - i * pitch_side - base_height - (0 if is_outer else cross_height)
                y_boundary = y_pos
            else:
                y_pos = base_height + i * pitch_side + (0 if not is_outer else cross_height)
                y_boundary = y_pos + h_side
            pad_placements.append({'label': label, 'position': (x_pos, y_pos), 'dimensions': (w_side, h_side), 'rotation': 90})
        return y_boundary

    y_end_left_bottom_outer = place_side_overflows_and_get_boundary(bg_overflow_left, 'left', True, False)
    y_end_left_bottom_inner = place_side_overflows_and_get_boundary(pg_overflow_left, 'left', False, False)
    y_end_right_bottom_outer = place_side_overflows_and_get_boundary(bg_overflow_right, 'right', True, False)
    y_end_right_bottom_inner = place_side_overflows_and_get_boundary(pg_overflow_right, 'right', False, False)
    place_side_overflows_and_get_boundary(t1_overflow_left, 'left', True, True)
    place_side_overflows_and_get_boundary(t2_overflow_left + t3_overflow_left, 'left', False, True)
    place_side_overflows_and_get_boundary(t1_overflow_right, 'right', True, True)
    place_side_overflows_and_get_boundary(t2_overflow_right + t3_overflow_right, 'right', False, True)

    sg_at_pos1, sg_at_pos2, sg_at_pos3 = sg_labels[2], sg_labels[0], sg_labels[1]
    qdd_label, qds_label = 'QD_D', 'QD_S'
    outer_left_manual_group = ['non', 'non', qdd_label, sg_at_pos2]
    inner_left_manual_group = ['GND', sg_at_pos1, 'GND']
    qds_group = [qds_label]
    sg2_group = [sg_at_pos3]
    x_outer_left, x_inner_left = 0, w_side + row_spacing
    x_outer_right, x_inner_right = layout_width - w_side, layout_width - 2*w_side - row_spacing
    y_start_pos_ol = y_end_left_bottom_outer
    if bg_overflow_left: y_start_pos_ol += pad_spacing
    for i, label in enumerate(outer_left_manual_group):
        pad_placements.append({'label': label, 'position': (x_outer_left, y_start_pos_ol + i * pitch_side), 'dimensions': (w_side, h_side), 'rotation': 90})
    y_start_pos_il = y_end_left_bottom_inner
    if pg_overflow_left: y_start_pos_il += pad_spacing
    for i, label in enumerate(inner_left_manual_group):
        pad_placements.append({'label': label, 'position': (x_inner_left, y_start_pos_il + i * pitch_side), 'dimensions': (w_side, h_side), 'rotation': 90})
    y_start_pos_or = y_end_right_bottom_outer
    if bg_overflow_right: y_start_pos_or += pad_spacing
    for i, label in enumerate(sg2_group):
        pad_placements.append({'label': label, 'position': (x_outer_right, y_start_pos_or + i * pitch_side), 'dimensions': (w_side, h_side), 'rotation': 90})
    y_start_pos_ir = y_end_right_bottom_inner
    if pg_overflow_right: y_start_pos_ir += pad_spacing
    for i, label in enumerate(qds_group):
        pad_placements.append({'label': label, 'position': (x_inner_right, y_start_pos_ir + i * pitch_side), 'dimensions': (w_side, h_side), 'rotation': 90})

    # 5. Visualization
    total_pads = len(pad_placements)
    fig, ax = plt.subplots(figsize=(14, 14))
    ax.set_xlim(-200, layout_width + 200); ax.set_ylim(-200, layout_height + 200)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f'Conceptual Layout for N={N} Quantum Dots ({total_pads} pads total)', fontsize=16)
    chip = patches.Rectangle((0, 0), layout_width, layout_height, lw=2, ec='black', fc='lightgrey')
    ax.add_patch(chip)
    active_area = patches.Rectangle((layout_width/4, layout_height/4), layout_width/2, layout_height/2, lw=1, ec='red', fc='lightblue')
    ax.add_patch(active_area)
    ax.text(layout_width/2, layout_height/2, 'Active Area', ha='center', va='center', fontsize=12, color='darkred')
    for pad in pad_placements:
        x, y = pad['position']
        w, h = pad['dimensions']
        rect = patches.Rectangle((x, y), w, h, linewidth=1, edgecolor='black', facecolor='gold')
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, pad['label'], ha='center', va='center', fontsize=5, rotation=pad['rotation'])
    plt.show()

generate_quantum_dot_layout_v22_modified(N=32, set_group_spacing=-40, set_horizontal_margin=100)