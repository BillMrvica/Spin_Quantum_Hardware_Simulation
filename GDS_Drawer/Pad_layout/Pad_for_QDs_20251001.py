import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import re

def generate_quantum_dot_layout_v16_modified(N):
    """
    Generates a layout with locally grouped triangular SETs, as per user's drawing.
    - Adds two GND pads around SG1 for microwave signals.
    - Repositions QD_D to the inner left column, between SG1 and SG2.
    - Repositions QD_S to the inner right column, below SG3.
    """
    # 1. Constants from user's code
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
    
    # QD_D and QD_S are removed from this list for manual placement
    qd_sd_labels = [] 
    
    set_d_labels = [f"SET{i+1}_D" for i in range(num_sets)]
    set_s_label = ["SET_S"]

    # New SET labels grouped by SET number
    set_groups = []
    for i in range(1, num_sets + 1):
        set_groups.append({
            'b1': f'SET{i}_B1',
            'b2': f'SET{i}_B2',
            'g': f'SET{i}_G'
        })
        
    top_row_3_labels = sorted(qd_sd_labels + set_d_labels + set_s_label, key=natural_sort_key)
    
    pad_placements = []
    pitch = pad_width + pad_spacing
    available_horizontal_width = layout_width - 2 * horizontal_margin

    # 4. Placement Logic

    def place_row_and_get_overflow(labels, y_pos, available_width, start_x_abs):
        num_can_fit = 0
        if labels:
            num_can_fit = math.floor((layout_width - horizontal_margin - start_x_abs + pad_spacing) / pitch)
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
    bg_overflow = place_row_and_get_overflow(bg_labels, 0, available_horizontal_width, start_x_bg_abs)

    start_x_pg_abs = start_x_bg_abs + pitch / 2
    pg_overflow = place_row_and_get_overflow(pg_labels, pad_height + row_spacing, layout_width, start_x_pg_abs)
    
    # --- New Top Row Logic (Grouped SETs) ---
    t1_overflow, t2_overflow = [], []
    group_width = 2 * pad_width + pad_spacing
    group_pitch = group_width + row_spacing

    two_group_width = 3 * pad_width + (pad_width+pad_spacing)/2
    
    num_groups_can_fit = math.floor((available_horizontal_width + row_spacing)*2 / two_group_width) if set_groups else 0
    groups_on_top = set_groups[:num_groups_can_fit]
    overflow_groups = set_groups[num_groups_can_fit:]
    
    total_fit_width = len(groups_on_top)/2 * two_group_width + (len(groups_on_top) - 1) * row_spacing /2 if groups_on_top else 0
    start_x_groups = horizontal_margin + (available_horizontal_width - total_fit_width) / 2
    
    y1_top = layout_height - pad_height
    y2_top = layout_height - 2 * pad_height - row_spacing
    
    for i, group in enumerate(groups_on_top):
        if (i+1)%2==1:
            group_start_x = start_x_groups + 3*(i) * (pad_width + pad_spacing)/2
            pad_placements.append({'label': group['b1'], 'position': (group_start_x, y1_top), 'dimensions': (pad_width, pad_height), 'rotation': 90})
            pad_placements.append({'label': group['b2'], 'position': (group_start_x + pad_width + pad_spacing, y1_top), 'dimensions': (pad_width, pad_height), 'rotation': 90})
            g_x = group_start_x + (group_width - pad_width) / 2
            pad_placements.append({'label': group['g'], 'position': (g_x, y2_top), 'dimensions': (pad_width, pad_height), 'rotation': 90})    
        else:
            group_start_x = start_x_groups + 2 * (pad_width + pad_spacing) + 3*(i-1) * (pad_width + pad_spacing)/2
            pad_placements.append({'label': group['g'], 'position': (group_start_x, y1_top), 'dimensions': (pad_width, pad_height), 'rotation': 90})
            g_x = group_start_x - (group_width - pad_width) / 2
            pad_placements.append({'label': group['b1'], 'position': (g_x, y2_top), 'dimensions': (pad_width, pad_height), 'rotation': 90})
            pad_placements.append({'label': group['b2'], 'position': (g_x + pad_width + pad_spacing, y2_top), 'dimensions': (pad_width, pad_height), 'rotation': 90})
            

    for group in overflow_groups:
        t1_overflow.extend([group['b1'], group['b2']])
        t2_overflow.append(group['g'])

    # --- Top Row 3 (D/S) ---
    num_t3_on_top = math.floor((available_horizontal_width + pad_spacing) / pitch) if top_row_3_labels else 0
    fit_width_t3 = num_t3_on_top * pad_width + (num_t3_on_top - 1) * pad_spacing if num_t3_on_top > 0 else 0
    start_x_t3_abs = horizontal_margin + (available_horizontal_width - fit_width_t3) / 2
    t3_overflow = place_row_and_get_overflow(top_row_3_labels, layout_height - 3*pad_height - 2*row_spacing, available_horizontal_width, start_x_t3_abs)

    # --- Side Placement ---
    w_side, h_side = pad_height, pad_width
    pitch_side = h_side + pad_spacing
    base_height = vertical_margin
    cross_height = pad_width*0.5+pad_spacing
    
    bg_overflow_left, bg_overflow_right = bg_overflow[:len(bg_overflow)//2], bg_overflow[len(bg_overflow)//2:]
    pg_overflow_left, pg_overflow_right = pg_overflow[:len(pg_overflow)//2], pg_overflow[len(pg_overflow)//2:]
    t1_overflow_left, t1_overflow_right = sorted(t1_overflow, key=natural_sort_key)[:len(t1_overflow)//2], sorted(t1_overflow, key=natural_sort_key)[len(t1_overflow)//2:]
    t2_overflow_left, t2_overflow_right = sorted(t2_overflow, key=natural_sort_key)[:len(t2_overflow)//2], sorted(t2_overflow, key=natural_sort_key)[len(t2_overflow)//2:]
    t3_overflow_left, t3_overflow_right = sorted(t3_overflow, key=natural_sort_key)[:len(t3_overflow)//2], sorted(t3_overflow, key=natural_sort_key)[len(t3_overflow)//2:]

    def place_side_overflows(overflows, side, is_outer, is_top):
        y_boundary = base_height if not is_top else layout_height - base_height
        x_pos = 0 if is_outer else w_side + row_spacing
        if side == 'right': x_pos = layout_width - w_side if is_outer else layout_width - 2*w_side - row_spacing
        
        for i, label in enumerate(overflows):
            if is_top:
                y_pos = layout_height - h_side - i * pitch_side - base_height - (0 if is_outer else cross_height)
                if i == len(overflows) - 1: y_boundary = min(y_boundary, y_pos)
            else:
                y_pos = base_height + i * pitch_side + (0 if not is_outer else cross_height)
                if i == len(overflows) - 1: y_boundary = max(y_boundary, y_pos + h_side)
            pad_placements.append({'label': label, 'position': (x_pos, y_pos), 'dimensions': (w_side, h_side), 'rotation': 90})
        return y_boundary

    y_end_left_bottom = max(place_side_overflows(bg_overflow_left, 'left', True, False), place_side_overflows(pg_overflow_left, 'left', False, False))
    y_end_right_bottom = max(place_side_overflows(bg_overflow_right, 'right', True, False), place_side_overflows(pg_overflow_right, 'right', False, False))
    y_start_left_top = min(place_side_overflows(t1_overflow_left, 'left', True, True), place_side_overflows(t2_overflow_left + t3_overflow_left, 'left', False, True))
    y_start_right_top = min(place_side_overflows(t1_overflow_right, 'right', True, True), place_side_overflows(t2_overflow_right + t3_overflow_right, 'right', False, True))

    # --- New Manual Side Pad Placement for SG, GND, QD_D, QD_S ---
    # Apply the SG swap from the previous request: SG1->SG3, SG2->SG1, SG3->SG2
    sg_at_pos1 = sg_labels[2] # SG3
    sg_at_pos2 = sg_labels[0] # SG1
    sg_at_pos3 = sg_labels[1] # SG2

    # Define labels for manual placement based on new request
    qdd_label = 'QD_D'
    qds_label = 'QD_S'
    left_outer_labels = [sg_at_pos2] # This is SG1
    left_inner_labels = ['GND', sg_at_pos1, 'GND', qdd_label] # This is the SG3 group + QD_D
    right_outer_label = sg_at_pos3 # This is SG2

    # Left side placement
    available_height_left = y_start_left_top - y_end_left_bottom
    if available_height_left > 0:
        x_outer_left = 0
        x_inner_left = w_side + row_spacing

        # Place left outer column (SG1)
        required_height_outer = len(left_outer_labels) * h_side + (len(left_outer_labels) - 1) * pad_spacing
        if required_height_outer <= available_height_left:
            y_offset_outer = y_end_left_bottom + (available_height_left - required_height_outer) / 2
            y_sg1 = y_offset_outer
            pad_placements.append({'label': left_outer_labels[0], 'position': (x_outer_left, y_sg1), 'dimensions': (w_side, h_side), 'rotation': 90})

        # Place left inner column (GND, SG3, GND, QD_D)
        required_height_inner = len(left_inner_labels) * h_side + (len(left_inner_labels) - 1) * pad_spacing
        if required_height_inner <= available_height_left:
            y_offset_inner = y_end_left_bottom + (available_height_left - required_height_inner) / 2
            for i, label in enumerate(left_inner_labels):
                y_pos = y_offset_inner + i * pitch_side
                pad_placements.append({'label': label, 'position': (x_inner_left, y_pos), 'dimensions': (w_side, h_side), 'rotation': 90})

    # Right side placement
    available_height_right = y_start_right_top - y_end_right_bottom
    if available_height_right > 0:
        # Place SG2 in the outer right column, centered in the available space
        y_sg2 = y_end_right_bottom + (available_height_right - h_side) / 2
        x_outer_right = layout_width - w_side
        pad_placements.append({'label': right_outer_label, 'position': (x_outer_right, y_sg2), 'dimensions': (w_side, h_side), 'rotation': 90})

        # Place QD_S in the inner right column, below SG2
        y_qds = y_sg2 - pitch_side
        x_inner_right = layout_width - 2*w_side - row_spacing
        pad_placements.append({'label': qds_label, 'position': (x_inner_right, y_qds), 'dimensions': (w_side, h_side), 'rotation': 90})

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

generate_quantum_dot_layout_v16_modified(N=20)