import gdstk
import matplotlib.pyplot as plt

# --- 1. Global Settings & Layer Definitions ---

# --- 1A. Pad and Extension Dimensions (µm) ---
pad_size = 100.0
sin_size = 120.0
pad_spacing = 200.0
# maskless litography
extension_1_width = 20.0
extension_1_height = 200.0
extension_2_width = 10.0
extension_2_height = 100.0
extension_3_width = 3.0
extension_3_height = 100.0
# EBL layers
extension_4_width = 1.0
extension_4_height = 40.0
extension_5_width = 0.5
extension_5_height = 5.0
extension_6_width = 0.1
extension_6_height = 5.2

# --- 1B. Label, Metal Bar, and Ground Pad Settings ---
text_size = 20.0
text_y_offset = -85.0
bottom_metal_height = 50.0
bottom_metal_y_bottom = -50.0
bottom_metal_left_margin = 50.0
bottom_metal_right_margin = 100.0
right_pad_label = "Bottom Metal"
right_pad_size = (100.0, 100.0)
right_pad_spacing_from_last = 200.0

# --- 1C. Via Array Settings ---
via_size = (10.0, 10.0)
via_array_dims = (4, 4)
via_array_spacing = (20.0, 20.0)

# --- 1D. Final Layout Assembly Settings ---
vertical_spacing = 1000.0

# --- 1E. GDS Layer Definitions ---
pad_layer = 1
sin_protection_layer = 2
text_layer = 3
bottom_metal_layer = 4
via_layer = 5
EBL_metal_layer_1_1 = 6
EBL_metal_layer_1_2 = 7
EBL_metal_layer_1_3 = 8
EBL_metal_layer_2_1 = 9
EBL_metal_layer_2_2 = 10
EBL_metal_layer_2_3 = 11

# --- 2. Core Reusable Function ---

def create_single_row_assembly(library, cell_name_prefix, pad_labels, EBL_layer1, EBL_layer2, EBL_layer3, mismatch_6):
    """
    Creates a GDSII cell containing a full row of extended pads and ground connections.
    All cells created within this function will be prefixed with cell_name_prefix to ensure they are unique.
    """
    num_pads = len(pad_labels)

    # --- Create a reusable cell for a single via with a UNIQUE name ---
    via_cell = library.new_cell(f"{cell_name_prefix}_VIA_CELL")
    via_cell.add(gdstk.rectangle((-via_size[0] / 2, -via_size[1] / 2),
                               (via_size[0] / 2, via_size[1] / 2), layer=via_layer))

    # --- Create the main cell for the entire row with a UNIQUE name ---
    row_cell = library.new_cell(f"{cell_name_prefix}_ROW_ASSEMBLY")

    # --- Place the extended pads and their labels ---
    for i in range(num_pads):
        # Create a unique pad cell for each pad
        pad_cell_name = f"{cell_name_prefix}_EXTENDED_PAD_{i}"
        extended_pad_cell = library.new_cell(pad_cell_name)
        
        # Add common parts
        extended_pad_cell.add(gdstk.rectangle((-sin_size / 2, -sin_size / 2), (sin_size / 2, sin_size / 2), layer=sin_protection_layer))
        extended_pad_cell.add(gdstk.rectangle((-pad_size / 2, -pad_size / 2), (pad_size / 2, pad_size / 2), layer=pad_layer))
        
        # Add extensions 1, 2, and 3 (fixed layer)
        ext1_y_end = pad_size / 2 + extension_1_height
        extended_pad_cell.add(gdstk.rectangle((-extension_1_width / 2, pad_size / 2), (extension_1_width / 2, ext1_y_end), layer=pad_layer))
        ext2_y_end = ext1_y_end + extension_2_height
        extended_pad_cell.add(gdstk.rectangle((-extension_2_width / 2, ext1_y_end), (extension_2_width / 2, ext2_y_end), layer=pad_layer))
        ext3_y_end = ext2_y_end + extension_3_height
        extended_pad_cell.add(gdstk.rectangle((-extension_3_width / 2, ext2_y_end), (extension_3_width / 2, ext3_y_end), layer=pad_layer))
        
        # Add EBL extensions 4, 5 and 6
        ext4_y_end = ext3_y_end + extension_4_height
        extended_pad_cell.add(gdstk.rectangle((-extension_4_width / 2, ext3_y_end), (extension_4_width / 2, ext4_y_end), layer=EBL_layer1))
        ext5_y_end = ext4_y_end + extension_5_height
        extended_pad_cell.add(gdstk.rectangle((-extension_5_width / 2, ext4_y_end), (extension_5_width / 2, ext5_y_end), layer=EBL_layer2))
        ext6_y_end = ext5_y_end + extension_6_height
        extended_pad_cell.add(gdstk.rectangle((-extension_6_width/2 + mismatch_6, ext5_y_end), (extension_6_width/2 + mismatch_6, ext6_y_end), layer=EBL_layer3))

        # Add reference and label
        pos_x = i * pad_spacing
        row_cell.add(gdstk.Reference(extended_pad_cell, (pos_x, 0)))
        row_cell.add(*gdstk.text(pad_labels[i], text_size, (pos_x, text_y_offset), layer=text_layer))

    # Add extra pad and ground connections...
    row_cell.add(gdstk.rectangle((-pad_size/2 + num_pads*pad_spacing, -pad_size/2), (pad_size/2 + num_pads*pad_spacing, pad_size/2), layer=pad_layer))
    last_pad_x = (num_pads - 1) * pad_spacing
    right_pad_center_x = last_pad_x + right_pad_spacing_from_last
    start_x = -sin_size / 2 - bottom_metal_left_margin
    end_x = (num_pads - 1) * pad_spacing + sin_size / 2 + bottom_metal_right_margin
    bottom_metal_bar = gdstk.rectangle((start_x, bottom_metal_y_bottom), (end_x, bottom_metal_y_bottom + bottom_metal_height))
    right_pad_body = gdstk.rectangle((right_pad_center_x - right_pad_size[0]/2, -right_pad_size[1]/2), (right_pad_center_x + right_pad_size[0]/2, right_pad_size[1]/2))
    combined_metal = gdstk.boolean(bottom_metal_bar, right_pad_body, 'or', layer=bottom_metal_layer)
    row_cell.add(*combined_metal)
    right_pad_sin_size = (right_pad_size[0] + 20, right_pad_size[1] + 20)
    row_cell.add(gdstk.rectangle((right_pad_center_x - right_pad_sin_size[0]/2, -right_pad_sin_size[1]/2), (right_pad_center_x + right_pad_sin_size[0]/2, right_pad_sin_size[1]/2), layer=sin_protection_layer))
    row_cell.add(*gdstk.text(right_pad_label, text_size, (right_pad_center_x, text_y_offset), layer=text_layer))
    cols, rows = via_array_dims
    spacing_x, spacing_y = via_array_spacing
    array_width, array_height = (cols - 1) * spacing_x, (rows - 1) * spacing_y
    via_array_origin = (right_pad_center_x - array_width/2, -array_height/2)
    repetition = gdstk.Repetition(columns=cols, rows=rows, spacing=(spacing_x, spacing_y))
    via_array_ref = gdstk.Reference(via_cell, origin=via_array_origin)
    via_array_ref.repetition = repetition
    row_cell.add(via_array_ref)
    
    return row_cell

# --- 3. Main Script Execution ---
if __name__ == "__main__":
    library = gdstk.Library()
    pad_labels = ["100 nm", "150 nm", "200 nm", "250 nm", "300 nm"]

    # --- Call the function with UNIQUE prefixes to create the single row cells ---
    single_row_cell_1 = create_single_row_assembly(library, "ROW_1", pad_labels, EBL_metal_layer_1_1, EBL_metal_layer_1_2, EBL_metal_layer_1_3, 0.05)
    single_row_cell_2 = create_single_row_assembly(library, "ROW_2", pad_labels, EBL_metal_layer_2_1, EBL_metal_layer_2_2, EBL_metal_layer_2_3, -0.05)

    # --- Create the FINAL TOP-LEVEL CELL with two rows ---
    main_cell = library.new_cell("FINAL_LAYOUT")
    main_cell.add(gdstk.Reference(single_row_cell_1, (0, 0)))
    main_cell.add(gdstk.Reference(single_row_cell_2, (0, vertical_spacing), x_reflection=True))

    # --- Visualization using Matplotlib ---
    print("Generating plot of the final layout...")
    fig, ax = plt.subplots(figsize=(16, 12))
    layer_colors = {
        pad_layer: ('#8A2BE2', 0.9), sin_protection_layer: ('#FF69B4', 0.6),
        text_layer: ('black', 1.0), bottom_metal_layer: ('#00BFFF', 0.9),
        via_layer: ('#808080', 0.95),
        EBL_metal_layer_1_1: ('#32CD32', 0.9), EBL_metal_layer_1_2: ('#228B22', 0.9),
        EBL_metal_layer_1_3: ('#98FB98', 0.9), EBL_metal_layer_2_1: ('#DC143C', 0.9),
        EBL_metal_layer_2_2: ('#B22222', 0.9), EBL_metal_layer_2_3: ('#F08080', 0.9)
    }
    all_polygons = main_cell.get_polygons()
    for poly in all_polygons:
        layer = poly.layer
        if layer in layer_colors:
            color, alpha = layer_colors[layer]
            patch = plt.Polygon(poly.points, facecolor=color, alpha=alpha, edgecolor='black', linewidth=0.5)
            ax.add_patch(patch)

    ax.set_title("Final Layout with Mirrored Top Copy")
    ax.set_xlabel("x-position (µm)"), ax.set_ylabel("y-position (µm)")
    ax.grid(True, which='both', linestyle='--'), ax.set_aspect('equal', 'box'), ax.autoscale_view()
    plt.show()
    print("Plot generation complete.")

    # --- Save the final layout to a GDS file ---
    # output_file = "New_Test_Stru_20250918_fixed.gds"
    # library.write_gds(output_file)
    # print(f"\nSuccessfully created the GDS file: '{output_file}'")