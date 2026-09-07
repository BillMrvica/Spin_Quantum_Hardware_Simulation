import gdstk
import matplotlib.pyplot as plt

# --- 1A. User-Customizable Pad Settings ---
pad_labels = ["100 nm", "150 nm", "200 nm", "250 nm", "300 nm"]
num_pads = len(pad_labels) 
pad_size = 100.0
sin_size = 120.0
extension_width = 30.0
extension_height = 80.0
pad_spacing = 200.0

# --- 1B. Label properties ---
text_size = 20.0
text_y_offset = -85.0

# --- 1C. Bottom Metal Bar Settings ---
bottom_metal_height = 50.0
bottom_metal_y_bottom = -50.0
bottom_metal_left_margin = 50.0
bottom_metal_right_margin = 100.0

# --- 1D. Right-Side Ground Pad Settings ---
right_pad_label = "Bottom Metal"
right_pad_size = (100.0, 100.0)
right_pad_spacing_from_last = 200.0

# --- 1E. Via Array Settings ---
via_size = (10.0, 10.0)
via_array_dims = (4, 4)
via_array_spacing = (20.0, 20.0)

# --- 2. Setup the GDS Library and Layers ---
library = gdstk.Library()
pad_layer = 1
sin_protection_layer = 2
text_layer = 3
bottom_metal_layer = 4
via_layer = 5

# --- 3. Create Reusable Cells ---
# Extended Pad Cell
extended_pad_cell = library.new_cell("EXTENDED_PAD")
extended_pad_cell.add(gdstk.rectangle((-sin_size / 2, -sin_size / 2), (sin_size / 2, sin_size / 2), layer=sin_protection_layer))
extended_pad_cell.add(gdstk.rectangle((-pad_size / 2, -pad_size / 2), (pad_size / 2, pad_size / 2), layer=pad_layer))
extended_pad_cell.add(gdstk.rectangle((-extension_width / 2, pad_size / 2), (extension_width / 2, pad_size / 2 + extension_height), layer=pad_layer))

# Single Via Cell
via_cell = library.new_cell("VIA_CELL")
via_cell.add(gdstk.rectangle((-via_size[0] / 2, -via_size[1] / 2), (via_size[0] / 2, via_size[1] / 2), layer=via_layer))

# --- 4. Create the Main Cell and place all components ---
main_cell = library.new_cell("FINAL_ASSEMBLY")

# Place the 6 extended pads and labels
for i in range(num_pads):
    pos_x = i * pad_spacing
    main_cell.add(gdstk.Reference(extended_pad_cell, (pos_x, 0)))
    main_cell.add(*gdstk.text(pad_labels[i], text_size, (pos_x, text_y_offset), layer=text_layer))

# main_cell.add(gdstk.Reference(extended_pad_cell, ((num_pads+1) * pad_spacing, 0)))
main_cell.add( gdstk.rectangle((-pad_size/2 + (num_pads)*pad_spacing, -pad_size / 2), (pad_size/2 + (num_pads)*pad_spacing, pad_size / 2), layer=pad_layer) )

# Create and place the Right-Side Ground Pad
last_pad_x = (num_pads - 1) * pad_spacing
right_pad_center_x = last_pad_x + right_pad_spacing_from_last
right_pad_center_y = 0

# Merge the bottom bar and pad body
start_x = -sin_size / 2 - bottom_metal_left_margin
end_x = (num_pads - 1) * pad_spacing + sin_size / 2 + bottom_metal_right_margin
bottom_metal_bar = gdstk.rectangle((start_x, bottom_metal_y_bottom), (end_x, bottom_metal_y_bottom + bottom_metal_height))
right_pad_body = gdstk.rectangle((right_pad_center_x - right_pad_size[0] / 2, right_pad_center_y - right_pad_size[1] / 2), (right_pad_center_x + right_pad_size[0] / 2, right_pad_center_y + right_pad_size[1] / 2))
combined_metal = gdstk.boolean(bottom_metal_bar, right_pad_body, 'or', layer=bottom_metal_layer)
main_cell.add(*combined_metal)

# Add SiN protection and label for the right pad
right_pad_sin_size = (right_pad_size[0] + 20, right_pad_size[1] + 20)
main_cell.add(gdstk.rectangle((right_pad_center_x - right_pad_sin_size[0] / 2, right_pad_center_y - right_pad_sin_size[1] / 2), (right_pad_center_x + right_pad_sin_size[0] / 2, right_pad_center_y + right_pad_sin_size[1] / 2), layer=sin_protection_layer))
main_cell.add(*gdstk.text(right_pad_label, text_size, (right_pad_center_x, text_y_offset), layer=text_layer))

# --- CORRECTED: Place the Via Array (Backwards-Compatible Method) ---
cols, rows = via_array_dims
spacing_x, spacing_y = via_array_spacing

# Calculate the starting position (origin) of the array to center it
array_width = (cols - 1) * spacing_x
array_height = (rows - 1) * spacing_y
via_array_origin_x = right_pad_center_x - array_width / 2
via_array_origin_y = right_pad_center_y - array_height / 2

# 1. Create a Repetition object
repetition = gdstk.Repetition(columns=cols, rows=rows, spacing=(spacing_x, spacing_y))

# 2. Create a Reference to a single via cell at the array's origin
via_array_ref = gdstk.Reference(via_cell, origin=(via_array_origin_x, via_array_origin_y))

# 3. Assign the repetition to the reference object
via_array_ref.repetition = repetition

# 4. Add the fully-configured reference to the main cell
main_cell.add(via_array_ref)


# --- 5. Visualization using Matplotlib ---
print("Generating plot of the final layout...")
fig, ax = plt.subplots()
layer_colors = {
    pad_layer: ('#8A2BE2', 0.9), sin_protection_layer: ('#FF69B4', 0.6),
    text_layer: ('black', 1.0), bottom_metal_layer: ('#00BFFF', 0.9),
    via_layer: ('#808080', 0.95)
}
all_polygons = main_cell.get_polygons()
for poly in all_polygons:
    layer = poly.layer
    if layer in layer_colors:
        color, alpha = layer_colors[layer]
        patch = plt.Polygon(poly.points, facecolor=color, alpha=alpha, edgecolor='black', linewidth=0.5)
        ax.add_patch(patch)

ax.set_title("Layout with Connected Ground Pad and Vias")
ax.set_xlabel("x-position (µm)"), ax.set_ylabel("y-position (µm)")
ax.grid(True, which='both', linestyle='--', linewidth=0.5)
ax.set_aspect('equal', adjustable='box'), ax.autoscale_view()
plt.show()
print("Plot generation complete.")

# --- 6. Optional: Save the final layout to a GDS file ---
output_file = "layout_with_metal_bar.gds"
library.write_gds(output_file)
print(f"\nSuccessfully created the GDS file: '{output_file}'")