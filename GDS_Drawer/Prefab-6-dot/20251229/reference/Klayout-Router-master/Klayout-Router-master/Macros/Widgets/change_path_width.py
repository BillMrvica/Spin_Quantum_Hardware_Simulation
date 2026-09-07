"""
This script changes the width of all path shapes on specified layers in the current cell.
This operation will change the original shapes and is irreversible, so please use it with caution and make backup if necessary.
"""
import pya

# Warning! This step is irreversable!
layers_to_change = [
                   "800/0",
                   ]

original_width = 30  # in microns
new_width = 20  # in microns

# standard code to get the current layout and cell: 
# application instance -&gt; main window -&gt; current view -&gt; active cellview -&gt; layout and active cell
active_cellview = pya.Application.instance().main_window().current_view().active_cellview()
layout = active_cellview.layout()
cell = active_cellview.cell

original_width = int(original_width / layout.dbu + 0.5)
new_width = int(new_width / layout.dbu + 0.5)

for layer in layers_to_change:
    layer_index = layout.find_layer(pya.LayerInfo.from_string(layer))
    shapes = cell.shapes(layer_index)
    if layer_index is not None:
        for shape in shapes.each():
            if shape.is_path():
                path = shape.path
                if path.width == original_width:
                    path.width = new_width
                    shapes.replace(shape, path)
        print(f"Paths on Layer {layer} are changed from width {original_width*layout.dbu} um to {new_width*layout.dbu} um.")
    else:
        print(f"Layer {layer} not found.")
