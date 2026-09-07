import pya

# standard code to get the current layout and cell: 
# application instance -&gt; main window -&gt; current view -&gt; active cellview -&gt; layout and active cell
active_cellview = pya.Application.instance().main_window().current_view().active_cellview()
layout = active_cellview.layout()
cell = active_cellview.cell

# Warning! This step is irreversable!
previous_layers = [
                   "102/1",
                   #"4/0",
                   ]

new_layer = "111/0"

keep_previous_layers = False

for previous_layer in previous_layers:
    layer_index = layout.find_layer(pya.LayerInfo.from_string(previous_layer))
    if layer_index is not None:
        new_layer_index = layout.layer(pya.LayerInfo.from_string(new_layer))
        cell.shapes(new_layer_index).insert(cell.shapes(layer_index))
        if not keep_previous_layers: cell.shapes(layer_index).clear()
        print(f"Shapes on Layer {previous_layer} are moved to Layer {new_layer}.")
    else:
        print(f"Layer {previous_layer} not found.")
