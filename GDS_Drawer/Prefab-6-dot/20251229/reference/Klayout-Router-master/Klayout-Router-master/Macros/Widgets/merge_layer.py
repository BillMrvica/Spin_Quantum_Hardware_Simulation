import pya

# standard code to get the current layout and cell: 
# application instance -&gt; main window -&gt; current view -&gt; active cellview -&gt; layout and active cell
active_cellview = pya.Application.instance().main_window().current_view().active_cellview()
layout = active_cellview.layout()
cell = active_cellview.cell

# Warning! This step is irreversable!
layers_to_merge = [
                   "12/0",
                   #"3/1",
                   ]

for layer in layers_to_merge:
    layer_index = layout.find_layer(pya.LayerInfo.from_string(layer))
    if layer_index is not None:
        region = pya.Region(cell.shapes(layer_index)).merged()
        cell.shapes(layer_index).clear()
        cell.shapes(layer_index).insert(region)
        print(f"Shapes on Layer {layer} are merged.")
    else:
        print(f"Layer {layer} not found.")
