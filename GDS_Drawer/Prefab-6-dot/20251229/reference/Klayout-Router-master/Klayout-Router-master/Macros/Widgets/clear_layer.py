import pya

# standard code to get the current layout and cell: 
# application instance -&gt; main window -&gt; current view -&gt; active cellview -&gt; layout and active cell
active_cellview = pya.Application.instance().main_window().current_view().active_cellview()
layout = active_cellview.layout()
cell = active_cellview.cell

# Warning! This step is irreversable!
layers_to_clear = [
                   # "1/2, test",
                   # "test",
                #    "102/0",
                   "5/1",
                   #"112/0"
                   #"Mesh_0",
                   #"Mesh_1",
                   #"Match"
                   ]

for layer in layers_to_clear:
    layer_index = layout.find_layer(pya.LayerInfo.from_string(layer))
    if layer_index is not None:
        cell.shapes(layer_index).clear()
        print(f"Layer {layer} cleared.")
    else:
        print(f"Layer {layer} not found.")
