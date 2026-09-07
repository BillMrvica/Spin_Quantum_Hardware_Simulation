import pya

# standard code to get the current layout and cell: 
# application instance -&gt; main window -&gt; current view -&gt; active cellview -&gt; layout and active cell
active_cellview = pya.Application.instance().main_window().current_view().active_cellview()
layout = active_cellview.layout()
cell = active_cellview.cell

# ===================== Warning! This step is irreversable! ===================== #

# This macro does this:
# Given a list of source_layers, a ref_layer and a dest_layer,
# for any shapes in source_layers overlapping with shapes in ref_layer,
# move the shapes from source_layers to dest_layer.

source_layers = [
                   "112/0",
                   #"4/0",
                   ]

ref_layer = "211/0"

dest_layer = "212/0"

keep_source_layers = False


# 预先获取目标图层的 Region 用于判断重叠
ref_layer_index = layout.find_layer(pya.LayerInfo.from_string(ref_layer))
if ref_layer_index is None:
    print(f"Reference Layer {ref_layer} not found!")
    exit()
ref_region = pya.Region(cell.shapes(ref_layer_index)).merge()

dest_layer_index = layout.layer(pya.LayerInfo.from_string(dest_layer))

for source_layer in source_layers:
    source_layer_index = layout.find_layer(pya.LayerInfo.from_string(source_layer))
    if source_layer_index is not None:
        source_shapes = cell.shapes(source_layer_index)
        # 对每个 shape 
        for source_shape in source_shapes.each():

            # 从 Shape 里取出具体几何对象
            if source_shape.is_box():
                geom = source_shape.box
            elif source_shape.is_polygon():
                geom = source_shape.polygon
            elif source_shape.is_path():
                geom = source_shape.path
            else:
                # 其他类型按需扩展，也可以先跳过
                continue

            # 用几何对象构造 Region
            source_region = pya.Region(geom)

            # 判断是否与目标图层重叠
            if not source_region.overlapping(ref_region).is_empty():
                # 重叠则移动
                cell.shapes(dest_layer_index).insert(source_shape)
                if not keep_source_layers:
                    cell.shapes(source_layer_index).erase(source_shape)
                    
        print(f"Shapes on Layer {source_layer} overlapping with {ref_layer} are moved to Layer {dest_layer}.")
    else:
        print(f"Layer {source_layer} not found.")

