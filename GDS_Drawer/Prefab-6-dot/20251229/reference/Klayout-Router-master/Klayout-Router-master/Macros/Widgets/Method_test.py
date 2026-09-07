import pya
# standard code to get the current layout and cell: 
# application instance -&amp;amp;gt; main window -&amp;amp;gt; current view -&amp;amp;gt; active cellview -&amp;amp;gt; layout and active cell
active_cellview = pya.Application.instance().main_window().current_view().active_cellview()
layout = active_cellview.layout()
cell = active_cellview.cell


# # test: What is the type of pya.Shapes.each()? 
# a = cell.shapes(layout.layer(pya.LayerInfo.from_string("1/0")))
# print(type(a.each()))  # <class 'pya._Iterator'>
# print(type(list(a.each()))) # <class 'list'>
# for i in list(a.each()):
#     print(type(i))  # <class 'pya.Shape'>

# # test: How does pya.Region(shapes, "*") work?
# # Answer: It does not work on pya.Shapes, only Recursive shape iterator.
# a = cell.shapes(layout.layer(pya.LayerInfo.from_string("13/0")))
# b = pya.Region(a, "*")
# b = pya.Region(a.each(), "*")
# print(b)
# # cell.shapes(layout.layer("55/0")).insert(b)



# test: can pya.Region() take a point? Answer: No!
# a = pya.Region(pya.Point(0, 0))
# a = pya.Region(pya.Box(-100, -100, 100, 100)).pull_interacting(pya.Point(0, 0))


# # test: pya.Region.inside()
# # Answer: It works as described, and even more, it can take a pya.Box as an argument.
# box = pya.Box(2050*1000, 1550*1000, 2150*1000, 1650*1000)
# layer_idx = layout.layer(pya.LayerInfo.from_string("1/0"))
# region = pya.Region(cell.shapes(layer_idx))
# region = region.inside(box)
# # cell.shapes(layout.layer("test")).insert(box)
# cell.shapes(layout.layer("test")).insert(region)




# # test pya.Region & pya.Box. Answer: It works.
# region = pya.Region(pya.Box(0, 0, 1000, 1000))
# box = pya.Box(0, 0, 500, 1000)
# region &= box
# cell.shapes(layout.layer("test")).insert(region)



# # test pya.Region.extents(pad)
# # Answer: It works as described.
# # Also, I tested pya.LayerInfo.from_string(), it can takes name like "test", or layer/datatype like"1/0", or "1/0, test"
# # and for "1/0, test", it will return the layer with name "test" and layer 1 and datatype 0, not the layer with name "1/0, test"
# layer_idx = layout.layer(pya.LayerInfo.from_string("1/0"))
# # region = pya.Region(cell.shapes(layer_idx)).extents()
# region = pya.Region(cell.shapes(layer_idx)).extents(5000)
# # layer_idx = layout.layer("test")
# layer_idx = layout.layer(1, 2, "test")
# # layer_idx = layout.layer("1/2, test")
# cell.shapes(layer_idx).insert(region)


# # test anonymous layer. Answer: Yes, it works.
# layer_mesh_pixel = layout.layer()
# print(layer_mesh_pixel)
# print(type(layer_mesh_pixel))
# print(layout.get_info(layer_mesh_pixel))


# # test: What if the LayerInfo for a layer with a name?
# # Answer: In general, a layer has either a layer and a datatype number (in GDS2), a name (for example in DXF or CIF) 
# # or both (in OASIS). In the latter case, the primary identification is through layer and datatype number and the name 
# # is some annotation attached to it. A LayerInfo object which specifies just a name returns true on is_named?. The LayerInfo 
# # object can also specify an anonymous layer (use LayerInfo#new without arguments). Such a layer will not be stored when 
# # saving the layout. They can be employed for temporary layers for example. Use LayerInfo#anonymous? to test whether a layer 
# # does not have a specification.
# layer_mesh_pixel = layout.layer(101,1,"Mesh Pixel")
# print(layer_mesh_pixel)
# print(type(layer_mesh_pixel))
# print(layout.get_info(layer_mesh_pixel))



# test: What does pya.Region.rasterize() do?
# a = pya.Region(cell.shapes(layout.find_layer(pya.LayerInfo.from_string("3/0"))))
# b = a.bbox()
# c = a.rasterize(pya.Point(b.left, b.bottom), pya.Vector(1000, 1000), 10, 1)
# print(c)
# print(type(c))
# cell.shapes(layout.layer(pya.LayerInfo.from_string("4/0"))).insert(c)


# # test: difference between |= and join()
# # Will |= merge? Yes
# # Will join() merge? No
# a = pya.Region(cell.shapes(layout.find_layer(pya.LayerInfo.from_string("3/0"))))
# b = pya.Region(cell.shapes(layout.find_layer(pya.LayerInfo.from_string("4/0"))))
# # a |= b
# # cell.shapes(layout.layer(pya.LayerInfo.from_string("5/0"))).insert(a)

# a.join_with(b)
# cell.shapes(layout.layer(pya.LayerInfo.from_string("5/0"))).insert(a)


# # test: if merged_semantics = False makes |= not merge 
# # No! Even if both a and b have merged_semantics = False, |= will still merge them.
# a = pya.Region(cell.shapes(layout.find_layer(pya.LayerInfo.from_string("3/0"))))
# b = pya.Region(cell.shapes(layout.find_layer(pya.LayerInfo.from_string("4/0"))))
# a.merged_semantics = False
# b.merged_semantics = False
# a |= b
# cell.shapes(layout.layer(pya.LayerInfo.from_string("5/0"))).insert(a)


# # test whether merged_semantics = False affects shapes inserted. Answer: No
# a = pya.Region(cell.shapes(layout.find_layer(pya.LayerInfo.from_string("3/0"))))
# cell.shapes(layout.layer(pya.LayerInfo.from_string("4/0"))).insert(a)
# print(a.merged_semantics)
# print(a.count())
# a.merged_semantics = False
# print(a.merged_semantics)
# print(a.count())
# cell.shapes(layout.layer(pya.LayerInfo.from_string("5/0"))).insert(a)

# # test whether merged_semantics = False affects edges generated. 
# # Answer: Yes
# a = pya.Region(cell.shapes(layout.find_layer(pya.LayerInfo.from_string("3/0"))))
# a.merged_semantics = False
# edges = a.edges()
# cell.shapes(layout.layer(pya.LayerInfo.from_string("4/0"))).insert(edges.extended(0, 0, 100, 100, True))

# points = []
# for edge in edges:
#     points.append(edge.p1)
#     points.append(edge.p2)
# cell.shapes(layout.find_layer(pya.LayerInfo.from_string("4/0"))).insert(pya.Path(points, width=500,))
