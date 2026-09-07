# To-do

1. 需要确认一下 region.left，box.width()等参数都是整数; yes
2. 设置自动删除grid网格选项
3. 提示grid网格必须是空图层
4. 设置脚本自动merge
5. 将脚本安装进Klayout menu
6. 运行完GUI后提示一系列结果，包括在什么图层上插入了什么内容等
7. 算法复杂度分析



# Auto-routing

### To-do

1. 检查警告机制：
   - grid太密
2. 或者也可以设置档位：密，中等，稀疏

#### New `Map`

为了应对后期越来越复杂的rasterize的需求，以往的`Map`方法太过繁琐，而且没有利用到Klayout的各种算法。

现在提出一种比较吃内存但很高效的算法：

- 首先，对于一张`Map`，不仅维护 `map` 矩阵，还维护一个 `mesh` (a pya.Region)，专门用来 `rasterize`，这个 `mesh` 就是一个真实的网格，由一个 array of square region (pixel) 组成，可以直接被导入 Klayout 做可视化，但一般情况下它待在内存里。

  > 生成这个网格的方法是 instance 阵列，用一个 pixel cell 生成阵列后再 flatten

- 注意 `mesh.merged_semantics=False`

- 对于任何需要 `rasterize `的 `region`（注意现在可以直接以一整个region作为input）：

  - 直接使用 `mesh.overlapping(region)` 来获取所有和`region`有交的 `pixels`

  - delete these pixels from mesh

  - 对于每个pixel，设置它在Map里的cost为对应值

    

### Outline

##### `Krouter`: Take layers, Get & insert shapes and regions; Take um, give dbu.

|------------------- `Map`

|------------------- `Graph`

|------------------- `PinMatcher`

1. 单一起终点的寻路  ***Done*** 

   > 1. 网格化layout，`create map`；需要`map coordinate <-> layout coordinate`   ***Done***
   >
   > 2. `create graph from map`  ***Done***
   >
   >    > - 将所有可走的坐标作为node放入graph，每个node和四周8个node有edge，所有edge的cost设为1
   >    > - 需要 `graph node <-> map coordinate`
   >
   > 3. `find_shortest_path`, A star path searching.  ***Done***
   >
   >    > - 起、终点坐标 -> map coor
   >    > - map coor -> graph node
   >    > - routing(graph nodes), 返回nodes
   >    > - nodes -> map coor -> 起终点坐标
   >    > - 生成path
   >
   > 4. 将寻得的路径转化回map上的坐标  ***Done***
   >
   > 5. 根据坐标生成 pya.Path()  ***Done***
   >
   > `layout coordinate <-> map coordinate <-> graph node`

2. 多个起终点的匹配问题  ***Done***

   > 1. 首先需要添加`obstacle_regulation`，对不同图层的障碍采取不同的避让规则，例如100wf布线时障碍物需要pad但不需要bbx，1000wf布线需要pad和bbx，对于外围的大金属电极则需要采取200um以上的pad   ***Done***
   >
   > 2. 自动匹配起终点 `PinMatcher`
   >
   >    - 读取坐标  ***Done***
   >
   >      - 根据图层内的box读取中心坐标
   >
   >      - 根据图层内的Path读取末端坐标
   >      - 直接读取手动标注的坐标
   >
   >    - 匹配起终点   
   >
   >      使用 Hungarian algorithm 处理 组合优化-线性分配问题

3. 多个起终点寻路问题  `Krouter` ，见后面的详细论述  ***Done with problem dangling***

   > 或许需要重新做一下`Map`的架构，使之精简化，不再与layout做交互，只处理抽象的对象，所有和layout做的交互都放在另一个类里完成，或许可以叫 `Krouter`.  ***Done***

4. 使用 `cell` + `instance` 方式获得更高效的 `visualize` 和 `rasterize`

5. 多个100um写场内的起终点匹配+寻路

   ```python
   for box in cell.shapes():
   	starts_ends = match_starts_ends()
   	map = create_map()
   	graph = create_graph(map)
   	path_coor = find_shortest_path(graph)
   	path = pya.Path(path_coor)
   ```

6. 大写场和100um写场分开处理，分开寻路

7. 做一个animation demo





#### 多个起终点的匹配问题 

1. 读取坐标  ***Done***

   - 根据图层内的box读取中心坐标
   - 根据图层内的Path读取末端坐标
   - 直接读取手动标注的坐标
2. 匹配起终点
   - 构造基于真实路径长度的成本矩阵（只需要把之前的直线距离改为真实距离即可）

     ```python
     [[shortest_path_length(pin_a, pin_b) for pin_b in Pin_B] for pin_a in Pin_A]
     ```

   - 使用 Hungarian algorithm 处理 组合优化-线性分配问题





## Auto-routing 2.0 Outline

抛弃 Klayout-Macro development 开发工作流，搭建自己的环境，结合 `klayout` 和 `scikit-image` 完成寻路。



Map:

- initialized as 1
- obstacles as -1
- Path is -2



- 使用 skimage.graph.MCP.find_costs(start_i, ends) 搭建 starts-ends real cost matrix
- hungarian_match 得到 start-end 配对
- mcp.find_costs(start, end) .traceback(end)得到 start-end 路径



### 多个起终点寻路问题 `Krouter`   ***Done***
### 多个起终点寻路问题 `Krouter`   ***Done***

**Overwriting with condition** + **Path soft safe distance** + **Obs & Pin soft safe distance**

Principle: 任何 rasterize 不要覆盖 负数 区域

但其实只要我每次寻路时保证：**1.obs仍然是obs  2.其他pins不可见  3.先前的路径不可见**

采用 *obs & pin soft safe distance* 之后情况变得有些复杂，思路是 1.取消还原，把 `hardness` 作为 `map` 的一部分，2. 每次寻路时清理 Pins 周围的 `hardness` 确保寻路通畅 3. `obs_safe_distance` 不用太大，主要是为了确保路径存在。4. obs >= pins > path

0然后进行如下寻路操作：

0. 对于所有的`Obs`，在 Map 上rasterize为 -1，将其周围 `obs_safe_distance` 范围设置为 `obs_hardness`。

1. 将所有 `padded_pins_region` 中的 `1-pixels` 覆盖为 `pin_hardness`, 再将 `pins_region` 覆盖为 `-2`

   > `rasterize( pins_region.sized(pin_safe_distance), pin_hardness, cost==0, merge=True )`
   >
   > `rasterize( pins_region, -2, always True, True )`

   至此寻路基底搭好

2. 取出 `i-th matched pins`:

   1. 将 `i-th pins` 自身恢复覆盖为 `0` 。**这一步本身违反了不覆盖负数区域的原则，但是由于在整个过程中除了对自身的寻路外，这一区域一直为负数，不可能被覆盖，仅在做寻路时恢复为0**

      > `rasterize( ith_pins_retion, 0, True, True )`

   2. 清理 pins 周围区域：将`i-th pins` 周围 `safe_distance` 区域 且不与其他 `padded pins region` 重叠的区域 中的 `>1-pixels` 覆盖为 `0`，确保 `Dijkstra` 寻路算法的通畅与高效

      > ```python
      > temp_region = sized_ith_pin_region - sized_pins_region.not_in(sized_ith_pin_region)
      > rasterize( temp_region, 0, cost>1, True )
      > ```

   3. 做 pin_a 到 pin_b 的寻路

   4. 将路径本身置为不可见 `-3`

      > `rasterize( path, -3, always true, merge=True )`

   5. 将得到的路径做 pad 后，将区域内 `1-pixels` 覆盖为 `path hardness`

      > `rasterize( path.sized(path_safe_distance), path_hardness, cost==1, merge=True )`

   6. 将自身的 `>=0-pixels` 覆盖为 `-2` , 将 `sized i-th pins region` 中的 `>=1 <obs -pixels` 覆盖为 `pin_hardness`，

      > ```python
      > rasterize( ith_pins_region, -2, cost>=0, True )
      > rasterize( ith_pins_region.sized(pin_safe_distance), pin_hardness, obs_hardness>cost>=1, True )
      > ```



### 重写Rasterize方法  ***Done***

基于 `kdb.Region.rasterize()` 的原生函数重写 `rasterize`。核心问题是：如何将其覆盖到 `map` 上满足指定 `condition`的 `pixel` 上？尝试以下思路：

1. `np.array(region.rasterize()).nonzero()` 获取相对栅格坐标
2. `np.where(map[:,:]=feature)` 获取 `map` 对应区域的特定 `condition` 的相对栅格坐标
3. 二者取交集得到需要修改的 `pixels` 相对坐标，之后取 `origin shift` 得到真实栅格坐标，然后赋值





### Animation  ***Done***

- 必须生成.mp4，否则很难做调试
- 边界处会出现bug？ **Done**
- 前几帧为什么不显示？ **Done** *but there is a bug in FuncAnimation I did not understand. It is compensated by yield the first frame twice and set blit=False*.
- 必须先配合visulize方法确定自己的map设置没错



#### Save all frames instead of yield them?



### Self-adaptive path density map

1. 正常寻路

2. 记录下上一次寻路的所有 `kdb.Path`，并初始化 `Map.map` （即回归到第一次运行 `prepare_map_for_routing()` 后的状态）

3. 将所有 `kdb.Path` 膨胀后再 `kdb.Region.rasterize()` （以 `Map.total_area_bounding_box` 作为参数），之后再把整个数组的正数部分覆盖到 `Map.map` 的 `0-pixels`  

   > ```python
   > path_density = np.array(+Path.size(path_safe_distance).rasterize(origin, resol, nx, ny) )
   > path_density = path_density // (resol*resol) - 1
   > path_density[path_density > 0] *= 10
   > map[np.where((path_density > 0) & (map == 0))] = path_density[np.where((path_density > 0) & (map == 1))]
   >```
   > 

4. 







## Tricks on parameters

- pad of obs should be small, because it may consume pins
- obs safe distance can be large or medium
- 





# Log

#### 2025年7月5日18:46:34

意识到浮点数计算的精度问题可能导致 um to dbu 单位转化时出现误差，因此制定了一些原则：

- only variables with `_um` suffix are float (self.__um_per_dbu is also float and is the only exception);
- all other numbers are integers;
- All caculations are done in database units (dbu), and so the numbers are integers.
- 每次在把 um 单位的浮点数转化为 dbu 的整数时，都采用 `int( a_um * self.__dbu_per_unit + 0.5) `的方式

并且规定了，整个`Map`类里的 `layout coordinate system` 应该都是以 `dbu` 为单位的整数，包括所有的 `(x, y)`, `(X, Y)`. 如果是 um 单位的浮点数，就在后面加 `_um`，例如 `(x_um, y_um)`

#### 2025年7月9日20:56:47

**规定：** 所有的 `layer` 变量都应该是一个字符串 `"LayerNum/DataType"`, 例如 `"1/0"`

#### 2025年7月14日20:39:02

我发现了关于 `pya.Region` 做布尔运算时是否会自动 merge 的一些奇怪的现象，目前还没有找到确切的结论，但将现象总结如下：

- `|` 操作一定会merge
-  `&` 有时会merge，`region1 & region2` 结果不会merge `region1/2`， 但 `region1 & region1` 会导致 `region1` 的一次merge

这个问题出现在 `Map.collect_shapes_from_layer()`，如果我逐层读取:

 `for layer in layers: region |= collect_shapes_from(layer)` 就会出现 padded bounding boxes 重叠并merge的问题，box就不再像是个box形状了。解决方法：`region = collect_shapes_from_layers(layers)`，在函数内部pad之前就做 `|=` 的操作，pad之后就不会merge了。且注意到  `region.sized(pad) & total_area_bounding_box` 这一步的 `&` 不会导致merge。

#### 2025年7月30日12:50:22  *Problem encountered*

<img src=".assets/image-20250730125140446.png" alt="image-20250730125140446" style="zoom:33%;" />

如图中被选中的match所示，以两点之间直线距离作为权重做 `hungarian linear sum assignment` ，最后的结果虽然在直线连接上不会出现交叉，但是在真实 routing 时的路径上会出现交叉，因为直线连接跨过了 obstacles.

**解决方案：** 不采用直线距离作为权重而采用真实的路径距离，这就要求对任意两个 start-end 之间都做寻路。

#### 2025年8月3日20:44:31

1. Mesh太粗导致很多obs的padding基本失效，很多跨100um写场小路找不到，加上path safe distance的限制，先来的路很容易直接占据了狭窄区域的中心，导致后面的路被挡住
   - Scale分区域规划，核心区域使用精细meshing，外围区域使用粗糙mesh
   - 在狭窄区域中心设置高成本权重，利用dijkstra算法促使程序优先贴边行走
   - 设置路径边缘的高成本权重，促使所有路径彼此远离
   - 使用scikit-image包
2. Visualize卡顿，难以迅速得知rasterize结果
   - 使用cell+instance的方法
3. 生成的路线难以调整
   - Path -> Polygon -> Points -> simple path

#### 2025年8月4日22:46:58

今天第一次在S499上尝试了Auto-routing，遇到了很多问题：

1. 首先计算速度太慢，这是最主要的问题，这导致尝试成本很高，以后必须是秒画
2. 需要手动一个个修改匹配也很有问题
3. 使用逻辑不够通畅，经常不知道下一步应该做什么：需要图形界面
4. 不同图层之间非常容易搞混，需要一个标准化的图层规范，什么图层放什么内容，这需要大量实践去总结
5. 后来尝试在主机上运行，发现networkx需要重新安装，主机的Github SSH也没有弄
6. 用主机发现Klayout Console能够正常工作，不会出现未响应的情况，而且最耗时的就是Initialize GraphM和寻路，其他都是秒算，果然GraphM的For 循环是最慢的。或许以后可以把GraphM存起来？不，以后一定要用scikit-image来做
7. 主机的CPU比华为的CPU单核算力竟然还慢很多倍，我的天哪
8. mesh不够细，path起点和pin接触不太够
9. 需要自动识别障碍物
10. 需要仔细设置规则，例如这次的锯齿状通道完全可以设计规则绕过的，不需要手动调试。要有耐心地把所有obsrule都考虑到
    - 划痕？设置rect obs阵列限制通道
    - 边缘锯齿？
    - 脏东西？

#### 2025年9月10日20:54:28

初步实现了 scikit-image 的使用，刚刚提交 Commit c9938d23a9e6603438acfd5cc0a3600208049d1c，完成了从`NetworkX`到`scikit-image`的转换，解决了上一次实践中出现的绝大多数致命问题。

实践中发现的问题：

1. 很多路径被后来的soft safe distance吞掉了，可能导致路径交叉

2. soft safe distance还会吞掉一些pin以及obs，非常不好

   > 以上两条可以通过优化 rasterize 方法解决

3. 计算速度还有待提升，目前是1s一条路

   > 通过分块寻路解决，预计速度优化5倍以上

使用soft safe distance的优点：

1. 只要匹配得当，一定能找到所有的路径

#### 2025年9月11日20:57:43

PinMatcher的问题还是没有解决，采用了 `MCP_Geometric` 寻路后的真实路径距离作为 `cost matrix`，然后再做 `Hungarian match`，结果是一团糟，各种交叉。

尝试加权平均直线距离和真实距离之后也不是很理想，权重比达到 1e-20:(1-1e-20) 之后真实距离的作用才开始显现，但远在修正直线距离的交叉之前，真实距离就已经开始导致交叉了。

或许后续通过分块可以把真实距离的作用放大。

#### 2025年9月18日16:30:39

解决：改进 `rasterization` 并重新设计 `Krouter.route()` 算法，把 3mm*5mm上 50 对 pins 的寻路时间从 50s 压缩到 5s，解决了上上次log中的路径交叉问题

问题：接下来尝试在 padded obs周围再添加一层 soft safe distance, 进一步避免无法找到路径的问题。

尝试 pre-routing self-adaptive path density map 解决路径靠得太近的问题。

#### 2025年9月19日17:36:46

解决：对 obs 和 pins 都采用了 **soft-safe-distance** 特征，彻底解决了路径被堵住的问题。

#### 2025年9月22日20:20:29

今天发现一个边界上的bug：如果 `total_area_bounding_box` 的尺寸恰好被 `resolution` 整除，那么我可以选择在最右侧和最上侧包一圈 `pixel`，也可以选择不包，但无论怎么选择，总是至少出现如下问题中的一个：

- 如果包，则 `rasterize` 无法触碰到最右、最上侧的一圈 `pixel`。（因为 clip）
- 如果不包，`layout_coor_to_map_coor` 对于右边界和上边界上的点会发生索引溢出，因为每个 `pixel` 的有效范围被认为是包括左下边而不包括右上边。

最好的解决方法是保证不被整除，可以选择修改分辨率或者修改 `total area bounding box` 的边界。这个后面再考虑，目前先不作处理，且选择包的情况，因为第一种bug出现的概率较小。

#### 2025年9月24日16:17:14

刚才尝试把 `Map.map` 的默认初始值改为 `0`，结果大失败，因为 `skimage.graph` 对全 0 的图像寻路似乎会出bug。最开始设计算法的时候我其实考虑到了这一点，幸运地选取了 `1` 作为初始值。Git 在这次尝试中起了重要的作用。做这个尝试主要是为了考虑 *self-adaptive path density map* 需要一个全 0 的 `map`，这可以用其他方法实现。

#### 2025年9月25日14:56:45

初次尝试 `self-adaptive routing` 发现了一些问题，目前性能还不如普通寻路

- 当前简单的 `path_safe_distance` 和 `path_hardness` 不够用，必须考虑 `distance-related hardness`，距离越近 `cost` 越高才行，否则会出现 `safe_distance` 失效问题，例如两条比较靠近的路的中间区域其实相当于没有 `safe_distance`
- 自适应寻路的算法需要改进
- 重新设计 `Animation` ？

#### 2025年9月25日23:03:06

- 尝试了 **damping soft path safe distance**，效果出奇得好！
- 接下来尝试对外接 Pad 做顺时针排序，我相信这基本能够解决 *path-path distance too close* 的问题，或许不需要做 *self-adaptive routing*.

#### 2025年9月26日21:08:02

- 尝试 **clockwise sorted pins**，调参后获得极好效果

#### 2025年9月28日17:42:59

- 尝试 **damping soft safe distance for obstacles and pins**，效果再次得到改进，而且意外地发现不再需要对 `pin_region` 做 recover 了：`damping raster` 的算法设计刚好可以把这种结构保留的很好，它非常优雅且鲁棒地实现了对所有其他 pins 都保持安全距离的这一目的。





# Ideas

### Scikit-image.graph

Abandon `NetworkX` but embrace `scikit-image.graph`, fully utilizing the power of `numpy.ndarray`. **1 hour -> 40s**

### Rasterization

`Klayout` native engine `klayout.db.Region.rasterize()` based `Map.rasterize_region()` method. Perfectly utilizing the power of `klayout` C++ engine and `numpy.ndarray` indexing to fulfill **conditionally overwriting** super super **fast**.  **40s -> 4s**

### Damping Soft path safe distance

Partly solving *no-path-found* problem. Damping is really useful.

### Pre-routed self-adaptive path density map

Fully dilute the density of paths and keep path away from each other as far as possible.

### Soft obs & pin safe distance

Absolute safe distance is too hard and always causes *no-path-found* problem. Now we add soft distance around it to alleviate this problem.



# Dictionary

- `layer`: A str,  `"1/0"`
- `map`: A numpy matrix, available cell is 1, inaccessable cell (obstacle cell) is 0. This cooresponding to a meshed grid from layout geometry.
- `pixel`
- `xxx_polygon`：A `pya.Region` that contains a single polygon
- `xxx_polygons`: A `list` of single polygons (which is a `pya.Region`)
- `region`: A `pya.Region` that may contain multiple polygons, merged or unmerged. This means I tend to use the original meaning of `pya.Region`.

Auto-routing/PinMatcher.py



Implemented a matcher for multiple pins in the layout.

  \- Get pins set A and B from shape bbox center, text or end of path.

  \- Match pins using a greedy approach or Hungarian algorithm.

​    \- Greedy approach does not guarantee no-crossing of pin-pin paths;

​    \- Hungarian algorithm finds the optimal match (minimum total distance) and guarantees no-crossing.

# Errors encountered during development

## pya.Region()

#### What can be "Region"ed?

```python
for shape in cell.shapes(box_layer_index):
    pya.Region(shape) # wrong! "shape" cannot be converted to "Region"
    pya.Region(shape.box()) # wrong! "box" is no callable.
    pya.Region(shape.box) # correct
```

1. shape cannot be Regioned
2. only box/polygon etc. objects can be Regioned.



## How to extract Edges?

```python
all_edges = pya.Edges() # empty object
for shape in cell.shapes(box_layer_index):
    all_edges.insert(shape.edges()) # wrong! shape has no attribute .edges()
    all_edges.insert(shape.edge) # wrong! shape.edge cannot be inserted into pya.Edges()
    
```

#### correct way:

```python
all_edges = pya.Edges() # empty object
for shape in cell.shapes(box_layer_index):
    temp_polygon = shape.polygon # turn any shapes to polygon
    shape_region = pya.Region(temp_polygon) # turn polygon into Region
    all_edges.insert(shape_region.edges()) # get edges from Region.edges()
```



## Tips

1. `Region.bbox()` 对整个region生成bbox，而不是region里的每个子region

2. `Region.sized()` 对整个region的子region扩大

3. `pya.Region.each()` 返回的列表里的元素不是 `pya.Region` 而是 `pya.PolygonWithProperties`，这个对象不能做 `pya.Region` 的各种布尔操作，因为有多种方法可以它转化为 `pya.Region` 而导致报错。

   > RuntimeError: There are multiple conversion constructors available to convert object of type PolygonWithProperties to type Region for argument #1 ('other') in Region.__ior__

4. `pya.Region` 的布尔操作可以应用于一般的对象，例如 `pya.Box`，运算符会自动把它变成 `pya.Region` 后再做 bool。但是对于上面那种情况，`pya.PolygonWithProperties` 有多种 conversion constructors 而导致自动转化过程报错。

5. `pya.Region.each()` 虽然不返回region，但是他返回的对象 `pya.PolygonWithProperties` 也有 `.bbox()` 

6. `|` 操作一定会merge. `&` 有时会merge，`region1 & region2` 结果不会merge `region1/2`， 但 `region1 & region1` 会导致 `region1` 的一次merge

7. `layout.layer(pya.LayerInfo.from_sting("1/0"))` 可以直接返回存在的 `layer_index` 或者创造不存在的 `layer`

8. **`pya.Region` 就是强化版的 `pya.Polygon` 或者是强化版的多个 `pya.Polygon`，所以我从来不用 `pya.Polygon` 这个类，而是全部用 `pya.Region`；当我要强调某个对象是一个单独的 `polygon`，我将它命名为`xxx_polygon` 而不是 `xxx_region`**





## Useful methods

|           | Return    | Method                                                       | Args                                                         | Description                                                  |
| --------- | --------- | ------------------------------------------------------------ | :----------------------------------------------------------- | ------------------------------------------------------------ |
| *[const]* | Region    | **[inside](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method124)** | (const Region other)                                         | Returns the polygons of this region which are completely inside polygons from the other region |
| *[const]* | EdgePairs | **[inside_check](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method125)** | (const Region other, unsigned int d, bool whole_edges = false, Metrics metrics = Euclidian, variant ignore_angle = default, variant min_projection = 0, variant max_projection = max, bool shielded = true, Region::OppositeFilter opposite_filter = NoOppositeFilter, Region::RectFilter rect_filter = NoRectFilter, bool negative = false, PropertyConstraint property_constraint = IgnoreProperties, ZeroDistanceMode zero_distance_mode = IncludeZeroDistanceWhenTouching) | Performs an inside check with options                        |
| *[const]* | Region    | in_                                                          | (const Region other)                                         | Returns all polygons which are members of the other region   |
| *[const]* | Region[]  | **[in_and_out](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method106)** | (const Region other)                                         | Returns all polygons which are members and not members of the other region |
| *[const]* | Region    | **[extents](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method90)** |                                                              | Returns a region with the bounding boxes of the polygons     |
| *[const]* | Region    | **[extents](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method91)** | (int d)                                                      | Returns a region with the enlarged bounding boxes of the polygons |
| *[const]* | Region    | **[extents](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method92)** | (int dx, int dy)                                             | Returns a region with the enlarged bounding boxes of the polygons |

| *[const]* | bool | **[is_box?](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method129)** |      | Returns true, if the region is a simple box |
| --------- | ---- | ------------------------------------------------------------ | ---- | ------------------------------------------- |
| *[const]* | bool | **[is_empty?](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method132)** |      | Returns true if the region is empty         |
| *[const]* | bool | **[is_merged?](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method133)** |      | Returns true if the region is merged        |



| *[const]* | Region | **[join](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method135)** | (const Region other) | Returns the combined region of self and the other region |
| --------- | ------ | ------------------------------------------------------------ | -------------------- | -------------------------------------------------------- |
|           | Region | **[join_with](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method136)** | (const Region other) | Adds the polygons of the other region to self            |

| *[const]* | Region | **[rectangles](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method213)** |      | Returns all polygons which are rectangles |
| --------- | ------ | ------------------------------------------------------------ | ---- | ----------------------------------------- |
| *[const]* | Region | **[squares](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method290)** |      | Returns all polygons which are squares    |

| *[const]* | Region | **[non_rectangles](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method180)** |      | Returns all polygons which are not rectangles  |
| --------- | ------ | ------------------------------------------------------------ | ---- | ---------------------------------------------- |
| *[const]* | Region | **[non_rectilinear](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method181)** |      | Returns all polygons which are not rectilinear |
| *[const]* | Region | **[non_squares](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method182)** |      | Returns all polygons which are not squares     |

|           | void   | **[merged_semantics=](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method145)** | (bool f)                                                     | Enables or disables merged semantics                         |
| --------- | ------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| *[const]* | bool   | **[merged_semantics?](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method146)** |                                                              | Gets a flag indicating whether merged semantics is enabled   |
| *[const]* | Region | **[overlapping](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method199)** | (const Region other, unsigned long min_count = 1, unsigned long max_count = unlimited) | Returns the polygons of this region which overlap polygons from the other region |

| *[const]* | Region     | **[pull_inside](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method206)** | (const Region other)                                         | Returns all polygons of "other" which are inside polygons of this region |
| --------- | ---------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| *[const]* | Region     | **[pull_interacting](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method207)** | (const Region other)                                         | Returns all polygons of "other" which are interacting with (overlapping, touching) polygons of this region |
| *[const]* | Edges      | **[pull_interacting](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method208)** | (const Edges other)                                          | Returns all edges of "other" which are interacting with polygons of this region |
| *[const]* | Texts      | **[pull_interacting](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method209)** | (const Texts other)                                          | Returns all texts of "other" which are interacting with polygons of this region |
| *[const]* | Region     | **[pull_overlapping](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method210)** | (const Region other)                                         | Returns all polygons of "other" which are overlapping polygons of this region |
| *[const]* | double[][] | **[rasterize](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method211)** | (const Point origin, const Vector pixel_size, unsigned int nx, unsigned int ny) | A grayscale rasterizer delivering the area covered per pixel |
| *[const]* | double[][] | **[rasterize](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method212)** | (const Point origin, const Vector pixel_distance, const Vector pixel_size, unsigned int nx, unsigned int ny) | A version of 'rasterize' that allows a pixel step distance which is larger than the pixel size |

| *[const]* | EdgePairs | **[grid_check](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method100)** | (int gx, int gy) | Returns a marker for all vertices not being on the given grid |
| :-------- | --------- | ------------------------------------------------------------ | ---------------- | ------------------------------------------------------------ |
|           |           |                                                              |                  |                                                              |

|           | void   | **[filter](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method96)** | (const PolygonFilterBase ptr filter) | Applies a generic filter in place (replacing the polygons from the Region) |
| --------- | ------ | ------------------------------------------------------------ | ------------------------------------ | ------------------------------------------------------------ |
|           | void   | **[filter_properties](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method97)** | (variant[] keys)                     | Filters properties by certain keys.                          |
| *[const]* | Region | **[filtered](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method98)** | (const PolygonFilterBase ptr filter) | Applies a generic filter and returns a filtered copy         |
|           | Region | **[flatten](https://www.klayout.org/downloads/master/doc-qt4/code/class_Region.html#method99)** |                                      | Explicitly flattens a region                                 |