<h1 align="center">Klayout-Router</h1>

<p align="center">
中文 | <a href="./Readme.md">English</a>
</p>
**Klayout-Router** 是一套用于在凝聚态物理（CMP）实验中自动化纳米图形设计工作的工具集。

> Klayout-Router 这个名字强调了布线这一主要功能，但本项目还提供了一套辅助工具（包括写场补丁拼接、图层操作等宏），都放在 `./Macros/` 目录下。

与工业级电路设计中高度标准化的版图不同，实验室级别的纳米器件设计往往高度客制化：不同芯片之间的器件几何形状、引线走线方式都可能差别很大。对于研究人员来说，在快速迭代工艺配方的过程中，如果每次都要手动画图、布线，将会极其耗时。**Klayout-Router** 旨在通过自动化布线、补丁生成、缺陷勾边等工作，为你节省大量重复劳动时间。

# 🛠 项目结构

```text
Klayout-Router
├── Auto-routing/           # 核心自动布线算法
│   ├── Krouter.py          
│   ├── Map.py              	
│   ├── PinMatcher.py      	
|   ├── routing.ipynb			# 布线功能用户界面
|   └── requirements.txt        # pip 依赖列表
├── Macros/                 # KLayout 宏（在 KLayout GUI 中运行的 Python 脚本）
│   ├── Auto-patching/     	# 自动补丁宏
│   └── Widgets/           	# 各类小工具（图层操作等）
└── Readme.md
```

# 目录

- [如何使用 Router？](#如何使用-router)
  - [环境要求](#环境要求)
  - [功能特性](#功能特性)
  - [🔥即刻启用](#🔥即刻启用)
- [如何使用其他 Klayout 宏？](#如何使用其他-klayout-宏)
  - [⚠️警告](#⚠️警告)
  - [环境要求](#环境要求-1)
  - [什么是 Klayout 宏？](#什么是-klayout-宏)
  - [功能一览](#功能一览)
- [许可证](#许可证)

# 如何使用 Router？

## 环境要求

#### 系统

项目是在 Windows 上开发和测试的，不保证 MacOS 和 Linux.

#### 运行环境

运行 Router 需要 **独立的 python/conda 环境**，与 *Klayout for Windows* 自带的 Python 环境相互独立。因此：

- 运行 Router 本身 **不需要** 安装 *Klayout for Windows*；
- 但某些中间过程仍然需要 *Klayout for Windows* 实时查看 `.gds` 文件进行辅助。

1. 创建 conda 环境

   ```bash
   conda create -n klayout_dev python=3.10 -y
   ```

2. 安装依赖

   ```bash
   conda activate klayout_dev
   cd Klayout-Router/Auto-routing/
   pip install -r requirements.txt
   ```

3. 注册 ipykernel（用于运行 `.ipynb`）

   ```bash
   python -m ipykernel install --user --name klayout_dev --display-name "klayout_dev"
   ```

#### `.mp4` 动画渲染

Router 支持将布线过程渲染到 `.mp4`。此功能依赖 `ffmpeg`（推荐 **conda** 安装）：

```bash
conda activate klayout_dev
conda install -c conda-forge ffmpeg
```

这不是必须的，也可以渲染到 `.gif`，这将不需要安装额外的依赖。

## ✨ 功能特性

Klayout-Router 可以在 KLayout 中自动生成 `Path` 对象，实现多对引脚的自动连线。

#### 动画示例

下面的芯片版图上，黑色边框为芯片边界，外围大约有 50 个红色大 pad，芯片核心区域分布着 50 个蓝色小 pad。绿色图形代表需要在布线时避开的障碍。示例中：芯片尺寸为 $5\text{mm} \times 5\text{mm}$；红色大 pad 为 $200\mu\text m \times 100\mu\text m$，蓝色小 pad 为 $4\mu\text m \times 4\mu\text m$。

<figure align="center">
   <img src="https://github.com/user-attachments/assets/78389716-95cb-4708-a7c5-be61fdaea4ee" alt="芯片概览动图，展示芯片轮廓、引脚和障碍物。" />
   <figcaption> 芯片整体示意图 </figcaption>
</figure>

如果手动画线，将每一个蓝色小 pad 分别连接到一个红色大 pad，即便是熟练的牛马，也往往需要 1 小时以上。耗时的原因在于布线存在大量物理约束：

- 连线之间不能相互交叉；
- 连线必须与所有绿色图形（障碍物）至少保持 1 µm 的安全距离；
- 连线之间一般需要保持 20 µm 的间距；
- 在红色大 pad 周围 200 µm 范围内不能出现任何连线——这是实验室配方中的一个自定义需求。

当 pad 数量达到 50 个、且存在如此多约束条件时，布线任务会变得相当复杂，人工作图不仅耗时，而且容易出错。

而本项目的 Router 可以在数秒内完成同样的布线任务：

<figure align="center">
   <img src="https://github.com/user-attachments/assets/69d262b3-166b-496d-9d3a-13920c79d71a" alt="布线过程的动画演示。" />
   <figcaption> 布线路径生成过程动画。由于帧率限制，动画速度看起来比真实布线要慢；在不渲染动画的情况下，实际布线通常在 3 秒内完成。 </figcaption>
</figure>

上面的动图只是算法过程的演示。最终生成的连线路径会实际写入 `.gds` 文件中：

<img src="./.assets/gds_after_routing.png" alt="image-20251121173423808" style="zoom: 33%;" />




#### Router 的行为

Router 在布线时会尽量满足以下约束：

1. 连线之间不发生交叉；
2. 连线与障碍物之间尽量保持不少于 $x$ µm 的安全距离；
3. 不同连线之间尽量保持不少于 $y$ µm 的间距；
4. 连线与其他非目标引脚之间尽量保持不少于 $z$ µm 的距离（除了其本身要连接的那一对 pad）。

其中，$x, y, z$ 都是你可以自由设定的参数。

需要注意的是：第 2、3、4 条被视为**软约束**。几何上如果不存在同时满足所有约束的完美解，Router 会采用“代价惩罚（cost penalty）”机制：

- 路径距离障碍物或其他路径越近，其代价（cost）越高；
- Router 会在代价场中寻找总代价尽可能小的解。

在可视化过程中，你可以看到所有路径、障碍物和引脚周围有一圈渐变“光晕”，它们就代表代价场。借助这一机制，即便完美解不存在，Router 也能给出一个“尽可能好”的方案，而不是简单报错失败。

#### 图像语义分割



## 🔥即刻启用

在完成 [环境要求](#环境要求) 中的环境配置后，你可以按照 `Auto-routing/routing.ipynb` 中的步骤，在示例文件 `./demo.gds` 上运行一次演示。更多技术细节请参考 `./Auto-routing/` 目录下的 [router 技术说明文档](./Auto-routing/manual.md)。在理解整个流程之后，你就可以把 Router 应用到自己的版图文件上了。

# 如何使用其他 Klayout 宏？

## ⚠️警告

通过 Klayout 宏对 `.gds` 文件做出的修改是 **不可逆的**！

在使用这些宏之前，请务必：

- 确认你清楚当前宏的功能；
- 或者先备份一份 `.gds` 文件，再进行操作。

由于未备份或误用本项目中的脚本/宏而导致的任何数据丢失或版图损坏，作者不承担任何责任。

## 环境要求

Klayout 宏完全在 *Klayout* 软件内部的 *Macro Development* 环境中运行。只需下载安装最新版本的 [KLayout Layout Viewer And Editor](https://www.klayout.de/build.html)，即可运行这些宏。

系统：推荐 Windows，不保证 MacOS 和 Linux.

## 什么是 Klayout 宏？

Klayout宏就是Ruby或Python脚本。

在 KLayout 中打开宏开发界面：

> 打开 Klayout GUI → 菜单栏 Macros → Macro development → 选择 Python → 在空白处右键选择 Add location → 选择 `/Klayout-Router/Macros/`

更多介绍请参考官方文档 [KLayout Macro Development](https://www.klayout.de/doc/about/macro_editor.html)。

<img src="./.assets/macro_development.png" alt="image-20251121173925408" style="zoom:33%;" />

## 功能一览

### Auto-patching（自动补丁）

此宏主要是为我们实验室的一个特定需求设计的，因此功能看上去会比较“客制化”。它可以在**版图图形与网格的每一个交点处自动放置一个 `Box` 补丁**。

> 更精确的描述：给定电极图层和网格图层后，脚本会在电极和网格的所有交点处，生成边长为 $d$ 的正方形补丁 `Box`，并放置在 `patch layer` 中。此外，它还会保留“扩展后的网格”（即用宽度为 $d$ 的网格替代原有网格）与电极的所有交叠区域。

注意：**网格必须是若干 box 或 polygon 的边缘。** 脚本是通过提取这些 polygon/box 的边缘来识别网格的，因此你需要在 `grid layer` 上画出这些 polygon 或 box。

使用方法：

1) 打开 KLayout → Macro Development；
2) 加载 `./Macros/Auto-patching/run_patch.py` 并运行（**点击带感叹号的绿色三角按钮**）；
3) 在弹出的对话框中设置 Electrode 图层、Writing Field 图层、Patch 图层（补丁放置图层）以及 Patch size（单位 µm）；
4) 点击 `Create Patches`，然后在 Patch 图层中检查生成结果。

<figure align="center">
   <img src="https://github.com/user-attachments/assets/59c04026-7672-4a00-8a12-f042b05b61b6" alt="展示如何在 Klayout 中使用 Auto-patching 宏的截图。" />
   <figcaption> 在 Klayout 中使用 Auto-patching 宏的示例 </figcaption>
</figure>

> **为什么需要这么“奇怪”的功能？**
>
> 这是一个高度定制的实验需求。在我们的 EBL 系统中，最大 writing field (WF) 只有 1mm×1mm，而整个芯片尺寸却是 5mm×5mm！这意味着整个芯片版图必须被分成 25 个 WF 分别曝光。这样在不同 WF 交界处可能会出现对位误差，导致实际曝光出的图形在 WF 边界处无法连上。
> 因此，我们需要在图形与 WF 边界的每个交点处添加补丁（patch），并用不同的写场设置单独曝光补丁图层，确保补丁本身不会再次与新的 WF 边界相交。这样可以修复 WF 交界处的断连问题。当然，与工业界的光刻方案相比，这种做法相当笨拙，也充满妥协。

### 其他 Widgets

在 `./Macros/Widgets/` 目录下还提供了一些可以显著减少鼠标点击的小工具：

- `change_layer.py`：将若干源图层上的所有图形移动到目标图层（可选是否清空源图层）；
- `change_layer_by_overlapping.py`：只在与参考图层发生几何重叠的情况下，将源图层上的图形移动到目标图层；
- `change_path_width.py`：修改指定图层上 `Path` 图形的线宽；
- `clear_layer.py`：清空当前 cell/layout 中指定图层上的所有图形；
- `merge_layer.py`：将多个图层上的图形合并到一起。

这些方法没有图形界面，但是使用方式与 Auto-patching 一样，参考上图。所有操作都是针对当前正在查看的版图进行。

# 许可证

本项目采用 MIT 许可证。

Copyright (c) 2025 Legendrexial

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
