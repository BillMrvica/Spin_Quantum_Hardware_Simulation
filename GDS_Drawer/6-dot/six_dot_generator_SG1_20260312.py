"""
6-Dot Quantum Device GDS Layout Generator (SG1 Custom)
======================================================
在 `six_dot_generator.py` 的基础上，仅修改 SG1 的几何外形，用于匹配手绘标注版本。
其余所有 gate/ohmic/引线/层号/参数保持与原脚本一致。

作者：Matrix Agent
日期：2026-03-12
"""

import os

import gdstk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon as MplPolygon


# ===================================================================
# 层定义
# ===================================================================
# L0: accumulation S/D（原来函数内的 S/D layer，语义改名但几何保持一致）
# L7: annealing S/D（新增橙色电极，用于 annealing S/D layer）
LAYER_SD_ANNEAL = 1      # Annealing Source/Drain - 退火源漏层（橙色）
LAYER_SD_ACCUM = 3       # Accumulation Source/Drain - 累积源漏层（浅蓝色）
LAYER_SG = 2             # Screen Gate - 屏蔽门层
LAYER_BG = 4             # Barrier Gate - 势垒门层
LAYER_PG = 3             # Plunger Gate - 柱栅层

LAYER_LABEL = 100        # 标签层

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False


def plot_gds(cell, title="量子点器件布局"):
    """
    可视化GDS布局
    
    Parameters:
    -----------
    cell : gdstk.Cell
        GDS单元
    title : str
        图表标题
    """
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.set_title(title, fontsize=18)
    
    # 层颜色配置
    # 说明：这里的 “第几层” 指的是绘制层级(视觉叠放顺序)，不是 GDS 的 layer number。
    # 目标绘制顺序（从底到顶 / 先画到后画）：
    #   第1层：SD annealing (L7)
    #   第2层：SG (L1)
    #   第3层：BG (L2)
    #   第4层：PG (L3) + SD accumulation (L0) （但用不同颜色区分）
    layer_config = {
        LAYER_SD_ANNEAL: {'color': '#FF8C00', 'alpha': 0.85, 'label': 'Annealing S/D (L7) - 退火源漏', 'zorder': 1},
        LAYER_SG: {'color': '#D3D3D3', 'alpha': 0.45, 'label': 'SG (L1) - 屏蔽门', 'hatch': '///', 'zorder': 2},
        LAYER_BG: {'color': '#FF69B4', 'alpha': 0.90, 'label': 'BG (L2) - 势垒门', 'zorder': 3},
        # 第4层：PG + Accum SD，颜色区分
        LAYER_PG: {'color': '#7B68EE', 'alpha': 0.90, 'label': 'PG (L3) - 柱栅', 'zorder': 4},
        LAYER_SD_ACCUM: {'color': '#4DB6FF', 'alpha': 0.65, 'label': 'Accumulation S/D (L0) - 累积源漏', 'zorder': 4},
    }
    
    drawn_labels = set()
    
    # 按绘制层级排序（先画底层，后画顶层）
    sorted_layers = [LAYER_SD_ANNEAL, LAYER_SG, LAYER_BG, LAYER_PG, LAYER_SD_ACCUM]
    
    for layer in sorted_layers:
        polygons = cell.get_polygons(layer=layer, datatype=0)
        if not polygons:
            continue
        
        cfg = layer_config.get(layer)
        label = cfg['label']
        
        for gds_poly in polygons:
            pts = gds_poly.points
            if pts.ndim != 2 or len(pts) < 3:
                continue
            
            current_label = label if label not in drawn_labels else None
            poly = MplPolygon(
                pts, closed=True,
                facecolor=cfg['color'],
                edgecolor='black',
                linewidth=0.5,
                alpha=cfg['alpha'],
                hatch=cfg.get('hatch', ''),
                label=current_label,
                zorder=cfg.get('zorder', 1)
            )
            ax.add_patch(poly)
            if current_label:
                drawn_labels.add(label)
    
    # 添加图例
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    ax.set_aspect('equal')
    
    # 自动调整视图
    bbox = cell.bounding_box()
    if bbox is not None:
        min_x, min_y = bbox[0]
        max_x, max_y = bbox[1]
        width = max_x - min_x
        height = max_y - min_y
        margin_x = width * 0.1 if width > 0 else 1
        margin_y = height * 0.1 if height > 0 else 1
        ax.set_xlim(min_x - margin_x, max_x + margin_x)
        ax.set_ylim(min_y - margin_y, max_y + margin_y)
        print(f"视图范围: X[{min_x:.2f}, {max_x:.2f}], Y[{min_y:.2f}, {max_y:.2f}]")
    else:
        ax.autoscale_view()
    
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')
    plt.grid(True, which='both', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()


def _add_label(cell, connection_points, name, xy):
    lbl = gdstk.Label(name, xy, layer=LAYER_LABEL)
    cell.add(lbl)
    connection_points[name] = lbl.origin
    return lbl


def _extend_sd_lead_path(start_x, start_y, x_sign, y_sign, s_extend, bend, end_extend):
    """
    统一四个 S/D 延长引线的折线（不改变原有几何公式）
    returns: path_pts list[tuple[float,float]]
    """
    return [
        (start_x, start_y),
        (start_x + x_sign * s_extend, start_y),
        (start_x + x_sign * (s_extend + bend), start_y + y_sign * bend),
        (start_x + x_sign * (s_extend + bend), start_y + y_sign * end_extend),
    ]


def create_6qd_device_core(
    # 量子点参数
    num_dots=6,
    
    # PG (柱栅) 参数
    pg_max_width=0.120,
    pg_vert_side_len=0.040,
    pg_chamfer_h=0.040,
    pg_top_flat_w=0.040,
    pg_bot_flat_w=0.042,
    
    # BG (势垒门) 参数
    bg_max_width=0.060,
    bg_vert_side_len=0.040,
    bg_top_flat_w=0.020,
    bg_bot_flat_w=0.042,
    
    # 门间距
    gap_pg_bg=0.002,
    d1_gap=0.020,  # PG到SG中线
    d2_gap=0.020,  # SG中线到上一行PG
    
    # SG (屏蔽门) 参数
    sg_mid_thick=0.100,
    sg_top_thick=0.300,
    sg_bot_thick=0.300,
    gap_gate_outer_sg=0.030,
    sg_extension=0.200,
    
    # Accumulation S/D 参数（原来 core 里的 S/D）
    sd_gap_to_gate=0.01,
    sd_height=0.10,
    sd_width=0.35,
    sd_lead_width=0.10,
    sd_outer_length=0.300,
    sd_lead_extend = 0.7,
    
    # Annealing S/D 参数（新增橙色电极，额外一层；5 个竖直长方形）
    # 规则（按你的描述）：
    # - SET_S anneal：x 与 SET_S accumulation 相同；y 在上方
    # - SET1_D/SET2_D anneal：相对各自 accumulation 的 (dx,dy) 分别为 (-0.5,+0.5)/(+0.5,+0.5)
    # - QD_D/QD_S anneal：相对各自 accumulation 的 (dx,dy) 分别为 (-0.5,-0.5)/(+0.5,-0.5)
    anneal_sd_enable=True,
    anneal_sd_gap_to_accum=0.50,   # 与对应 accumulation S/D 的边缘间距（500nm）- 用于计算 lead 延长/重叠
    anneal_sd_rect_w=0.10,         # 竖直长方形 x 宽度（100nm）
    anneal_sd_rect_h=0.40,         # 竖直长方形 y 高度（400nm）
    anneal_sd_dx=0.50,             # 相对 accumulation 的 x 偏移量（±500nm）
    anneal_sd_dy=0.30,  
    anneal_central_set_dx = 0.4,           # 相对 accumulation 的 y 偏移量（±500nm）
    anneal_central_set_dy= 0.8,

    anneal_sd_gap_override=False,  # True 时强制用“边缘间距”而不是 (dx,dy) 放置
    accum_sd_lead_extend_to_anneal=0.20,  # 蓝色 accumulation 的 lead 额外延长（朝向 anneal），用于与 anneal 重叠（um）
    anneal_set_s_enable=True,
    
    # 引线参数
    lead_width=0.042,
    lead_length_bot=0.4,
    lead_length_top=0.5,
    lead_overlap=0.025
):
    """
    创建6量子点器件核心布局
    
    基于图片中的设计，创建包含：
    - 下方主阵列：6个量子点
    - 上方传感阵列：2个量子点
    
    Parameters:
    -----------
    num_dots : int
        主阵列量子点数量 (默认6)
    其他参数 : float
        门和引线的几何参数
    
    Returns:
    --------
    tuple : (lib, cell, connection_points)
    """
    
    lib = gdstk.Library()
    cell = lib.new_cell('6QD_DEVICE_CORE')
    connection_points = {}
    
    # ========== 辅助函数 ==========
    
    def get_pg_points(c):
        """生成PG (柱栅) 形状点"""
        cx, cy = c
        yvt = pg_vert_side_len / 2
        yvb = -pg_vert_side_len / 2
        yt = yvt + pg_chamfer_h
        yb = yvb - pg_chamfer_h
        xm = pg_max_width / 2
        xtf = pg_top_flat_w / 2
        xbf = pg_bot_flat_w / 2
        return [(cx + xtf, cy + yt), (cx + xm, cy + yvt), (cx + xm, cy + yvb),
                (cx + xbf, cy + yb), (cx - xbf, cy + yb), (cx - xm, cy + yvb),
                (cx - xm, cy + yvt), (cx - xtf, cy + yt)], yt, yb
    
    def get_bg_points(c):
        """生成BG (势垒门) 形状点"""
        cx, cy = c
        yvt = bg_vert_side_len / 2
        yvb = -bg_vert_side_len / 2
        bcht = (bg_max_width - bg_top_flat_w) / 2
        bchb = (bg_max_width - bg_bot_flat_w) / 2
        yt = yvt + bcht
        yb = yvb - bchb
        xm = bg_max_width / 2
        xtf = bg_top_flat_w / 2
        xbf = bg_bot_flat_w / 2
        return [(cx + xtf, cy + yt), (cx + xm, cy + yvt), (cx + xm, cy + yvb),
                (cx + xbf, cy + yb), (cx - xbf, cy + yb), (cx - xm, cy + yvb),
                (cx - xm, cy + yvt), (cx - xtf, cy + yt)], yt, yb
    
    def create_rounded_rect(xmin, ymin, xmax, ymax, corners, radius=0.020, layer=LAYER_SD_ACCUM):
        """创建圆角矩形"""
        poly = gdstk.Polygon(
            [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)],
            layer=layer
        )
        poly.fillet([radius if c else 0 for c in corners], tolerance=0.001)
        return poly
    
    # ========== 计算几何参数 ==========
    
    gate_half_height = pg_vert_side_len / 2 + pg_chamfer_h
    row_pitch = gate_half_height + d1_gap + sg_mid_thick + d2_gap + gate_half_height
    
    # 引线端点
    pglb = -(pg_vert_side_len / 2 + pg_chamfer_h)
    bglb = -(bg_vert_side_len / 2 + (bg_max_width - bg_bot_flat_w) / 2)
    fbe = min(pglb, bglb) - lead_length_bot  # 底部引线终点
    
    yct = row_pitch  # 第二行（传感器阵列）y坐标
    
    pglt = pg_vert_side_len / 2 + pg_chamfer_h
    bglt = bg_vert_side_len / 2 + (bg_max_width - bg_top_flat_w) / 2
    fte = yct + max(pglt, bglt) + lead_length_top  # 顶部引线终点
    
    # ========== 创建形状 ==========
    all_shapes = []
    
    # ========== 底部主阵列 (6个量子点) ==========
    ycb = 0  # 底部行y坐标
    curr_x = 0
    bgn = []
    [bgn.extend([f'QD_B{i}', f'QD_PG{i}']) for i in range(1, num_dots + 1)]
    bgn.append(f'QD_B{num_dots + 1}')
    bot_tips_y = []
    
    for i, n in enumerate(bgn):
        gt = 'PG' if 'PG' in n else 'BG'
        w = pg_max_width if gt == 'PG' else bg_max_width
        lyr = LAYER_PG if gt == 'PG' else LAYER_BG
        
        if i > 0:
            pw = pg_max_width if 'PG' in bgn[i - 1] else bg_max_width
            curr_x += (pw / 2) + gap_pg_bg + (w / 2)
        
        cen = (curr_x, ycb)
        pts, tip_y, by = get_pg_points(cen) if gt == 'PG' else get_bg_points(cen)
        bot_tips_y.append(tip_y)
        
        poly = gdstk.Polygon(pts, layer=lyr)
        
        # 添加引线
        lys = by + lead_overlap
        lye = fbe
        lead = gdstk.rectangle((curr_x - lead_width / 2, lys),
                               (curr_x + lead_width / 2, lye), layer=lyr)
        
        all_shapes.extend(gdstk.boolean([poly, lead], [], 'or', layer=lyr))
        
        # 添加标签
        _add_label(cell, connection_points, n, (curr_x, lye))
    
    # 计算中心区域
    bmx = -bg_max_width / 2
    bax = curr_x + bg_max_width / 2
    cxg = (bmx + bax) / 2
    
    # 传感器阵列位置（上方2个量子点）
    igr = 2 * sd_gap_to_gate + sd_width
    off_cen = igr / 2 + bg_max_width / 2 + gap_pg_bg + pg_max_width / 2
    lscpgx = cxg - off_cen   # 左侧传感器中心
    rscpgx = cxg + off_cen   # 右侧传感器中心
    
    # ========== 传感器阵列 (2个量子点 - SET1和SET2) ==========
    
    def create_sensor_set(cpgx, side):
        """创建传感器门阵列"""
        shps = []
        offs = [-(pg_max_width / 2 + gap_pg_bg + bg_max_width / 2), 0,
                (pg_max_width / 2 + gap_pg_bg + bg_max_width / 2)]
        
        if side == 'L':
            sgn = ['SET1_B1', 'SET1_G', 'SET1_B2']
        else:
            sgn = ['SET2_B1', 'SET2_G', 'SET2_B2']
        
        for off, n in zip(offs, sgn):
            cx = cpgx + off
            gt = 'PG' if '_G' in n else 'BG'
            lyr = LAYER_PG if gt == 'PG' else LAYER_BG
            
            pts, _, _ = get_pg_points((cx, yct)) if gt == 'PG' else get_bg_points((cx, yct))
            poly = gdstk.Polygon(pts, layer=lyr)
            
            # 引线
            Wty = np.max([pt[1] for pt in pts])
            lys = Wty - lead_overlap
            lye = fte
            lead = gdstk.rectangle((cx - lead_width / 2, lys),
                                   (cx + lead_width / 2, lye), layer=lyr)
            
            shps.extend(gdstk.boolean([poly, lead], [], 'or', layer=lyr))
            
            # 标签
            _add_label(cell, connection_points, n, (cx, lye))
        
        return shps
    
    all_shapes.extend(create_sensor_set(lscpgx, side='L'))  # SET1
    all_shapes.extend(create_sensor_set(rscpgx, side='R'))  # SET2
    
    # ========== 侧边门 (SG1, SG2, SG3) ==========
    
    sgxs = bmx - sg_extension
    sgxe = bax + sg_extension
    
    # SG2 - 中间屏蔽门
    sg2xs = bmx - sd_gap_to_gate - sd_outer_length - anneal_sd_dx # SG2 x 坐标，向左偏移使其分离QD_D和SET1_D的accumulation
    sg2xe = bax + sd_gap_to_gate + sd_outer_length + anneal_sd_dx # SG2 x 坐标，向左偏移使其分离QD_D和SET1_D的accumulation
    ysgmb = max(bot_tips_y) + d1_gap
    ysgmt = ysgmb + sg_mid_thick
    sg_mid = gdstk.rectangle((sg2xs, ysgmb), (sg2xe, ysgmt), layer=LAYER_SG)
    
    # SG1 - 顶部屏蔽门（按手绘轮廓改形状；其余结构保持不变）
    # 说明：原始实现为两个矩形叠加（sg_top + sg_top2）。这里改为“带台阶/凹口”的单一多边形轮廓。
    # 基准线：
    #   - x 左右边界：sgxs / sgxe（保持不变）
    #   - 中间台阶左右边界：sg2xs / sg2xe（保持不变）
    #   - y 底边：ysgtb（保持不变）
    #   - y 第一层顶部：ysgtb + (anneal_sd_dy + anneal_sd_rect_w/2)（保持与原 sg_top 顶边一致）
    #   - y 第二层顶部：上述 + sg_top_thick（保持与原 sg_top2 顶边一致）
    ysgtb = yct + gate_half_height + gap_gate_outer_sg
    ysgtt = ysgtb + anneal_sd_dy + anneal_sd_rect_w / 2
    ysgtt2 = ysgtt + sg_top_thick
    ysgtt3 = ysgtt2 + anneal_central_set_dy + sd_height / 2
    # 轮廓点按逆时针给出：参考手绘红色阴影，顶边为“中间下凹”的形状（左右高、中间低）。
    # 结构：
    #   - 底边在 ysgtb，全宽 sgxs~sgxe
    #   - 左右两端上到 ysgtt2
    #   - 中间区域（sg2xs~sg2xe）只上到 ysgtt
    shw = bg_max_width + gap_pg_bg + pg_max_width / 2
    lsre = lscpgx + shw
    smxm = lsre + sd_gap_to_gate
    rsle = rscpgx - shw
    smxa = rsle - sd_gap_to_gate
    msdcx = (smxm + smxa) / 2

    sg1_pts = [
        (sgxs, ysgtt2),
        (sgxs, ysgtt3),
        (msdcx - sd_lead_width, ysgtt3),
        (msdcx + sd_lead_width, ysgtt3),
        (sgxe, ysgtt3),
        (sgxe, ysgtt2),
        (sg2xe, ysgtt2),
        (sg2xe, ysgtt),
        (sgxe, ysgtt),
        (sgxe, ysgtb),
        (msdcx + sd_lead_width, ysgtb),
        (msdcx + sd_lead_width, ysgtt3 - sg_top_thick),
        (msdcx - sd_lead_width, ysgtt3 - sg_top_thick),
        (msdcx - sd_lead_width, ysgtb),
        (sgxs, ysgtb),
        (sgxs, ysgtt),
        (sg2xs, ysgtt),
        (sg2xs, ysgtt2),
    ]
    sg_top = gdstk.Polygon(sg1_pts, layer=LAYER_SG)
    
    # SG3 - 底部屏蔽门
    ysgbt = ycb - gate_half_height - gap_gate_outer_sg
    # ysgbb = ysgbt - sg_bot_thick
    ysgbb = ysgbt - anneal_sd_dy - anneal_sd_rect_w / 2
    sg_bot1 = gdstk.rectangle((sgxs, ysgbb), (sgxe, ysgbt), layer=LAYER_SG)
    sg_bot2 = gdstk.rectangle((sg2xs, ysgbb), (sg2xe, ysgbb - sg_bot_thick), layer=LAYER_SG)
    
    sg3_pts = [
        (sg2xs, ysgbb),
        (sg2xs, ysgbb - sg_bot_thick),
        (sg2xe, ysgbb - sg_bot_thick),
        (sg2xe, ysgbb),
        (sgxe, ysgbb),
        (sgxe, ysgbt),
        (sgxs, ysgbt),
        (sgxs, ysgbb),
    ]
    sg_bot = gdstk.Polygon(sg3_pts, layer=LAYER_SG)

    # all_shapes.extend([sg_top, sg_mid, sg_bot1, sg_bot2])
    all_shapes.extend([sg_top, sg_mid, sg_bot])
    
    # SG标签
    _add_label(cell, connection_points, "SG1", (sgxs, (ysgtb + ysgtt) / 2))
    _add_label(cell, connection_points, "SG2", (sgxs, (ysgmb + ysgmt) / 2))
    _add_label(cell, connection_points, "SG3", (sgxs, (ysgbb + ysgbt) / 2))
    
    # ========== Accumulation S/D 区域（原 S/D） ==========

    sd_shapes = []
    
    # 主阵列 S/D (底部左右两端)
    # 左侧
    blbe = bmx
    blsix = blbe - sd_gap_to_gate
    golx = blsix - sd_outer_length
    
    sd_shapes.append(create_rounded_rect(golx, ycb - sd_height / 2,
                                          blsix, ycb + sd_height / 2,
                                          [False, True, True, False], layer=LAYER_SD_ACCUM))


    # QD_D extend lead - 为了与 anneal 电极有一部分重叠，向左额外延长 extend_to_anneal
    # 说明：anneal 对 QD_D 的中心相对 QD_D accumulation 的中心偏移为 (-anneal_sd_dx, -anneal_sd_dy)；
    #      QD_DA 在 add_anneal_rect 里被改成“水平长条”，其 y 尺寸为 anneal_sd_rect_w。
    # 目标：让 accumulation 的延长引线末端进入 anneal 形状内部一点，保证有重叠。
    s_extend = anneal_sd_dx - anneal_sd_rect_h / 2
    bend = 0.05
    path_pts = _extend_sd_lead_path(
        golx, ycb,
        x_sign=-1, y_sign=-1,
        s_extend=s_extend, bend=bend,
        end_extend=sd_lead_extend,
    )
    sd_shapes.append(gdstk.FlexPath(path_pts, width=0.1, layer=LAYER_SD_ACCUM, ends="flush"))
    
    _add_label(cell, connection_points, "QD_D", (golx - s_extend - bend, ycb - sd_lead_extend))

    # 右侧
    brbe = bax
    brsix = brbe + sd_gap_to_gate
    gorx = brsix + sd_outer_length
    
    sd_shapes.append(create_rounded_rect(brsix, ycb - sd_height / 2,
                                          gorx, ycb + sd_height / 2,
                                          [True, False, False, True], layer=LAYER_SD_ACCUM))
    


    # QD_S extend lead - 向右下额外延长，使其与 QD_SA (anneal) 有部分重叠
    # anneal(QD_SA) 中心在 (gorx + anneal_sd_dx, ycb - anneal_sd_dy)，且 QD_SA 为水平长条，y 尺寸为 anneal_sd_rect_w
    s_extend = anneal_sd_dx - anneal_sd_rect_h / 2
    bend = 0.05
    path_pts = _extend_sd_lead_path(
        gorx, ycb,
        x_sign=+1, y_sign=-1,
        s_extend=s_extend, bend=bend,
        end_extend=sd_lead_extend,
    )
    sd_shapes.append(gdstk.FlexPath(path_pts, width=0.1, layer=LAYER_SD_ACCUM, ends="flush"))
    _add_label(cell, connection_points, "QD_S", (gorx + s_extend + bend, ycb - sd_lead_extend))
    
    # 传感器阵列 S/D
    shw = bg_max_width + gap_pg_bg + pg_max_width / 2
    
    # SET1 S/D (左侧传感器)
    lsle = lscpgx - shw
    tlsix = lsle - sd_gap_to_gate
    
    sd_shapes.append(create_rounded_rect(golx, yct - sd_height / 2,
                                          tlsix, yct + sd_height / 2,
                                          [False, True, True, False], layer=LAYER_SD_ACCUM))
    

    # SET1_D extend lead - 向左上额外延长，使其与 SET1_DA (anneal) 有部分重叠
    # anneal(SET1_DA) 中心在 (golx - anneal_sd_dx, yct + anneal_sd_dy)，且 SET1_DA 为水平长条，y 尺寸为 anneal_sd_rect_w
    s_extend = anneal_sd_dx - anneal_sd_rect_h / 2
    bend = 0.05
    path_pts = _extend_sd_lead_path(
        golx, yct,
        x_sign=-1, y_sign=+1,
        s_extend=s_extend, bend=bend,
        end_extend=sd_lead_extend,
    )
    sd_shapes.append(gdstk.FlexPath(path_pts, width=0.1, layer=LAYER_SD_ACCUM, ends="flush"))

    _add_label(cell, connection_points, "SET1_D", (golx - s_extend - bend, yct + sd_lead_extend))
    
    # SET2 S/D (右侧传感器)
    rsre = rscpgx + shw
    trsix = rsre + sd_gap_to_gate # DO NOT MODIFY THIS - 这个点是 SET2_D 的边缘点，不是中心点，create_rounded_rect 里会根据这个点和 sd_width 自动计算中心点位置 
    
    sd_shapes.append(create_rounded_rect(trsix, yct - sd_height / 2,
                                          gorx, yct + sd_height / 2,
                                          [True, False, False, True], layer=LAYER_SD_ACCUM))
    

    # SET2_D extend lead - 向右上额外延长，使其与 SET2_DA (anneal) 有部分重叠
    # anneal(SET2_DA) 中心在 (gorx + anneal_sd_dx, yct + anneal_sd_dy)，且 SET2_DA 为水平长条，y 尺寸为 anneal_sd_rect_w
    s_extend = anneal_sd_dx - anneal_sd_rect_h / 2
    bend = 0.05
    path_pts = _extend_sd_lead_path(
        gorx, yct,
        x_sign=+1, y_sign=+1,
        s_extend=s_extend, bend=bend,
        end_extend=sd_lead_extend,
    )
    sd_shapes.append(gdstk.FlexPath(path_pts, width=0.1, layer=LAYER_SD_ACCUM, ends="flush"))

    _add_label(cell, connection_points, "SET2_D", (gorx + s_extend + bend, yct + sd_lead_extend))
    
    # 公共 S/D (中间传感器)
    lsre = lscpgx + shw
    smxm = lsre + sd_gap_to_gate
    rsle = rscpgx - shw
    smxa = rsle - sd_gap_to_gate
    
    msdp = create_rounded_rect(smxm, yct - sd_height / 2,
                                  smxa, yct + sd_height / 2,
                                  [True, True, True, True], layer=LAYER_SD_ACCUM)
    
    msdcx = (smxm + smxa) / 2
    mlye = yct + anneal_central_set_dy + anneal_sd_dy - anneal_sd_rect_h / 2 - sd_height / 2
    # 中间S/D引线（蓝色 accumulation）：
    # 为了与 annealing 电极有一部分重叠，向上额外延长 accum_sd_lead_extend_to_anneal
    msdl = gdstk.rectangle((msdcx - sd_lead_width / 2, yct),
                             (msdcx + sd_lead_width / 2, mlye + accum_sd_lead_extend_to_anneal),
                             layer=LAYER_SD_ACCUM)

    sd_shapes.extend(gdstk.boolean([msdp, msdl], [], 'or', layer=LAYER_SD_ACCUM))
    
    # 中间 S/D 增加覆盖anneal电极的引线
    bend = sd_lead_width*3/4  # 弯曲半径，确保引线宽度足够覆盖 anneal 电极
    extend = yct + anneal_central_set_dy + anneal_sd_dy - mlye - accum_sd_lead_extend_to_anneal
    print(f"计算中间 S/D 引线延长部分: bend={bend:.3f}, extend={extend:.3f}")
    path_pts_set_s = [
        (msdcx, mlye + accum_sd_lead_extend_to_anneal),
        # set_s_accum_cx+anneal_central_set_dx, set_accum_cy + anneal_central_set_dy + dy
        (msdcx, mlye + accum_sd_lead_extend_to_anneal + extend),
        (msdcx - bend, mlye + accum_sd_lead_extend_to_anneal + extend + bend),
        (msdcx - bend, ysgtt3 + anneal_sd_rect_h/2 )
    ]
    sd_shapes.append(gdstk.FlexPath(path_pts_set_s, width=0.1, layer=LAYER_SD_ACCUM, ends="flush"))

    _add_label(cell, connection_points, "SET_S", path_pts_set_s[-1])
    
    all_shapes.extend(sd_shapes)
    
    # ========== Annealing S/D 区域（新增橙色电极：5 个竖直长方形） ==========
    if anneal_sd_enable:
        anneal_shapes = []
        
        rect_w = anneal_sd_rect_w
        rect_h = anneal_sd_rect_h
        dx = anneal_sd_dx
        dy = anneal_sd_dy
        
        # 参考点：accumulation S/D 的“中心点”
        # QD：用左右 S/D 的中心 (x_center, ycb)
        
        # qd_d_accum_cx = (golx + blsix) / 2
        # qd_s_accum_cx = (brsix + gorx) / 2
        qd_d_accum_cx = golx
        qd_s_accum_cx = gorx
        
        qd_accum_cy = ycb
        
        # SET：SET1_D / SET2_D 用左右 S/D 的中心；SET_S 用中间 S/D 的中心
        
        # set1_d_accum_cx = (golx + tlsix) / 2
        # set2_d_accum_cx = (trsix + gorx) / 2
        set1_d_accum_cx = golx
        set2_d_accum_cx = gorx
        
        set_s_accum_cx = msdcx
        set_accum_cy = yct
        
        def add_anneal_rect(name, xc, yc):
            # 默认：竖直长方形 (rect_w x rect_h)
            w = rect_w
            h = rect_h
            # 将 QD_S/D 和 SET 左右两端的 D 这四个 annealing gate 改为“水平方向”
            # 即：交换宽高，让长边沿 x 方向
            if name in {"QD_DA", "QD_SA", "SET1_DA", "SET2_DA"}:
                w, h = rect_h, rect_w

            r = gdstk.rectangle(
                (xc - w / 2, yc - h / 2),
                (xc + w / 2, yc + h / 2),
                layer=LAYER_SD_ANNEAL
            )
            anneal_shapes.append(r)
            # _add_label(cell, connection_points, name, (xc, yc))
        
        # QD：2 个，在 accumulation 下方（dy 为负），并左右分别 ±dx
        add_anneal_rect("QD_DA", qd_d_accum_cx - dx, qd_accum_cy - dy)
        _add_label(cell, connection_points, "QD_DA", (qd_d_accum_cx - dx - rect_h/2, qd_accum_cy - dy))
        add_anneal_rect("QD_SA", qd_s_accum_cx + dx, qd_accum_cy - dy)
        _add_label(cell, connection_points, "QD_SA", (qd_s_accum_cx + dx + rect_h/2, qd_accum_cy - dy))
        
        # SET：3 个，在 accumulation 上方（dy 为正）
        # SET1_DA / SET2_DA：相对各自 accumulation 分别 (-dx,+dy) / (+dx,+dy)
        add_anneal_rect("SET1_DA", set1_d_accum_cx - dx, set_accum_cy + dy)
        _add_label(cell, connection_points, "SET1_DA", (set1_d_accum_cx - dx - rect_h/2, set_accum_cy + dy))
        add_anneal_rect("SET2_DA", set2_d_accum_cx + dx, set_accum_cy + dy)
        _add_label(cell, connection_points, "SET2_DA", (set2_d_accum_cx + dx + rect_h/2, set_accum_cy + dy))
        # SET_SA：x 坐标相同（不加 dx）
        # add_anneal_rect("SET_SA", set_s_accum_cx+anneal_central_set_dx, set_accum_cy + anneal_central_set_dy)  # 为了区分，SET_SA 在 y 方向上稍微偏移一点（dy*sqrt(2)/2）
        path_pts_set_sa = [
            (set_s_accum_cx+anneal_central_set_dx, set_accum_cy + anneal_central_set_dy + dy - rect_h / 2),
            (set_s_accum_cx+anneal_central_set_dx, set_accum_cy + anneal_central_set_dy + dy),
            (set_s_accum_cx+anneal_central_set_dx + rect_w*3/4, set_accum_cy + anneal_central_set_dy + dy + rect_w*3/4),
            (set_s_accum_cx+anneal_central_set_dx + rect_w*3/4, ysgtt3 + rect_h/2),
        ]
        anneal_shapes.append(gdstk.FlexPath(path_pts_set_sa, width=rect_w, layer=LAYER_SD_ANNEAL, ends="flush"))
        _add_label(cell, connection_points, "SET_SA", path_pts_set_sa[-1])
        all_shapes.extend(anneal_shapes)
    
    # 添加所有形状到cell
    cell.add(*all_shapes)
    
    print(f"6量子点器件核心已创建:")
    print(f"  - 主阵列量子点: {num_dots} 个")
    print(f"  - 传感器阵列量子点: 2 个")
    print(f"  - 总连接点数: {len(connection_points)}")
    
    return lib, cell, connection_points


def save_gds(lib, filepath):
    """保存GDS文件"""
    lib.write_gds(filepath)
    print(f"GDS已保存: {filepath}")


# ============================================================================
# 主函数
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("6量子点GDS布局生成器")
    print("基于图片 image_1773139486280.png 设计")
    print("=" * 70)
    
    # 创建输出目录（固定到当前脚本目录下，避免从项目根运行时写到根目录 ./output_6dot）
    output_dir = os.path.join(os.path.dirname(__file__), "output_6dot")
    os.makedirs(output_dir, exist_ok=True)
    
    # 生成6量子点器件核心
    print("\n生成6量子点器件核心...")
    lib, cell, points = create_6qd_device_core(
        num_dots=6,
        pg_max_width=0.120,
        pg_vert_side_len=0.040,
        pg_chamfer_h=0.040,
        bg_max_width=0.060,
        bg_vert_side_len=0.040,
        gap_pg_bg=0.002,
        d1_gap=0.020,
        d2_gap=0.020,
        sg_mid_thick=0.100,
        sg_top_thick=0.100,
        sg_bot_thick=0.100,
        
        anneal_sd_enable=True,
        anneal_sd_gap_to_accum=0.50,   # 与对应 accumulation S/D 的边缘间距（500nm）- 用于计算 lead 延长/重叠
        anneal_sd_rect_w=0.10,         # 竖直长方形 x 宽度（100nm）
        anneal_sd_rect_h=0.40,         # 竖直长方形 y 高度（400nm）
        anneal_sd_dx=0.50,             # 相对 accumulation 的 x 偏移量（±500nm）
        anneal_sd_dy=0.20,  
        anneal_central_set_dx = 0.0,           # 相对 accumulation 的 y 偏移量（±500nm）
        anneal_central_set_dy= 0.5,

        anneal_sd_gap_override=False,  # True 时强制用“边缘间距”而不是 (dx,dy) 放置
        accum_sd_lead_extend_to_anneal=0.20,  # 蓝色 accumulation 的 lead 额外延长（朝向 anneal），用于与 anneal 重叠（um）
        anneal_set_s_enable=True,
        
        lead_width=0.042,
        lead_length_bot=0.6,
        lead_length_top=1
    )
    
    # 可视化
    print("\n可视化布局...")
    plot_gds(cell, title="6量子点器件核心 (基于图片设计)")
    
    # 保存GDS
    print("\n保存GDS文件...")
    save_gds(lib, os.path.join(output_dir, "6qd_device_core_short.gds"))
