#!/usr/bin/env python3
"""
直接执行six_dot_generator_demo.ipynb中的代码
"""

import os
import sys

# 设置工作目录
os.chdir(r"C:\Users\Administrator\Desktop\Spin_Quantum_Hardware_Simulation\GDS_Drawer\Prefab-6-dot")

# 添加当前目录到路径
sys.path.insert(0, os.getcwd())

print("=" * 70)
print("6-Dot Quantum Device GDS Layout Generator - 演示")
print("=" * 70)
print()

# 导入必要的库
import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

# 导入自定义类
from six_dot_generator import SixDotQuantumDeviceGenerator, plot_gds

# 设置输出目录
output_dir = r"C:\Users\Administrator\Desktop\Spin_Quantum_Hardware_Simulation\GDS_Drawer\Prefab-6-dot\output_gds"
os.makedirs(output_dir, exist_ok=True)
print(f"输出目录: {os.path.abspath(output_dir)}")
print()

# ========================================
# 步骤1: 创建并可视化量子点器件核心
# ========================================
print("=" * 70)
print("步骤1: 创建量子点器件核心")
print("=" * 70)

generator = SixDotQuantumDeviceGenerator(output_dir=output_dir)
lib_device, cell_device, device_points = generator.create_quantum_dot_device()

print("量子点器件已创建!")
print(f"  器件连接点数: {len(device_points)}")
print(f"  连接点标签: {list(device_points.keys())[:5]}... (共{len(device_points)}个)")
print()

# 保存步骤1的GDS
step1_path = os.path.join(output_dir, "step1_device_with_leads.gds")
lib_device.write_gds(step1_path)
print(f"✓ 步骤1 GDS已保存: {step1_path}")
print()

# 可视化 (不显示图形，只打印信息)
print("图形可视化已跳过 (需要图形环境)")
print()

# ========================================
# 步骤2: 创建并可视化Pad框架
# ========================================
print("=" * 70)
print("步骤2: 创建Pad框架")
print("=" * 70)

ACTIVE_SIZE = 200

lib_pads, cell_pads, pad_points, active_center, all_pads_info = generator.create_pad_frame(
    active_size=ACTIVE_SIZE,
    layout_width=1400,
    layout_height=1400,
    pad_width=100,
    pad_height=100,
    pad_spacing=30,
    edge_margin=50,
    trace_width=10,
    trace_spacing=10
)

print("Pad框架已创建!")
print(f"  有源区中心: {active_center}")
print(f"  Pad数量: {len(all_pads_info)}")
print(f"  Pad标签: {[p['label'] for p in all_pads_info[:5]]}... (共{len(all_pads_info)}个)")
print()

# 保存步骤2的GDS
step2_path = os.path.join(output_dir, "step2_pad_frame.gds")
lib_pads.write_gds(step2_path)
print(f"✓ 步骤2 GDS已保存: {step2_path}")
print()

# ========================================
# 步骤3: 组装器件和Pad框架
# ========================================
print("=" * 70)
print("步骤3: 组装器件和Pad框架")
print("=" * 70)

# 重新初始化生成器用于完整流程
generator2 = SixDotQuantumDeviceGenerator(output_dir=output_dir)

# 生成完整器件
lib_assembly, cell_assembly = generator2.generate_full_device(active_size=ACTIVE_SIZE)

print("器件组装完成!")
print(f"  活动中心: {generator2.active_center}")
print(f"  器件连接点数: {len(generator2.device_points)}")
print(f"  Pad连接点数: {len(generator2.pad_points)}")
print()

# 保存步骤3的GDS
step3_path = os.path.join(output_dir, "step3_placed_unrouted.gds")
lib_assembly.write_gds(step3_path)
print(f"✓ 步骤3 GDS已保存: {step3_path}")
print()

# ========================================
# 步骤4: 单线布线演示
# ========================================
print("=" * 70)
print("步骤4: 单线布线演示")
print("=" * 70)

label_to_connect = 'QD_PG4'
print(f"正在布设单条线: {label_to_connect}")

generator2.route_all_connections(single_label=label_to_connect)
print(f"✓ 单条线 {label_to_connect} 已布设完成!")
print()

# 保存步骤4的GDS
step4_path = os.path.join(output_dir, "step4_single_route.gds")
generator2.lib.write_gds(step4_path)
print(f"✓ 步骤4 GDS已保存: {step4_path}")
print()

# ========================================
# 步骤5: 完整布线
# ========================================
print("=" * 70)
print("步骤5: 完整布线")
print("=" * 70)

print("正在使用分层Z型扇出算法连接所有线...")
generator2.route_all_connections()
print("✓ 所有线已布设完成!")
print()

# 保存步骤5的GDS
step5_path = os.path.join(output_dir, "step5_fully_routed.gds")
generator2.lib.write_gds(step5_path)
print(f"✓ 步骤5 GDS已保存: {step5_path}")
print()

# ========================================
# 总结与输出文件列表
# ========================================
print("=" * 70)
print("生成的GDS文件列表:")
print("=" * 70)

gds_files = [
    ("step1_device_with_leads.gds", "量子点器件核心 (带长引线)"),
    ("step2_pad_frame.gds", "Pad框架 (带引线)"),
    ("step3_placed_unrouted.gds", "器件已放置 (未布线)"),
    ("step4_single_route.gds", "单线布线演示"),
    ("step5_fully_routed.gds", "最终完整布线")
]

total_size = 0
for filename, description in gds_files:
    filepath = os.path.join(output_dir, filename)
    if os.path.exists(filepath):
        size = os.path.getsize(filepath)
        total_size += size
        print(f"✓ {filename}")
        print(f"    描述: {description}")
        print(f"    大小: {size:,} bytes")
        print(f"    路径: {os.path.abspath(filepath)}")
        print()

print("-" * 70)
print(f"总大小: {total_size:,} bytes")
print("=" * 70)
print()
print("全部分步流程执行完毕!")
print(f"所有GDS文件已保存到: {os.path.abspath(output_dir)}")
print("可以在KLayout等EDA软件中打开查看。")
