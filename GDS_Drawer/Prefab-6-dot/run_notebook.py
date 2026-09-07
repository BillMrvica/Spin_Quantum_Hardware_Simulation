#!/usr/bin/env python3
# 运行Jupyter Notebook的脚本

import subprocess
import sys
import os

# 切换到正确的工作目录
notebook_path = r"C:\Users\Administrator\Desktop\Spin_Quantum_Hardware_Simulation\GDS_Drawer\Prefab-6-dot\six_dot_generator_demo.ipynb"
output_dir = r"C:\Users\Administrator\Desktop\Spin_Quantum_Hardware_Simulation\GDS_Drawer\Prefab-6-dot\output_gds"

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

# 使用jupyter nbconvert执行notebook
print("正在执行Jupyter Notebook...")
print(f"Notebook路径: {notebook_path}")
print("-" * 60)

# 执行notebook并生成输出
cmd = [
    "jupyter", "nbconvert",
    "--to", "notebook",
    "--execute",
    "--ExecutePreprocessor.timeout=300",
    "--output", notebook_path,
    notebook_path
]

result = subprocess.run(cmd, capture_output=True, text=True)

print("STDOUT:")
print(result.stdout)

if result.stderr:
    print("STDERR:")
    print(result.stderr)

print("-" * 60)
print(f"执行完成! 返回码: {result.returncode}")

if result.returncode == 0:
    print("\n所有步骤执行成功!")
    print(f"输出文件位置: {output_dir}")
    
    # 列出生成的文件
    print("\n生成的文件:")
    for f in os.listdir(output_dir):
        filepath = os.path.join(output_dir, f)
        size = os.path.getsize(filepath)
        print(f"  - {f} ({size:,} bytes)")
