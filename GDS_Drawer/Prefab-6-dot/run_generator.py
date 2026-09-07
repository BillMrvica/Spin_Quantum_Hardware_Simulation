# -*- coding: utf-8 -*-
from six_dot_generator import SixDotQuantumDeviceGenerator
import os

# 创建输出目录
output_dir = './output_gds'
os.makedirs(output_dir, exist_ok=True)

# 初始化生成器
generator = SixDotQuantumDeviceGenerator(output_dir=output_dir)

# 步骤1: 创建量子点器件
print('步骤1: 创建量子点器件...')
lib_device, cell_device, device_points = generator.create_quantum_dot_device()
lib_device.write_gds(os.path.join(output_dir, 'step1_device_with_leads.gds'))
print(f'  器件点数: {len(device_points)}')

# 步骤2: 创建Pad框架
print('步骤2: 创建Pad框架...')
ACTIVE_SIZE = 200
lib_pads, cell_pads, pad_points, active_center, all_pads_info = generator.create_pad_frame(active_size=ACTIVE_SIZE)
lib_pads.write_gds(os.path.join(output_dir, 'step2_pad_frame.gds'))
print(f'  Pad数量: {len(all_pads_info)}')

# 步骤3: 组装
print('步骤3: 组装器件...')
generator2 = SixDotQuantumDeviceGenerator(output_dir=output_dir)
lib_assembly, cell_assembly = generator2.generate_full_device(active_size=ACTIVE_SIZE)
lib_assembly.write_gds(os.path.join(output_dir, 'step3_placed_unrouted.gds'))

# 步骤4: 单线布线
print('步骤4: 单线布线...')
generator2.route_all_connections(single_label='QD_PG4')
lib_assembly.write_gds(os.path.join(output_dir, 'step4_single_route.gds'))

# 步骤5: 完整布线
print('步骤5: 完整布线...')
generator2.route_all_connections()
lib_assembly.write_gds(os.path.join(output_dir, 'step5_fully_routed.gds'))

print('')
print('所有GDS文件生成完成!')
print(f'输出目录: {os.path.abspath(output_dir)}')
