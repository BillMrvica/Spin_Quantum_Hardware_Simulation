import pandas as pd
import numpy as np

def parse_resistance(r_str):
    """将字符串格式的电阻值（如 '1.118 M' 或 '5.09 k'）转换为欧姆数值。"""
    if not isinstance(r_str, str) or r_str.strip() == "":
        return 0.0
    
    r_str = r_str.strip().lower()
    val = float(r_str.split()[0])
    
    if 'm' in r_str:
        return val * 1_000_000
    if 'k' in r_str:
        return val * 1_000
    return val

def format_resistance(value):
    """将欧姆数值格式化为带有单位的字符串。"""
    if value == 0:
        return "0"
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.4f} M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:.4f} k"
    return f"{value:.4f}"

# 步骤 1: 从第一张图片手动提取数据
# 创建一个从原始电极名称到电阻字符串的映射
raw_data = {
    "SET_B2": "1.118 M",
    "SET_D": "5.09 k",
    "QD_D": "5.1 k",
    "QD_B1": "1.108 M",
    "SG3": "1.09 M",
    "QD_PG1": "1.095 M",
    "QD_B2": "1.096 M",
    "QD_PG2": "1.107 M",
    "QD_B3": "1.113 M",
    # 没有电阻值的电极默认为0
    "DQD1": "", "SET_S": "", "SET_G": "", "SET_B1": "",
    "SG1": "", "SG2": "", "QD_S": ""
}

# 步骤 2: 定义矩阵的行和列名称 (来自第二张图片)
# 注意：我将根据图片中的名称进行匹配，例如将 V_DQD_PG1 映射到 QD_PG1
matrix_electrodes = [
    "V_SET_G", "V_SET_B1", "V_SET_B2", 
    "V_DQD_PG1", "V_DQD_PG2", 
    "V_DQD_B1", "V_DQD_B2", "V_DQD_B3",
    "V_SG1", "V_SG2", "V_SG3"
]

# 步骤 3: 创建从矩阵电极名称到其电阻值（欧姆）的最终映射
# 这里处理名称不一致的问题 (例如 V_DQD_PG1 -> QD_PG1)
resistance_map_ohms = {}
for electrode in matrix_electrodes:
    # 移除前缀 'V_' 并将 'DQD' 替换为 'QD' 以匹配 raw_data 中的键
    key = electrode.replace("V_", "").replace("DQD_", "QD_")
    resistance_map_ohms[electrode] = parse_resistance(raw_data.get(key, "0"))

# 步骤 4: 计算寄生电阻矩阵
# 初始化一个空的矩阵
matrix_size = len(matrix_electrodes)
resistance_matrix_values = np.zeros((matrix_size, matrix_size))

# 填充矩阵
for i in range(matrix_size):
    for j in range(matrix_size):
        electrode_i = matrix_electrodes[i]
        electrode_j = matrix_electrodes[j]
        
        r_i = resistance_map_ohms[electrode_i]
        r_j = resistance_map_ohms[electrode_j]
        
        if i == j:
            # 对角线元素为电极自身的电阻
            resistance_matrix_values[i, j] = r_i
        else:
            # 非对角线元素为两者之和
            resistance_matrix_values[i, j] = r_i + r_j

# 步骤 5: 将计算出的数值矩阵格式化为字符串矩阵
formatted_matrix = [[format_resistance(val) for val in row] for row in resistance_matrix_values]

# 步骤 6: 使用 pandas 创建 DataFrame 并保存到 Excel
df = pd.DataFrame(formatted_matrix, index=matrix_electrodes, columns=matrix_electrodes)

output_filename = "parasitic_resistance_matrix.xlsx"
df.to_excel(output_filename)

print(f"成功生成Excel文件: '{output_filename}'")
print("\n矩阵预览:")
print(df)
