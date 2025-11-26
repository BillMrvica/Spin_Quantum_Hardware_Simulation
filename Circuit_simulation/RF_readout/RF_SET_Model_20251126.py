import numpy as np
import matplotlib.pyplot as plt

def calculate_s11_quantum_dot(frequencies, R_load_state, line_length_mm=10):
    """
    计算特定负载状态下的 S11 反射系数。
    
    参数:
    - frequencies: 频率数组 (Hz)
    - R_load_state: 量子点等效电阻 (Ohm)
    - line_length_mm: 电感与量子点之间的传输线长度 (毫米)
    """
    
    # --- 1. 常数与组件定义 ---
    Z0_sys = 50.0       # 系统阻抗
    
    # Bias-T 参数
    C_bias = 100e-9     # 100 nF (Series)
    R_bias = 5100.0     # 5.1 kOhm (Shunt)
    
    # 匹配网络/谐振器参数
    L_c = 1e-6          # 1 uH (Series Inductor)
    
    # 传输线 (CPW) 参数
    Z0_line = 50.0      # 传输线阻抗
    # 假设有效介电常数 (Silicon/PCB), 影响波速
    # v_p = c / sqrt(epsilon_eff). 假设 epsilon_eff approx 6 (FR4/Si mix) -> v_p approx 1.2e8
    # 或者简化设为 2/3 c
    v_phase = 3e8 * 0.66 
    
    # 量子点/器件参数
    C_p = 0.6e-12       # 0.6 pF (Shunt)
    
    # 角频率
    w = 2 * np.pi * frequencies
    w[0] = 1e-9 # 避免除以零
    
    # --- 2. 阻抗计算 (从右向左推导) ---
    
    # A. 负载点阻抗 Z_load (最右端: Cp 并联 Rs)
    # Z_load = 1 / (1/Rs + jwCp)
    Y_load = (1/R_load_state) + (1j * w * C_p)
    Z_load = 1 / Y_load
    
    # B. 传输线变换 (CPW)
    # 将 Z_load 通过传输线变换到电感右侧
    # beta = w / v_phase
    beta = w / v_phase
    l_meters = line_length_mm / 1000.0
    
    # 传输线输入阻抗公式 (Lossless)
    # Z_in(l) = Z0 * (ZL + j*Z0*tan(bl)) / (Z0 + j*ZL*tan(bl))
    tan_bl = np.tan(beta * l_meters)
    num = Z_load + 1j * Z0_line * tan_bl
    den = Z0_line + 1j * Z_load * tan_bl
    Z_after_line = Z0_line * (num / den)
    
    # C. 串联电感 Lc
    Z_after_Lc = Z_after_line + (1j * w * L_c)
    
    # D. Bias-T 分流电阻 (PC point)
    # 此处 Bias-T 电阻是并联到地的
    Y_after_Lc = 1 / Z_after_Lc
    Y_bias = 1 / R_bias
    Y_node_bias = Y_after_Lc + Y_bias
    Z_node_bias = 1 / Y_node_bias
    
    # E. Bias-T 串联电容 (Si-Cap) -> 输入阻抗 Z_in
    Z_in_total = Z_node_bias + (1 / (1j * w * C_bias))
    
    # --- 3. 计算反射系数 Gamma (S11) ---
    gamma = (Z_in_total - Z0_sys) / (Z_in_total + Z0_sys)
    
    return gamma, Z_in_total

# --- 主程序：模拟与绘图 ---

# 1. 设置频率范围
# 谐振频率估算 f = 1/(2*pi*sqrt(L*C)) = 1/(2*pi*sqrt(1uH * 0.6pF)) ≈ 205 MHz
# 我们扫描 10MHz 到 260MHz
freqs = np.linspace(10e6, 260e6, 2000)

# 2. 设置传输线长度 (假设芯片上的走线长度为 20mm)
cpw_length = 0.1 # mm

# 3. 计算两种状态
# 状态 A: 自旋阻塞/空态 (电阻无穷大)
gamma_inf, Z_inf = calculate_s11_quantum_dot(freqs, R_load_state=10e6, line_length_mm=cpw_length)
s11_inf_db = 20 * np.log10(np.abs(gamma_inf))
phase_inf = np.angle(gamma_inf, deg=True)

# 状态 B: 自旋隧穿态 (电阻 50 kOhm)
gamma_50k, Z_50k = calculate_s11_quantum_dot(freqs, R_load_state=50000.0, line_length_mm=cpw_length)
s11_50k_db = 20 * np.log10(np.abs(gamma_50k))
phase_50k = np.angle(gamma_50k, deg=True)

# 4. 绘图结果
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# 幅度响应 (Magnitude)
ax1.set_title(f"RF Reflectometry Response (LC Resonator + {cpw_length}mm CPW)")
ax1.plot(freqs/1e6, s11_inf_db, 'b-', label='State 0: R_load = $\infty$ (High Q)')
ax1.plot(freqs/1e6, s11_50k_db, 'r--', label='State 1: R_load = 50 k$\Omega$ (Damped)')
ax1.set_ylabel("Reflection Magnitude |S11| (dB)")
ax1.grid(True, which='both', alpha=0.6)
ax1.legend()

# 相位响应 (Phase)
ax2.plot(freqs/1e6, phase_inf, 'b-', label='State 0 Phase')
ax2.plot(freqs/1e6, phase_50k, 'r--', label='State 1 Phase')
ax2.set_ylabel("Phase (Degrees)")
ax2.set_xlabel("Frequency (MHz)")
ax2.grid(True, which='both', alpha=0.6)

# 标注谐振点
f_res_idx = np.argmin(s11_inf_db)
f_res = freqs[f_res_idx]
ax1.axvline(f_res/1e6, color='k', linestyle=':', alpha=0.5)
ax2.axvline(f_res/1e6, color='k', linestyle=':', alpha=0.5)

plt.tight_layout()
plt.show()

# 5. 计算灵敏度 (差异)
# RF 读出通常测量特定频率点的反射差异 (Magnitude or Phase contrast)
diff_mag = np.abs(s11_inf_db - s11_50k_db)
diff_phase = np.abs(phase_inf - phase_50k)

max_diff_mag = np.max(diff_mag)
max_diff_phase = np.max(diff_phase)

print(f"Resonance Frequency: {f_res/1e6:.2f} MHz")
print(f"Max Magnitude Contrast: {max_diff_mag:.2f} dB")
print(f"Max Phase Contrast: {max_diff_phase:.2f} degrees")