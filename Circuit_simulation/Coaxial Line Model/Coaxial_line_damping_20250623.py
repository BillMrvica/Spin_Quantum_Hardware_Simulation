import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# --- 设置中文字体 ---
# 尝试使用系统自带的黑体或宋体，如果找不到，Matplotlib会回退到默认字体
try:
    font = FontProperties(fname='C:/Windows/Fonts/simhei.ttf', size=14)
except IOError:
    font = FontProperties(size=14)

# --- 主程序 ---
if __name__ == "__main__":
    
    # 1. 定义电路元件 (L和C的值固定，它们决定了系统的固有振荡频率)
    L = 1e-6  # 电感: 1 uH
    C = 100e-12 # 电容: 100 pF

    # 2. 定义三种阻尼状态的电阻R值和标签
    # 理论上的临界阻尼电阻 R_crit = 2 * sqrt(L/C) = 2 * sqrt(1e-6 / 100e-12) = 200 Ω
    cases = [
        {'R': 30,    'label': '欠阻尼 (Underdamped) - 产生振铃'},
        {'R': 200,   'label': '临界阻尼 (Critically Damped) - 理想状态'},
        {'R': 800,   'label': '过阻尼 (Overdamped) - 边沿迟缓'}
    ]

    # 3. 创建绘图窗口
    plt.figure(figsize=(12, 7))
    
    # 4. 循环处理每一种情况
    for case in cases:
        R = case['R']
        label = case['label']
        
        # 定义系统的传递函数 H(s) = V_out(s) / V_in(s)
        # 对于串联RLC电路，输出电压在电容两端，传递函数为:
        # H(s) = (1/LC) / (s^2 + (R/L)s + 1/LC)
        numerator = [1 / (L * C)]
        denominator = [1, R / L, 1 / (L * C)]
        
        # 创建一个LTI(线性时不变)系统对象
        system = signal.lti(numerator, denominator)
        
        # 计算并绘制系统对阶跃信号(从0到1的跳变)的响应
        # 设置一个足够长的时间轴来观察所有状态
        t = np.linspace(0, 1.5e-6, 2000) # 1.5 us
        t, yout = signal.step(system, T=t)
        
        # 绘制曲线
        plt.plot(t * 1e9, yout, label=label, linewidth=2.5)

        # 打印阻尼比zeta，用于验证
        zeta = (R / 2) * np.sqrt(C / L)
        print(f"当 R = {R} Ω, 阻尼比 ζ ≈ {zeta:.2f}")

    # 5. 美化图表
    plt.axhline(1.0, color='gray', linestyle='--', label='目标电压 (V_final)')
    plt.title('不同阻尼状态下的二阶系统阶跃响应', fontproperties=font, fontsize=18)
    plt.xlabel('时间 (ns)', fontproperties=font)
    plt.ylabel('输出电压 (V)', fontproperties=font)
    plt.grid(True, which='both', linestyle=':')
    plt.legend(prop=font)
    plt.ylim(-0.2, 1.8)
    plt.tight_layout()
    
    # 显示图像
    plt.show()