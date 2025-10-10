import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Parâmetros iniciais
f_out = 60  # Frequência de saída desejada (Hz)
V_dc = 400  # Tensão de entrada DC (Volts)
T_pwm_initial = 0.0001  # Período do PWM inicial (segundos)

# Função para gerar sinais e plotar
def plot_signals(T_pwm):
    # Calculando a frequência de chaveamento
    f_chaveamento = 1 / T_pwm  # Frequência de chaveamento (Hz)

    # Tempo de simulação
    t_final = 0.1  # Tempo final da simulação (segundos)
    t = np.arange(0, t_final, T_pwm)  # Vetor de tempo

    # Geração da onda senoidal de referência
    V_out_peak = V_dc / 2  # Amplitude de pico da tensão de saída CA
    omega = 2 * np.pi * f_out  # Frequência angular
    V_ref = V_out_peak * np.sin(omega * t)  # Onda senoidal de referência

    # Geração do sinal PWM
    V_carrier = V_out_peak * np.sign(np.sin(2 * np.pi * f_chaveamento * t))

    # Sinal PWM
    PWM_signal = np.where(V_ref >= V_carrier, V_dc, 0)

    # Plotando os sinais
    ax1.clear()
    ax1.plot(t, V_ref, label='Onda Senoidal de Referência')
    ax1.set_title('Onda Senoidal de Referência')
    ax1.set_xlabel('Tempo (s)')
    ax1.set_ylabel('Tensão (V)')
    ax1.legend()

    ax2.clear()
    ax2.plot(t, V_carrier, label='Onda Portadora')
    ax2.set_title('Onda Portadora')
    ax2.set_xlabel('Tempo (s)')
    ax2.set_ylabel('Tensão (V)')
    ax2.legend()

    ax3.clear()
    ax3.plot(t, PWM_signal, label='Sinal PWM')
    ax3.set_title('Sinal PWM')
    ax3.set_xlabel('Tempo (s)')
    ax3.set_ylabel('Tensão (V)')
    ax3.legend()

    fig.canvas.draw_idle()

# Configuração inicial do gráfico
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8))
plt.subplots_adjust(left=0.1, bottom=0.25)

# Plot inicial
plot_signals(T_pwm_initial)

# Eixo para o slider de controle de T_pwm
ax_slider = plt.axes([0.1, 0.1, 0.8, 0.03], facecolor='lightgoldenrodyellow')
T_pwm_slider = Slider(ax_slider, 'T_pwm (s)', 0.00001, 0.1, valinit=T_pwm_initial)

# Função de atualização do slider
def update(val):
    T_pwm = T_pwm_slider.val
    plot_signals(T_pwm)

T_pwm_slider.on_changed(update)

plt.show()
