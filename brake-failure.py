import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
# =====================
# Parâmetros do veículo
# =====================
W = 15000  # peso do veículo [N]
R = 0.35   # raio efetivo do pneu [m]
r = 0.15   # raio do disco/tambor [m]
Awc = 0.0025  # área do cilindro de roda [m²]
BF = 1.8      # fator de freio (adimensional)
eta_c = 0.95  # eficiência do cilindro de roda
n_wheels = 4  # nº de rodas com freio
po = 0.0      # pressão de retorno [Pa]

# Parâmetros do cilindro mestre / pedal
Amc = 0.0015   # área do cilindro mestre [m²]
lp = 4.0       # relação da alavanca do pedal
eta_p = 0.9    # eficiência do pedal
B = 3.5        # fator de assistência
a = 0          # aceleração (analítico)
t = 0          # tempo (analítico)
# =====================
# Funções principais
# =====================

# Força de frenagem por roda
def brake_force(pt, po=0.0):
    return 2 * (pt - po) * Awc * eta_c * BF * (r / R)

# Desaceleração total do veículo
def vehicle_deceleration(pt, failed_wheels=0):
    active_wheels = n_wheels - failed_wheels
    Fx_total = active_wheels * brake_force(pt, po)
    return Fx_total / W  # em g-units podemos dividir por 9.81 depois
# Força no pedal (para gerar pt desejado)
def pedal_force(pt):
    return (pt * Amc * B) / (lp * eta_p)

#Velocidade em função do tempo (analítico)

def speed(t,v):
    return - vehicle_deceleration((5e6)* t)

t_eval = np.linspace(0, 10, 1000)
sol = solve_ivp(speed, [0, 10], [30.0], t_eval = t_eval)

# =====================
# Simulação
# =====================

pt_range = np.linspace(0, 6e6, 200)  # pressão de linha [Pa]

# Cenários
a_no_failure = [vehicle_deceleration(pt, failed_wheels=0) for pt in pt_range]
a_half_failure = [vehicle_deceleration(pt, failed_wheels=2) for pt in pt_range]
Fp_curve = [pedal_force(pt) for pt in pt_range]

# =====================
# Gráficos
# =====================

plt.figure(figsize=(12,5))

# Desaceleração
plt.subplot(1,2,1)
plt.plot(pt_range/1e6, np.array(a_no_failure)*9.81, label="Sem falhas")
plt.plot(pt_range/1e6, np.array(a_half_failure)*9.81, label="Falha em 2 rodas")
plt.xlabel("Pressão de linha [MPa]")
plt.ylabel("Desaceleração [m/s²]")
plt.legend()
plt.title("Desaceleração do veículo")

# Força no pedal
plt.subplot(1,2,2)
plt.plot(pt_range/1e6, Fp_curve, color="red")
plt.xlabel("Pressão de linha [MPa]")
plt.ylabel("Força no pedal [N]")
plt.title("Força do pedal vs pressão")

plt.tight_layout()
plt.show()

#Velocidade em função do tempo
plt.plot(sol.t, sol.y[0])
plt.xlabel("t [s]")
plt.ylabel("v(t) [m/s]")
plt.title("Velocidade em função do tempo (solve_ivp)")

plt.grid()
plt.show()
