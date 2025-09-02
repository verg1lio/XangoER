# -*- coding: utf-8 -*-
"""
Modelo de frenagem veicular em questão do tempo:
-Desaceleração
-Tensão nos discos
-Aquecimento dos discos
-Coeficiente de atrito efetivo
"""

import math
import numpy as np
import matplotlib.pyplot as plt

class BrakeSystem:
    def __init__(self, params):
        # Parâmetros do sistema de freio e veículo
        self.params = params
        self.RedP = params['RedP']                  # Redução mecânica do pedal
        self.μl = params['μl']                      # coef. atrito pastilha/disco nominal
        self.pi = params['pi']                      # π (valor mantido do código original)
        self.μ = params['μ']                        # coef. atrito pneu/solo
        self.Rdp = params['Rdp']                    # raio dinâmico do pneu [m]
        self.Dcm = params['Dcm']                    # diâmetro cilindro mestre [m]
        self.Dwc = params['Dwc']                    # diâmetro cilindro roda [m]
        self.Npast = params['Npast']                # nº de pastilhas
        self.atrito_coeficiente = params['atrito_coeficiente']  # μ nominal
        self.red = params['red']                    # raio efetivo do disco (raio médio de atrito) [m]
        self.Mt = params['Mt']                      # massa do veículo [kg]
        self.Ro = params.get('R_ext_disc', self.red)  # raio externo do disco [m]
        self.Ri = params.get('R_int_disc', 0.5*self.red)  # raio interno do disco [m]
        self.n_wheels = params.get('n_wheels', 4)   # nº de rodas freando

        # Parâmetros térmicos do disco
        self.m_disc = params.get('m_disc', 8.0)     # massa disco [kg]
        self.c_disc = params.get('c_disc', 500.0)   # calor específico [J/kgK]
        self.hA = params.get('hA', 60.0)            # dissipação térmica [W/K]
        self.T_amb = params.get('T_amb', 25.0)      # temperatura ambiente [°C]
        self.T_limit = params.get('T_limit', 400.0) # limite de fade [°C]
        self.fade_k = params.get('fade_k', 0.8)     # severidade da queda de μ com T

    def effective_mu(self, T):
        """
        Coeficiente de atrito efetivo em função da temperatura do disco.
        """
        if T <= self.T_limit:
            return self.atrito_coeficiente
        else:
            drop = self.fade_k * (T - self.T_limit) / self.T_limit
            return max(0.05, self.atrito_coeficiente * (1 - drop))

    def hydraulics_and_torque_per_wheel(self, pedal_force, mu_eff):
        """
        Converte a força no pedal em torque por disco:
        """
        Acm = (self.pi * (self.Dcm ** 2)) / 4
        Awc = (self.pi * (self.Dwc ** 2)) / 4
        p_mestre = pedal_force / Acm
        F_aperto = p_mestre * Awc
        F_atrito = F_aperto * mu_eff
        T_disc = 2 * F_atrito * self.red
        return T_disc

    def total_brake_force_from_pedal(self, pedal_force, mu_eff):
        """
        Torque por disco → força de frenagem na pista.
        """
        T_one = self.hydraulics_and_torque_per_wheel(pedal_force, mu_eff)
        F_one = T_one / self.Rdp
        F_total = self.n_wheels * F_one
        return F_total, T_one

    def disc_shear_stress(self, torque_disc):
        """
        Calcula a tensão de cisalhamento máxima no disco pelo modelo de torção.
        """
        Ro, Ri = self.Ro, self.Ri
        J = (math.pi / 2.0) * (Ro**4 - Ri**4)
        tau_max = torque_disc * Ro / J
        return tau_max, J


# ======================
# Parâmetros do sistema
# ======================
params = {
    'RedP': 4,
    'μl': 0.45,
    'pi': 3.14,
    'μ': 0.60,
    'Rdp': 0.30,
    'Dcm': 0.02,
    'Dwc': 0.032,
    'Npast': 2,
    'atrito_coeficiente': 0.35,
    'red': 0.15,
    'Mt': 250,
    'R_ext_disc': 0.16,
    'R_int_disc': 0.06,
    'n_wheels': 2,
    # térmicos
    'm_disc': 8.0,
    'c_disc': 500.0,
    'hA': 60.0,
    'T_amb': 25.0,
    'T_limit': 400.0,
    'fade_k': 0.8,
}

brake = BrakeSystem(params)

# ======================
# Simulação
# ======================
v0 = 30.0       # velocidade inicial (m/s)
t_max = 15.0    # tempo máximo simulação (s)
dt = 0.01
Fmax = 600.0    # força máxima no pedal
ramp_time = 0.8 # tempo para atingir Fmax

def pedal_force_of_t(t):
    if t <= 0: return 0.0
    elif t < ramp_time: return Fmax * (t / ramp_time)
    else: return Fmax

# Estados
t = 0.0
v = v0
T_disc = params['T_amb']

# Histórico
ts, vs, taus, Ts, mus = [], [], [], [], []

while t < t_max and v > 0:
    mu_eff = brake.effective_mu(T_disc)
    F_total, T_one = brake.total_brake_force_from_pedal(pedal_force_of_t(t), mu_eff)

    # Limite por atrito pneu/solo
    F_max_tire = params['μ'] * params['Mt'] * 9.81
    if F_total > F_max_tire:
        F_total = F_max_tire
        T_one = (F_total / brake.n_wheels) * brake.Rdp

    # Tensão no disco
    tau_max, _ = brake.disc_shear_stress(T_one)

    # Dinâmica do veículo
    a = F_total / params['Mt']
    v = max(0, v - a*dt)

    # Balanço térmico
    P_fric = F_total * v
    dTdt = (P_fric - brake.hA*(T_disc - brake.T_amb)) / (brake.m_disc * brake.c_disc)
    T_disc += dTdt * dt

    # Armazenar histórico
    ts.append(t); vs.append(v); taus.append(tau_max); Ts.append(T_disc); mus.append(mu_eff)
    t += dt

ts, vs, taus, Ts, mus = map(np.array, [ts, vs, taus, Ts, mus])

# ======================
# Gráficos
# ======================

# Velocidade e Tensão lado a lado
fig, ax = plt.subplots(1, 2, figsize=(12,5))

ax[0].plot(ts, vs, color='blue')
ax[0].set_xlabel("Tempo [s]")
ax[0].set_ylabel("Velocidade [m/s]", color='blue')
ax[0].set_title("Velocidade do veículo até a parada")
ax[0].grid(True)

ax[1].plot(ts, taus, color='red')
ax[1].set_xlabel("Tempo [s]")
ax[1].set_ylabel("Tensão no disco [Pa]", color='red')
ax[1].set_title("Tensão de cisalhamento no disco")
ax[1].grid(True)

plt.tight_layout()
plt.show()

# Coeficiente de atrito efetivo
fig, ax = plt.subplots(1, 2, figsize=(12,5))
ax[0].plot(ts, mus, color='green')
ax[0].set_xlabel("Tempo [s]")
ax[0].set_ylabel("Coef. de atrito efetivo μ")
ax[0].set_title("Variação do coeficiente de atrito (fade)")
ax[0].grid(True)

# Temperatura do disco
ax[1].plot(ts, Ts, label="Temperatura do disco")
ax[1].axhline(params['T_limit'], color='r', linestyle='--', label="Limite de fade")
ax[1].set_xlabel("Tempo [s]")
ax[1].set_ylabel("Temperatura [°C]")
ax[1].set_title("Aquecimento do disco de freio")
ax[1].legend()
ax[1].grid(True)
plt.show()
