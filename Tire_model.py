import numpy as np
import matplotlib.pyplot as plt

class Dynamics:
    
    def __init__(self, spring_type=None, spring_k=None, spring_F=None, spring_non_lin_coef=None,tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, damper_type=None, damper_V=None, damper_F_viscous=None, damper_K_friction=None, damper_F_static=None):
        # Modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # rigidez da mola [N/m]
        self.spring_F = spring_F  # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # coeficiente de ganho não-linear
        # Modelo de pneu
        self.tire_Fz = tire_Fz  # carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # slip angle do pneu [rad]
        self.tire_Ls = tire_Ls  # longitudinal slip do pneu [Admensional]
        self.tire_type = 'Default'
        self.tire_friction_coef = tire_friction_coef # coeficiente de fricção entre o pneu e a pista
        # Modelo de amortecedor
        self.damper_type = damper_type # Coulumb, Integrated, Stribeck
        self.damper_V = damper_V # velocidade relativa amortecedor [m/s]
        self.damper_F_viscous = damper_F_viscous # força viscosa do fluído [N]
        self.damper_F_static = damper_F_static # coeficiente de fricção estática de coulumb [N]
        self.damper_K_friction = damper_K_friction # rigidez de fricção [N/m]

    def Tire(self, params):
        E, Cy, Cx, Cz, c1, c2 = params
        Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))
        D = self.tire_friction_coef * self.tire_Fz
        Bz = Cs / (Cz * D)
        Bx = Cs / (Cx * D)
        By = Cs / (Cy * D)
        tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
        tire_longitudinal_force = D * np.sin(Cx * np.arctan(9 * Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
        tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))


        return tire_lateral_force, (10 + (tire_auto_align_moment/55)), tire_longitudinal_force

# Valores a serem usados na instância da classe Dynamics
slip_ratio = np.linspace(-1, 1, 1000)
slip_angles = np.linspace(-9, 9, 1000)

# Dados experimentais
ratio = np.linspace(-1, 1, 19)
angles = np.array([-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
tire_lateral_forces_1 = np.array([-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890, 2000, 2090])
tire_auto_align_moment_1 = np.array([-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09, 50.90, 50.10, 50.81, 48.12, 48.83])

# Instanciando a classe Dynamics
dynamics_instance = Dynamics(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45)

# Parâmetros de Pacejka 
result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]

# Imprimindo os parâmetros otimizados
print("Parâmetros otimizados:")
print("E:", result[0])
print("Cy:", result[1])
print("Cx:", result[2])
print("Cz:", result[3])
print("c1:", result[4])
print("c2:", result[5])

# Plotagem da curva otimizada com os dados experimentais
predicted_tire_lateral_forces, predicted_tire_auto_align_moment, predicted_tire_longitudinal_forces = dynamics_instance.Tire(result)

# Definindo um tamanho para a figura
plt.figure(figsize=(20, 7))

# Plotagem força lateral
plt.subplot(1, 3, 1)
plt.plot(slip_angles, predicted_tire_lateral_forces, label='Curva Otimizada')
plt.scatter(angles, tire_lateral_forces_1, color='red', label='Dados Experimentais')
plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
plt.ylabel('Força Lateral do Pneu (N)')
plt.title('Curva Otimizada com os Dados Experimentais')
plt.legend()
plt.grid(True)

# Plotagem torque auto-alinhante
plt.subplot(1, 3, 2)
plt.plot(slip_angles, predicted_tire_auto_align_moment, label='Curva Otimizada')
plt.scatter(angles, tire_auto_align_moment_1, color='blue', label='Dados Experimentais')
plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
plt.ylabel('Torque auto-alinhante (N.m)')
plt.title('Curva Otimizada com os Dados Experimentais')
plt.legend()
plt.grid(True)

# Plotagem força longitudinal
plt.subplot(1, 3, 3)
plt.plot(slip_ratio, predicted_tire_longitudinal_forces, label='Curva Sem Otimizar')
plt.xlabel('Taxa de Escorregamento Longitudinal (Admensional)')
plt.ylabel('Força Longitudinal (N)')
plt.title('Força Longitudinal - Sem Dados Experimentais')
plt.legend()
plt.grid(True)

plt.tight_layout(pad=3.0)  # Aumentando a distância entre os subplots
plt.show()
