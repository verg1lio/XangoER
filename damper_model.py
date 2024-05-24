import numpy as np
import matplotlib.pyplot as plt

class Dynamics:
    def __init__(self, spring_type=None, spring_k=None, spring_F=None, spring_non_lin_coef=None,tire_Fz=None, tire_Sa=None, tire_Ls=None,damper_type=None, damper_V=None, damper_F_viscous=None, damper_K_friction=None, damper_F_static=None):
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
        # Modelo de amortecedor
        self.damper_type = damper_type # Coulumb, Integrated, Stribeck
        self.damper_V = damper_V # velocidade relativa amortecedor [m/s]
        self.damper_F_viscous = damper_F_viscous # força viscosa do fluído [N]
        self.damper_F_static = damper_F_static # coeficiente de fricção estática de coulumb [N]
        self.damper_K_friction = damper_K_friction # rigidez de fricção [N/m]



    def Spring(self):

        if self.spring_type == 'Hooke':
            spring_x = self.spring_F / self.spring_k
        if self.spring_type == 'Softening':
            spring_x = self.spring_F / (self.spring_non_lin_coef * (self.spring_k) ** 2)

        return spring_x

    def Damper(self):

        if self.damper_type == 'Coulumb':
            damper_F = self.damper_F_static * np.tanh(self.damper_K_friction * self.damper_V)
        if self.damper_type == 'Integrated':
            damper_F = self.damper_F_static * np.tanh(self.damper_K_friction * self.damper_V) + np.sign(self.damper_V) * self.damper_F_viscous * (self.damper_V ** 2)

        return damper_F

#Teste Damper

damper_V_values = np.linspace(-5, 5, 1000)
Teste = Dynamics(damper_type='Integrated', damper_V=damper_V_values, damper_F_static=50, damper_K_friction=10, damper_F_viscous=10)
Plot = Teste.Damper()
plt.figure(figsize=(8, 6))
plt.plot(damper_V_values, Plot, label='Velocidade')
plt.title('Força em função da velocidade')
plt.xlabel('Velocidade[m/s]')
plt.ylabel('Força[N]')
plt.grid(True)
plt.legend()
plt.show()

Teste_2 = Dynamics(damper_type='Coulumb', damper_V= damper_V_values, damper_F_static=50, damper_K_friction=10)
Plot = Teste_2.Damper()
plt.figure(figsize=(8, 6))
plt.plot(damper_V_values, Plot, label='Velocidade')
plt.title('Força em função da velocidade')
plt.xlabel('Velocidade[m/s]')
plt.ylabel('Força[N]')
plt.grid(True)
plt.legend()
plt.show()
