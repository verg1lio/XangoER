import numpy as np
import matplotlib.pyplot as plt


class Dynamics:
    def __init__(self, spring_type, spring_k, spring_F, spring_non_lin_coef, tire_Fz, tire_Sa, tire_Ls, tire_type):
        # Modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # rigidez da mola [N/m]
        self.spring_F = spring_F  # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # coeficiente de ganho não-linear
        # Modelo de pneu
        self.tire_Fz = tire_Fz  # carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # slip angle do pneu [rad]
        self.tire_Ls = tire_Ls  # longitudinal slip do pneu [Admensional]
        self.tire_type = tire_type  # Default, Admensional

    def Tire(self):
        # Pacejka parâmetros
        E = -3
        Cy = 1.30  # C para força lateral
        Cx = 1.65  # C para força longitudinal
        Cz = 2.40  # C para momento de torque auto-alinhante
        c1 = 60000
        c2 = 4000
        Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))

        if self.tire_type == 'Default':
            D = 0.85 * self.tire_Fz
            Bz = Cs / (Cz * D)
            Bx = Cs / (Cx * D)
            By = Cs / (Cy * D)
            tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
            tire_longitudinal_force = D * np.sin(Cx * np.arctan(Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
            tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))

        elif self.tire_type == 'Admensional':
            D = 1
            Bz = Cs / 10000*(Cz * D)
            Bx = Cs / 10000*(Cx * D)
            By = Cs / 10000*(Cy * D)
            tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
            tire_longitudinal_force = D * np.sin(Cx * np.arctan(Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
            tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))


        return tire_lateral_force, tire_longitudinal_force, tire_auto_align_moment


# Valores de tire_Fz para plotagem
tire_Fz_values = [1000, 1500, 2000, 2500]

# Instanciando a classe Dynamics para cada valor de tire_Fz
dynamics_list = [Dynamics(spring_type="Hooke", spring_k=1000, spring_F=500, spring_non_lin_coef=0.1, tire_Fz=tire_Fz, tire_Sa=5, tire_Ls=0.1, tire_type="Default") for tire_Fz in tire_Fz_values]

# Criando arrays de valores para os parâmetros relevantes
tire_Sa_values = np.linspace(-np.pi / 10, np.pi / 10, 1000)  # Valores de ângulo de deslizamento lateral
tire_Ls_values = np.linspace(-1, 1, 1000)  # Valores de slip longitudinal

# Calculando as variáveis retornadas pela função Tire() para cada conjunto de parâmetros e para cada valor de tire_Fz
tire_lateral_force_values_list = [[dynamics.Tire()[0] for dynamics.tire_Sa in tire_Sa_values] for dynamics in dynamics_list]
tire_longitudinal_force_values_list = [[dynamics.Tire()[1] for dynamics.tire_Ls in tire_Ls_values] for dynamics in dynamics_list]
tire_auto_align_moment_values_list = [[dynamics.Tire()[2] for dynamics.tire_Sa in tire_Sa_values] for dynamics in dynamics_list]

# Plotagem dos resultados
plt.figure(figsize=(15, 5))

# Plotagem da Força Lateral em função do Ângulo de Deslizamento Lateral
plt.subplot(1, 3, 1)
for i, tire_lateral_force_values in enumerate(tire_lateral_force_values_list):
    plt.plot(np.degrees(tire_Sa_values), tire_lateral_force_values, label=f'tire_Fz = {tire_Fz_values[i]}')
plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
plt.ylabel('Força Lateral do Pneu (N)')
plt.title('Força Lateral X Slip Angle')
plt.legend()

# Plotagem da Força Longitudinal em função do Slip Longitudinal
plt.subplot(1, 3, 2)
for i, tire_longitudinal_force_values in enumerate(tire_longitudinal_force_values_list):
    plt.plot(tire_Ls_values, tire_longitudinal_force_values, label=f'tire_Fz = {tire_Fz_values[i]}')
plt.xlabel('Slip Longitudinal')
plt.ylabel('Força Longitudinal do Pneu (N)')
plt.title('Força Longitudinal X Longitudinal Slip')
plt.legend()

# Plotagem do Momento de Auto-Alinhamento em função do Ângulo de Deslizamento Lateral
plt.subplot(1, 3, 3)
for i, tire_auto_align_moment_values in enumerate(tire_auto_align_moment_values_list):
    plt.plot(np.degrees(tire_Sa_values), tire_auto_align_moment_values, label=f'tire_Fz = {tire_Fz_values[i]}')
plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
plt.ylabel('Momento de Auto-Alinhamento do Pneu (N)')
plt.title('Torque Auto-Alinhante X Slip Angle')
plt.legend()

plt.tight_layout()
plt.show()
