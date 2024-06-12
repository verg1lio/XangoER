import math


class BrakeSystem:
    def __init__(self, params):
        self.params = params
        self.RedP = params['RedP']  # Redução do pedal
        self.a = params['a']  # Desaceleração [g]
        self.psi = params['psi']  # Distribuição de carga estática por eixo [adm]
        self.μl = params['μl']  # Coeficiente de atrito do contato pastilha/disco
        self.pi = params['pi']  # Valor de pi
        self.HCG = params['HCG']  # Altura do centro de gravidade [m]
        self.μ = params['μ']  # Coeficiente de atrito do pneu/solo
        self.FzF = params['FzF']  # Força de reação estática na dianteira [N]
        self.FzR = params['FzR']  # Força de reação estática na traseira [N]
        self.Rdp = params['Rdp']  # Raio do pneu [m]
        self.Rdr = params['Rdr'] # Raio da roda [m]
        self.Dcm = params['Dcm']  # Diâmetro do cilindro mestre em metros
        self.Dwc = params['Dwc']  # Diâmetro do cilindro da roda em metros
        self.Npast = params['Npast']  # Número de pastilhas por disco
        self.atrito_coeficiente = params['atrito_coeficiente']  # Coeficiente de atrito da pastilha
        self.red = params['red']  # Raio efetivo do disco [m]
        self.Mt = params['Mt']  # Massa total do veículo [kg]
        self.L = params['L']  # Distância entre eixos [m]
        self.m_wheel = params['m_wheel']  # Massa da roda [kg]

    def calculate_params(self, pedal_force):
        BF = 2 * self.μl  # Fator de freio
        χ = self.HCG / self.L  # Razão entre HCG e L
        W = self.Mt * 9.81  # Peso do veículo
        FzF_dyn = (1 - self.psi + self.a * χ) * W  # Força de reação dinâmica dianteira
        FzR_dyn = (self.psi - self.a * χ) * W  # Força de reação dinâmica traseira
        τF = FzF_dyn * self.μ * self.Rdp  # Torque na dianteira
        τR = FzR_dyn * self.μ * self.Rdp  # Torque na traseira
        FnF = τF / self.Npast * self.RedP * self.red  # Força normal das pastilhas dianteira
        FnR = τR / self.Npast * self.RedP * self.red  # Força normal das pastilhas traseira
        Awc = (self.pi * (self.Dwc ** 2)) / 4  # Área do cilindro de roda
        Acm = (self.pi * (self.Dcm ** 2)) / 4  # Área do cilindro mestre
        PF = FnF / Awc  # Pressão necessária na dianteira
        PR = FnR / Awc  # Pressão necessária na traseira
        FaCM = PF * Acm  # Força necessária para acionar o cilindro mestre
        lF = self.psi * self.L  # Distância do eixo ao centro de gravidade em relação a distribuição de carga dianteiro
        lR = (
                         1 - self.psi) * self.L  # Distância do eixo ao centro de gravidade em relação a distribuição de carga traseiro
        return BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR

    def apply_brake(self, pedal_force, initial_speed):
        resultados = self.calculate_params(pedal_force)
        BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR = resultados

        # Calculando a pressão do cilindro mestre
        pressao_cilindro_mestre = pedal_force / Acm

        # Calculando a pressurização do fluido
        pressao_fluido = pressao_cilindro_mestre

        # Calculando a transmissão de pressão da linha de freio
        transmissao_pressao = pressao_fluido * Awc

        # Calculando a força de aperto da pinça
        forca_aperto_pinca = transmissao_pressao

        # Calculando a força de atrito da pastilha de freio
        forca_atrito_pastilha = forca_aperto_pinca * self.atrito_coeficiente

        # Calculando o torque do disco de freio
        torque_disco_freio = forca_atrito_pastilha * self.red


        # Calculando a força de frenagem considerando todos os fatores
        forca_frenagem = (FnF + FnR) / self.Npast  # Força normal total
        forca_frenagem *= self.atrito_coeficiente  # Força de atrito total
        forca_frenagem *= self.red  # Torque total
        forca_frenagem /= self.Rdp  # Força de frenagem total

        # Calculando a velocidade angular do pneu a partir do torque de freio
        inertia_wheel = 0.5 * self.m_wheel * (self.Rdp ** 2 + self.Rdr ** 2) # Inércia do pneu
        angular_deceleration = torque_disco_freio / inertia_wheel  # Desaceleração angular
        angular_velocity = initial_speed / (self.Rdp + self.Rdr) # Velocidade angular inicial
        new_angular_velocity = angular_velocity - angular_deceleration  # Nova velocidade angular


        # Retornando os resultados calculados
        return resultados, forca_frenagem, new_angular_velocity


# Definindo os parâmetros
params = {
    'RedP': 4,  # Redução do pedal
    'a': 0.8,  # Desaceleração [g]
    'psi': 0.40,  # Distribuição de carga estática por eixo [adm]
    'μl': 0.45,  # Coeficiente de atrito do contato pastilha/disco
    'pi': 3.14,  # Valor de pi
    'HCG': 0.5,  # Altura do centro de gravidade [m]
    'μ': 0.60,  # Coeficiente de atrito do pneu/solo
    'FzF': 1471.5,  # Força de reação estática na dianteira [N]
    'FzR': 981.0,  # Força de reação estática na traseira [N]
    'Rdp': 0.30,  # Raio do pneu [m]
    'Rdr': 0.1651, # Raio da roda [m]
    'Dcm': 0.02,  # Diâmetro do cilindro mestre em metros
    'Dwc': 0.032,  # Diâmetro do cilindro da roda em metros
    'Npast': 2,  # Número de pastilhas por disco
    'atrito_coeficiente': 0.35,  # Coeficiente de atrito da pastilha
    'red': 0.75,  # Raio efetivo do disco [m]
    'Mt': 250,  # Massa total do veículo [kg]
    'm_wheel': 88,  # Massa da roda [kg] <-- estimativa
    'L': 1.5,  # Distância entre eixos [m]
}

# Criação da instância do sistema de freios e aplicação da força no pedal de freio
BrakeSystem = BrakeSystem(params)
pedal_force = 445  # N
initial_speed = 20 # m/s (Velocidade inicial do veículo)
resultados, forca_frenagem, new_angular_velocity = BrakeSystem.apply_brake(pedal_force, initial_speed)

# Exibindo os resultados calculados
print("Resultados Calculados:")
for i, result in enumerate(resultados):
    print(f"Resultado {i + 1}: {result}")
print("Força de frenagem:", forca_frenagem, "N")
print("Velocidade Angular:", new_angular_velocity, "rad/s")

import numpy as np
import matplotlib.pyplot as plt

class Dynamics:
    def __init__(self, spring_type, spring_k, spring_F, spring_non_lin_coef, tire_Fz, tire_Sa, tire_Ls):
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

    def Tire(self):
        # Pacejka parâmetros
        E = -2
        Cy = 1.4  # C para força lateral
        Cx = 1.65  # C para força longitudinal
        Cz = 2.4  # C para momento de torque auto-alinhante
        c1 = 54000
        c2 = 6600
        Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))

        if self.tire_type == 'Default':
            D = 1.45 * self.tire_Fz
            Bz = Cs / (Cz * D)
            Bx = Cs / (Cx * D)
            By = Cs / (Cy * D)
            tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
            tire_longitudinal_force = D * np.sin(Cx * np.arctan(Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
            tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))




        return tire_lateral_force, tire_longitudinal_force, (10 + (tire_auto_align_moment/55))


# Valores de tire_Fz para plotagem
tire_Fz_values = [1500]

# Instanciando a classe Dynamics para cada valor de tire_Fz
dynamics_list = [Dynamics(spring_type="Hooke", spring_k=1000, spring_F=500, spring_non_lin_coef=0.1, tire_Fz=tire_Fz, tire_Sa=5, tire_Ls=0.1) for tire_Fz in tire_Fz_values]

# Criando arrays de valores para os parâmetros relevantes
tire_Sa_values = np.linspace(-np.pi / 20, np.pi / 20, 1000)  # Valores de ângulo de deslizamento lateral
tire_Ls_values = np.linspace(-0.4, 0.4, 1000)  # Valores de slip longitudinal

# Calculando as variáveis retornadas pela função Tire() para cada conjunto de parâmetros e para cada valor de tire_Fz
tire_lateral_force_values_list = [[dynamics.Tire()[0] for dynamics.tire_Sa in tire_Sa_values] for dynamics in dynamics_list]
tire_longitudinal_force_values_list = [[dynamics.Tire()[1] for dynamics.tire_Ls in tire_Ls_values] for dynamics in dynamics_list]
tire_auto_align_moment_values_list = [[dynamics.Tire()[2] for dynamics.tire_Sa in tire_Sa_values] for dynamics in dynamics_list]

# Dados experimentais
angles = [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

forces_1 = [-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890, 2000, 2090]
torque_1 = [-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09, 50.90, 50.10, 50.81, 48.12, 48.83]

forces_2 = [-2930, -2750, -2530, -2310, -2050, -1760, -1410, -990, -480, 90, 620, 1080, 1450, 1760, 2010, 2230, 2420, 2590, 2750]
torque_2 = [-48.55, -49.34, -47.81, -45.55, -46.18, -43.39, -40.66, -28.43, -9.02, 17.80, 42.97, 59.03, 64.82, 67.67, 68.95, 67.69, 67.29, 64.64, 64.16]

tire_Sa = angles  # Ângulos de deslizamento lateral
tire_Ls = (new_angular_velocity * params['Rdp'] / initial_speed) - 1 # Slip ratio baseado na nova velocidade angular
print("Slip Ratio Longitudinal:", tire_Ls)

# Plotagem dos resultados
plt.figure(figsize=(15, 5))

# Plotagem da Força Lateral em função do Ângulo de Deslizamento Lateral
plt.subplot(1, 3, 1)
for i, tire_lateral_force_values in enumerate(tire_lateral_force_values_list):
    plt.plot(np.degrees(tire_Sa_values), tire_lateral_force_values, label=f'tire_Fz = {tire_Fz_values[i]}')
plt.scatter(angles, forces_1, color='black', label='Dados Experimentais 1500N')
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
plt.scatter(angles, torque_1, color='black', label='Dados Experimentais 1500N')
plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
plt.ylabel('Momento de Auto-Alinhamento do Pneu (N.m)')
plt.title('Torque Auto-Alinhante X Slip Angle')
plt.legend()

plt.tight_layout()
plt.show()
