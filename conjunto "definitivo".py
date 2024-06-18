import numpy as np
import matplotlib.pyplot as plt

class BrakeSystem:
    def __init__(self, params):
        self.params = params
        self.RedP = params['RedP']
        self.a = params['a']
        self.psi = params['psi']
        self.μl = params['μl']
        self.pi = params['pi']
        self.HCG = params['HCG']
        self.μ = params['μ']
        self.FzF = params['FzF']
        self.FzR = params['FzR']
        self.Rdp = params['Rdp']
        self.Rdr = params['Rdr']
        self.Dcm = params['Dcm']
        self.Dwc = params['Dwc']
        self.Npast = params['Npast']
        self.atrito_coeficiente = params['atrito_coeficiente']
        self.red = params['red']
        self.Mt = params['Mt']
        self.L = params['L']
        self.m_wheel = params['m_wheel']
        self.m_tire = params['m_tire']
        self.c_rr = params['c_rr']

    def calculate_params(self, pedal_force):
        BF = 2 * self.μl
        χ = self.HCG / self.L
        W = self.Mt * 9.81
        FzF_dyn = (1 - self.psi + self.a * χ) * W
        FzR_dyn = (self.psi - self.a * χ) * W
        τF = FzF_dyn * self.μ * self.Rdp
        τR = FzR_dyn * self.μ * self.Rdp
        FnF = τF / self.Npast * self.RedP * self.red
        FnR = τR / self.Npast * self.RedP * self.red
        Awc = (self.pi * (self.Dwc ** 2)) / 4
        Acm = (self.pi * (self.Dcm ** 2)) / 4
        PF = FnF / Awc
        PR = FnR / Awc
        FaCM = PF * Acm
        lF = self.psi * self.L
        lR = (1 - self.psi) * self.L
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

        # Calculando a resistência ao rolamento
        resistencia_rolamento = self.c_rr * W

        # Calculando o torque gerado pela resistência ao rolamento
        torque_resistencia_rolamento = resistencia_rolamento * (self.Rdp + self.Rdr)

        # Calculando o torque de freio ajustado considerando a resistência ao rolamento
        torque_ajustado = torque_disco_freio - torque_resistencia_rolamento

        # Calculando a força gerada pelo disco de freio ajustado
        forca_f = torque_ajustado / self.Rdp

        # Calculando a força de frenagem considerando todos os fatores
        forca_frenagem = (FnF + FnR) / self.Npast
        forca_frenagem *= self.atrito_coeficiente
        forca_frenagem *= self.red
        forca_frenagem /= self.Rdp
        forca_frenagem -= resistencia_rolamento

        # Calculando a desaceleração linear
        desaceleracao_linear = forca_frenagem / self.Mt

        # Calculando a velocidade angular do pneu a partir do torque de freio ajustado
        inertia_wheel = 0.5 * (self.m_tire + self.m_wheel) * (self.Rdr * 2 + self.Rdp * 2)
        angular_deceleration = (torque_ajustado / 4) / inertia_wheel
        initial_angular_velocity = initial_speed / (self.Rdp + self.Rdr)
        time_step = time_intervals[1] - time_intervals[0]
        angular_velocity = initial_angular_velocity
        angular_velocities = [angular_velocity]

        for t in time_intervals:
            angular_velocity -= angular_deceleration * time_step
            angular_velocities.append(angular_velocity)

        angular_velocity = initial_angular_velocity + angular_velocity

        # Retornando os resultados calculados
        return resultados, forca_f, inertia_wheel, torque_ajustado, angular_deceleration, angular_velocity, angular_velocities, torque_resistencia_rolamento, desaceleracao_linear, torque_disco_freio, resistencia_rolamento

# Definindo os parâmetros
params = {
    'RedP': 4,  # Redução do pedal
    'a': 0.8,  # Desaceleração [g]
    'psi': 0.48,  # Distribuição de carga estática por eixo [adm]
    'μl': 0.45,  # Coeficiente de atrito do contato pastilha/disco
    'pi': 3.14,  # Valor de pi
    'HCG': 0.5,  # Altura do centro de gravidade [m]
    'μ': 1.5,  # Coeficiente de atrito do pneu/solo
    'FzF': 1471.5,  # Força de reação estática na dianteira [N]
    'FzR': 981.0,  # Força de reação estática na traseira [N]
    'Rdp': 0.15,  # Raio do pneu [m]
    'Rdr': 0.1651,  # Raio da roda [m]
    'Dcm': 0.02,  # Diâmetro do cilindro mestre em metros
    'Dwc': 0.032,  # Diâmetro do cilindro da roda em metros
    'Npast': 2,  # Número de pastilhas por disco
    'atrito_coeficiente': 0.35,  # Coeficiente de atrito da pastilha
    'red': 0.75,  # Raio efetivo do disco [m]
    'Mt': 250,  # Massa total do veículo [kg]
    'm_wheel': 40,  # Massa da roda [kg] <-- estimativa
    'm_tire': 40,  # Massa do pneu [kg] <-- estimativa
    'L': 1.5,  # Distância entre eixos [m]]
    'c_rr': 0.015  # Coeficiente de resistência ao rolamento
}


# Criação da instância do sistema de freios e aplicação da força no pedal de freio
BrakeSystem = BrakeSystem(params)
pedal_force = 445  # N
initial_speed = 20  # m/s (Velocidade inicial do veículo)
tempo = 2  # s (Tempo total de acionamento do pedal)
time_intervals = np.linspace(0, tempo, 100)
resultados, forca_f, inertia_wheel, torque_ajustado, angular_deceleration, angular_velocity, angular_velocities, torque_resistencia_rolamento, desaceleracao_linear, torque_dicos_freio, resistencia_rolamento = BrakeSystem.apply_brake(pedal_force, initial_speed)

# Exibindo os resultados calculados
print("Resultados Calculados:")
for i, result in enumerate(resultados):
    print(f"Resultado {i + 1}: {result}")
print("Força de frenagem:", forca_f, "N")
print("Aceleração tangencial:", BrakeSystem.apply_brake(pedal_force, initial_speed)[4]*(params['Rdr'] + params['Rdp']),
      "m/s*2")
print("Inércia da roda:", BrakeSystem.apply_brake(pedal_force, initial_speed)[2], "kg.m^2")
print("Torque do disco de freio:", BrakeSystem.apply_brake(pedal_force, initial_speed)[9], "N.m")
print("aceleração angular:", BrakeSystem.apply_brake(pedal_force, initial_speed)[4], "rad/s^2")
print("velocidade tangente:", BrakeSystem.apply_brake(pedal_force, initial_speed)[5] * (params['Rdr'] + params['Rdp']), "m/s")
print("Velocidade angular final:", angular_velocity, "rad/s")
print("Resistência ao rolamento:", resistencia_rolamento, "N")
print("Torque de resistência ao rolamento:", torque_resistencia_rolamento, "N.m")
print("Torque ajustado:", torque_ajustado, "N.m")

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
            tire_lateral_force = D * np.sin(
                Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
            tire_longitudinal_force = D * np.sin(
                Cx * np.arctan(Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
            tire_auto_align_moment = D * np.sin(
                Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))

        return tire_lateral_force, tire_longitudinal_force, (10 + (tire_auto_align_moment / 55))


# Valores de tire_Fz para plotagem
tire_Fz_values = [1500]

# Instanciando a classe Dynamics para cada valor de tire_Fz
dynamics_list = [
    Dynamics(spring_type="Hooke", spring_k=1000, spring_F=500, spring_non_lin_coef=0.1, tire_Fz=tire_Fz, tire_Sa=5,
             tire_Ls=0.1) for tire_Fz in tire_Fz_values]

# Criando arrays de valores para os parâmetros relevantes
tire_Sa_values = np.linspace(-np.pi / 20, np.pi / 20, 1000)  # Valores de ângulo de deslizamento lateral
tire_Ls_values = np.linspace(-0.4, 0.4, 1000)  # Valores de slip longitudinal

# Calculando as variáveis retornadas pela função Tire() para cada conjunto de parâmetros e para cada valor de tire_Fz
tire_lateral_force_values_list = [[dynamics.Tire()[0] for dynamics.tire_Sa in tire_Sa_values] for dynamics in
                                  dynamics_list]
tire_longitudinal_force_values_list = [[dynamics.Tire()[1] for dynamics.tire_Ls in tire_Ls_values] for dynamics in
                                       dynamics_list]
tire_auto_align_moment_values_list = [[dynamics.Tire()[2] for dynamics.tire_Sa in tire_Sa_values] for dynamics in
                                      dynamics_list]

# Dados experimentais
angles = [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

forces_1 = [-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890,
            2000, 2090]
torque_1 = [-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09,
            50.90, 50.10, 50.81, 48.12, 48.83]

forces_2 = [-2930, -2750, -2530, -2310, -2050, -1760, -1410, -990, -480, 90, 620, 1080, 1450, 1760, 2010, 2230, 2420,
            2590, 2750]
torque_2 = [-48.55, -49.34, -47.81, -45.55, -46.18, -43.39, -40.66, -28.43, -9.02, 17.80, 42.97, 59.03, 64.82, 67.67,
            68.95, 67.69, 67.29, 64.64, 64.16]

tire_Sa = angles  # Ângulos de deslizamento lateral
tire_Ls = (angular_velocity * params['Rdp'] / initial_speed) - 1  # Slip ratio baseado na nova velocidade angular
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

plt.plot(time_intervals, angular_velocities[:-1])
plt.xlabel('Tempo (s)')
plt.ylabel('Velocidade Angular (rad/s)')
plt.title('Velocidade Angular ao Longo do Tempo')

plt.tight_layout()
plt.show()
