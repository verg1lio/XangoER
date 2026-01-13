import numpy as np
import matplotlib.pyplot as plt

initial_speed = 20  # m/s (Velocidade inicial do veículo)
if initial_speed != 0:
    class BrakeSystem:
        def __init__(self, params=None):
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
            # print("Área do pistão no cilindro mestre:",Acm)
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
        'm_wheel': 10,  # Massa da roda [kg] <-- estimativa
        'm_tire': 10,  # Massa do pneu [kg] <-- estimativa
        'L': 1.5,  # Distância entre eixos [m]]
        'c_rr': 0.015  # Coeficiente de resistência ao rolamento
    }

    # Criação da instância do sistema de freios e aplicação da força no pedal de freio
    BrakeSystem = BrakeSystem(params)
    pedal_force = 885 # N

    # Faixa de forças no pedal
    pedal_forces = np.linspace(0, pedal_force, 1000)

    tempo = 2  # s (Tempo total de acionamento do pedal)
    time_intervals = np.linspace(0, tempo, 100)

    slip_ratios = []
    forces = []

    for pedal_force in pedal_forces:
        resultados, forca_f, inertia_wheel, torque_ajustado, angular_deceleration, angular_velocity, angular_velocities, torque_resistencia_rolamento, desaceleracao_linear, torque_disco_freio, resistencia_rolamento = BrakeSystem.apply_brake(
            pedal_force, initial_speed)

        # Slip ratio
        slip_ratio = (angular_velocity * params['Rdp'] / initial_speed) - 1
        slip_ratios.append(slip_ratio)
        forces.append(pedal_force)

    # Exibindo os resultados calculados
    print("Resultados Calculados:")
    for i, result in enumerate(resultados):
        print(f"Resultado {i + 1}: {result}")
    print("Força de frenagem:", forca_f, "N")
    print("Aceleração tangencial:",
          BrakeSystem.apply_brake(pedal_force, initial_speed)[4] * (params['Rdr'] + params['Rdp']),
          "m/s*2")
    print("Inércia da roda:", BrakeSystem.apply_brake(pedal_force, initial_speed)[2], "kg.m^2")
    print("Torque do disco de freio:", BrakeSystem.apply_brake(pedal_force, initial_speed)[9], "N.m")
    print("aceleração angular:", BrakeSystem.apply_brake(pedal_force, initial_speed)[4], "rad/s^2")
    print("velocidade tangente:",
          BrakeSystem.apply_brake(pedal_force, initial_speed)[5] * (params['Rdr'] + params['Rdp']), "m/s")
    print("Velocidade angular final:", angular_velocity, "rad/s")
    print("Resistência ao rolamento:", resistencia_rolamento, "N")
    print("Torque de resistência ao rolamento:", torque_resistencia_rolamento, "N.m")
    print("Torque ajustado:", torque_ajustado, "N.m")


    class Dynamics:

        def __init__(self, spring_type=None, spring_k=None, spring_F=None, spring_non_lin_coef=None, tire_Fz=None,
                     tire_Sa=None, tire_Ls=None, tire_friction_coef=None, damper_type=None, damper_V=None,
                     damper_F_viscous=None, damper_K_friction=None, damper_F_static=None):
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
            self.tire_friction_coef = tire_friction_coef  # coeficiente de fricção entre o pneu e a pista
            # Modelo de amortecedor
            self.damper_type = damper_type  # Coulumb, Integrated, Stribeck
            self.damper_V = damper_V  # velocidade relativa amortecedor [m/s]
            self.damper_F_viscous = damper_F_viscous  # força viscosa do fluído [N]
            self.damper_F_static = damper_F_static  # coeficiente de fricção estática de coulumb [N]
            self.damper_K_friction = damper_K_friction  # rigidez de fricção [N/m]

        def Tire(self, params):
            E, Cy, Cx, Cz, c1, c2 = params
            Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))
            D = self.tire_friction_coef * self.tire_Fz
            Bz = Cs / (Cz * D)
            Bx = Cs / (Cx * D)
            By = Cs / (Cy * D)
            tire_lateral_force = D * np.sin(
                Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
            tire_longitudinal_force = D * np.sin(
                Cx * np.arctan(9 * Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
            tire_auto_align_moment = D * np.sin(
                Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))

            return tire_lateral_force, (10 + (tire_auto_align_moment / 55)), tire_longitudinal_force


    # Valores a serem usados na instância da classe Dynamics
    slip_ratio = np.linspace(-1, 1, 1000)
    slip_angles = np.linspace(-9, 9, 1000)

    # Dados experimentais
    ratio = np.linspace(-1, 1, 19)
    angles = np.array(
        [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
    tire_lateral_forces_1 = np.array(
        [-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890, 2000,
         2090])
    tire_auto_align_moment_1 = np.array(
        [-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09,
         50.90, 50.10, 50.81, 48.12, 48.83])

    # Instanciando a classe Dynamics
    dynamics_instance = Dynamics(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45)

    # Parâmetros de Pacejka
    result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057),
              (366.4936818126405)]
    tire_Sa = angles  # Ângulos de deslizamento lateral
    tire_Ls = (angular_velocity * params['Rdp'] / initial_speed) - 1  # Slip ratio baseado na nova velocidade angular
    print("Slip Ratio Longitudinal:", tire_Ls)
    # Imprimindo os parâmetros otimizados
    print("Parâmetros otimizados:")
    print("E:", result[0])
    print("Cy:", result[1])
    print("Cx:", result[2])
    print("Cz:", result[3])
    print("c1:", result[4])
    print("c2:", result[5])

    # Plotagem da curva otimizada com os dados experimentais
    predicted_tire_lateral_forces, predicted_tire_auto_align_moment, predicted_tire_longitudinal_forces = dynamics_instance.Tire(
        result)

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

    # Plotando a curva

    plt.figure(figsize=(10, 6))
    plt.plot(forces, slip_ratios, label="Slip Ratio vs. Força do Pedal")

    # Invertendo o eixo x
    plt.gca().invert_yaxis()
    plt.xlabel("Força no Pedal (N)")
    plt.ylabel("Slip Ratio")
    plt.title("Força no Pedal em Relação ao Slip Ratio de Frenagem")
    plt.grid(True)
    plt.legend()
    plt.show()


else:
    print("A velocidade inicial do veículo não pode ser zero")
