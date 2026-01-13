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

            # Calculando o slip ratio longitudinal
            slip_ratio_longitudinal = 1 - (angular_velocity * self.Rdp / initial_speed)

            class Dynamics:
                def __init__(self, spring_type=None, spring_k=None, spring_F=None, spring_non_lin_coef=None, tire_Fz=None, tire_Sa=None, tire_Ls=None):
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

            # Calculando a força longitudinal gerada pelo pneu utilizando Pacejka
            dynamics = Dynamics(tire_Fz=(FzF_dyn + FzR_dyn) / 2, tire_Sa=0, tire_Ls=slip_ratio_longitudinal)
            tire_lateral_force, tire_longitudinal_force, tire_auto_align_moment = dynamics.Tire()

            # Retornando os resultados calculados
            return resultados, forca_f, inertia_wheel, torque_ajustado, angular_deceleration, angular_velocity, angular_velocities, torque_resistencia_rolamento, desaceleracao_linear, torque_disco_freio, resistencia_rolamento, tire_longitudinal_force


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

    # Criação da instância do sistema de freios e aplicação da força do pedal
    brake_system = BrakeSystem(params)
    pedal_forces = np.linspace(0, 885, 100)
    time_intervals = np.linspace(0, 2, 100)

    # Listas para armazenar os resultados
    tire_longitudinal_forces = []
    forcas_frenagem = []

    for pedal_force in pedal_forces:
        resultados, forca_f, inertia_wheel, torque_ajustado, angular_deceleration, angular_velocity, angular_velocities, torque_resistencia_rolamento, desaceleracao_linear, torque_disco_freio, resistencia_rolamento, tire_longitudinal_force = brake_system.apply_brake(pedal_force, initial_speed)
        tire_longitudinal_forces.append(tire_longitudinal_force)
        forcas_frenagem.append(forca_f)

    # Plotar o gráfico da força do pedal em relação à força longitudinal do pneu
    plt.plot(pedal_forces, tire_longitudinal_forces, label='Força Longitudinal do Pneu')
    plt.plot(pedal_forces, forcas_frenagem, label='Força de Frenagem')
    plt.xlabel('Força do Pedal (N)')
    plt.ylabel('Força Longitudinal (N)')
    plt.title('Relação entre a Força do Pedal e a Força Longitudinal')
    plt.legend()
    plt.grid(True)
    plt.show()
else:
    print("A velocidade inicial do veículo não pode ser zero.")
