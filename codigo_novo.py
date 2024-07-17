import math
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
    
    def slip_ratio_1(self, velocidade_angular, velocidade_linear, raio_pneu):
        
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1
        

        return value

    def show_results(self, slip_ratio):
        print("Valores do Slip Ratio: ")

        for dado in slip_ratio:
            print(dado)
    
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
        self.Dcm = params['Dcm']  # Diâmetro do cilindro mestre em metros
        self.Dwc = params['Dwc']  # Diâmetro do cilindro da roda em metros
        self.Npast = params['Npast']  # Número de pastilhas por disco
        self.atrito_coeficiente = params['atrito_coeficiente']  # Coeficiente de atrito da pastilha
        self.red = params['red']  # Raio efetivo do disco [m]
        self.Mt = params['Mt']  # Massa total do veículo [kg]
        self.L = params['L']  # Distância entre eixos [m]
        self.c_rr = params['c_rr'] # Coeficiente de resistência ao rolamento
        self.Rdr = params['Rdr']
        self.m_tire = params['m_tire']
        self.m_wheel = params['m_wheel']

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
        lR = (1 - self.psi) * self.L  # Distância do eixo ao centro de gravidade em relação a distribuição de carga traseiro
        return BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR
        
    def apply_brake(self, pedal_force):
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
        torque_resistencia_rolamento = resistencia_rolamento * (self.Rdp)

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

        return resultados, forca_frenagem, torque_ajustado

    def calculate_angular_velocity(self, torque_ajustado, tempo, initial_speed):
        time_intervals = np.linspace(0, tempo, 100)
        inertia_wheel = 0.5 * (self.m_tire + self.m_wheel) * (self.Rdp * 2)
        angular_deceleration = (torque_ajustado / 4) / inertia_wheel
        initial_angular_velocity = initial_speed / (self.Rdp)
        time_step = time_intervals[1] - time_intervals[0]
        angular_velocity = initial_angular_velocity
        angular_velocities = []

        for i in time_intervals:
            angular_velocity -= angular_deceleration * time_step
            angular_velocities.append(angular_velocity)
            angular_velocity = initial_angular_velocity + angular_velocity


        # Retornando os resultados calculados
        return angular_deceleration, angular_velocity, angular_velocities
    
    def show_graph(self, longitudinal_force, pedal_force):
        plt.figure(figsize=(10, 6))
        plt.plot(pedal_force, longitudinal_force, label='Longitudinal Force', color='blue')
        plt.xlabel('Pedal Force (N)')
        plt.ylabel('Longitudinal Force (N)')
        plt.title('Longitudinal Force vs. Pedal Force')
        plt.legend()
        plt.grid(True)
        plt.show()

        
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
        'Rdp': 0.30,  # Raio do pneu [m]
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
        'c_rr': 0  # Coeficiente de resistência ao rolamento
}

tempo = 2
initial_speed = 20

# Criação da instância do sistema de freios e aplicação da força no pedal de freio
BrakeSystem = BrakeSystem(params)
pedal_force = np.linspace(0, 0, 1)  # N
resultados, forca_frenagem, torque_ajustado = BrakeSystem.apply_brake(pedal_force) 
desaceleracao_angular, velocidade_angular, lista_velocidade = BrakeSystem.calculate_angular_velocity(torque_ajustado, initial_speed, tempo) 
print(f'Velocidade angular {velocidade_angular}, {torque_ajustado}' )  

# Valores a serem usados na instância da classe Dynamics
slip_ratio = np.linspace(-1, 1, 1000)
slip_angles = np.linspace(-9, 9, 1000)

result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]

#Instância dynamics
dynamics_instance = Dynamics()

calculo_slip_ratio = dynamics_instance.slip_ratio_1(velocidade_angular, initial_speed, params['Rdp'])


dynamics_instance_2 = Dynamics(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=calculo_slip_ratio, tire_friction_coef=1.45)
forca_longitudinal = dynamics_instance_2.Tire(result)[2]
dynamics_instance_2.show_results(calculo_slip_ratio)



# Plotting the longitudinal force against pedal force
pedal_force = np.linspace(0, 823, 1000)

# Since forca_longitudinal is an array, we need to ensure it's the same length as pedal_force
forca_longitudinal = np.array(forca_longitudinal)

# BrakeSystem.show_graph(forca_longitudinal, pedal_force)