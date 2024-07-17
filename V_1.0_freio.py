import numpy as np
import matplotlib.pyplot as plt

'''
Classe que encapsula modelos e métodos para simular o comportamento dinâmico de um sistema veicular,
especificamente focando em molas, pneus e amortecedores.
'''
class Dynamics:

    # Função para definir todos os parâmetros necessários
    def __init__(self, spring_type=None, spring_k=None, spring_F=None, spring_non_lin_coef=None, tire_Fz=None,
               tire_Sa=None, tire_Ls=None, tire_friction_coef=None, damper_type=None, damper_V=None,
               damper_F_viscous=None, damper_K_friction=None, damper_F_static=None):
        # Inicializa o modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # rigidez da mola [N/m]
        self.spring_F = spring_F  # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # coeficiente de ganho não-linear
        # Inicializa o modelo de pneu
        self.tire_Fz = tire_Fz  # carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # slip angle do pneu [rad]
        self.tire_Ls = tire_Ls  # longitudinal slip do pneu [Admensional]
        self.tire_type = 'Default'
        self.tire_friction_coef = tire_friction_coef  # coeficiente de fricção entre o pneu e a pista
        # Inicializa o modelo de amortecedor
        self.damper_type = damper_type  # Coulumb, Integrated, Stribeck
        self.damper_V = damper_V  # velocidade relativa amortecedor [m/s]
        self.damper_F_viscous = damper_F_viscous  # força viscosa do fluído [N]
        self.damper_F_static = damper_F_static  # coeficiente de fricção estática de coulumb [N]
        self.damper_K_friction = damper_K_friction  # rigidez de fricção [N/m]

    # Função que calcula as forças no pneu (lateral, longitudinal e torque auto-alinhante)
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

    # Função que calcula a razão de escorregamento (slip ratio) com base na velocidade angular e no raio do pneu
    def slip_ratio_1(self, velocidade_angular, raio_pneu):
        velocidade_linear = initial_speed
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1

        return value

    # Função que exibe os valores do slip ratio
    def show_results(self, value):
        print("Valores do Slip Ratio: ")

        for dado in value:
            print(dado)

    # Função que plota o gráfico do slip ratio em relação à força do pedal
    def plot_graph_slip_ratio(self, value=None):
        plt.figure(figsize=(10, 6))
        plt.plot(pedal_forces, value, label="Slip Ratio vs. Força do Pedal")

        # Invertendo o eixo y
        plt.gca().invert_yaxis()
        plt.xlabel("Força no Pedal (N)")
        plt.ylabel("Slip Ratio")
        plt.title("Força no Pedal em Relação ao Slip Ratio de Frenagem")
        plt.grid(True)
        plt.legend()
        plt.show()
        
    # Função que plota os gráficos das forças laterais, longitudinais e torque auto-alinhante do pneu
    def plot_graph(self, predicted_tire_lateral_forces, predicted_tire_auto_align_moment,
                   predicted_tire_longitudinal_forces, tire_lateral_experimental=None,
                   tire_auto_align_experimental=None, angles=None, ratio=None, pedal_forces=None, value=None):
        plt.figure(figsize=(20, 7))

        # Plotagem força lateral
        plt.subplot(1, 3, 1)
        plt.plot(self.tire_Sa, predicted_tire_lateral_forces, label='Curva Otimizada')
        plt.scatter(angles, tire_lateral_experimental, color='red', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Força Lateral do Pneu (N)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Plotagem torque auto-alinhante
        plt.subplot(1, 3, 2)
        plt.plot(self.tire_Sa, predicted_tire_auto_align_moment, label='Curva Otimizada')
        plt.scatter(angles, tire_auto_align_experimental, color='blue', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Torque auto-alinhante (N.m)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Plotagem força longitudinal
        plt.subplot(1, 3, 3)
        plt.plot(self.tire_Ls, predicted_tire_longitudinal_forces, label='Curva Sem Otimizar')
        plt.xlabel('Taxa de Escorregamento Longitudinal (Admensional)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal - Sem Dados Experimentais')
        plt.legend()
        plt.grid(True)

        plt.tight_layout(pad=3.0)  # Aumentando a distância entre os subplots
        plt.show()


'''
Classe que modela o sistema de frenagem de um veículo e calcula diversos parâmetros relacionados à frenagem, 
incluindo forças, torques, desaceleração, pressões necessárias e forças de atrito.
'''
class BrakeSystem:

    # Inicializa a classe com os parâmetros necessários
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
        self.c_rr = params['c_rr']  # Coeficiente de resistência ao rolamento
        self.m_tire = params['m_tire']  # Massa do pneu
        self.m_wheel = params['m_wheel']  # Massa da roda

    # Calcula diversos parâmetros relacionados à frenagem com base na força aplicada ao pedal
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
        lF = self.psi * self.L  # Distância do eixo ao centro de gravidade em relação à distribuição de carga dianteira
        lR = (1 - self.psi) * self.L  # Distância do eixo ao centro de gravidade em relação à distribuição de carga traseira
        return BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR

    # Aplica o freio e calcula os resultados com base na força aplicada ao pedal
    def apply_brake(self, pedal_force):
        resultados = self.calculate_params(pedal_force)
        BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR = resultados

        # Calcula a pressão do cilindro mestre
        pressao_cilindro_mestre = pedal_force / Acm

        # Calcula a pressurização do fluido
        pressao_fluido = pressao_cilindro_mestre

        # Calcula a transmissão de pressão da linha de freio
        transmissao_pressao = pressao_fluido * Awc

        # Calcula a força de aperto da pinça
        forca_aperto_pinca = transmissao_pressao

        # Calcula a força de atrito da pastilha de freio
        forca_atrito_pastilha = forca_aperto_pinca * self.atrito_coeficiente

        # Calcula o torque do disco de freio
        torque_disco_freio = forca_atrito_pastilha * self.red

        # Calcula a resistência ao rolamento
        resistencia_rolamento = self.c_rr * W

        # Calcula o torque gerado pela resistência ao rolamento
        torque_resistencia_rolamento = resistencia_rolamento * self.Rdp

        # Calcula o torque de freio ajustado considerando a resistência ao rolamento
        torque_ajustado = torque_disco_freio - torque_resistencia_rolamento

        # Calcula a força gerada pelo disco de freio ajustado
        forca_f = torque_ajustado / self.Rdp

        # Calcula a força de frenagem considerando todos os fatores
        forca_frenagem = (FnF + FnR) / self.Npast
        forca_frenagem *= self.atrito_coeficiente
        forca_frenagem *= self.red
        forca_frenagem /= self.Rdp
        forca_frenagem -= resistencia_rolamento

        return resultados, forca_frenagem, torque_ajustado, forca_f, torque_disco_freio, resistencia_rolamento, torque_resistencia_rolamento

    # Calcula a velocidade angular das rodas durante a frenagem
    def calculate_angular_velocity(self, torque_ajustado):

        # Cria uma série de intervalos de tempo igualmente espaçados de 0 até 'tempo' com 100 pontos
        time_intervals = np.linspace(0, tempo, 100)
        
        # Calcula a inércia da roda considerando a massa do pneu e da roda e o raio do pneu
        inertia_wheel = 0.5 * (self.m_tire + self.m_wheel) * (self.Rdp * 2)
        
        # Calcula a desaceleração angular com base no torque ajustado e na inércia da roda
        angular_desaceleration = (torque_ajustado / 4) / inertia_wheel
        
        # Calcula a velocidade angular inicial das rodas com base na velocidade inicial e no raio do pneu
        initial_angular_velocity = initial_speed / self.Rdp
        
        # Calcula o intervalo de tempo entre cada ponto de 'time_intervals'
        time_step = time_intervals[1] - time_intervals[0]
        
        # Define a velocidade angular inicial para o cálculo
        angular_velocity = initial_angular_velocity
        
        # Inicializa uma lista para armazenar as velocidades angulares ao longo do tempo
        angular_velocities = []

        # Itera sobre cada intervalo de tempo
        for i in time_intervals:
            # Atualiza a velocidade angular subtraindo a desaceleração angular multiplicada pelo intervalo de tempo
            angular_velocity -= angular_desaceleration * time_step
            # Adiciona a velocidade angular atual à lista de velocidades angulares
            angular_velocities.append(angular_velocity)

        return angular_desaceleration, angular_velocity, angular_velocities, inertia_wheel


    # Imprime os resultados calculados
    def print_resultados(self):
        print("Resultados Calculados:")
        print("Força de frenagem teórica:", BrakeSystem.apply_brake(pedal_force=823)[3], 'N')
        print("Inércia da roda:", inertia_wheel, "kg.m^2")
        print("Resistência ao rolamento:", resistencia_rolamento, "N")
        print("Torque do disco de freio:", BrakeSystem.apply_brake(pedal_force=823)[4], "N.m")
        print("Torque de resistência ao rolamento:", torque_resistencia_rolamento, "N.m")
        print("Torque ajustado:", BrakeSystem.apply_brake(pedal_force=823)[2], "N.m")

    # Plota o gráfico da força longitudinal em relação à força aplicada ao pedal
    def show_graph(self, longitudinal_force, pedal_forces):
        plt.figure(figsize=(10, 6))
        plt.plot(pedal_forces, longitudinal_force, label='Longitudinal Force', color='blue')
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

# Define a velocidade inicial e o tempo para o cálculo
initial_speed = 20  # Velocidade inicial do veículo [km/h]
tempo = 2  # Tempo para o cálculo [s]

# Parâmetros de Pacejka
result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]

# Criação de arrays com valores para plotagem e análise
slip_ratio = np.linspace(-1, 1, 1000)  # Razão de escorregamento variando de -1 a 1
slip_angles = np.linspace(-9, 9, 1000)  # Ângulos de deslizamento variando de -9 a 9 graus
pedal_forces = np.linspace(0, 823, 1000)  # Forças do pedal variando de 0 a 823 N

# Dados experimentais dos ângulos de deslizamento e forças do pneu
angles = np.array([-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
tire_lateral_forces_1 = np.array([-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890, 2000, 2090])
tire_auto_align_moment_1 = np.array([-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09, 50.90, 50.10, 50.81, 48.12, 48.83])

