import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.optimize import minimize

class Drivetrain:
    
    #Construtores com as variáveis mais recorrentes como entrada do código
    def __init__(self, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, matriz_dados):
        self.cgx = cgx
        self.cgy = cgy
        self. massa = massa
        self.entre_eixos =entre_eixos
        self.coeficiente_atrito = coeficiente_atrito
        self.raio_pneu = raio_pneu
        self.aceleracao_ideal = aceleracao_ideal
        self.reducao_primaria = reducao_primaria
        self. reducao_unica = reducao_unica
        self.matriz_dados = matriz_dados
        self.cp = 1
        self.objective_function = 0

    def CalculateOutputs(self):
        peso = self.massa * 9.81
        reacao_traseira = (peso * (self.cgx * 0.001)) / (self.entre_eixos * 0.001)
        reacao_dianteira = self.massa * 9.81 - reacao_traseira
        forca_trativa = (reacao_traseira * self.coeficiente_atrito) / (1 - ((self.cgy * 0.001) * self.coeficiente_atrito) / (self.entre_eixos * 0.001))
        transferencia_longitudinal = (forca_trativa * self.cgy * 0.001) / (self.entre_eixos * 0.001)
        carga_traseira = reacao_traseira + transferencia_longitudinal
        pico_torque_traseiro = carga_traseira * self.raio_pneu * 0.001 * self.coeficiente_atrito
        carga_pneu = carga_traseira / 2
        torque_pneu = forca_trativa * self.raio_pneu * 0.001
        reducao_final = self.reducao_primaria * self.reducao_unica
        torque_necessario_motor = pico_torque_traseiro / (self.reducao_unica * self.reducao_primaria * self.cp)
        aceleracao_primaria_real = (forca_trativa / self.massa) / 9.81
        aceleracao_primaria_ideal = self.aceleracao_ideal * 9.81
        aceleraco_real_final = aceleracao_primaria_real * 9.81
        forca_trativa_ideal = self.massa * aceleracao_primaria_ideal
        torque_pneu_ideal = forca_trativa_ideal * self.raio_pneu * 0.001
        torque_motor_ideal = torque_pneu_ideal / (self.reducao_unica * self.reducao_primaria * self.cp)
        transferencia_carga_ideal = (forca_trativa_ideal * self.cgy) / self.entre_eixos
        transferencia_carga_real = (forca_trativa * self.cgy) / self.entre_eixos
        carga_traseira_ideal = reacao_traseira + transferencia_carga_ideal
        return peso, reacao_traseira, reacao_dianteira, forca_trativa, transferencia_longitudinal, carga_traseira, pico_torque_traseiro, carga_pneu, torque_pneu, reducao_final, torque_necessario_motor, aceleracao_primaria_real, aceleracao_primaria_ideal, aceleraco_real_final, forca_trativa_ideal, torque_pneu_ideal, torque_motor_ideal, transferencia_carga_ideal, transferencia_carga_real, carga_traseira_ideal

    def show_results(self):
        peso, reacao_traseira, reacao_dianteira, forca_trativa, transferencia_longitudinal, carga_traseira, pico_torque_traseiro, carga_pneu, torque_pneu, reducao_final, torque_necessario_motor, aceleracao_primaria_real, aceleracao_primaria_ideal, aceleraco_real_final, forca_trativa_ideal, torque_pneu_ideal, torque_motor_ideal, transferencia_carga_ideal, transferencia_carga_real, carga_traseira_ideal = Drivetrain.CalculateOutputs(self)
        rpm_values, torque_values, power_values = Drivetrain.CurveTorquePower(self)
        
        # Print dos resultados obtidos
        print(f'''Resultados:
            peso: {peso}N
            Reação no eixo traseiro: {reacao_traseira}N
            Reação no eixo dianteiro: {reacao_dianteira}N
            Força trativa: {forca_trativa}N
            Transferência de carga longitudinal: {transferencia_longitudinal}N
            Carga no eixo traseiro: {carga_traseira}N
            Pico de torque no eixo traseiro: {pico_torque_traseiro}Nm
            Carga no pneu: {carga_pneu}N
            Torque no pneu: {torque_pneu}Nm
            Redução final: {reducao_final}
            Torque necessário no motor: {torque_necessario_motor}Nm
            Aceleração primária real (g): {aceleracao_primaria_real}
            Aceleração final ideal: {aceleracao_primaria_ideal}m/s²
            Aceleração final real: {aceleraco_real_final}m/s²
            Força trativa ideal: {forca_trativa_ideal}N
            Torque no pneu ideal: {torque_pneu_ideal}Nm
            Torque no motor ideal: {torque_motor_ideal}Nm
            Transferência de carga ideal: {transferencia_carga_ideal}N
            Transferência de carga real: {transferencia_carga_real}N
            Carga total no eixo traseiro ideal: {carga_traseira_ideal}N\n
            ''')
        
        print("\nMatriz de RPM, Torque e Potência:")
        print("RPM\t\tTorque (Nm)\tPotência (kW)")
        
        for data in self.matriz_dados:
            rpm = data["rpm"]
            torque = data["trq"]
            potencia = data["ptc"]
            print("{:.2f}\t\t{:.2f}\t\t{:.2f}".format(rpm, torque, potencia))

        # Plotando o gráfico
        plt.figure(figsize=(10, 6))
        plt.plot(rpm_values, torque_values, label='Torque [Nm]', color='blue')
        plt.plot(rpm_values, power_values, label='Power [kW]', color='orange')
        plt.title("Curva de Torque e Potência")
        plt.xlabel('RPM')
        plt.ylabel('Torque [Nm] / Power [kW]')
        plt.legend()
        plt.grid(True)
        plt.show()

    # Definindo a função objetivo para otimização de cp
    def objective_function(cp, self): #Por conta da funcão de otimização, o self é passado como segundo argumento
        
        self.cp = cp
        
        # Calcula os outputs
        outputs = Drivetrain.CalculateOutputs(self)
        torque_necessario_motor = outputs[10]  # Torque necessário no motor
        target_torque_necessario_motor = 80.0  # Alvo para o torque necessário no motor
        
        # Calcular o desempenho do carro com o cp atual
        performance = Drivetrain.CarPerformance(self)
        
        # Extrair a velocidade linear máxima e a menor força final
        max_velocidade_linear = max([p['velocidade_linear'] for p in performance])
        min_forca_final = min([p['forca_final'] for p in performance])
        
        # Penalidades para as novas condições
        penalty = 0
        if max_velocidade_linear < 100:
            penalty += (100 - max_velocidade_linear) ** 2
        if min_forca_final <= 0:
            penalty += abs(min_forca_final) * 1000
        
        
        
        # A função objetivo será a diferença absoluta entre torque_necessario_motor e o valor alvo com penalidade
        return abs(torque_necessario_motor - target_torque_necessario_motor) + penalty

    # Função que realiza a otimização do parâmetro cp
    def optimize_cp(self):
        
        initial_guess = 1.5  # Chute inicial para cp
        result = minimize(Drivetrain.objective_function, initial_guess, args = (self), method='BFGS')
        
        self.cp = result.x[0]
        self.objective_function = result.fun
        
    # Função de desempenho do carro
    def CarPerformance(self):
        peso = self.massa * 9.81
        rendimento_transmissao = 0.9
        transmissao_motor_roda = 0.9
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        basic_f = 0.015
        speed_f = 0.012

        parametros = []

        for dado in self.matriz_dados:
            # Cálculo da força trativa (N)
            forca_trativa = ((dado["trq"] * self.reducao_primaria * self.reducao_unica * self.cp) / (self.raio_pneu * 0.001)) * rendimento_transmissao
        
            # Cálculo da velocidade angular (rad/s)
            velocidade_angular = (dado["rpm"] * 2 * math.pi) / (60 * self.reducao_primaria * self.reducao_unica * self.cp)
            
            # Cálculo da velocidade linear (km/h)
            velocidade_linear = ((velocidade_angular * (self.raio_pneu * 0.001)) * transmissao_motor_roda) * 3.6

            # Cálculo da força de arrasto (N)
            fa = (densidade_ar * velocidade_linear ** 2 * coeficiente_arrasto * area_frontal) / 2

            # Cálculo da resistência de rolamento (N)
            rr = (basic_f + (3.24 * speed_f * ((velocidade_linear / 100 * 0.44704) ** 2.5))) * peso

            # Cálculo da força final (N)
            forca_final = forca_trativa - fa - rr

            # Armazenar os parâmetros calculados em um dicionário
            parametro = {"forca_trativa": forca_trativa,"va": velocidade_angular,"velocidade_linear": velocidade_linear,"fa": fa,"rr": rr,"forca_final": forca_final}

            parametros.append(parametro)

        # Retornar os parâmetros calculados
        return parametros 
    
    def print_car_performance(self):
        performance = Drivetrain.CarPerformance(self)
        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]")
        
        for param in performance:
            print(f"{param['forca_trativa']}\t{param['va']}\t{param['velocidade_linear']}")
        
        print("\nForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for param in performance:
            print(f"{param['fa']}\t{param['rr']}\t{param['forca_final']}")

    def HalfShaftsSizing(self, fsi=1.25, tet=786, tec=471.6, dif=1):
        # Obtendo o maior torque do motor a partir dos dados experimentais 
        torque_max_motor = max(data["trq"] for data in self.matriz_dados)
        
        # Calculando o torque máximo nos semieixos
        torque_max_semieixo = torque_max_motor * self.reducao_primaria * self.reducao_unica * self.cp * dif
        
        # Calculando o torque máximo de projeto
        torque_max_projeto = torque_max_semieixo * fsi
        
        # Calculando o diâmetro dos semieixos (mm)
        diametro_semieixo = (((2 * torque_max_projeto) / (math.pi * tec * 10 ** 6)) ** (1 / 3)) * 2000
        
        # Calculando o fator de segurança obtido
        fator_seguranca_obtido = (math.pi * (((diametro_semieixo / 1000) / 2) ** 3) * tec * (10 ** 6)) / (2 * torque_max_semieixo)
        
        # Calculando o fator de segurança para 1 polegada
        fs1p = (math.pi * ((0.0254 / 2) ** 3) * tec * (10 ** 6)) / (2 * torque_max_semieixo)
        
        # Print dos resultados obtidos
        print("Dimensionamento dos Semieixos:")
        print("Torque máximo do motor:", torque_max_motor, "Nm")
        print("Torque máximo nos semieixos:", torque_max_semieixo, "Nm")
        print("Torque máximo de projeto:", torque_max_projeto, "Nm")
        print("Diâmetro dos semieixos:", diametro_semieixo, "mm")
        print("Fator de segurança ideal:", fsi)
        print("Fator de segurança obtido:", fator_seguranca_obtido)
        print("Fator de segurança para 1 polegada:", fs1p)
        
        return torque_max_motor, torque_max_semieixo, torque_max_projeto, diametro_semieixo, fator_seguranca_obtido, fs1p

    def CurveTorquePower(self):
        matriz_dados = self.matriz_dados
        rpm_values = [data["rpm"] for data in matriz_dados]
        torque_values = [data["trq"] for data in matriz_dados]
        power_values = [data["ptc"] for data in matriz_dados]
        return rpm_values, torque_values, power_values

    @classmethod
    def generate_model(cls, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, matriz_dados):
        model = cls(cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, matriz_dados)
        model.optimize_cp()
        return model
    
    def __str__(self):
        return f'Relação Coroa-Pinhão: {self.cp}\nFunção Objetivo: {self.objective_function}'

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

# Definindo uma classe chamada Tire
class Tire:
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=0, B0=0, B1=None, 
                B2=0, B3=None, omega=315, slip_angle_start=-10, slip_angle_end=10, angle_step=0.5, WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
        # Inicializando parâmetros para o modelo de pneu
        self.tire_Fz = tire_Fz  # Carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # Ângulo de deslizamento lateral do pneu [rad]
        self.tire_Ls = tire_Ls  # Escorregamento longitudinal do pneu [Adimensional]
        self.tire_type = 'Default'  # Tipo de pneu
        self.tire_friction_coef = tire_friction_coef  # Coeficiente de fricção entre pneu e pista
        self.tire_Ca = tire_Ca  # Ângulo de camber do pneu

        # Comprimentos das barras do mecanismo de quatro barras
        self.L0 = B0  # Comprimento da barra de direção
        self.L1 = B1  # Comprimento do braço de direção
        self.L2 = B2  # Comprimento da bitola
        self.L3 = B3  # Comprimento do braço de direção

        # Ângulo de orientação da barra longitudinal em relação ao eixo horizontal
        self.alpha = np.radians(omega)

        # Array de ângulos de deslizamento lateral do pneu em graus
        self.theta2 = np.arange(slip_angle_start + 90, slip_angle_end + 91, angle_step)

        # Convertendo os ângulos de deslizamento lateral do pneu para radianos
        self.theta2_rad = np.radians(self.theta2)

        # Inicialização das listas de resultados para armazenar dados ao longo do cálculo
        self.AC = []  # Lista para armazenar AC
        self.beta = []  # Lista para armazenar beta
        self.psi = []  # Lista para armazenar psi
        self.lamda = []  # Lista para armazenar lambda
        self.theta3 = []  # Lista para armazenar theta3
        self.theta4 = []  # Lista para armazenar theta4
        self.Ox, self.Oy = [0] * len(self.theta2), [0] * len(self.theta2)  # Lista para armazenar coordenadas de O
        self.Ax, self.Ay = [], []  # Listas para armazenar coordenadas de A
        self.Bx, self.By = [], []  # Listas para armazenar coordenadas de B
        self.Cx, self.Cy = [B0] * len(self.theta2), [0] * len(self.theta2)  # Listas para armazenar coordenadas de C
        self.w = []  # Lista para armazenar w
        self.om2, self.om4 = [], []  # Listas para armazenar om2 e om4
        self.alpha_dot = []  # Lista para armazenar alpha_dot
        self.outer_slip = []  # Lista para armazenar ângulos de deslizamento lateral externo
        self.inner_slip = []  # Lista para armazenar ângulos de deslizamento lateral interno
        self.static_slip_angle = None  # Variável para armazenar o ângulo de deslizamento estático (inicialmente não definido)

        # Coordenadas da barra longitudinal (entre-eixos)
        self.WB = WB  # Entre-eixos (wheelbase) fixo
        self.rear_axle_length = rear_axle_length  # Comprimento do eixo traseiro fixo
        self.rear_axle_center = (B0 / 2, 0)

        # Entradas de rigidez do pneu
        self.track_y = track_y
        self.tire_k = tire_k


    def calculate_kinematics(self):
        for i in range(len(self.theta2)):
            # Cálculos intermediários
            AC_i = np.sqrt(self.L0**2 + self.L1**2 - 2 * self.L0 * self.L1 * np.cos(self.theta2_rad[i]))
            beta_i = np.arccos((self.L0**2 + AC_i**2 - self.L1**2) / (2 * self.L0 * AC_i))
            psi_i = np.arccos((self.L2**2 + AC_i**2 - self.L3**2) / (2 * self.L2 * AC_i))
            lamda_i = np.arccos((self.L3**2 + AC_i**2 - self.L2**2) / (2 * self.L3 * AC_i))

            theta3_i = psi_i - beta_i
            theta4_i = np.pi - lamda_i - beta_i

            if self.theta2[i] > 180:
                theta3_i = psi_i + beta_i
                theta4_i = np.pi - lamda_i + beta_i

            # Armazenamento dos resultados
            self.AC.append(AC_i)
            self.beta.append(beta_i)
            self.psi.append(psi_i)
            self.lamda.append(lamda_i)
            self.theta3.append(theta3_i)
            self.theta4.append(theta4_i)

            # Definição das posições das juntas
            Ax_i = self.L1 * np.cos(self.theta2_rad[i])
            Ay_i = self.L1 * np.sin(self.theta2_rad[i])
            Bx_i = Ax_i + self.L2 * np.cos(theta3_i)
            By_i = Ay_i + self.L2 * np.sin(theta3_i)

            self.Ax.append(Ax_i)
            self.Ay.append(Ay_i)
            self.Bx.append(Bx_i)
            self.By.append(By_i)

            r = np.array([
                [self.L2 * np.cos(theta3_i), -self.L3 * np.cos(theta4_i)],
                [self.L2 * np.sin(theta3_i), -self.L3 * np.sin(theta4_i)]
            ])
            v = np.array([-self.L1 * self.alpha * np.cos(self.theta2_rad[i]), -self.L1 * self.alpha * np.sin(self.theta2_rad[i])])

            w_i = np.linalg.solve(r, v)
            self.w.append(w_i)
            self.om2.append(w_i[0])
            self.om4.append(w_i[1])


            # Cálculo dos ângulos outer_slip e inner_slip
            outer_slip_i = -(np.pi / 2 - theta4_i)
            inner_slip_i = -(np.pi / 2 - self.theta2_rad[i])

            self.outer_slip.append(np.degrees(outer_slip_i))
            self.inner_slip.append(np.degrees(inner_slip_i))

            # Verificação do caso estático
            if np.isclose(np.abs(outer_slip_i), np.abs(inner_slip_i), atol=1e-2):
                self.static_slip_angle = self.inner_slip[i]

    def Tire_forces(self, params):
        # Desembalando os parâmetros de Pacejka
        E, Cy, Cx, Cz, c1, c2 = params

        # Calculando parâmetros intermediários
        Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))
        D = self.tire_friction_coef * self.tire_Fz
        Bz = Cs / (Cz * D)
        Bx = Cs / (Cx * D)
        By = Cs / (Cy * D)

        # Calculando a força lateral do pneu
        tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))

        # Calculando a força longitudinal do pneu
        tire_longitudinal_force = D * np.sin(Cx * np.arctan(9 * Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))

        # Calculando o momento de auto-alinhamento do pneu
        tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))

        # Calculando a força de camber
        camber_thrust = D * np.sin(Cy * np.arctan(By * self.tire_Ca))

        # Retornando as forças calculadas e o momento de auto-alinhamento
        return tire_lateral_force + 0.5 * camber_thrust, (10 + (tire_auto_align_moment / 55)), tire_longitudinal_force
    
    def calcular_forca(self):
        """Calcula a força com base na deformção do pneu"""
        force = self.track_y * self.tire_k
        return force
    
    # Função que calcula a razão de escorregamento (slip ratio) com base na velocidade angular e no raio do pneu
    def slip_ratio_1(self, velocidade_angular, raio_pneu):
        velocidade_linear = initial_speed
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1

        return value
    
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

    # Função que mostra valores de slip ratio em relação ao rpm
    def show_slip_ratio(self, rpm_values, slip_ratio, velocidade_angular):
        
        print("Valores do Slip Ratio: ")
        for dado in slip_ratio:
            print(dado)
            
        plt.figure(figsize=(15, 5))
        
        plt.subplot(1, 2, 1)
        plt.plot(rpm_values, slip_ratio, label = 'Slip Ratio', color = 'blue')
        plt.xlabel('RPM')
        plt.ylabel('Slip Ratio')
        plt.title('Slip Ratio x RPM') 
        plt.legend()
        
        plt.subplot(1, 2, 2)
        plt.plot(rpm_values, velocidade_angular, label = 'Velocidade Angular', color = 'red')
        plt.xlabel('RPM')
        plt.ylabel('Velocidade Angular (rad/s)')
        plt.title('Velocidade Angular x RPM')
        plt.legend()
        plt.tight_layout()
        plt.show()

    # Função que plota valores de força longitudinal em função do rpm
    def show_longitudinal_rpm(rpm_values, tire_longitudinal_force):
        plt.figure(figsize=(10, 6))
        plt.plot(rpm_values, tire_longitudinal_force, label='Longitudinal Force', color='blue')
        plt.xlabel('RPM')
        plt.ylabel('Longitudinal Force (N)')
        plt.title('Longitudinal Force vs. RPM')
        plt.legend()
        plt.grid(True)
        plt.show()
    
    def plotar_deformacao(self, track_y_values):
        """Plota o gráfico com base em uma lista de valores fornecidos que representam a variação da pista"""
        force_values = []

        for y in track_y_values:
            self.track_y = y
            force_values.append(self.calcular_forca())

        # Criando o gráfico
        plt.figure(figsize=(8, 6))
        plt.plot(track_y_values, force_values, color='b')
        plt.title('Carregamento Aplicado vs Deformação')
        plt.xlabel('Deformação (mm)')
        plt.ylabel('Carregamento (N)')
        plt.grid(True)
        plt.show()
    
    def plot_camber(self, predicted_tire_lateral_forces, predicted_tire_lateral_forces_1, predicted_tire_lateral_forces_2, tire_lateral_experimental=None, tire_lateral_experimental_1=None, tire_lateral_experimental_2=None, angles=None, ratio=None):
        # Configuração da figura e tamanho
        plt.figure(figsize=(20, 7))

        # Subplot para força lateral do pneu
        plt.subplot(1, 3, 1)
        plt.plot(self.tire_Sa, predicted_tire_lateral_forces, label='Curva com Camber')
        plt.scatter(angles, tire_lateral_experimental, color='red', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Força Lateral do Pneu (N)')
        plt.title('Curva com Camber e Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Subplot para torque auto-alinhante
        plt.subplot(1, 3, 2)
        plt.plot(self.tire_Sa, predicted_tire_lateral_forces_1, label='Curva com Camber')
        plt.scatter(angles, tire_lateral_experimental_1, color='red', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Força Lateral do Pneu (N)')
        plt.title('Curva com Camber e Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Subplot para força longitudinal do pneu
        plt.subplot(1, 3, 3)
        plt.plot(self.tire_Sa, predicted_tire_lateral_forces_2, label='Curva com Camber')
        plt.scatter(angles, tire_lateral_experimental_2, color='red', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Força Lateral do Pneu (N)')
        plt.title('Curva com Camber e Dados Experimentais')
        plt.legend()
        plt.grid(True)
    
    def plot_graph(self, predicted_tire_lateral_forces, predicted_tire_auto_align_moment, predicted_tire_longitudinal_forces, tire_lateral_experimental=None, tire_auto_align_experimental=None, angles=None, ratio=None):
        # Configuração da figura e tamanho
        plt.figure(figsize=(20, 7))

        # Subplot para força lateral do pneu
        plt.subplot(1, 3, 1)
        plt.plot(self.tire_Sa, predicted_tire_lateral_forces, label='Curva Otimizada')
        plt.scatter(angles, tire_lateral_experimental, color='red', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Força Lateral do Pneu (N)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Subplot para torque auto-alinhante
        plt.subplot(1, 3, 2)
        plt.plot(self.tire_Sa, predicted_tire_auto_align_moment, label='Curva Otimizada')
        plt.scatter(angles, tire_auto_align_experimental, color='blue', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Torque auto-alinhante (N.m)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Subplot para força longitudinal do pneu
        plt.subplot(1, 3, 3)
        plt.plot(self.tire_Ls, predicted_tire_longitudinal_forces, label='Curva Sem Otimizar')
        plt.xlabel('Taxa de Escorregamento Longitudinal (Admensional)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal - Sem Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Ajusta a disposição dos subplots para evitar sobreposição
        plt.tight_layout(pad=3.0)
        plt.show()

    def plot_mechanism(self):
        

        for i in range(len(self.theta2)):
            plt.figure(1)
            plt.plot([self.Ox[i], self.Ax[i]], [self.Oy[i], self.Ay[i]], 'r', linewidth=1)
            plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'k', linewidth=2)
            plt.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i]], 'r', linewidth=1)
            plt.plot([self.Ox[i], self.Cx[i]], [self.Oy[i], self.Cy[i]], 'r', linewidth=1, linestyle = 'dotted')

            # Desenho da barra do entre-eixos
            midpoint_x = (self.Ox[i] + self.Cx[i]) / 2
            rear_x1 = midpoint_x - (self.rear_axle_length / 2)
            rear_x2 = midpoint_x + (self.rear_axle_length / 2)
            

            plt.plot([midpoint_x, midpoint_x], [self.L3, (self.L3 - self.WB)], 'g', linewidth=1)  # Barra do entre-eixos
            plt.plot([rear_x1, rear_x2], [(self.L3 - self.WB), (self.L3 - self.WB)], 'g', linewidth=1)  # Barra do eixo traseiro

            # Adicionando as coordenadas dos pontos
            plt.text(self.Ox[i], self.Oy[i], f'({self.Ox[i]:.2f}, {self.Oy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Ax[i], self.Ay[i], f'Manga de eixo ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Bx[i], self.By[i], f'Manga de eixo ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Cx[i], self.Cy[i], f'({self.Cx[i]:.2f}, {self.Cy[i]:.2f})', fontsize=8, ha='right')

            # Adicionando os ângulos de slip
            plt.text((self.Cx[i] + self.Bx[i]) / 2, (self.Cy[i] + self.By[i]) / 2, f'{self.outer_slip[i]:.2f}°', fontsize=10, ha='center')
            plt.text((self.Ox[i] + self.Ax[i]) / 2, (self.Oy[i] + self.Ay[i]) / 2, f'{self.inner_slip[i]:.2f}°', fontsize=10, ha='center')

            # Achando primeiro ponto de Ackerman

            plt.plot([(self.Ax[i] + self.Ox[i])/2, (self.Ax[i] + self.Ox[i])/2 - (self.L3/2 - self.WB)/np.sin(self.theta2_rad[i] + np.pi/2)], [(self.Ay[i] + self.Oy[i])/2, (self.L3 - self.WB)], 'g', linewidth=1, linestyle = 'dotted')  # Barra do entre-eixos
            plt.plot([(self.Bx[i] + self.Cx[i])/2, (self.Bx[i] + self.Cx[i])/2 - (self.L3/2 - self.WB)/np.sin(self.theta4[i] + np.pi/2)], [(self.By[i] + self.Cy[i])/2, (self.L3 - self.WB)], 'r', linewidth=1, linestyle = 'dotted')


            if i == 0 or i == (len(self.theta2) - 1):
                plt.show()

            plt.grid()
            plt.axis('equal')
            plt.axis([-300, 1800, -1000, 1000])
            plt.draw()
            plt.pause(0.01)
            plt.clf()

        for i in range(len(self.theta2)-1, -1, -1):
            plt.plot([self.Ox[i], self.Ax[i]], [self.Oy[i], self.Ay[i]], 'r', linewidth=1)
            plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'k', linewidth=2)
            plt.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i]], 'r', linewidth=1)
            plt.plot([self.Ox[i], self.Cx[i]], [self.Oy[i], self.Cy[i]], 'r', linewidth=1, linestyle = 'dotted')

            # Desenho da barra do entre-eixos
            midpoint_x = (self.Ox[i] + self.Cx[i]) / 2
            rear_x1 = midpoint_x - (self.rear_axle_length / 2) * np.sin(self.alpha)
            rear_x2 = midpoint_x + (self.rear_axle_length / 2) * np.sin(self.alpha)
            

            plt.plot([midpoint_x, midpoint_x], [self.L3, (self.L3 - self.WB)], 'g', linewidth=1)  # Barra do entre-eixos
            plt.plot([rear_x1, rear_x2], [(self.L3 - self.WB), (self.L3 - self.WB)], 'g', linewidth=1)  # Barra do eixo traseiro

            # Adicionando as coordenadas dos pontos
            plt.text(self.Ox[i], self.Oy[i], f'({self.Ox[i]:.2f}, {self.Oy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Ax[i], self.Ay[i], f'Manga de eixo ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Bx[i], self.By[i], f'Manga de eixo ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='center')

            # Adicionando os ângulos de slip
            plt.text((self.Cx[i] + self.Bx[i]) / 2, (self.Cy[i] + self.By[i]) / 2, f'{self.outer_slip[i]:.2f}°', fontsize=10, ha='center')
            plt.text((self.Ox[i] + self.Ax[i]) / 2, (self.Oy[i] + self.Ay[i]) / 2, f'{self.inner_slip[i]:.2f}°', fontsize=10, ha='center')

            # Achando primeiro ponto de Ackerman

            plt.plot([(self.Ax[i] + self.Ox[i])/2, (self.Ax[i] + self.Ox[i])/2 - (self.L3/2 - self.WB)/np.sin(self.theta2_rad[i] + np.pi/2)], [(self.Ay[i] + self.Oy[i])/2, (self.L3 - self.WB)], 'g', linewidth=1, linestyle = 'dotted')  # Barra do entre-eixos
            plt.plot([(self.Bx[i] + self.Cx[i])/2, (self.Bx[i] + self.Cx[i])/2 - (self.L3/2 - self.WB)/np.sin(self.theta4[i] + np.pi/2)], [(self.By[i] + self.Cy[i])/2, (self.L3 - self.WB)], 'r', linewidth=1, linestyle = 'dotted')

            plt.grid()
            plt.axis('equal')
            plt.axis([-300, 1800, -1000, 1100])
            plt.draw()
            plt.pause(0.01)
            plt.clf()

    def run(self):
        self.calculate_kinematics()
        self.plot_mechanism()

        if self.static_slip_angle is not None:
            print(f"O ângulo de toe é : {self.static_slip_angle:.2f}°")
        else:
            print("Não foi possível determinar um ângulo estático para as rodas dentro do intervalo fornecido.")

class Kinematics:
    def __init__(self, L0=None, L1=None, L2=0, L3=None, alpha=60, Py_start=100, Py_end=400, Py_step=5,
                spring_type=None, spring_k=None, spring_x=None, spring_y=None, spring_non_lin_coef=None, spring_angle = 45,
                damper_type=None, damper_V=None, damper_F_viscous=None, damper_K_friction=None, damper_F_static=None):
        # Inicialização dos comprimentos das barras e outros parâmetros do mecanismo
        self.L0 = L0  # Comprimento da barra fixa
        self.L1 = L1  # Comprimento da barra de entrada
        self.L2 = L2  # Comprimento da barra de acoplamento
        self.L3 = L3  # Comprimento da barra de saída
        self.alpha = np.radians(alpha)  # Ângulo inicial da barra de entrada em radianos
        self.L_AP = 0.5 * L2  # Metade do comprimento da barra de acoplamento (ponto médio)
        self.Py_step = Py_step  # Passo do movimento vertical de P
        self.Py = np.arange(Py_start, Py_end, Py_step)  # Intervalo de posições verticais de P
        self.num_points = int((Py_end - Py_start)/Py_step) # Número de pontos no intervalo

        # Inicialização das listas de resultados
        self.AC = []  # Lista para armazenar os comprimentos AC
        self.beta = []  # Lista para armazenar os ângulos beta
        self.psi = []  # Lista para armazenar os ângulos psi
        self.lamda = []  # Lista para armazenar os ângulos lambda
        self.theta2 = []  # Lista para armazenar os ângulos theta2 (entrada)
        self.theta3 = []  # Lista para armazenar os ângulos theta3 (acoplamento)
        self.theta4 = []  # Lista para armazenar os ângulos theta4 (saída)
        self.Ox, self.Oy = [0] * len(self.Py), [0] * len(self.Py)  # Posições das âncoras fixas no chassi
        self.Ax, self.Ay = [], []  # Listas para armazenar as posições de A
        self.Bx, self.By = [], []  # Listas para armazenar as posições de B
        self.Cy, self.Cx = [L0] * len(self.Py), [0] * len(self.Py)  # Posições das âncoras fixas no chassi
        self.Px = []  # Lista para armazenar as posições horizontais de P
        self.w = []  # Lista para armazenar as velocidades angulares
        self.om2, self.om4 = [], []  # Listas para armazenar as velocidades angulares de theta2 e theta4
        self.alpha_dot = []  # Lista para armazenar as acelerações angulares
        self.V_Px, self.V_Py = [], []  # Listas para armazenar as velocidades de P
        self.V_Ax, self.V_Ay = [], []  # Listas para armazenar as velocidades de A
        self.V_A, self.P_A = [], [] # # Listas para armazenar a magnitude das velocidades e posições de A

        # Inicialização dos parâmetros do modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # rigidez da mola [N/m]
        self.spring_x = spring_x  # força que a mola recebe [N]
        self.spring_y = spring_y  # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # coeficiente de ganho não-linear
        self.spring_angle = np.radians(spring_angle) # Ângulo da mola em relação a vertical

        # Inicialização dos parâmetros do modelo de amortecedor
        self.damper_type = damper_type # Coulumb, Integrated, Stribeck
        self.damper_V = damper_V # velocidade relativa amortecedor [m/s]
        self.damper_F_viscous = damper_F_viscous # força viscosa do fluído [N]
        self.damper_F_static = damper_F_static # coeficiente de fricção estática de coulumb [N]
        self.damper_K_friction = damper_K_friction # rigidez de fricção [N/m]

    def Spring(self):
        """Calcula a força da mola com base na deformação e no tipo de mola."""
        if self.spring_type == 'Hooke':
            spring_F = self.spring_x * self.spring_k
        if self.spring_type == 'Softening':
            spring_F = self.spring_k * (self.spring_x ** self.spring_non_lin_coef)
        return spring_F

    def Damper(self):
        """Calcula a força do amortecedor com base no tipo de amortecedor."""
        if self.damper_type == 'Coulumb':
            damper_F = self.damper_F_static * np.tanh(self.damper_K_friction * self.damper_V)
        if self.damper_type == 'Integrated':
            damper_F = self.damper_F_static * np.tanh(self.damper_K_friction * self.damper_V) + np.sign(self.damper_V) * self.damper_F_viscous * (self.damper_V ** 2)
        return damper_F
    
    def calcular_theta2(self, Py_i):
        """Calcula theta2 a partir de Py."""
        theta2_rad_i = np.arccos((Py_i - self.L1 * np.cos(self.alpha)) / self.L1)
        return theta2_rad_i

    def calcular_camber(self, Ax, Ay, Bx, By):
        """Calcula o ângulo que a barra de acoplamento faz com a vertical."""
        delta_x = Bx - Ax
        delta_y = By - Ay
        angulo = np.degrees(np.arctan2(delta_x, delta_y))
        return angulo

    def calcular_cinematica(self):
        """Calcula as posições das juntas, ângulos e outras propriedades cinemáticas do mecanismo."""
        for i in range(len(self.Py)):
            # Calcula theta2 para a posição atual de Py
            theta2_rad_i = self.calcular_theta2(self.Py[i])
            self.theta2.append(np.degrees(theta2_rad_i))

            # Cálculos intermediários
            AC_i = np.sqrt(self.L0**2 + self.L1**2 - 2 * self.L0 * self.L1 * np.cos(theta2_rad_i))
            beta_i = np.arccos((self.L0**2 + AC_i**2 - self.L1**2) / (2 * self.L0 * AC_i))
            psi_i = np.arccos((self.L2**2 + AC_i**2 - self.L3**2) / (2 * self.L2 * AC_i))
            lamda_i = np.arccos((self.L3**2 + AC_i**2 - self.L2**2) / (2 * self.L3 * AC_i))

            # Calcula os ângulos theta3 e theta4
            theta3_i = psi_i - beta_i
            theta4_i = np.pi - lamda_i - beta_i

            # Ajuste de ângulos para diferentes quadrantes
            if self.theta2[i] > 180:
                theta3_i = psi_i + beta_i
                theta4_i = np.pi - lamda_i + beta_i

            # Armazenamento dos resultados intermediários
            self.AC.append(AC_i)
            self.beta.append(beta_i)
            self.psi.append(psi_i)
            self.lamda.append(lamda_i)
            self.theta3.append(theta3_i)
            self.theta4.append(theta4_i)

            # Definição das posições das juntas A e B
            Ax_i = self.L1 * np.sin(theta2_rad_i)
            Ay_i = self.L1 * np.cos(theta2_rad_i)
            Bx_i = Ax_i + self.L2 * np.sin(theta3_i)
            By_i = Ay_i + self.L2 * np.cos(theta3_i)
            Px_i = Ax_i + self.L_AP * np.sin(self.alpha + theta3_i)
            Py_i = Ay_i + self.L_AP * np.cos(self.alpha + theta3_i)

            # Armazenamento das posições das juntas
            self.Ax.append(Ax_i)
            self.Ay.append(Ay_i)
            self.Bx.append(Bx_i)
            self.By.append(By_i)
            self.Px.append(Px_i)

            # Sistema de equações para encontrar velocidades angulares
            r = np.array([
                [self.L2 * np.cos(theta3_i), -self.L3 * np.cos(theta4_i)],
                [self.L2 * np.sin(theta3_i), -self.L3 * np.sin(theta4_i)]
            ])
            v = np.array([-self.L1 * self.alpha * np.cos(theta2_rad_i), -self.L1 * self.alpha * np.sin(theta2_rad_i)])

            # Solução do sistema para encontrar as velocidades angulares
            w_i = np.linalg.solve(r, v)
            self.w.append(w_i)
            self.om2.append(w_i[0])
            self.om4.append(w_i[1])

            # Cálculo da aceleração angular
            alpha_dot_i = - (self.L1 * w_i[0] * np.cos(theta2_rad_i) - self.L3 * w_i[1] * np.cos(theta4_i)) / \
                          (self.L2 * np.cos(self.alpha) + self.L_AP * np.cos(self.alpha))
            self.alpha_dot.append(alpha_dot_i)

            # Cálculo das velocidadees de A
            V_Ax_i = self.om2[i] * Ax_i
            V_Ay_i = self.om2[i] * Ay_i

            self.V_Ax.append(V_Ax_i)
            self.V_Ay.append(V_Ay_i)

            # Cálculo das velocidades de P
            V_Px_i = alpha_dot_i * Px_i
            V_Py_i = alpha_dot_i * Py_i

            self.V_Px.append(V_Px_i)
            self.V_Py.append(V_Py_i)

            # Cálculo da magnitude da velocidade e da posição de A
            V_A_i = np.sqrt((V_Ax_i**2) + (V_Ay_i**2))
            P_A_i = np.sqrt((Ax_i**2) + (Ay_i**2))

            self.V_A.append(V_A_i)
            self.P_A.append(P_A_i)

        # Calculo da variação de Ax em relação ao caso estático
        Ax_d = []
        for i in range(len(self.Py)):
            Ax_d_i = self.Ax[int(self.num_points/2)] - self.Ax[i]
            Ax_d.append(Ax_d_i)

        # Calculando ângulo de câmber
        angulo_camber = self.calcular_camber(self.Ax[int(self.num_points - 1)], self.Ay[int(self.num_points - 1)], self.Bx[int(self.num_points - 1)], self.By[int(self.num_points - 1)])
        return angulo_camber, Ax_d, self.V_Ax, self.Ay, self.V_Ay
    
    def plotar_cinematica(self):
        """Plota a cinemática do mecanismo."""
        for i in range(len(self.Py)):
            angulo_camber = self.calcular_camber(self.Ax[i], self.Ay[i], self.Bx[i], self.By[i])
            
            plt.figure(1)
            plt.plot([self.Ox[i], self.Ax[i]], [self.Oy[i], self.Ay[i]], 'b', linewidth=2, label='Barra de Entrada')
            plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'b', linewidth=2, label='Barra de Acoplamento', linestyle='dotted')
            plt.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i]], 'b', linewidth=2, label='Barra de Saída')
            plt.plot([self.Ax[i], self.Px[i]], [self.Ay[i], self.Py[i]], 'r', linewidth=1, label='Barra Adicional')
            plt.plot([self.Bx[i], self.Px[i]], [self.By[i], self.Py[i]], 'g', linewidth=1, label='Barra Adicional')
            
            # Adicionando as coordenadas dos pontos
            plt.text(self.Ox[i], self.Oy[i], f'Âncoragem no chassis ({self.Ox[i]:.2f}, {self.Oy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Ax[i], self.Ay[i], f'A ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='left')
            plt.text(self.Bx[i], self.By[i], f'B ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='left')
            plt.text(self.Cx[i], self.Cy[i], f'Âncoragem no chassis ({self.Cx[i]:.2f}, {self.Cy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Px[i], self.Py[i], f'Cubo de Roda ({self.Px[i]:.2f}, {self.Py[i]:.2f})', fontsize=8, ha='left')
            
            # Adicionando o ângulo da barra de acoplamento
            plt.text((self.Ax[i] + self.Bx[i]) / 2, (self.Ay[i] + self.By[i]) / 2, f'{angulo_camber:.2f}°', fontsize=15, ha='center')

            if i == 0 or i == self.num_points/2 or i == (self.num_points - 1):
                plt.show()

            plt.grid()
            plt.axis('equal')
            plt.axis([-500, 1000, -350, 1100])
            plt.draw()
            plt.pause(0.01)
            plt.clf()

        for i in range(len(self.Py)-1, -1, -1):
            angulo_camber = self.calcular_camber(self.Ax[i], self.Ay[i], self.Bx[i], self.By[i])
            
            plt.plot([self.Ox[i], self.Ax[i]], [self.Oy[i], self.Ay[i]], 'b', linewidth=2, label='Barra de Entrada')
            plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'b', linewidth=2, label='Barra de Acoplamento', linestyle='dotted')
            plt.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i]], 'b', linewidth=2, label='Barra de Saída')
            plt.plot([self.Ax[i], self.Px[i]], [self.Ay[i], self.Py[i]], 'r', linewidth=1, label='Barra Adicional')
            plt.plot([self.Bx[i], self.Px[i]], [self.By[i], self.Py[i]], 'g', linewidth=1, label='Barra Adicional')

            # Adicionando as coordenadas dos pontos
            plt.text(self.Ox[i], self.Oy[i], f'Âncoragem no chassis ({self.Ox[i]:.2f}, {self.Oy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Ax[i], self.Ay[i], f'A ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='left')
            plt.text(self.Bx[i], self.By[i], f'B ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='left')
            plt.text(self.Cx[i], self.Cy[i], f'Âncoragem no chassis ({self.Cx[i]:.2f}, {self.Cy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Px[i], self.Py[i], f'Cubo de Roda ({self.Px[i]:.2f}, {self.Py[i]:.2f})', fontsize=8, ha='left')
            
            # Adicionando o ângulo da barra de acoplamento
            plt.text((self.Ax[i] + self.Bx[i]) / 2, (self.Ay[i] + self.By[i]) / 2, f'{angulo_camber:.2f}°', fontsize=15, ha='center')

            plt.grid()
            plt.axis('equal')
            plt.axis([-500, 1000, -350, 1100])
            plt.draw()
            plt.pause(0.01)
            plt.clf()

        # Imprimindo os ângulos de câmber inicial, estático, final e a taxa de variação
        camber_inicial = self.calcular_camber(self.Ax[0], self.Ay[0], self.Bx[0], self.By[0])
        camber_estatico = self.calcular_camber(self.Ax[int(self.num_points/2)], self.Ay[int(self.num_points/2)], self.Bx[int(self.num_points/2)], self.By[int(self.num_points/2)])
        camber_final = self.calcular_camber(self.Ax[int(self.num_points - 1)], self.Ay[int(self.num_points - 1)], self.Bx[int(self.num_points - 1)], self.By[int(self.num_points - 1)])
        taxa_variacao = (camber_final - camber_inicial) / (self.Py_step * len(self.Py))
        
        print(f"Ângulo de câmber rebound: {camber_inicial:.2f}°")
        print(f"Ângulo de câmber estático: {camber_estatico:.2f}°")
        print(f"Ângulo de câmber bump: {camber_final:.2f}°")
        print(f"Taxa de variação de câmber: {taxa_variacao:.2f}°/ mm")

    def plot_damper(self, damper_Vx_values, damper_Vy_values):
        """Plota a força do amortecedor em função da velocidade."""
        damper_Fy_values = []
        damper_Fx_values = []

        for vx in damper_Vx_values:
            self.damper_V = vx
            damper_Fx_values.append(self.Damper())

        for vy in damper_Vy_values:
            self.damper_V = vy
            damper_Fy_values.append(self.Damper())

        plt.figure(figsize=(8, 6))
        plt.plot(damper_Vx_values, damper_Fx_values, label='Velocidade')
        plt.title('Força em função da velocidade em X')
        plt.xlabel('Velocidade [mm/s]')
        plt.ylabel('Força [N]')
        plt.grid(True)
        plt.legend()
        plt.show()

        plt.figure(figsize=(8, 6))
        plt.plot(damper_Vy_values, damper_Fy_values, label='Velocidade')
        plt.title('Força em função da velocidade em Y')
        plt.xlabel('Velocidade [mm/s]')
        plt.ylabel('Força [N]')
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot_spring(self, spring_x_values, spring_y_values):
        """Plota a força do amortecedor em função da velocidade."""
        spring_Fx_values = []
        spring_Fy_values = []

        for x in spring_x_values:
            self.spring_x = x
            spring_Fx_values.append(self.Spring())
        for y in spring_y_values:
            self.spring_x = y
            spring_Fy_values.append(self.Spring())

        plt.figure(figsize=(8, 6))
        plt.plot(spring_x_values, spring_Fx_values, label='Deformação')
        plt.title('Força em função da deformação em X')
        plt.xlabel('Deformação [mm]')
        plt.ylabel('Força [N]')
        plt.grid(True)
        plt.legend()
        plt.show()

        plt.figure(figsize=(8, 6))
        plt.plot(spring_y_values, spring_Fy_values, label='Deformação')
        plt.title('Força em função da deformação em Y')
        plt.xlabel('Deformação [mm]')
        plt.ylabel('Força [N]')
        plt.grid(True)
        plt.legend()
        plt.show()

# Variáveis globais

# Carregue os dados da matriz e os parâmetros de uma planilha Excel
arquivo_excel = r'C:\Users\joao_\Downloads\Datasheet - Dinâmica (1).xlsx'#ADICIONE O CAMINHO DO ARQUIVO EXCEL

# Planilha com Matriz de Torque e Potência"
df_matriz = pd.read_excel(arquivo_excel, sheet_name='Matriz Torque x Potência')
matriz_dados = df_matriz.dropna().to_dict('records')  # Remove linhas com células vazias

# Parâmetros do sistema
df_parametros = pd.read_excel(arquivo_excel, sheet_name='Transmissão')

# Parâmetros gerais
df_gerais = pd.read_excel(arquivo_excel, sheet_name='Variáveis Globais')

# Extraia os parâmetros 
cgx = df_gerais.at[1, 'Valor']
cgy = df_gerais.at[2, 'Valor']
massa = df_gerais.at[0, 'Valor']
etex = df_gerais.at[8, 'Valor']
cfat = df_parametros.at[0, 'Valor']
rpneu = df_gerais.at[6, 'Valor']
acpi = df_parametros.at[1, 'Valor']
redp = df_parametros.at[2, 'Valor']
red1 = df_parametros.at[3, 'Valor']


# Definindo os parâmetros de freio
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

# Dados experimentais
ratio = np.linspace(-1, 1, 19)
angles = np.array([-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
tire_lateral_forces_1 = np.array([-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890, 2000, 2090])
tire_auto_align_moment_1 = np.array([-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09, 50.90, 50.10, 50.81, 48.12, 48.83])

camber_experimental = [-2060.0, -1950.0, -1840.0, -1700.0, -1540.0, -1350.0, -1130.0, -860.0, -480.0, -30.0, 460.0, 880.0, 1230.0, 1490.0, 1720.0, 1910.0, 2090.0, 2230.0, 2310.0]
camber_experimental_1 = [-1940.0, -1860.0, -1750.0, -1610.0, -1450.0, -1260.0, -1050.0, -760.0, -400.0, 60.0, 540.0, 970.0, 1290.0, 1550.0, 1780.0, 1980.0, 2150.0, 2280.0, 2370.0]
camber_experimental_15 = [-1840.0, -1750.0, -1670.0, -1520.0, -1370.0, -1180.0, -960.0, -680.0, -310.0, 130.0, 610.0, 1020.0, 1360.0, 1630.0, 1850.0, 2040.0, 2220.0, 2360.0, 2430.0]

# Instâncias gerais
dynamics_instance = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45, tire_Ca=0)

#Instância de Drivetrain

# Crie uma instância da classe Drivetrain
drivetrain = Drivetrain.generate_model(cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, matriz_dados)

# Calcule a transmissão
drivetrain.show_results()

# Calcule o desempenho do carro
drivetrain.print_car_performance()

# Dimensione os semi-eixos
drivetrain.HalfShaftsSizing()

####### Cálculo de Slip Ratio ############

#Recebendo os dados da performance do carro. Aqui que se encontra dados da velocidade angular
performance_veiculo = drivetrain.CarPerformance()

#Filtrando a velocidade angular
velocidade_angular = [dado["va"] for dado in performance_veiculo]

#Chama a função de slip ratio e salva seus valores numa lista utilizando uma velocidade arbitrária e fixa para o veículo
slip_ratio = dynamics_instance.slip_ratio_1(velocidade_angular=np.array(velocidade_angular), raio_pneu=params['Rdp'])

rpm = drivetrain.CurveTorquePower()[0]

#Plotagem do gráfico de slip ratio e velcidade angular X RPM
dynamics_instance.show_slip_ratio(rpm, slip_ratio, velocidade_angular)

#Transformando num arranjo para permitir a realização de contas dentro do código de Tire
slip_ratio_array = np.array(slip_ratio)

#Recebendo slip angles
slip_angles = np.linspace(-9, 9, 1000)  # Ângulos de deslizamento variando de -9 a 9 graus

#Criando uma instância específica de Drivetrain
dynamics_instance_d = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio_array, tire_friction_coef=1.45, tire_Ca=0)

#Parâmetros Pacejka
result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]

#Recebendo valores de força lateral, momento auto alinhante e força lateral
tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = dynamics_instance_d.Tire_forces(result)

#Plotagem dos gráficos
dynamics_instance_d.plot_graph(tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces)

#Plotagem Força Longitudinal x RPM
Tire.show_longitudinal_rpm(rpm, tire_longitudinal_forces)

######################

#Instâncias de Freio

# Cria uma instância da classe Dynamics e da classe BrakeSystem
BrakeSystem = BrakeSystem(params)

# Aplica o sistema de freio às forças do pedal e obtém resultados
resultados, forca_frenagem, torque_ajustado, forca_f, torque_disco_freio, resistencia_rolamento, torque_resistencia_rolamento = BrakeSystem.apply_brake(pedal_forces)

# Calcula a velocidade angular, desaceleração angular e lista de velocidades
desaceleracao_angular, velocidade_angular, lista_velocidade, inertia_wheel = BrakeSystem.calculate_angular_velocity(torque_ajustado)

# Calcula a razão de escorregamento com base na velocidade angular, raio do pneu e forças do pedal
calculo_slip_ratio = dynamics_instance.slip_ratio_1(velocidade_angular, params['Rdp'])

# Cria instâncias da classe Dynamics com diferentes parâmetros de pneus e coeficiente de atrito
dynamics_instance_b = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=calculo_slip_ratio, tire_friction_coef=1.45, tire_Ca=0)

# Calcula a força longitudinal do pneu
forca_longitudinal = dynamics_instance_b.Tire_forces(result)[2]

# Plotagem dos gráficos e print dos resultados de freio
dynamics_instance_b.plot_graph_slip_ratio(calculo_slip_ratio)
BrakeSystem.show_graph(forca_longitudinal, pedal_forces)
BrakeSystem.print_resultados()

######################

# Instâncias de Cinemática

# Criação do objeto Kinematics e execução dos cálculos
kinematics = Kinematics(L0=500, L1=500, L2=450, L3=500)
camber_angle = kinematics.calcular_cinematica()[0]
Ax, Ay = kinematics.calcular_cinematica()[1], kinematics.calcular_cinematica()[3]
Vx, Vy = kinematics.calcular_cinematica()[2], kinematics.calcular_cinematica()[4]
kinematics.plotar_cinematica()

# Instância para o amortecedor do tipo 'Integrated'
teste_integrated = Kinematics(damper_type='Integrated', damper_F_static=50, damper_K_friction=10, damper_F_viscous=10)
teste_integrated.plot_damper(damper_Vx_values=np.linspace(1, 10, 100), damper_Vy_values=Vy)

# Instância para o amortecedor do tipo 'Coulumb'
teste_coulumb = Kinematics(damper_type='Coulumb', damper_F_static=50, damper_K_friction=10)
teste_coulumb.plot_damper(damper_Vx_values=np.linspace(1, 10, 100), damper_Vy_values=Vy)

# Instância para a mola do tipo 'Hooke'
teste_hooke = Kinematics(spring_type= 'Hooke', spring_k=40, spring_x=0, spring_non_lin_coef=None)
teste_hooke.plot_spring(spring_x_values=Ax , spring_y_values=Ay)

# Instância para a mola do tipo 'Softening'
teste_hooke = Kinematics(spring_type= 'Softening', spring_k=40, spring_x=0, spring_non_lin_coef=0.1)
teste_hooke.plot_spring(spring_x_values=Ax , spring_y_values=Ay)

# Plotar forças do pneu com câmber
dynamics_instance_k = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45, tire_Ca=camber_angle)
dynamics_instance_k.plot_graph(dynamics_instance_k.Tire_forces(result)[0], dynamics_instance_k.Tire_forces(result)[1], dynamics_instance_k.Tire_forces(result)[2])


######################

# Instâncias de Pneu

# Criação do mecanismo
mechanism = Tire(B0 = 1393, B1 = 300, B2 = 1400, B3 = 300)
# Cálculos
mechanism.calculate_kinematics()
# Plots
mechanism.run()

# Calculando as forças previstas e o momento de auto-alinhamento
predicted_tire_lateral_forces, predicted_tire_auto_align_moment, predicted_tire_longitudinal_forces = dynamics_instance.Tire_forces(result)

# Instanciando a classe Tire com diferentes valores de camber
dynamics_camber = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45, tire_Ca=0.5)
dynamics_camber_1 = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45, tire_Ca=1)
dynamics_camber_15 = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45, tire_Ca=1.5)

# Calculando as forças laterais previstas para diferentes valores de camber
forca_camber = dynamics_camber.Tire_forces(result)[0]
forca1_camber_1 = dynamics_camber_1.Tire_forces(result)[0]
forca2_camber_15 = dynamics_camber_15.Tire_forces(result)[0]

# Plotando as curvas com camber comparadas com dados experimentais
dynamics_camber.plot_camber(forca_camber, forca1_camber_1, forca2_camber_15, camber_experimental, camber_experimental_1, camber_experimental_15, angles)

# Plotando as curvas previstas com os dados experimentais
dynamics_instance.plot_graph(predicted_tire_lateral_forces, predicted_tire_auto_align_moment, predicted_tire_longitudinal_forces,
                            tire_lateral_experimental= tire_lateral_forces_1, tire_auto_align_experimental= tire_auto_align_moment_1, angles=angles, ratio=ratio)

# Criando uma instância da classe e executando os métodos
rigidez_deformacao = Tire(track_y=0, tire_k =120)
rigidez_deformacao.calcular_forca()
rigidez_deformacao.plotar_deformacao(track_y_values=np.linspace(0, 25, 1000))
