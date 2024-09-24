import math
import matplotlib.pyplot as plt
import numpy as np

class Drivetrain:

    '''
    Drivetrain Model

    Representação do modelo de Drivetrain de um veículo

    Parameters
    ---------
    cgx: centro de massa no eixo x
    cgy: centro de massa no eixo y
    massa: massa do veículo
    entre_eixos: distância entre eixos
    coeficiente_atrito: coeficiente de atrito
    reio_pneu: raio do pneu do veículo
    aceleracao_ideal:
    reducao_primaria:
    reducao_unica: 
    rpm: rpm de entrada(dado do motor)
    torque: torque de entrada(dado do motor)
    cp: relação coroa-pinhão
    '''

    def __init__(self, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, rpm, torque, cp, potencia):
        self.cgx = cgx
        self.cgy = cgy
        self. massa = massa
        self.entre_eixos =entre_eixos
        self.coeficiente_atrito = coeficiente_atrito
        self.raio_pneu = raio_pneu
        self.aceleracao_ideal = aceleracao_ideal
        self.reducao_primaria = reducao_primaria
        self. reducao_unica = reducao_unica
        self.rpm = rpm
        self.torque = torque
        self.cp = cp
        self.new_rpm = 0
        self.potencia = potencia

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

    def showResults(self):
        
        peso, reacao_traseira, reacao_dianteira, forca_trativa, transferencia_longitudinal, carga_traseira, pico_torque_traseiro, carga_pneu, torque_pneu, reducao_final, torque_necessario_motor, aceleracao_primaria_real, aceleracao_primaria_ideal, aceleraco_real_final, forca_trativa_ideal, torque_pneu_ideal, torque_motor_ideal, transferencia_carga_ideal, transferencia_carga_real, carga_traseira_ideal = Drivetrain.CalculateOutputs(self)
                
        # Print dos resultados obtidos
        print(f'''
            Resultados:
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

    def interpolatedTorque(self):
            
            if self.new_rpm:
                variacao = range(self.rpm, self.new_rpm, 30)
            else:
                variacao = [self.rpm]
            
            torques = []
            
            for rpm in variacao:
                if rpm > 0:
                    torque_interpolado = self.potencia / (2 * math.pi * rpm)
                else:
                    torque_interpolado = 0 
                torques.append(torque_interpolado)
            
            return variacao, torques

    def CarPerformance(self):
        
        peso = self.massa * 9.81
        rendimento_transmissao = 0.9
        transmissao_motor_roda = 0.9
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        basic_f = 0.015
        speed_f = 0.012
        
        variacao = []
        
        if self.new_rpm:
            variacao = range(self.rpm, self.new_rpm, 30)
            self.rpm = self.new_rpm
            self.new_rpm = 0
            
        else:
            variacao = [self.rpm]
        
        parametros = []
        

        _, torque_interpolados = self.interpolatedTorque()

        for rpm, torque in zip(variacao, torque_interpolados):
            # Cálculo da força trativa (N) usando o torque interpolado
            forca_trativa = ((torque * self.reducao_primaria * self.reducao_unica * self.cp) / (self.raio_pneu * 0.001)) * rendimento_transmissao

            # Cálculo da velocidade angular (rad/s)
            velocidade_angular = (rpm * 2 * math.pi) / (60 * self.reducao_primaria * self.reducao_unica * self.cp)

            # Cálculo da velocidade linear (km/h)
            velocidade_linear = ((velocidade_angular * (self.raio_pneu * 0.001)) * transmissao_motor_roda) * 3.6

            # Cálculo da força de arrasto (N)
            fa = (densidade_ar * velocidade_linear ** 2 * coeficiente_arrasto * area_frontal) / 2

            # Cálculo da resistência de rolamento (N)
            rr = (basic_f + (3.24 * speed_f * ((velocidade_linear / 100 * 0.44704) ** 2.5))) * peso

            # Cálculo da força final (N)
            forca_final = forca_trativa - fa - rr

            # Armazenar os parâmetros calculados em um dicionário
            parametro = {"forca_trativa": forca_trativa, "va": velocidade_angular, "velocidade_linear": velocidade_linear,
                         "fa": fa, "rr": rr, "forca_final": forca_final}

            parametros.append(parametro)

        # Retornar os parâmetros calculados
        return parametros, variacao
    
    def printCarPerformance(self):
            
            rpm = self.rpm
            new_rpm = self.new_rpm
            
            performance, rpm_faixa = Drivetrain.CarPerformance(self)
            
            print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]")
            
            for param in performance:
                print(f"{param['forca_trativa']}\t{param['va']}\t{param['velocidade_linear']}")
            
            print("\nForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
            for param in performance:
                print(f"{param['fa']}\t{param['rr']}\t{param['forca_final']}")

            self.rpm = rpm
            self.new_rpm = new_rpm

    def HalfShaftsSizing(self, fsi=1.25, tet=786, tec=471.6, dif=1):
        
        _, torque_interpolados = self.interpolatedTorque()

        # Obtendo o maior torque interpolado
        torque_max_motor = max(torque_interpolados)

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

class Tire:
    
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=0, B0=0, B1=0, 
                B2=1, B3=1, omega=315, slip_angle_start=-9, slip_angle_end=9, angle_step=0.5, WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
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
        self.angle = self.theta2_rad[0]

        # Inicialização das listas de resultados para armazenar dados ao longo do cálculo
        self.AC = []  # Lista para armazenar AC
        self.beta = []  # Lista para armazenar beta
        self.psi = []  # Lista para armazenar psi
        self.lamda = []  # Lista para armazenar lambda
        self.theta3 = []  # Lista para armazenar theta3
        self.theta4 = []  # Lista para armazenar theta4
        self.Ox, self.Oy = 0, 0  # Lista para armazenar coordenadas de O
        self.Ax, self.Ay = [], []  # Listas para armazenar coordenadas de A
        self.Bx, self.By = [], []  # Listas para armazenar coordenadas de B
        self.Cx, self.Cy = B0, 0 # Listas para armazenar coordenadas de C
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

    @staticmethod
    def slip_ratio_1(velocidade_angular, raio_pneu):
       
        velocidade_linear = initial_speed
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1

        return value

    @staticmethod
    def show_slip_ratio(rpm_values, slip_ratio, velocidade_angular):
       
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
        
    def plot_graph(self, tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces, tire_lateral_experimental=None, tire_auto_align_experimental=None, angles=None, ratio=None):

        # Definindo um tamanho para a figura
        plt.figure(figsize=(20, 7))

        # Plotagem força lateral
        plt.subplot(1, 3, 1)
        plt.plot(self.tire_Sa, tire_lateral_forces, label='Curva Otimizada')
        plt.scatter(angles, tire_lateral_experimental, color='red', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Força Lateral do Pneu (N)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Plotagem torque auto-alinhante
        plt.subplot(1, 3, 2)
        plt.plot(self.tire_Sa, tire_auto_align_moment, label='Curva Otimizada')
        plt.scatter(angles, tire_auto_align_experimental, color='blue', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Torque auto-alinhante (N.m)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Plotagem força longitudinal
        plt.subplot(1, 3, 3)
        plt.plot(self.tire_Ls, tire_longitudinal_forces, label='Curva Sem Otimizar')
        plt.xlabel('Taxa de Escorregamento Longitudinal (Admensional)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal - Sem Dados Experimentais')
        plt.legend()
        plt.grid(True)

        plt.tight_layout(pad=3.0)  # Aumentando a distância entre os subplots
        plt.show()


def uniteExample():
    dt_model = Drivetrain(
        cgx = 853,  # mm
        cgy = 294,  # mm
        massa = 347,  # kg
        entre_eixos = 1567,  # mm
        coeficiente_atrito = 0.9 , # coeficiente de atrito
        raio_pneu = 259,  # mm
        aceleracao_ideal = 1.2,  # g
        reducao_primaria = 2.12,  # redução primária
        reducao_unica = 2.76,
        rpm = 0,
        torque= 244.14,
        cp = 2.22,
        potencia = 23000)
    
    dt_model.showResults()
    
    dt_model.interpolatedTorque()

    dt_model.HalfShaftsSizing()
    
    dt_model.new_rpm = 6000

    #Exemplo em lista
    dt_model.printCarPerformance() 

    #Recebendo os dados da performance do carro. Aqui que se encontra dados da velocidade angular
    performance_veiculo, rpm_faixa = dt_model.CarPerformance()

    #Filtrando a velocidade angular
    velocidade_angular = [dado['va'] for dado in performance_veiculo]
       
    #Transformando num array    
    velocidade_angular = np.array(velocidade_angular)

    #Transformando num array  
    rpm_faixa = np.array(rpm_faixa)

    #Calcular o slip ratio
    slip_ratio = Tire.slip_ratio_1(velocidade_angular, 0.259)

    #Plotagem de gráfico do slip ratio e saídas de seus dados no terminal
    Tire.show_slip_ratio(rpm_faixa, slip_ratio, velocidade_angular)

    #Salvando os dados como array para cálculo de força longitudinal
    slip_ratios = np.array(slip_ratio)
    
    #Criando instância da classe Tire
    Slip_model = Tire(tire_Fz=1500, tire_Sa=0, tire_Ls=slip_ratios, tire_friction_coef=1.45, tire_Ca=0)
    
    #Dados experimentais para instância em Tire
    result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]
    
    #Recebendo valores de força lateral, torque auto allinhante e força longitudinal
    tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = Slip_model.Tire_forces(result)
    
    #Plotagem de gráficos
    Slip_model.plot_graph( tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces)

initial_speed = 1

uniteExample()

Tire()
