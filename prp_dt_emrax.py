import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

class Torque:

    def __init__(self):
        self.rpm_interp = None
        self.torque_continuo_interp = None
        self.torque_pico_interp = None
        self.potencia_continua_interp = None
        self.potencia_pico_interp = None
        self._setup_data()

    def _setup_data(self): 
        rpm = np.array([0, 1000, 2000, 3000, 3500, 4000, 4500, 5000, 5500])
        potencia_continua = np.array([0, 20, 40, 60, 75, 85, 90, 95, 97])
        potencia_pico = np.array([0, 50, 100, 150, 180, 200, 210, 215, 215])
        torque_continuo = np.array([200, 200, 200, 200, 190, 170, 150, 130, 110])
        torque_pico = np.array([500, 500, 500, 500, 480, 450, 400, 370, 350])

        self.rpm_interp = np.linspace(0, 5500, 1000)
        self.torque_continuo_interp = interp1d(rpm, torque_continuo, kind='cubic')(self.rpm_interp)
        self.torque_pico_interp = interp1d(rpm, torque_pico, kind='cubic')(self.rpm_interp)
        self.potencia_continua_interp = interp1d(rpm, potencia_continua, kind='cubic')(self.rpm_interp)
        self.potencia_pico_interp = interp1d(rpm, potencia_pico, kind='cubic')(self.rpm_interp)

    def plot(self):
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.set_xlabel('Motor speed (RPM)', fontsize=12)
        ax1.set_ylabel('Power (kW)', color='tab:green', fontsize=12)
        ax1.plot(self.rpm_interp, self.potencia_continua_interp, label='Continuous power', color='tab:green', linestyle='--', linewidth=2)
        ax1.plot(self.rpm_interp, self.potencia_pico_interp, label='Peak power', color='tab:green', linewidth=2)
        ax1.set_ylim(0, 220)
        ax1.tick_params(axis='y', labelcolor='tab:green')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Torque (Nm)', color='tab:blue', fontsize=12)
        ax2.plot(self.rpm_interp, self.torque_continuo_interp, label='Continuous torque', color='tab:blue', linestyle='--', linewidth=2)
        ax2.plot(self.rpm_interp, self.torque_pico_interp, label='Peak torque', color='tab:blue', linewidth=2)
        ax2.set_ylim(0, 550)
        ax2.tick_params(axis='y', labelcolor='tab:blue')

        fig.suptitle('EMRAX 268MV LC - Torque and Power Curves', fontsize=14)
        ax1.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98))
        ax2.legend(loc="upper right", bbox_to_anchor=(0.98, 0.98))
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def example():
        motor = Torque()
        return motor.torque_pico_interp[:]
        

# Teste de execução
valores = Torque.example()
print(valores)



class Drivetrain:

    def __init__(self, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, rpm, cp, tempo_i):
        self.cgx = cgx
        self.cgy = cgy
        self. massa = massa
        self.entre_eixos =entre_eixos
        self.coeficiente_atrito = coeficiente_atrito
        self.raio_pneu = raio_pneu
        self.aceleracao_ideal = aceleracao_ideal
        self.reducao_primaria = reducao_primaria
        self.reducao_unica = reducao_unica
        self.rpm = rpm
        self.cp = cp
        self.new_rpm = 0
        self.tempo_i = tempo_i
        self.tempo_f = 0

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
        
    def CarPerformance(self):
        # Definindo parâmetros do veículo
        peso = self.massa * 9.81  # Força peso
        rendimento_transmissao = 0.9
        transmissao_motor_roda = 0.9
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        c_r = 0.025  # Coeficiente de resistência ao rolamento 
        f_v = 0.012  
        Iw = (self.massa * (self.raio_pneu * 0.001)) / 2  # Momento de inércia da roda

        # Listas para armazenar os dados
        variacao_rpm = []
        variacao_tempo = []
        variacao_torque = []

        if self.new_rpm:
            variacao_rpm = range(self.rpm, self.new_rpm, 30)
            variacao_tempo = np.linspace(self.tempo_i, self.tempo_f, len(variacao_rpm))
            variacao_torque = Torque.example()
            self.tempo_i = self.tempo_f
            self.rpm = self.new_rpm
            self.new_rpm = 0
            self.tempo_f = 0
           
        else:
            variacao_rpm = [self.rpm]
            variacao_tempo = [self.tempo_i]

        parametros = []

        for rpm, tempo, torque in zip(variacao_rpm, variacao_tempo, variacao_torque):

            # Inicializa a velocidade angular (rad/s) a partir da RPM
            velocidade_angular = (rpm * 2 * math.pi) / (60 * self.reducao_primaria * self.reducao_unica)

            # Cálculo da força trativa (N)
            forca_trativa = ((torque * self.reducao_primaria * self.reducao_unica) / (self.raio_pneu * 0.001)) * rendimento_transmissao

            # Cálculo da resistência ao rolamento (N)
            resistencia_rolamento = (c_r + (3.24 * f_v * (((velocidade_angular) * (self.raio_pneu * 0.001)) ** 2.5))) * peso

            # Cálculo do momento My (Nm)
            momento_y = c_r * peso * self.raio_pneu * 0.001  

            # Cálculo da aceleração angular (rad/s^2)
            aceleracao_angular = (torque - momento_y) / Iw

            # Atualizando a velocidade angular (rad/s)
            velocidade_angular += aceleracao_angular * tempo  # Atualiza a velocidade angular com base no tempo

            # Cálculo da velocidade tangencial (km/h)
            velocidade_tangencial = ((velocidade_angular * (self.raio_pneu * 0.001)) * transmissao_motor_roda) * 3.6

            # Cálculo da força de arrasto (N)
            forca_arrasto = (densidade_ar * (velocidade_tangencial / 3.6)** 2 * coeficiente_arrasto * area_frontal) / 2

            # Cálculo da força final (N)
            forca_final = max(0, forca_trativa - forca_arrasto - resistencia_rolamento)
            
            # Armazenar os parâmetros calculados em um dicionário
            parametro = {
                "ft": forca_trativa,
                "va": velocidade_angular,
                "vt": velocidade_tangencial,
                "fa": forca_arrasto,
                "rr": resistencia_rolamento,
                "my": momento_y,
                "ff": forca_final
            }

            parametros.append(parametro)
      
        return parametros, variacao_rpm, variacao_tempo

    def printCarPerformance(self):
                  
        rpm = self.rpm
        new_rpm = self.new_rpm
        
        performance, variacao_rpm, variacao_tempo = self.CarPerformance()
        
        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Tangencial [km/h]")
        
        for param in performance:
            print(f"{param['ft']}\t{param['va']}\t{param['vt']}")
        
        print("\nForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for param in performance:
            print(f"{param['fa']}\t{param['rr']}\t{param['ff']}")

        self.rpm = rpm
        self.new_rpm = new_rpm

    def HalfShaftsSizing(self, fsi=1.25, tet=786, tec=471.6, dif=1):
        
        # Obtendo o maior torque do motor a partir dos dados experimentais 
        torque_max_motor = max(Torque.example())
        
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

    def calcular_velocidade_linear(self, forcas_finais, massa, tempos):
        tempos = np.array(tempos)[::-1]  # Garante que os tempos estão em ordem crescente
        forcas_finais = np.array(forcas_finais)  # Ajusta as forças para bater com os tempos
        
        aceleracoes = forcas_finais / massa  # Calcula aceleração em cada instante
        velocidades = np.zeros_like(tempos)  # Inicializa vetor de velocidades
        
        print("\nValores de Aceleração e Tempo:")
        for i in range(len(tempos)):
            print(f"t = {tempos[i]:.2f}s, a = {aceleracoes[i]:.2f} m/s²")

        # Integração numérica usando método de Euler
        print("\nCálculo de Velocidade:")
        for i in range(1, len(tempos)):
            dt = tempos[i] - tempos[i - 1]  # Intervalo de tempo
            velocidades[i] = velocidades[i - 1] + aceleracoes[i] * dt
            print(f"t = {tempos[i]:.2f}s, dt = {dt:.2f}s, v = {velocidades[i]:.2f} m/s")

        print("\nVelocidades Calculadas:", velocidades)
        return velocidades

    @staticmethod
    def slip_ratio_1(velocidade_angular, raio_pneu, velocidade_linear):
        v_roda = velocidade_angular * raio_pneu
        numerador = v_roda - velocidade_linear
        denominador = np.maximum(v_roda, velocidade_linear)
        denominador = np.maximum(denominador, 1e-3)  # evita divisão por zero
        return numerador / denominador

        
    @staticmethod
    def show_slip_ratio(rpm_values, slip_ratio, velocidade_angular):
         
         # Ajuste para garantir tamanhos compatíveis
        min_length = min(len(rpm_values), len(slip_ratio))
        rpm_values = np.array(rpm_values[:min_length])
        slip_ratio = np.array(slip_ratio[:min_length])
        velocidade_angular = np.array(velocidade_angular[:min_length])
       

        print("Valores do Slip Ratio: ")
        for dado in slip_ratio:
            print(dado)
            
        plt.figure(figsize=(10, 8))
        
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
        
    def plot_graph(self, tire_longitudinal_forces):

        # Definindo um tamanho para a figura
        plt.figure(figsize=(10, 8))

        # Plotagem força longitudinal
        plt.subplot(1, 1, 1)
        plt.plot(self.tire_Ls, tire_longitudinal_forces, label='Curva Sem Otimizar')
        plt.xlabel('Slip Ratio (Admensional)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal - Sem Dados Experimentais')
        plt.legend()
        plt.grid(True)

        plt.tight_layout(pad=3.0)  # Aumentando a distância entre os subplots
        plt.show()


def unite_example():
    # Instância do modelo Drivetrain com os parâmetros fornecidos
    dt_model = Drivetrain(
        cgx=853,  # Centro de gravidade no eixo X [mm]
        cgy=294,  # Centro de gravidade no eixo Y [mm]
        massa=347,  # Massa do veículo [kg]
        entre_eixos=1567,  # Entre-eixos [mm]
        coeficiente_atrito=0.9,  # Coeficiente de atrito
        raio_pneu=259,  # Raio do pneu [mm]
        aceleracao_ideal=1.2,  # Aceleração ideal [g]
        reducao_primaria=2.12,  # Redução primária
        reducao_unica=2.76,  # Redução única
        rpm=0,  # RPM inicial
        cp=2.22,  # Constante de potência
        tempo_i=0  # Tempo inicial [s]
    )

    motor = Torque()
    motor.plot()

    # Exibindo resultados iniciais do Drivetrain
    dt_model.showResults()

    # Dimensionamento dos semi-eixos
    dt_model.HalfShaftsSizing()

    # Atualização de parâmetros dinâmicos
    dt_model.new_rpm = 5500  # Novo RPM final
    dt_model.tempo_f = 10  # Tempo final [s]

    # Exibição da performance do carro em forma de lista
    dt_model.printCarPerformance()

    # Coleta dos dados de performance do carro
    performance_veiculo, variacao_rpm, variacao_tempo = dt_model.CarPerformance()

    # Extração das forças finais dos dados de performance
    forcas_finais = [dado["ff"] for dado in performance_veiculo]

    # Cálculo das velocidades lineares
    velocidades_lineares = Tire().calcular_velocidade_linear(
        forcas_finais, dt_model.massa, variacao_tempo
    )

    # Cálculo do slip ratio para cada ponto
    velocidade_angular = [dado["va"] for dado in performance_veiculo]
    velocidade_angular = np.array(velocidade_angular)
    raio_pneu_metros = dt_model.raio_pneu / 1000  # Conversão do raio para metros
    slip_ratios = Tire.slip_ratio_1(velocidade_angular, raio_pneu_metros, velocidades_lineares)

    # Plotagem dos resultados de slip ratio
    Tire.show_slip_ratio(variacao_rpm, slip_ratios, velocidade_angular)

    # Reobtendo os dados de performance do carro para análise adicional
    performance_veiculo, variacao_rpm, variacao_tempo = dt_model.CarPerformance()

    # Filtrando a velocidade angular dos dados de performance
    velocidade_angular = [dado["va"] for dado in performance_veiculo]
    velocidade_angular = np.array(velocidade_angular)
    variacao_rpm = np.array(variacao_rpm)

    # Cálculo do slip ratio em forma de array
    slip_ratio = Tire.slip_ratio_1(velocidade_angular, 0.259, velocidades_lineares)
    slip_ratios = np.array(slip_ratio)

    # Instância da classe Tire com os parâmetros fornecidos
    slip_model = Tire(
        tire_Fz=1500,  # Carga vertical no pneu [N]
        tire_Sa=0,  # Ângulo de escorregamento lateral [rad]
        tire_Ls=slip_ratios,  # Escorregamento longitudinal
        tire_friction_coef=1.45,  # Coeficiente de fricção
        tire_Ca=0  # Ângulo de camber
    )

    # Parâmetros experimentais fornecidos
    result = [
        0.3336564873588197,
        1.6271741344929977,
        1,
        4.3961693695846655,
        931.4055775279057,
        366.4936818126405,
    ]

    # Cálculo das forças e do momento de auto-alinhamento
    tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = slip_model.Tire_forces(result)

    # Plotagem do gráfico das forças longitudinais
    slip_model.plot_graph(tire_longitudinal_forces)


# Chamada da função principal
unite_example()

