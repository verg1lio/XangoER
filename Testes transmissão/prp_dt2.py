import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

class Motor:

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
    def get_torque():
        motor = Motor()
        return motor.torque_continuo_interp[:] * 2

    @staticmethod
    def get_potencia():
        motor = Motor()
        return motor.potencia_continua_interp[:] 
    @staticmethod
    def get_rpm():
        motor = Motor()
        return motor.rpm_interp[:]
    
# Teste de execução

# valores_t = Motor.get_torque()
# valores_p = Motor.get_potencia()
# print(f'torque:{valores_t}')
# print(f'potencia:{valores_p}')

class Drivetrain:

    def __init__(self, cgx, cgy, massa, entre_eixos, raio_pneu, reducao_primaria, reducao_final, cp, tempo_i):
        self.cgx = cgx
        self.cgy = cgy
        self.massa = massa
        self.entre_eixos =entre_eixos
        self.raio_pneu = raio_pneu
        self.reducao_primaria = reducao_primaria
        self.reducao_final = reducao_final
        self.cp = cp
        self.tempo_i = tempo_i
        self.tempo_f = 0

    def Reduções(self, eficiencia_corrente=0.96):
     
        # Dados do motor
        rpm_motor = Motor.get_rpm()
        torque_motor = Motor.get_torque()
        potencia_motor = Motor.get_potencia()


        # Etapa 1: Redução primária (Redutor)
        rpm_pos_redutor = rpm_motor / self.reducao_primaria
        torque_pos_redutor = torque_motor * self.reducao_primaria
        potencia_pos_redutor = potencia_motor  # desprezando perdas aqui

        # Etapa 2: Perdas por corrente de rolos (efeito poligonal e atrito)
        rpm_pos_corrente = rpm_pos_redutor  # corrente não altera rotação, só transmite
        torque_pos_corrente = torque_pos_redutor * eficiencia_corrente
        potencia_pos_corrente = potencia_pos_redutor * eficiencia_corrente
        
        # Etapa 3: Redução final (Diferencial)
        rpm_rodas = rpm_pos_corrente / self.reducao_final
        torque_rodas = torque_pos_corrente * self.reducao_final
        potencia_rodas = potencia_pos_corrente  # perdas desprezadas nesta etapa

        rendimento_transmissao = np.nanmean((potencia_rodas[1:] / potencia_motor[1:]))

        return  [torque_rodas, potencia_rodas, rpm_rodas,rendimento_transmissao]
    
    def CarPerformance(self):
        # Parâmetros do veículo
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        c_r = 0.025
        f_v = 0.012

        # Dados já reduzidos pelo sistema de transmissão
        torque_reduzido, potencia_reduzida, rpm_reduzido, rendimento_transmissao = self.Reduções()

        # Tempo de simulação
        self.tempo_f = len(rpm_reduzido) * 0.01  # tempo final automático baseado na resolução
        variacao_tempo = np.linspace(self.tempo_i, self.tempo_f, len(rpm_reduzido))

        parametros = []
        velocidade_veiculo = 0  # m/s inicial

        for rpm, torque, potencia, tempo in zip(rpm_reduzido, torque_reduzido, potencia_reduzida, variacao_tempo):
            
            # Velocidade angular da roda (rad/s)
            velocidade_angular = (rpm * 2 * np.pi) / 60

            # Velocidade linear na roda (m/s)
            velocidade_linear_roda = velocidade_angular * (self.raio_pneu * 0.001)

            # Força trativa 
            forca_trativa = torque / (self.raio_pneu * 0.001) * rendimento_transmissao

            # Resistência ao rolamento
            resistencia_rolamento = (c_r + 3.24 * f_v * velocidade_linear_roda ** 2.5) * peso

            # Força de arrasto
            forca_arrasto = 0.5 * densidade_ar * coeficiente_arrasto * area_frontal * velocidade_veiculo ** 2

            # Força líquida resultante
            forca_final = np.maximum(forca_trativa - forca_arrasto - resistencia_rolamento, 0)

            # Aceleração (F = ma)
            aceleracao = forca_final / self.massa

            # Integração da velocidade com método de Euler
            if len(parametros) == 0:
                dt = 0
            else:
                dt = tempo - parametros[-1]["tempo"]

            velocidade_veiculo += aceleracao * dt  # m/s
            velocidade_linear_carro = velocidade_veiculo * 3.6  # km/h

            parametro = {
                "ft": forca_trativa,
                "va": velocidade_angular,
                "vlc": velocidade_linear_carro,
                "fa": forca_arrasto,
                "rr": resistencia_rolamento,
                "ff": forca_final,
                "tempo": tempo
            }

            parametros.append(parametro)

        return parametros, rpm_reduzido, variacao_tempo


    def printCarPerformance(self, performance):

        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear do Carro [km/h]")
        for i in range(0, len(performance), 10):  
            param = performance[i]
            print(f"{param['ft']:.2f}\t\t\t{param['va']:.2f}\t\t\t\t{param['vlc']:.2f}")

        print("\nForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for i in range(0, len(performance), 10): 
            param = performance[i]
            print(f"{param['fa']:.2f}\t\t\t{param['rr']:.2f}\t\t\t{param['ff']:.2f}")

        rpm_motor = Motor.get_rpm()

       
        velocidade_angular = [param["va"] for param in performance]
        velocidade_linear = [param["vlc"] for param in performance]

        # Plot dos gráficos
        plt.figure(figsize=(12, 6))

        # Plot Velocidade Angular (rad/s) x RPM
        plt.subplot(1, 2, 1)
        plt.plot(rpm_motor, velocidade_angular, color='tab:red')
        plt.xlabel("RPM Motor")
        plt.ylabel("Velocidade Angular [rad/s]")
        plt.title("Velocidade Angular x RPM")
        plt.grid(True)

        # Plot Velocidade Linear (km/h) x RPM
        plt.subplot(1, 2, 2)
        plt.plot(rpm_motor, velocidade_linear, color='tab:purple')
        plt.xlabel("RPM Motor")
        plt.ylabel("Velocidade Linear [km/h]")
        plt.title("Velocidade Linear x RPM")
        plt.grid(True)

        plt.tight_layout()
        plt.show()

class Tire:
    
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=0, B0=0, B1=0, B2=1, B3=1, omega=315, slip_angle_start=-9, slip_angle_end=9, angle_step=0.5, WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
        
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
    def SlipRatio(velocidade_angular, raio_pneu, velocidade_linear):

            value = (1-(velocidade_angular * raio_pneu / velocidade_linear))*100
            return value
       
    
    def printSlipRatio(self, rpm_values, slip_ratio, velocidade_angular, tire_longitudinal_forces):
        # Ajuste para garantir tamanhos compatíveis
        min_length = min(len(rpm_values), len(slip_ratio), len(velocidade_angular), len(tire_longitudinal_forces))
        rpm_values = np.array(rpm_values[:min_length])
        slip_ratio = np.array(slip_ratio[:min_length])
        velocidade_angular = np.array(velocidade_angular[:min_length])
        tire_longitudinal_forces = np.array(tire_longitudinal_forces[:min_length])

        passo = 10  # Intervalo para impressão
        print("Valores do Slip Ratio:")
        for i in range(0, len(slip_ratio), passo):
            print(f"{slip_ratio[i]:.2f}")

        # Figura com dois subgráficos
        plt.figure(figsize=(14, 6))

        # Gráfico 1: Slip Ratio x RPM
        plt.subplot(1, 2, 1)
        plt.plot(rpm_values, slip_ratio, label='Slip Ratio', color='blue')
        plt.xlabel('RPM')
        plt.ylabel('Slip Ratio (%)')
        plt.title('Slip Ratio (%) x RPM')
        plt.grid(True)
        plt.legend()

        # Gráfico 2: Força Longitudinal x Slip Ratio
        plt.subplot(1, 2, 2)
        plt.plot(self.tire_Ls[:min_length], tire_longitudinal_forces, label='Força Longitudinal', color='green')
        plt.xlabel('Slip Ratio (%)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal x Slip Ratio')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()

            
    
def Instancias():
    # Instância do modelo Drivetrain com os parâmetros fornecidos
    dt_model = Drivetrain(
        cgx=853,             # Centro de gravidade no eixo X [mm]
        cgy=294,             # Centro de gravidade no eixo Y [mm]
        massa=347,           # Massa do veículo [kg]
        entre_eixos=1567,    # Entre-eixos [mm]
        raio_pneu=259,       # Raio do pneu [mm]
        reducao_primaria=2.12,  # Redução primária
        reducao_final=2.76,     # Redução final (diferencial)
        cp=2.22,             # Constante de potência (pode ser usada futuramente)
        tempo_i=0            # Tempo inicial [s]
    )

    # Instância do motor (para plot)
    motor = Motor()
    motor.plot()

    # Cálculo da performance do veículo com os valores reduzidos
    performance_veiculo, variacao_rpm, variacao_tempo = dt_model.CarPerformance()

    # Impressão dos resultados
    dt_model.printCarPerformance(performance_veiculo)

    # Extração das velocidadesdos dados de performance
    velocidade_angular = np.array([dado["va"] for dado in performance_veiculo])
    velocidade_linear_carro = np.array([dado["vlc"] for dado in performance_veiculo])
  

    # Cálculo do slip ratio
    slip_ratio = Tire.SlipRatio(velocidade_angular, dt_model.raio_pneu * 0.001, velocidade_linear_carro)
    
    # Instância da classe Tire com os slip ratios calculados
    slip_model = Tire(
        tire_Fz=250,               # Carga vertical no pneu [N]
        tire_Sa=0,                 # Ângulo de escorregamento lateral [rad]
        tire_Ls=slip_ratio,       # Slip ratio já calculado
        tire_friction_coef=1.45,   # Coeficiente de fricção
        tire_Ca=0                  # Ângulo de camber
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

    # Plotagem dos gráficos de slip ratio
    slip_model.printSlipRatio(Motor.get_rpm(), slip_ratio, velocidade_angular, tire_longitudinal_forces)
    

Instancias()
