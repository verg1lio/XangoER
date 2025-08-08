# ============================================================================
# Classe Tire: modelo simplificado de pneu
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt

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

            value = ((velocidade_angular * raio_pneu / velocidade_linear)- 1) 
            return value
       
    def printSlipRatio(self, tempo, slip_ratio, tire_longitudinal_forces):

        passo = 10  # Intervalo para impressão
        print("Valores do Slip Ratio:")
        for i in range(0, len(slip_ratio), passo):
            print(f"{slip_ratio[i]:.2f}")

        # Figura com dois subgráficos
        plt.figure(figsize=(14, 6))

        # Gráfico 1: Slip Ratio x RPM
        plt.subplot(1, 2, 1)
        plt.plot(tempo, slip_ratio, label='Slip Ratio', color='blue')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Slip Ratio')
        plt.title('Slip Ratio x Tempo')
        plt.grid(True)
        plt.legend()

        # Gráfico 2: Força Longitudinal x Slip Ratio
        plt.subplot(1, 2, 2)
        plt.plot(slip_ratio, tire_longitudinal_forces, label='Força Longitudinal', color='green')
        plt.xlabel('Slip Ratio')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal x Slip Ratio')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()
