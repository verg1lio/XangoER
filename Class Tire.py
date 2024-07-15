import numpy as np
import matplotlib.pyplot as plt

# Definindo uma classe chamada Tire
class Tire:
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=None, B0=0, B1=None, 
                B2=0, B3=None, omega=315, slip_angle_start=-10, slip_angle_end=10, angle_step=0.5, WB=1500, rear_axle_length=1600):
        # Inicializando parâmetros para o modelo de pneu
        self.tire_Fz = tire_Fz  # Carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # Ângulo de deslizamento lateral do pneu [rad]
        self.tire_Ls = tire_Ls  # Escorregamento longitudinal do pneu [Adimensional]
        self.tire_type = 'Default'  # Tipo de pneu
        self.tire_friction_coef = tire_friction_coef  # Coeficiente de fricção entre pneu e pista
        self.tire_Ca = tire_Ca  # Ângulo de camber do pneu

        # Comprimentos das barras do mecanismo de quatro barras
        self.L0 = B0  # Comprimento da primeira barra
        self.L1 = B1  # Comprimento do braço de direção
        self.L2 = B2  # Comprimento da barra de direção
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

# Criação do mecanismo
mechanism = Tire(B0 = 1393, B1 = 300, B2 = 1400, B3 = 300)
# Cálculos
mechanism.calculate_kinematics()
# Plots
mechanism.run()

# Valores para a instância da classe Tire
slip_ratio = np.linspace(-1, 1, 1000)
slip_angles = np.linspace(-9, 9, 1000)

# Dados experimentais
ratio = np.linspace(-1, 1, 19)
angles = np.array([-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
tire_lateral_forces_1 = np.array([-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890, 2000, 2090])
tire_auto_align_moment_1 = np.array([-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09, 50.90, 50.10, 50.81, 48.12, 48.83])

camber_experimental = [-2060.0, -1950.0, -1840.0, -1700.0, -1540.0, -1350.0, -1130.0, -860.0, -480.0, -30.0, 460.0, 880.0, 1230.0, 1490.0, 1720.0, 1910.0, 2090.0, 2230.0, 2310.0]
camber_experimental_1 = [-1940.0, -1860.0, -1750.0, -1610.0, -1450.0, -1260.0, -1050.0, -760.0, -400.0, 60.0, 540.0, 970.0, 1290.0, 1550.0, 1780.0, 1980.0, 2150.0, 2280.0, 2370.0]
camber_experimental_15 = [-1840.0, -1750.0, -1670.0, -1520.0, -1370.0, -1180.0, -960.0, -680.0, -310.0, 130.0, 610.0, 1020.0, 1360.0, 1630.0, 1850.0, 2040.0, 2220.0, 2360.0, 2430.0]

# Parâmetros de Pacejka Atrelados ao pneu
result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]

# Instanciando a classe Tire sem camber
dynamics_instance = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45, tire_Ca=0)

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
dynamics_instance.plot_graph(predicted_tire_lateral_forces, predicted_tire_auto_align_moment, predicted_tire_longitudinal_forces)