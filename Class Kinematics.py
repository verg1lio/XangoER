import numpy as np
import matplotlib.pyplot as plt

class Kinematics:
    def __init__(self, L0=None, L1=None, L2=0, L3=None, alpha=60, Py_start=100, Py_end=400, Py_step=5,
                spring_type=None, spring_k=None, spring_F=None, spring_non_lin_coef=None,
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

        # Inicialização dos parâmetros do modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # rigidez da mola [N/m]
        self.spring_F = spring_F  # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # coeficiente de ganho não-linear

        # Inicialização dos parâmetros do modelo de amortecedor
        self.damper_type = damper_type # Coulumb, Integrated, Stribeck
        self.damper_V = damper_V # velocidade relativa amortecedor [m/s]
        self.damper_F_viscous = damper_F_viscous # força viscosa do fluído [N]
        self.damper_F_static = damper_F_static # coeficiente de fricção estática de coulumb [N]
        self.damper_K_friction = damper_K_friction # rigidez de fricção [N/m]

    def Spring(self):
        """Calcula a deformação da mola com base no tipo de mola."""
        if self.spring_type == 'Hooke':
            spring_x = self.spring_F / self.spring_k
        if self.spring_type == 'Softening':
            spring_x = self.spring_F / (self.spring_non_lin_coef * (self.spring_k) ** 2)
        return spring_x

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

            # Cálculo das velocidades de P
            V_Px_i = alpha_dot_i * Px_i
            V_Py_i = alpha_dot_i * Py_i

            self.V_Px.append(V_Px_i)
            self.V_Py.append(V_Py_i)

        # Calculando ângulo de câmber
        angulo_camber = self.calcular_camber(self.Ax[int(self.num_points - 1)], self.Ay[int(self.num_points - 1)], self.Bx[int(self.num_points - 1)], self.By[int(self.num_points - 1)])
        return angulo_camber
    
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

    def plot_damper(self):
        """Plota a força do amortecedor em função da velocidade."""
        damper_V_values = np.linspace(-5, 5, 1000)
        self.damper_V = damper_V_values
        damper_F_values = self.Damper()

        plt.figure(figsize=(8, 6))
        plt.plot(damper_V_values, damper_F_values, label='Velocidade')
        plt.title('Força em função da velocidade')
        plt.xlabel('Velocidade [m/s]')
        plt.ylabel('Força [N]')
        plt.grid(True)
        plt.legend()
        plt.show()

# Instância para o amortecedor do tipo 'Integrated'
teste_integrated = Kinematics(damper_type='Integrated', damper_F_static=50, damper_K_friction=10, damper_F_viscous=10)
teste_integrated.plot_damper()

# Instância para o amortecedor do tipo 'Coulumb'
teste_coulumb = Kinematics(damper_type='Coulumb', damper_F_static=50, damper_K_friction=10)
teste_coulumb.plot_damper()

# Criação do objeto Kinematics e execução dos cálculos
kinematics = Kinematics(L0=500, L1=500, L2=450, L3=500)
kinematics.calcular_cinematica()
kinematics.plotar_cinematica()
