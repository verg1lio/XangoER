import numpy as np
import matplotlib.pyplot as plt

class Kinematics:
    def __init__(self, L0=None, L1=None, L2=None, L3=None, alpha=np.radians(45), theta2_start=75, theta2_end=105, theta2_step=0.6):
        self.L0 = L0
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.alpha = alpha
        self.L_AP = L2 / np.sqrt(2)
        self.theta2 = np.arange(theta2_start, theta2_end, theta2_step)
        self.theta2_rad = np.radians(self.theta2)
        
        # Inicialização das listas de resultados
        self.AC = []
        self.beta = []
        self.psi = []
        self.lamda = []
        self.theta3 = []
        self.theta4 = []
        self.Ox, self.Oy = [0] * len(self.theta2), [0] * len(self.theta2)
        self.Ax, self.Ay = [], []
        self.Bx, self.By = [], []
        self.Cy, self.Cx = [L0] * len(self.theta2), [0] * len(self.theta2)
        self.Px, self.Py = [], []
        self.w = []
        self.om2, self.om4 = [], []
        self.alpha_dot = []
        self.V_Px, self.V_Py = [], []

    def calcular_camber(self, Ax, Ay, Bx, By):
        """Calcula o ângulo que a barra de acoplamento faz com a vertical."""
        delta_x = Bx - Ax
        delta_y = By - Ay
        angulo = np.degrees(np.arctan2(delta_x, delta_y))
        return angulo

    def calcular_cinematica(self):
        for i in range(len(self.theta2)):
            # Cálculos intermediários
            AC_i = np.sqrt(self.L0**2 + self.L1**2 - 2*self.L0*self.L1*np.cos(self.theta2_rad[i]))
            beta_i = np.arccos((self.L0**2 + AC_i**2 - self.L1**2) / (2*self.L0*AC_i))
            psi_i = np.arccos((self.L2**2 + AC_i**2 - self.L3**2) / (2*self.L2*AC_i))
            lamda_i = np.arccos((self.L3**2 + AC_i**2 - self.L2**2) / (2*self.L3*AC_i))

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
            Ax_i = self.L1 * np.sin(self.theta2_rad[i])
            Ay_i = self.L1 * np.cos(self.theta2_rad[i])
            Bx_i = Ax_i + self.L2 * np.sin(theta3_i)
            By_i = Ay_i + self.L2 * np.cos(theta3_i)
            Px_i = Ax_i + self.L_AP * np.sin(self.alpha + theta3_i)
            Py_i = Ay_i + self.L_AP * np.cos(self.alpha + theta3_i)

            self.Ax.append(Ax_i)
            self.Ay.append(Ay_i)
            self.Bx.append(Bx_i)
            self.By.append(By_i)
            self.Px.append(Px_i)
            self.Py.append(Py_i)

            r = np.array([
                [self.L2 * np.cos(theta3_i), -self.L3 * np.cos(theta4_i)],
                [self.L2 * np.sin(theta3_i), -self.L3 * np.sin(theta4_i)]
            ])
            v = np.array([-self.L1 * self.alpha * np.cos(self.theta2_rad[i]), -self.L1 * self.alpha * np.sin(self.theta2_rad[i])])

            w_i = np.linalg.solve(r, v)
            self.w.append(w_i)
            self.om2.append(w_i[0])
            self.om4.append(w_i[1])

            alpha_dot_i = - (self.L1 * w_i[0] * np.cos(self.theta2_rad[i]) - self.L3 * w_i[1] * np.cos(theta4_i)) / \
                          (self.L2 * np.cos(self.alpha) + self.L_AP * np.cos(self.alpha))
            self.alpha_dot.append(alpha_dot_i)

            V_Px_i = alpha_dot_i * Px_i
            V_Py_i = alpha_dot_i * Py_i

            self.V_Px.append(V_Px_i)
            self.V_Py.append(V_Py_i)
    
    def plotar_cinematica(self):
        for i in range(len(self.theta2)):
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

            if i == 0 or i == 24 or i == 49:
                plt.show()

            plt.grid()
            plt.axis('equal')
            plt.axis([-500, 1000, -350, 1100])
            plt.draw()
            plt.pause(0.01)
            plt.clf()

        for i in range(len(self.theta2)-1, -1, -1):
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
        camber_estatico = self.calcular_camber(self.Ax[24], self.Ay[24], self.Bx[24], self.By[24])
        camber_final = self.calcular_camber(self.Ax[49], self.Ay[49], self.Bx[49], self.By[49])
        taxa_variacao = (camber_final - camber_inicial) / len(self.theta2)
        
        print(f"Ângulo de câmber bump: {camber_inicial:.2f}°")
        print(f"Ângulo de câmber estático: {camber_estatico:.2f}°")
        print(f"Ângulo de câmber rebound: {camber_final:.2f}°")
        print(f"Taxa de variação de câmber: {taxa_variacao:.2f}°/passo")

# Instanciando a classe e executando os cálculos e plots
kinematics = Kinematics(L0=500, L1=500, L2=450, L3=500)
kinematics.calcular_cinematica()
kinematics.plotar_cinematica()
