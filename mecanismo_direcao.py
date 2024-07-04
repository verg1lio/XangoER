import numpy as np
import matplotlib.pyplot as plt

class FourBarMechanism:
    def __init__(self, L0=None, L1=None, L2=None, L3=None, alpha=None, slip_angle_start=-9, slip_angle_end=9, theta2_step=0.5):
        self.L0 = L0
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.alpha = alpha
        self.L_AP = L2 / np.sqrt(2)
        self.theta2 = np.arange(slip_angle_start + 90 , slip_angle_end + 91, theta2_step)
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
        self.Cx, self.Cy = [L0] * len(self.theta2), [0] * len(self.theta2)
        self.Px, self.Py = [], []
        self.w = []
        self.om2, self.om4 = [], []
        self.alpha_dot = []
        self.V_Px, self.V_Py = [], []
        self.outer_slip = []
        self.inner_slip = []

        self.static_slip_angle = None

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
            Px_i = Ax_i + self.L_AP * np.cos(self.alpha + theta3_i)
            Py_i = Ay_i + self.L_AP * np.sin(self.alpha + theta3_i)

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

            # Cálculo dos ângulos outer_slip e inner_slip
            outer_slip_i = -(np.pi / 2 - theta4_i)
            inner_slip_i = -(np.pi / 2 - self.theta2_rad[i])

            self.outer_slip.append(np.degrees(outer_slip_i))
            self.inner_slip.append(np.degrees(inner_slip_i))

            # Verificação do caso estático
            if np.isclose(np.abs(outer_slip_i), np.abs(inner_slip_i), atol=1e-2):
                self.static_slip_angle = self.inner_slip[i]

    def plot_mechanism(self):
        for i in range(len(self.theta2)):
            plt.figure(1)
            plt.plot([self.Ox[i], self.Ax[i]], [self.Oy[i], self.Ay[i]], 'r', linewidth=1)
            plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'k', linewidth=2)
            plt.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i]], 'r', linewidth=1)
            plt.plot([self.Ax[i], self.Px[i]], [self.Ay[i], self.Py[i]], 'r', linewidth=1, linestyle='dotted')
            plt.plot([self.Bx[i], self.Px[i]], [self.By[i], self.Py[i]], 'b', linewidth=1, linestyle='dotted')

            # Adicionando as coordenadas dos pontos
            plt.text(self.Ox[i], self.Oy[i], f'({self.Ox[i]:.2f}, {self.Oy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Ax[i], self.Ay[i], f'Manga de eixo ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Bx[i], self.By[i], f'Manga de eixo ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Cx[i], self.Cy[i], f'({self.Cx[i]:.2f}, {self.Cy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Px[i], self.Py[i], f'Ponto de Ackerman ({self.Px[i]:.2f}, {self.Py[i]:.2f})', fontsize=8, ha='right')

            # Adicionando os ângulos de slip
            plt.text((self.Cx[i] + self.Bx[i]) / 2, (self.Cy[i] + self.By[i]) / 2, f'{self.outer_slip[i]:.2f}°', fontsize=10, ha='center')
            plt.text((self.Ox[i] + self.Ax[i]) / 2, (self.Oy[i] + self.Ay[i]) / 2, f'{self.inner_slip[i]:.2f}°', fontsize=10, ha='center')

            if i == 0 or i == 18 or i == 36:
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
            plt.plot([self.Ax[i], self.Px[i]], [self.Ay[i], self.Py[i]], 'r', linewidth=1, linestyle='dotted')
            plt.plot([self.Bx[i], self.Px[i]], [self.By[i], self.Py[i]], 'b', linewidth=1, linestyle='dotted')

            # Adicionando as coordenadas dos pontos
            plt.text(self.Ox[i], self.Oy[i], f'({self.Ox[i]:.2f}, {self.Oy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Ax[i], self.Ay[i], f'Manga de eixo ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Bx[i], self.By[i], f'Manga de eixo ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Cx[i], self.Cy[i], f'({self.Cx[i]:.2f}, {self.Cy[i]:.2f})', fontsize=8, ha='right')
            plt.text(self.Px[i], self.Py[i], f'Ponto de Ackerman ({self.Px[i]:.2f}, {self.Py[i]:.2f})', fontsize=8, ha='right')

            # Adicionando os ângulos de slip
            plt.text((self.Cx[i] + self.Bx[i]) / 2, (self.Cy[i] + self.By[i]) / 2, f'{self.outer_slip[i]:.2f}°', fontsize=10, ha='center')
            plt.text((self.Ox[i] + self.Ax[i]) / 2, (self.Oy[i] + self.Ay[i]) / 2, f'{self.inner_slip[i]:.2f}°', fontsize=10, ha='center')

            plt.grid()
            plt.axis('equal')
            plt.axis([-300, 1800, -1000, 1000])
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
mechanism = FourBarMechanism(L0 = 1375, L1 = 300, L2 = 1400, L3 = 300, alpha = np.radians(315))

# Cálculos
mechanism.calculate_kinematics()

# Plots
mechanism.run()

