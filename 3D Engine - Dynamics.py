import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Kinematics:
    def __init__(self, L0=None, L1=None, L2=None, L3=None, alpha=np.radians(45), theta2_start=75, theta2_end=105, theta2_step=0.6, space=200, mirror=230):
        self.L0 = L0
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.alpha = alpha
        self.L_AP = L2 / np.sqrt(2)
        self.theta2 = np.arange(theta2_start, theta2_end, theta2_step)
        self.theta2_rad = np.radians(self.theta2)
        self.space = space # Espaçamento entre os braços
        self.mirror = mirror # Distância do espelhamento (ou largura do chassis)

        # Inicialização das listas de resultados
        self.AC = []
        self.beta = []
        self.psi = []
        self.lamda = []
        self.theta3 = []
        self.theta4 = []
        self.Ox, self.Oz = [0] * len(self.theta2), [0] * len(self.theta2)
        self.Ax, self.Az = [], []
        self.Bx, self.Bz = [], []
        self.Cz, self.Cx = [L0] * len(self.theta2), [0] * len(self.theta2)
        self.Px, self.Pz = [], []
        self.w = []
        self.om2, self.om4 = [], []
        self.alpha_dot = []
        self.V_Px, self.V_Pz = [], []

    def calcular_camber(self, Ax, Az, Bx, Bz):
        """Calcula o ângulo que a barra de acoplamento faz com a vertical."""
        delta_x = Bx - Ax
        delta_z = Bz - Az
        angulo = np.degrees(np.arctan2(delta_x, delta_z))
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
            Az_i = self.L1 * np.cos(self.theta2_rad[i])
            Bx_i = Ax_i + self.L2 * np.sin(theta3_i)
            Bz_i = Az_i + self.L2 * np.cos(theta3_i)
            Px_i = Ax_i + self.L_AP * np.sin(self.alpha + theta3_i)
            Pz_i = Az_i + self.L_AP * np.cos(self.alpha + theta3_i)

            self.Ax.append(Ax_i)
            self.Az.append(Az_i)
            self.Bx.append(Bx_i)
            self.Bz.append(Bz_i)
            self.Px.append(Px_i)
            self.Pz.append(Pz_i)

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
            V_Pz_i = alpha_dot_i * Pz_i

            self.V_Px.append(V_Px_i)
            self.V_Pz.append(V_Pz_i)

    def plotar_cinematica(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        self.Oy = [0] * len(self.theta2)
        self.Ay = [0] * len(self.theta2)
        self.By = [0] * len(self.theta2)
        self.Cy = [0] * len(self.theta2)
        self.Py = [0] * len(self.theta2)
        self.Ay.append(0)
        self.By.append(0)
        self.Py.append(0)


        for i in range(len(self.theta2)):
                ax.cla()  # Limpa o gráfico a cada frame

                angulo_camber = self.calcular_camber(self.Ax[i], self.Az[i], self.Bx[i], self.Bz[i])

                # Desenha as barras
                ax.plot([self.Ox[i], self.Ax[i]], [self.Oy[i] - self.space/2, self.Ay[i]], [self.Oz[i], self.Az[i]], 'b', linewidth=2)  # Entrada
                ax.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], [self.Az[i], self.Bz[i]], 'b--', linewidth=2)  # Acoplamento
                ax.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i] - self.space/2], [self.Bz[i], self.Cz[i]], 'b', linewidth=2)  # Saída
                ax.plot([self.Ax[i], self.Px[i]], [self.Ay[i], self.Py[i]], [self.Az[i], self.Pz[i]], 'r', linewidth=1)  # Braço adicional
                ax.plot([self.Bx[i], self.Px[i]], [self.By[i], self.Py[i]], [self.Bz[i], self.Pz[i]], 'g', linewidth=1)  # Braço adicional
                ax.plot([self.Ox[i], self.Ax[i]], [self.Oy[i] + self.space/2 , self.Ay[i] ], [self.Oz[i], self.Az[i]], 'b', linewidth=2) # Entrada com espaçamento
                ax.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i] + self.space/2] , [self.Bz[i], self.Cz[i]], 'b', linewidth=2) # Saída com espaçamento

                #### ESPELHAMENTO

                ax.plot([(self.Ox[i] * -1) - self.mirror, (self.Ax[i]* -1) - self.mirror], [self.Oy[i] - self.space/2, self.Ay[i]], [self.Oz[i], self.Az[i]], 'b', linewidth=2) # Entrada espelhada
                ax.plot([(self.Ax[i] * -1) - self.mirror, (self.Bx[i]* -1) - self.mirror], [self.Ay[i], self.By[i]], [self.Az[i], self.Bz[i]], 'b--', linewidth=2) # Acoplamento espelhado
                ax.plot([(self.Bx[i] * -1) - self.mirror, (self.Cx[i]* -1) - self.mirror], [self.By[i], self.Cy[i] - self.space/2], [self.Bz[i], self.Cz[i]], 'b', linewidth=2) # Saída espelhada
                ax.plot([(self.Ax[i] * -1) - self.mirror, (self.Px[i]* -1) - self.mirror], [self.Ay[i], self.Py[i]], [self.Az[i], self.Pz[i]], 'r', linewidth=1) # Braço adicional espelhada
                ax.plot([(self.Bx[i]* -1) - self.mirror, (self.Px[i]* -1) - self.mirror], [self.By[i], self.Py[i]], [self.Bz[i], self.Pz[i]], 'g', linewidth=1) # Braço adicional espelhada
                ax.plot([(self.Ox[i]* -1) - self.mirror, (self.Ax[i]* -1) - self.mirror], [self.Oy[i] + self.space/2 , self.Ay[i] ], [self.Oz[i], self.Az[i]], 'b', linewidth=2) # Entrada com espaçamento espelhada
                ax.plot([(self.Bx[i]* -1) - self.mirror, (self.Cx[i]* -1) - self.mirror], [self.By[i], self.Cy[i] + self.space/2] , [self.Bz[i], self.Cz[i]], 'b', linewidth=2) # Saída com espaçamento espelhada




                # Texto do ângulo de câmber em 3D
                meio_x = (self.Ax[i] + self.Bx[i]) / 2
                meio_z = (self.Az[i] + self.Bz[i]) / 2
                meio_y = (self.Ay[i] + self.By[i]) / 2
                ax.text(meio_x, meio_y, meio_z + 10, f'{angulo_camber:.2f}°', fontsize=10)

                # Texto do ângulo de câmber em 3D espelhada
                meio_x = (self.Ax[i]*-1 + self.Bx[i]*-1) / 2
                meio_z = (self.Az[i] + self.Bz[i]) / 2
                meio_y = (self.Ay[i] + self.By[i]) / 2
                ax.text(meio_x, meio_y, meio_z + 10, f'{angulo_camber:.2f}°', fontsize=10)




                ax.set_xlim(-1000, 1000)
                ax.set_ylim(-400, 400)
                ax.set_zlim(-1000, 1100)
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Z')

                plt.pause(0.01)

        plt.show()

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')

        for i in range(len(self.theta2)-1, -1, -1):
            ax2.cla()

            angulo_camber = self.calcular_camber(self.Ax[i], self.Az[i], self.Bx[i], self.Bz[i])

            ax2.plot([self.Ox[i], self.Ax[i]], [self.Oy[i] - self.space/2, self.Ay[i]], [self.Oz[i], self.Az[i]], 'b', linewidth=2)  # Entrada
            ax2.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], [self.Az[i], self.Bz[i]], 'b--', linewidth=2)  # Acoplamento
            ax2.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i] - self.space/2], [self.Bz[i], self.Cz[i]], 'b', linewidth=2)  # Saída
            ax2.plot([self.Ax[i], self.Px[i]], [self.Ay[i], self.Py[i]], [self.Az[i], self.Pz[i]], 'r', linewidth=1)  # Braço adicional
            ax2.plot([self.Bx[i], self.Px[i]], [self.By[i], self.Py[i]], [self.Bz[i], self.Pz[i]], 'g', linewidth=1)  # Braço adicional
            ax2.plot([self.Ox[i], self.Ax[i]], [self.Oy[i] + self.space/2 , self.Ay[i] ], [self.Oz[i], self.Az[i]], 'b', linewidth=2) # Entrada com espaçamento
            ax2.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i] + self.space/2] , [self.Bz[i], self.Cz[i]], 'b', linewidth=2) # Saída com espaçamento

            #### ESPELHAMENTO

            ax2.plot([(self.Ox[i] * -1) - self.mirror, (self.Ax[i]* -1) - self.mirror], [self.Oy[i] - self.space/2, self.Ay[i]], [self.Oz[i], self.Az[i]], 'b', linewidth=2) # Entrada espelhada
            ax2.plot([(self.Ax[i] * -1) - self.mirror, (self.Bx[i]* -1) - self.mirror], [self.Ay[i], self.By[i]], [self.Az[i], self.Bz[i]], 'b--', linewidth=2) # Acoplamento espelhado
            ax2.plot([(self.Bx[i] * -1) - self.mirror, (self.Cx[i]* -1) - self.mirror], [self.By[i], self.Cy[i] - self.space/2], [self.Bz[i], self.Cz[i]], 'b', linewidth=2) # Saída espelhada
            ax2.plot([(self.Ax[i] * -1) - self.mirror, (self.Px[i]* -1) - self.mirror], [self.Ay[i], self.Py[i]], [self.Az[i], self.Pz[i]], 'r', linewidth=1) # Braço adicional espelhada
            ax2.plot([(self.Bx[i]* -1) - self.mirror, (self.Px[i]* -1) - self.mirror], [self.By[i], self.Py[i]], [self.Bz[i], self.Pz[i]], 'g', linewidth=1) # Braço adicional espelhada
            ax2.plot([(self.Ox[i]* -1) - self.mirror, (self.Ax[i]* -1) - self.mirror], [self.Oy[i] + self.space/2 , self.Ay[i] ], [self.Oz[i], self.Az[i]], 'b', linewidth=2) # Entrada com espaçamento espelhada
            ax2.plot([(self.Bx[i]* -1) - self.mirror, (self.Cx[i]* -1) - self.mirror], [self.By[i], self.Cy[i] + self.space/2] , [self.Bz[i], self.Cz[i]], 'b', linewidth=2) # Saída com espaçamento espelhada




            # Texto do ângulo de câmber em 3D
            meio_x = (self.Ax[i] + self.Bx[i]) / 2
            meio_z = (self.Az[i] + self.Bz[i]) / 2
            meio_y = (self.Ay[i] + self.By[i]) / 2
            ax2.text(meio_x, meio_y, meio_z + 10, f'{angulo_camber:.2f}°', fontsize=10)

            # Texto do ângulo de câmber em 3D espelhada
            meio_x = (self.Ax[i]*-1 + self.Bx[i]*-1) / 2
            meio_z = (self.Az[i] + self.Bz[i]) / 2
            meio_y = (self.Ay[i] + self.By[i]) / 2
            ax2.text(meio_x, meio_y, meio_z + 10, f'{angulo_camber:.2f}°', fontsize=10)



            ax2.set_xlim(-1000, 1000)
            ax2.set_ylim(-400, 400)
            ax2.set_zlim(-1000, 1100)
            ax2.set_xlabel('X')
            ax2.set_ylabel('Y')
            ax2.set_zlabel('Z')
            plt.pause(0.01)

        # plt.show()


        # Imprimindo os ângulos de câmber inicial, estático, final e a taxa de variação
        camber_inicial = self.calcular_camber(self.Ax[0], self.Az[0], self.Bx[0], self.Bz[0])
        camber_estatico = self.calcular_camber(self.Ax[24], self.Az[24], self.Bx[24], self.Bz[24])
        camber_final = self.calcular_camber(self.Ax[49], self.Az[49], self.Bx[49], self.Bz[49])
        taxa_variacao = (camber_final - camber_inicial) / len(self.theta2)

        print(f"Ângulo de câmber bump: {camber_inicial:.2f}°")
        print(f"Ângulo de câmber estático: {camber_estatico:.2f}°")
        print(f"Ângulo de câmber rebound: {camber_final:.2f}°")
        print(f"Taxa de variação de câmber: {taxa_variacao:.2f}°/passo")

# Instanciando a classe e executando os cálculos e plots
kinematics = Kinematics(L0=500, L1=500, L2=450, L3=500)
kinematics.calcular_cinematica()
kinematics.plotar_cinematica()
