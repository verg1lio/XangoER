import numpy as np
import matplotlib.pyplot as plt

class Kinematics:
    """
    Classe que modela a cinemática de um mecanismo de quatro barras com inclusão de elementos como mola e 
    amortecedor, permitindo a análise do comportamento dinâmico ao longo do movimento.
    """
    def __init__(self, L0=None, L1=None, L2=0, L3=None, alpha=60, Py_start=100, Py_end=400, Py_step=5,
                spring_type=None, spring_k=None, spring_x=None, spring_y=None, spring_non_lin_coef=None, spring_angle=45,
                damper_type=None, damper_V=None, damper_F_viscous=None, damper_K_friction=None, damper_F_static=None):
        # Inicialização dos comprimentos das barras e outros parâmetros do mecanismo
        self.L0 = L0  # Comprimento da barra fixa
        self.L1 = L1  # Comprimento da barra de entrada
        self.L2 = L2  # Comprimento da barra de acoplamento
        self.L3 = L3  # Comprimento da barra de saída
        self.alpha = np.radians(alpha)  # Ângulo inicial da barra de entrada em radianos
        self.L_AP = 0.5 * L2  # Metade do comprimento da barra de acoplamento (ponto médio)
        self.Py_step = Py_step  # Passo do movimento vertical de P
        self.Py_start = Py_start 
        self.Py = np.arange(Py_start, Py_end, Py_step)  # Intervalo de posições verticais de P
        self.num_points = int((Py_end - Py_start) / Py_step)  # Número de pontos no intervalo

        # Inicialização das listas de resultados
        self.AC = []  # Lista para armazenar os comprimentos AC
        self.beta = []  # Lista para armazenar os ângulos beta
        self.psi = []  # Lista para armazenar os ângulos psi
        self.lamda = []  # Lista para armazenar os ângulos lambda
        self.theta2 = []  # Lista para armazenar os ângulos theta2 (entrada)
        self.theta3 = []  # Lista para armazenar os ângulos theta3 (acoplamento)
        self.theta4 = []  # Lista para armazenar os ângulos theta4 (saída)
        self.Ox, self.Oy = 0, 0  # Posições das âncoras fixas no chassi
        self.Ax, self.Ay = 0, 0  # Posições de A
        self.Bx, self.By = 0, 0  # Posições de B
        self.Cy, self.Cx = L0, 0  # Posições das âncoras fixas no chassi
        self.Px = 0  # Posições horizontais de P
        self.w = [0, 0]  # Velocidades angulares
        self.om2, self.om4 = 0, 0  # Velocidades angulares de theta2 e theta4
        self.alpha_dot = 0  # Aceleração angular
        self.V_Px, self.V_Py = 0, 0  # Velocidades de P
        self.V_Ax, self.V_Ay = 0, 0  # Velocidades de A
        self.V_A, self.P_A = 0, 0  # Magnitude das velocidades e posições de A
        

        # Inicialização dos parâmetros do modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # Rigidez da mola [N/m]
        self.spring_x = spring_x  # Força que a mola recebe [N]
        self.spring_y = spring_y  # Força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # Coeficiente de ganho não-linear
        self.spring_angle = np.radians(spring_angle)  # Ângulo da mola em relação à vertical

        # Inicialização dos parâmetros do modelo de amortecedor
        self.damper_type = damper_type  # Coulomb, Integrated, Stribeck
        self.damper_V = damper_V  # Velocidade relativa do amortecedor [m/s]
        self.damper_F_viscous = damper_F_viscous  # Força viscosa do fluído [N]
        self.damper_F_static = damper_F_static  # Coeficiente de fricção estática de Coulomb [N]
        self.damper_K_friction = damper_K_friction  # Rigidez de fricção [N/m]

    def calcular_theta2(self, Py_i):
        """Calcula o ângulo theta2 a partir da posição vertical de P (Py_i).

        Parameters
        ----------
        Py_i : float
        >>> Posição vertical do ponto P no instante i [mm].

        Returns
        -------
        theta2_rad_i : float
        >>> Ângulo theta2 em radianos.

        Examples
        --------
        >>> mechanism = Kinematics(L1=150, alpha=60)
        >>> theta2 = mechanism.calcular_theta2(Py_i=200)
        >>> print(f"Ângulo theta2: {theta2:.2f} rad")

        Resultado Calculado:
        Ângulo theta2: 1.23 rad

        Notes
        -----
        - A função utiliza a relação geométrica do triângulo formado entre as barras e a posição de P para calcular o ângulo.
        - É importante garantir que Py_i esteja dentro do alcance físico permitido pela configuração do mecanismo para evitar erros de cálculo.
        """

        theta2_rad_i = np.arccos((Py_i - self.L1 * np.cos(self.alpha)) / self.L1)
        return theta2_rad_i

    def calcular_camber(self, Ax, Ay, Bx, By):
        """Calcula o ângulo de câmber da barra de acoplamento em relação à vertical.

        Parameters
        ----------
        Ax : float
            Posição x do ponto A [mm].
        Ay : float
            Posição y do ponto A [mm].
        Bx : float
            Posição x do ponto B [mm].
        By : float
            Posição y do ponto B [mm].

        Returns
        -------
        angulo : float
            Ângulo de câmber em graus.

        Examples
        --------
        >>> mechanism = Kinematics()
        >>> camber_angle = mechanism.calcular_camber(Ax=0, Ay=100, Bx=50, By=150)
        >>> print(f"Ângulo de Câmber: {camber_angle:.2f} graus")

        Resultado Calculado:
        Ângulo de Câmber: 26.57 graus

        Notes
        -----
        - O ângulo de câmber representa a inclinação da barra de acoplamento em relação à vertical, o que é importante para a dinâmica e estabilidade do mecanismo.
        - A função usa a função `np.arctan2` para calcular o ângulo, o que garante a correção de sinal e o ajuste do quadrante.
        """

        delta_x = Bx - Ax
        delta_y = By - Ay
        angulo = np.degrees(np.arctan2(delta_x, delta_y))
        return angulo

    def calcular_cinematica(self):
        """Calcula as posições das juntas, ângulos e outras propriedades cinemáticas do mecanismo para um único ponto Py.

        Returns
        -------
        camber_angle : float
            Ângulo de câmber em graus.
        Ax_i : float
            Posição x da junta A [mm].
        Ay_i : float
            Posição y da junta A [mm].
        Bx_i : float
            Posição x da junta B [mm].
        By_i : float
            Posição y da junta B [mm].
        Px_i : float
            Posição x do ponto P [mm].

        Examples
        --------
        >>> mechanism = Kinematics(L0=100, L1=150, L2=200, L3=250, alpha=60, Py_start=100)
        >>> camber_angle, Ax, Ay, Bx, By, Px = mechanism.calcular_cinematica()
        >>> print(f"Ângulo de Câmber: {camber_angle:.2f} graus")
        >>> print(f"Posição de A: ({Ax:.2f}, {Ay:.2f})")
        >>> print(f"Posição de B: ({Bx:.2f}, {By:.2f})")
        >>> print(f"Posição de P: ({Px:.2f})")

        Resultado Calculado:
        Ângulo de Câmber: 26.57 graus
        Posição de A: (129.90, 75.00)
        Posição de B: (208.85, 220.33)
        Posição de P: (158.82)

        References
        ----------
        title={FourBarLinkage-MATLAB},
        author={Hussein Shata},
        year={2020},
        url={https://www.youtube.com/watch?v=4O-XPJ7flLU},
        url={https://github.com/hshata/FourBarLinkage-MATLAB}

        Notes
        -----
        - Esta função assume que as dimensões do mecanismo são especificadas de maneira que o cálculo de `theta2` seja válido.
        - Os cálculos são feitos para uma única posição de Py, definida por `Py_start`.
        - Certifique-se de que as variáveis de entrada, como os comprimentos das barras e o ângulo inicial, sejam consistentes com a configuração física do mecanismo.
        """

        Py = self.Py_start

        # Calcula theta2 para a posição atual de Py
        theta2_rad_i = self.calcular_theta2(Py)
        self.theta2 = np.degrees(theta2_rad_i)

        # Cálculos intermediários
        AC_i = np.sqrt(self.L0 ** 2 + self.L1 ** 2 - 2 * self.L0 * self.L1 * np.cos(theta2_rad_i))
        beta_i = np.arccos((self.L0 ** 2 + AC_i ** 2 - self.L1 ** 2) / (2 * self.L0 * AC_i))
        psi_i = np.arccos((self.L2 ** 2 + AC_i ** 2 - self.L3 ** 2) / (2 * self.L2 * AC_i))
        lamda_i = np.arccos((self.L3 ** 2 + AC_i ** 2 - self.L2 ** 2) / (2 * self.L3 * AC_i))

        # Calcula os ângulos theta3 e theta4
        theta3_i = psi_i - beta_i
        theta4_i = np.pi - lamda_i - beta_i

        # Ajuste de ângulos para diferentes quadrantes
        if self.theta2 > 180:
            theta3_i = psi_i + beta_i
            theta4_i = np.pi - lamda_i + beta_i

        # Armazenamento dos resultados intermediários
        self.AC = AC_i
        self.beta = beta_i
        self.psi = psi_i
        self.lamda = lamda_i
        self.theta3 = np.degrees(theta3_i)
        self.theta4 = np.degrees(theta4_i)

        # Definição das posições das juntas A e B
        Ax_i = self.L1 * np.sin(theta2_rad_i)
        Ay_i = self.L1 * np.cos(theta2_rad_i)
        Bx_i = Ax_i + self.L2 * np.sin(theta3_i)
        By_i = Ay_i + self.L2 * np.cos(theta3_i)
        Px_i = Ax_i + self.L_AP * np.sin(self.alpha + theta3_i)
        Py_i = Ay_i + self.L_AP * np.cos(self.alpha + theta3_i)

        # Sistema de equações para encontrar velocidades angulares
        r = np.array([
            [self.L2 * np.cos(theta3_i), -self.L3 * np.cos(theta4_i)],
            [self.L2 * np.sin(theta3_i), -self.L3 * np.sin(theta4_i)]
        ])
        v = np.array([-self.L1 * self.alpha * np.cos(theta2_rad_i), -self.L1 * self.alpha * np.sin(theta2_rad_i)])

        # Solução do sistema para encontrar as velocidades angulares
        w_i = np.linalg.solve(r, v)
        self.w = w_i
        self.om2 = w_i[0]
        self.om4 = w_i[1]

        # Cálculo da aceleração angular
        alpha_dot_i = -(self.L1 * w_i[0] * np.cos(theta2_rad_i) - self.L3 * w_i[1] * np.cos(theta4_i)) / \
                      (self.L2 * np.cos(self.alpha) + self.L_AP * np.cos(self.alpha))
        self.alpha_dot = alpha_dot_i

        # Cálculo das velocidades de A
        V_Ax_i = self.om2 * Ax_i
        V_Ay_i = self.om2 * Ay_i
        self.V_Ax = V_Ax_i
        self.V_Ay = V_Ay_i

        # Cálculo da magnitude das velocidades e posições de A
        V_A_i = np.sqrt(V_Ax_i ** 2 + V_Ay_i ** 2)
        self.V_A = V_A_i
        P_A_i = np.sqrt(Ax_i ** 2 + Ay_i ** 2)
        self.P_A = P_A_i

        # Cálculo das velocidades de P
        V_Px_i = self.om2 * Px_i
        V_Py_i = self.om2 * Py_i
        self.V_Px = V_Px_i
        self.V_Py = V_Py_i

        # Cálculo do ângulo de câmber
        camber_angle = self.calcular_camber(Ax_i, Ay_i, Bx_i, By_i)

        Ax_d = self.L1 - Ax_i 

        return camber_angle, Ax_i, Ay_i, Bx_i, By_i, Px_i

    def Spring(self):
        """Calcula a força da mola com base na deformação e no tipo de mola.

        Returns
        -------
        spring_F : float
            Força aplicada pela mola [N].

        Examples
        --------
        >>> mechanism = Kinematics(spring_type='Hooke', spring_k=1000, spring_x=0.02)
        >>> spring_force = mechanism.Spring()
        >>> print(f"Força da Mola: {spring_force} N")
        
        Instanciando:
        >>> Mechanism = Kinematics(params)
        >>> Mechanism.Spring()
        >>> Resultado Calculado:
        Força da Mola: 20.0 N

        References
        ----------
        url={https://www.sciencebuddies.org/science-fair-projects/references/linear-nonlinear-springs-tutorial}
        """

        if self.spring_type == 'Hooke':
            spring_F = self.spring_x * self.spring_k
        if self.spring_type == 'Softening':
            if self.spring_x < 0:
                spring_F = self.spring_k * (np.absolute(self.spring_x) ** self.spring_non_lin_coef)
                spring_F = -spring_F
            else:
                spring_F = self.spring_k * (self.spring_x ** self.spring_non_lin_coef)

        return spring_F

    def Damper(self):
        """Calcula a força do amortecedor com base no tipo de amortecedor.

        Returns
        -------
        damper_F : float
            Força aplicada pelo amortecedor [N].

        Examples
        --------
        >>> mechanism = Kinematics(damper_type='Coulomb', damper_F_static=50, damper_K_friction=0.3, damper_V=1.5)
        >>> damper_force = mechanism.Damper()
        >>> print(f"Força do Amortecedor: {damper_force} N")

        Instanciando:
        >>> Mechanism = Kinematics(params)
        >>> Mechanism.Damper()
        >>> Resultado Calculado:
        Força do Amortecedor: 14.83 N

        References
        ----------
        title={Modelling of Automotive Suspension Damper},
        author={VENKATA DINESH RAJU JONNALAGADDA; SAURABH VYAS},
        year={2020},
        pages={10-11},
        url={https://www.diva-portal.org/smash/get/diva2:1547559/FULLTEXT01.pdf}
        """

        if self.damper_type == 'Coulumb':
            damper_F = self.damper_F_static * np.tanh(self.damper_K_friction * self.damper_V)
        if self.damper_type == 'Integrated':
            damper_F = self.damper_F_static * np.tanh(self.damper_K_friction * self.damper_V) + np.sign(self.damper_V) * self.damper_F_viscous * (self.damper_V ** 2)
        return damper_F
    
    def plotar_cinematica(self):
        """Plota a cinemática do mecanismo.

        Esta função gera uma animação do movimento do mecanismo de quatro barras,
        mostrando as posições das juntas e a variação do ângulo de câmber ao longo
        do movimento. Ela também calcula e imprime os ângulos de câmber inicial,
        estático, final e a taxa de variação.

        Parameters
        ----------
        self : object
        >>> Referência à instância atual do objeto que contém os atributos e métodos
            necessários para o cálculo e plotagem da cinemática.

        Attributes Utilizados
        ---------------------
        self.Py : list
        >>> Lista de posições do ponto Py ao longo do movimento.
        self.Ox, self.Oy : float
        >>> Coordenadas do ponto de âncoragem O no chassis.
        self.Cx, self.Cy : float
        >>> Coordenadas do ponto de âncoragem C no chassis.
        self.num_points : int
        >>> Número total de pontos no movimento, utilizado para determinar quando
            mostrar os plots durante a animação.
        self.Py_step : float
        >>> Passo de movimento do ponto Py, utilizado para calcular a taxa de variação
            do ângulo de câmber.

        Returns
        -------
        None
        >>> Esta função não retorna valores, mas exibe plots e imprime informações
            sobre o ângulo de câmber.

        Examples
        --------
        >>> mechanism = Kinematics(L0=100, L1=150, L2=200, L3=250, alpha=60, Py_start=100)
        >>> mechanism.plotar_cinematica()
        
        Notas
        -----
        - Certifique-se de que todos os atributos necessários estão devidamente
        inicializados antes de chamar esta função.
        - O gráfico é atualizado em cada iteração para criar uma animação suave do
        movimento do mecanismo.
        - A função utiliza `plt.pause(0.01)` para criar a animação, o que pode ser
        ajustado para modificar a velocidade da animação.
        """

        self.Ax = []
        self.Ay = []
        self.Bx = []
        self.By = []
        self.Px = []

        for i in range(len(self.Py)):

                self.Py_start = self.Py[i]
                angulo_camber, Ax_i, Ay_i, Bx_i, By_i, Px_i = self.calcular_cinematica()
                self.Ax.append(Ax_i)
                self.Ay.append(Ay_i)
                self.Bx.append(Bx_i)
                self.By.append(By_i)
                self.Px.append(Px_i)
                
                angulo_camber = self.calcular_camber(self.Ax[i], self.Ay[i], self.Bx[i], self.By[i])
                
                plt.figure(1)
                plt.plot([self.Ox, self.Ax[i]], [self.Oy, self.Ay[i]], 'b', linewidth=2, label='Barra de Entrada')
                plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'b', linewidth=2, label='Barra de Acoplamento', linestyle='dotted')
                plt.plot([self.Bx[i], self.Cx], [self.By[i], self.Cy], 'b', linewidth=2, label='Barra de Saída')
                plt.plot([self.Ax[i], self.Px[i]], [self.Ay[i], self.Py[i]], 'r', linewidth=1, label='Barra Adicional')
                plt.plot([self.Bx[i], self.Px[i]], [self.By[i], self.Py[i]], 'g', linewidth=1, label='Barra Adicional')
                
                # Adicionando as coordenadas dos pontos
                plt.text(self.Ox, self.Oy, f'Âncoragem no chassis ({self.Ox:.2f}, {self.Oy:.2f})', fontsize=8, ha='right')
                plt.text(self.Ax[i], self.Ay[i], f'A ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='left')
                plt.text(self.Bx[i], self.By[i], f'B ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='left')
                plt.text(self.Cx, self.Cy, f'Âncoragem no chassis ({self.Cx:.2f}, {self.Cy:.2f})', fontsize=8, ha='right')
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
                
                plt.plot([self.Ox, self.Ax[i]], [self.Oy, self.Ay[i]], 'b', linewidth=2, label='Barra de Entrada')
                plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'b', linewidth=2, label='Barra de Acoplamento', linestyle='dotted')
                plt.plot([self.Bx[i], self.Cx], [self.By[i], self.Cy], 'b', linewidth=2, label='Barra de Saída')
                plt.plot([self.Ax[i], self.Px[i]], [self.Ay[i], self.Py[i]], 'r', linewidth=1, label='Barra Adicional')
                plt.plot([self.Bx[i], self.Px[i]], [self.By[i], self.Py[i]], 'g', linewidth=1, label='Barra Adicional')

                # Adicionando as coordenadas dos pontos
                plt.text(self.Ox, self.Oy, f'Âncoragem no chassis ({self.Ox:.2f}, {self.Oy:.2f})', fontsize=8, ha='right')
                plt.text(self.Ax[i], self.Ay[i], f'A ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='left')
                plt.text(self.Bx[i], self.By[i], f'B ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='left')
                plt.text(self.Cx, self.Cy, f'Âncoragem no chassis ({self.Cx:.2f}, {self.Cy:.2f})', fontsize=8, ha='right')
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
        """Plota a força do amortecedor em função da velocidade.

        Esta função calcula a força aplicada pelo amortecedor para uma série de
        valores de velocidade, tanto no eixo X quanto no eixo Y, e plota os
        resultados em gráficos separados.

        Parameters
        ----------
        damper_Vx_values : list
            Lista de valores de velocidade em mm/s para o eixo X.
        damper_Vy_values : list
            Lista de valores de velocidade em mm/s para o eixo Y.

        Attributes Utilizados
        ---------------------
        self.damper_V : float
            Velocidade atual do amortecedor, utilizada no cálculo da força.
        self.Damper() : method
            Método que calcula e retorna a força do amortecedor com base na
            velocidade atual `self.damper_V`.

        Returns
        -------
        None
            Esta função não retorna valores, mas exibe gráficos.

        Examples
        --------
        >>> damper_Vx = [0, 50, 100, 150, 200]
        >>> damper_Vy = [0, 50, 100, 150, 200]
        >>> suspension_system.plot_damper(damper_Vx, damper_Vy)
        """
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
        """Plota a força da mola em função da deformação.

        Esta função calcula a força aplicada pela mola para uma série de valores de
        deformação, tanto no eixo X quanto no eixo Y, e plota os resultados em
        gráficos separados.

        Parameters
        ----------
        spring_x_values : list
            Lista de valores de deformação em mm para o eixo X.
        spring_y_values : list
            Lista de valores de deformação em mm para o eixo Y.

        Attributes Utilizados
        ---------------------
        self.spring_x : float
            Deformação atual da mola, utilizada no cálculo da força.
        self.Spring() : method
            Método que calcula e retorna a força da mola com base na deformação
            atual `self.spring_x`.

        Returns
        -------
        None
            Esta função não retorna valores, mas exibe gráficos.

        Examples
        --------
        >>> spring_x = [0, 10, 20, 30, 40]
        >>> spring_y = [0, 10, 20, 30, 40]
        >>> suspension_system.plot_spring(spring_x, spring_y)
        """
        
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

    def kinematics_example():
        """Exemplo de uso da classe Kinematics.

        Esta função demonstra como criar instâncias da classe Kinematics para diferentes
        tipos de amortecedores e molas e como usar seus métodos de plotagem.

        Instâncias Criadas
        ------------------
        teste_integrated : Kinematics
        >>> Instância representando um amortecedor do tipo 'Integrated', que utiliza
            forças estáticas, de fricção e viscosas.

        teste_coulumb : Kinematics
        >>> Instância representando um amortecedor do tipo 'Coulumb', focando
            em forças de fricção.

        teste_hooke : Kinematics
        >>> Instância representando uma mola linear do tipo 'Hooke', que segue
            a lei de Hooke para molas.

        teste_softening : Kinematics
        >>> Instância representando uma mola do tipo 'Softening', que possui
            características de amolecimento não linear.

        kinematics : Kinematics
        >>> Instância para a análise cinemática de um mecanismo de quatro barras.

        Métodos Chamados
        ----------------
        kinematics.plotar_cinematica() : method
        >>> Plota a cinemática do mecanismo de quatro barras.

        teste_integrated.plot_damper() : method
        >>> Plota as características do amortecedor do tipo 'Integrated', variando as
            velocidades Vx e Vy.

        teste_coulumb.plot_damper() : method
        >>> Plota as características do amortecedor do tipo 'Coulumb', variando as
            velocidades Vx e Vy.

        teste_hooke.plot_spring() : method
        >>> Plota as características da mola do tipo 'Hooke', variando os deslocamentos
            x e y.

        teste_softening.plot_spring() : method
        >>> Plota as características da mola do tipo 'Softening', variando os deslocamentos
            x e y.

        Returns
        -------
        None
        >>> Esta função não retorna valores, mas gera plots das características dos
            amortecedores e molas, além da cinemática do mecanismo.

        Examples
        --------
        >>> kinematics_example()
        """
        
        teste_integrated = Kinematics(damper_type='Integrated', damper_F_static=50, damper_K_friction=1, damper_F_viscous=1)
        teste_coulumb = Kinematics(damper_type='Coulumb', damper_F_static=50, damper_K_friction=1)
        teste_hooke = Kinematics(spring_type= 'Hooke', spring_k=40, spring_x=0, spring_non_lin_coef=None)
        teste_softening = Kinematics(spring_type= 'Softening', spring_k=4, spring_x=0, spring_non_lin_coef=3)
        kinematics = Kinematics(L0=500, L1=500, L2=450, L3=500)


        kinematics.plotar_cinematica()
        teste_integrated.plot_damper(damper_Vx_values=np.linspace(-10, 10, 100), damper_Vy_values=np.linspace(-10, 10, 100))
        teste_coulumb.plot_damper(damper_Vx_values=np.linspace(-10, 10, 100), damper_Vy_values=np.linspace(-10, 10, 100))
        teste_hooke.plot_spring(spring_x_values=np.linspace(-10, 10, 100) , spring_y_values=np.linspace(-10, 10, 100))
        teste_softening.plot_spring(spring_x_values=np.linspace(-10, 10, 100) , spring_y_values=np.linspace(-10, 10, 100))
