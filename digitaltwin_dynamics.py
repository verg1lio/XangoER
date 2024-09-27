import numpy as np
import matplotlib.pyplot as plt

class Tire:

    """
    Modelo de Pneu.

    Esta classe representa um modelo de pneu com vários parâmetros para calcular
    forças de pneus, cinemática e plotar gráficos relevantes. O modelo se baseia
    nas equações descritas em "Tire and Vehicle Dynamics" de Hans B. Pacejka, um 
    dos principais recursos sobre a dinâmica de pneus.
    """
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

    def calculate_kinematics(self):
        """
        Calcula a cinemática para um único valor de theta2_rad, considerando a posição das barras e os ângulos de deslizamento.

        Returns
        --------
        Ax_i : float
            Coordenada X do ponto A [m].
        Ay_i : float
            Coordenada Y do ponto A [m].
        Bx_i : float
            Coordenada X do ponto B [m].
        By_i : float
            Coordenada Y do ponto B [m].
        outer_slip : float
            Ângulo de deslizamento lateral externo [graus].

        Examples
        --------
        >>> tire = Tire(tire_Fz=3500, tire_Sa=0.05, tire_Ls=0.1, tire_friction_coef=0.9, B0=1.2, B1=0.8, B2=0.5, B3=0.9)
        >>> Ax, Ay, Bx, By, outer_slip = tire.calculate_kinematics()
        >>> print(f"Coordenadas de A: ({Ax}, {Ay}), Coordenadas de B: ({Bx}, {By}), Outer Slip Angle: {outer_slip}")

        Instanciando:        
        >>> Tire = Tire(params)
        >>> Tire.calculate_kinematics()
        >>> Coordenadas Calculadas:
        (Coordenadas de A: (0.6930124858148588, 0.7201704500352653), Coordenadas de B: (0.933739272127332, 1.1295564707481776), Outer Slip Angle: -81.61429844626039)

        References
        ----------
        title={FourBarLinkage-MATLAB},
        author={Hussein Shata},
        year={2020},
        url={https://www.youtube.com/watch?v=4O-XPJ7flLU},
        url={https://github.com/hshata/FourBarLinkage-MATLAB}
        """

        theta2_rad = self.angle

        # Cálculos intermediários
        AC_i = np.sqrt(self.L0**2 + self.L1**2 - 2 * self.L0 * self.L1 * np.cos(theta2_rad))
        
        # Verifique e ajuste os valores antes de passar para arccos
        beta_value = (self.L0**2 + AC_i**2 - self.L1**2) / (2 * self.L0 * AC_i)
        beta_value = np.clip(beta_value, -1, 1)
        beta_i = np.arccos(beta_value)

        psi_value = (self.L2**2 + AC_i**2 - self.L3**2) / (2 * self.L2 * AC_i)
        psi_value = np.clip(psi_value, -1, 1)
        psi_i = np.arccos(psi_value)

        lamda_value = (self.L3**2 + AC_i**2 - self.L2**2) / (2 * self.L3 * AC_i)
        lamda_value = np.clip(lamda_value, -1, 1)
        lamda_i = np.arccos(lamda_value)

        theta3_i = psi_i - beta_i
        theta4_i = np.pi - lamda_i - beta_i

        if np.degrees(theta2_rad) > 180:
            theta3_i = psi_i + beta_i
            theta4_i = np.pi - lamda_i + beta_i

        # Armazenamento dos resultados
        self.AC = AC_i
        self.beta = beta_i
        self.psi = psi_i
        self.lamda = lamda_i
        self.theta3 = theta3_i
        self.theta4 = theta4_i

        # Definição das posições das juntas
        Ax_i = self.L1 * np.cos(theta2_rad)
        Ay_i = self.L1 * np.sin(theta2_rad)
        Bx_i = Ax_i + self.L2 * np.cos(theta3_i)
        By_i = Ay_i + self.L2 * np.sin(theta3_i)

        # Verifica se L2 ou L3 é zero para evitar matriz singular
        if self.L2 == 0 or self.L3 == 0:
            raise ValueError("L2 and L3 must be non-zero and positive values.")
        
        r = np.array([
            [self.L2 * np.cos(theta3_i), -self.L3 * np.cos(theta4_i)],
            [self.L2 * np.sin(theta3_i), -self.L3 * np.sin(theta4_i)]
        ])
        v = np.array([-self.L1 * self.alpha * np.cos(theta2_rad), -self.L1 * self.alpha * np.sin(theta2_rad)])

        # Verifica se a matriz é singular e tenta ajustar se necessário
        try:
            w_i = np.linalg.solve(r, v)
        except np.linalg.LinAlgError:
            print("Matriz singular detectada. Ajustando matriz para resolução.")
            epsilon = 1e-6
            w_i = np.linalg.solve(r + np.eye(2) * epsilon, v)

        self.w = w_i
        self.om2 = w_i[0]
        self.om4 = w_i[1]

        # Cálculo dos ângulos outer_slip e inner_slip
        outer_slip_i = -(np.pi / 2 - theta4_i)
        inner_slip_i = -(np.pi / 2 - theta2_rad)

        self.outer_slip = np.degrees(outer_slip_i)
        self.inner_slip = np.degrees(inner_slip_i)

        # Verificação do caso estático
        if np.isclose(np.abs(outer_slip_i), np.abs(inner_slip_i), atol=1e-2):
            self.static_slip_angle = self.inner_slip

        return Ax_i, Ay_i, Bx_i, By_i, self.outer_slip

    def Tire_forces(self, params):
        """Calcular Forças do Pneu.

        Calcula as forças do pneu e o momento de auto-alinhamento usando os parâmetros do modelo Pacejka.

        Parâmetros
        ----------
        params : tuple
        >>> Os parâmetros para o modelo de Pacejka.

        Returns:
        -------
        tuple : (float, float, float)
        >>> Força lateral, momento de auto-alinhamento e força longitudinal do pneu.

        Examples:
        --------
        >>> tire = Tire(tire_Fz=5000, tire_Sa=0.1, tire_Ls=0.05, tire_friction_coef=1.0, tire_Ca=0.02)
        >>> tire.Tire_forces((1.5, 1.3, 1.1, 1.0, 2000, 3000))
        (400.0, 10.3, 150.0)
        
        References
        ----------
        Pacejka, H. B. (2002). Tire and Vehicle Dynamics. Elsevier.
        """

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

    def calcular_forca(self):
        """Calcular Força.

        Calcula a força com base na deformação do pneu.

        Returns:
        -------
        force : float
        >>> A força calculada com base em track_y e tire_k.

        Examples:
        --------
        >>> tire = Tire(track_y=10, tire_k=100)
        >>> tire.calcular_forca()
        1000
        """

        force = self.track_y * self.tire_k
        return force

    def slip_ratio_1(velocidade_angular, raio_pneu):
        """Calcular Razão de Escorregamento.

        Calcula a razão de escorregamento com base na velocidade angular e no raio do pneu.

        Parâmetros
        ----------
        velocidade_angular : float
        >>> A velocidade angular do pneu.
        raio_pneu : float
        >>> O raio do pneu.

        Returns:
        -------
        slip_ratio : float
        >>> A razão de escorregamento calculada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.slip_ratio_1(100, 0.3)
        0.33333333333333337
        
        References
        ----------
        Curso do Bob, Módulo 2, 2.6.
        """

        velocidade_linear = initial_speed
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1

        return value

    def plot_graph_slip_ratio(self, pedal_forces, value=None):
        """Mostrar Razão de Escorregamento.

        Exibe os valores da razão de escorregamento em relação ao RPM e à velocidade angular.

        Parâmetros
        ----------
        rpm_values : array-like
        >>> Os valores de RPM.
        slip_ratio : array-like
        >>> Os valores da razão de escorregamento.
        velocidade_angular : array-like
        >>> Os valores da velocidade angular.

        Returns:
        -------
        Nenhum : NoneType
        >>> O método plota o gráfico e imprime os valores da razão de escorregamento.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.show_slip_ratio([1000, 2000], [0.1, 0.2], [10, 20])
        """

        plt.figure(figsize=(10, 6))
        plt.plot(pedal_forces, value, label="Slip Ratio vs. Força do Pedal")

        # Invertendo o eixo y
        plt.gca().invert_yaxis()
        plt.xlabel("Força no Pedal (N)")
        plt.ylabel("Slip Ratio")
        plt.title("Força no Pedal em Relação ao Slip Ratio de Frenagem")
        plt.grid(True)
        plt.legend()
        plt.show()

    def show_slip_ratio(self, rpm_values, slip_ratio, velocidade_angular):
        """Mostrar Força Longitudinal vs. RPM.

        Plota a força longitudinal em relação ao RPM.

        Parâmetros
        ----------
        rpm_values : array-like
        >>> Os valores de RPM.
        tire_longitudinal_force : array-like
        >>> Os valores da força longitudinal.

        Returns:
        -------
        None
        >>> O método plota o gráfico e não retorna nada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.show_longitudinal_rpm([1000, 2000], [1500, 1600])
        """

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

    def show_longitudinal_rpm(rpm_values, tire_longitudinal_force):
        """Mostrar Força Longitudinal vs. RPM.

        Plota a força longitudinal em relação ao RPM.

        Parâmetros
        ----------
        rpm_values : array-like
        >>> Os valores de RPM.
        tire_longitudinal_force : array-like
        >>> Os valores da força longitudinal.

        Returns:
        -------
        None
        >>> O método plota o gráfico e não retorna nada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.show_longitudinal_rpm([1000, 2000], [1500, 1600])
        """

        plt.figure(figsize=(10, 6))
        plt.plot(rpm_values, tire_longitudinal_force, label='Longitudinal Force', color='blue')
        plt.xlabel('RPM')
        plt.ylabel('Longitudinal Force (N)')
        plt.title('Longitudinal Force vs. RPM')
        plt.legend()
        plt.grid(True)
        plt.show()

    def plotar_deformacao(self, track_y_values):
        """Plota a Deformação.

        Plota o gráfico com base em uma lista de valores fornecidos que representam a variação da pista,
        inicializa uma lista para receber os valores de força e adiciona para cada ponto fornecido.

        Parâmetros
        ----------
        track_y_values : array-like
        >>> Os valores que representam a deformação da pista.

        Returns:
        -------
        None
        >>> O método plota o gráfico e não retorna nada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.plotar_deformacao([0, 5, 10])
        """

        force_values = []

        for y in track_y_values:
            self.track_y = y
            force_values.append(self.calcular_forca())

        # Criando o gráfico
        plt.figure(figsize=(8, 6))
        plt.plot(track_y_values, force_values, color='b')
        plt.title('Carregamento Aplicado vs Deformação')
        plt.xlabel('Deformação (mm)')
        plt.ylabel('Carregamento (N)')
        plt.grid(True)
        plt.show()

    def plot_camber(self, predicted_tire_lateral_forces, predicted_tire_lateral_forces_1, predicted_tire_lateral_forces_2, tire_lateral_experimental=None, tire_lateral_experimental_1=None, tire_lateral_experimental_2=None, angles=None, ratio=None):
        """Plota o Câmber.

        Plota as curvas de camber e dados experimentais.

        Parâmetros
        ----------
        predicted_tire_lateral_forces : array-like
            As forças laterais previstas.
        predicted_tire_lateral_forces_1 : array-like
            As forças laterais previstas para a primeira curva.
        predicted_tire_lateral_forces_2 : array-like
            As forças laterais previstas para a segunda curva.
        tire_lateral_experimental : array-like, opcional
            As forças laterais experimentais (default é None).
        tire_lateral_experimental_1 : array-like, opcional
            As forças laterais experimentais para a primeira curva (default é None).
        tire_lateral_experimental_2 : array-like, opcional
            As forças laterais experimentais para a segunda curva (default é None).
        angles : array-like, opcional
            Os ângulos de deslizamento (default é None).
        ratio : float, opcional
            A proporção para escalonamento (default é None).

        Returns:
        -------
        None
        >>> Esta função não retorna valores, mas plota o gráfico com câmber.

        Examples:
        --------
        >>> tire = Tire(tire_Sa=np.linspace(-0.1, 0.1, 100))
        >>> tire.plot_camber([0, 5, 10], [1, 6, 11], [2, 7, 12])
        """

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
        """Plotar Gráfico de Forças do Pneu.

        Plota as forças do pneu e o momento de auto-alinhamento em relação aos ângulos de deslizamento.

        Parâmetros
        ----------
        predicted_tire_lateral_forces : array-like
            As forças laterais previstas.
        predicted_tire_auto_align_moment : array-like
            Os momentos de auto-alinhamento previstos.
        predicted_tire_longitudinal_forces : array-like
            As forças longitudinais previstas.
        tire_lateral_experimental : array-like, opcional
            As forças laterais experimentais (default é None).
        tire_auto_align_experimental : array-like, opcional
            Os momentos de auto-alinhamento experimentais (default é None).
        angles : array-like, opcional
            Os ângulos de deslizamento (default é None).
        ratio : float, opcional
            A proporção para escalonamento (default é None).

        Returns:
        -------
        None
        >>> Esta função não retorna valores, mas exibe um gráfico com dados experimentais.

        Examples:
        --------
        >>> tire = Tire(tire_Sa=np.linspace(-0.1, 0.1, 100), tire_Ls=np.linspace(0, 0.1, 100))
        >>> tire.plot_graph([0, 5, 10], [1, 6, 11], [2, 7, 12])
        """

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
        """Anima o mecanismo de quatro barras.

        Esta função desenha e anima o movimento do mecanismo de quatro barras,
        atualizando a posição dos pontos e as barras em cada iteração.

        Parâmetros
        ---------------------
        self.theta2 : list
            Lista de ângulos em graus para a entrada do mecanismo (barra de entrada).
        self.Ox, self.Oy : float
            Coordenadas do ponto fixo O no plano.
        self.Cx, self.Cy : float
            Coordenadas do ponto fixo C no plano.
        self.WB : float
            Distância entre o eixo dianteiro e traseiro do veículo (entre-eixos).
        self.rear_axle_length : float
            Comprimento do eixo traseiro do veículo.
        self.L3 : float
            Comprimento da barra L3 do mecanismo.
        self.theta4 : float
            Ângulo da barra de saída em radianos.
        self.outer_slip, self.inner_slip : float
            Ângulos de deslizamento calculados para as barras.

        Methods Utilizados
        ------------------
        self.calculate_kinematics() : method
        >>> Calcula e retorna as coordenadas dos pontos do mecanismo para o ângulo
        atual de entrada `self.angle`.

        Returns
        -------
        None
        >>> Esta função não retorna valores, mas exibe a animação do mecanismo.

        Examples
        --------
        >>> suspension_system.plot_mechanism()
        """

        self.Ax = []
        self.Ay = []
        self.Bx = []
        self.By = []
        self.Px = []

        for i in range(len(self.theta2)):
            self.angle = np.radians(self.theta2[i])
            Ax_i, Ay_i, Bx_i, By_i, outer_slip = self.calculate_kinematics()

            self.Ax.append(Ax_i)
            self.Ay.append(Ay_i)
            self.Bx.append(Bx_i)
            self.By.append(By_i)

            plt.figure(1)
            plt.plot([self.Ox, self.Ax[i]], [self.Oy, self.Ay[i]], 'r', linewidth=1)
            plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'k', linewidth=2)
            plt.plot([self.Bx[i], self.Cx], [self.By[i], self.Cy], 'r', linewidth=1)
            plt.plot([self.Ox, self.Cx], [self.Oy, self.Cy], 'r', linewidth=1, linestyle = 'dotted')

            # Desenho da barra do entre-eixos
            midpoint_x = (self.Ox + self.Cx) / 2
            rear_x1 = midpoint_x - (self.rear_axle_length / 2)
            rear_x2 = midpoint_x + (self.rear_axle_length / 2)
            

            plt.plot([midpoint_x, midpoint_x], [self.L3, (self.L3 - self.WB)], 'g', linewidth=1)  # Barra do entre-eixos
            plt.plot([rear_x1, rear_x2], [(self.L3 - self.WB), (self.L3 - self.WB)], 'g', linewidth=1)  # Barra do eixo traseiro

            # Adicionando as coordenadas dos pontos
            plt.text(self.Ox, self.Oy, f'({self.Ox:.2f}, {self.Oy:.2f})', fontsize=8, ha='right')
            plt.text(self.Ax[i], self.Ay[i], f'Manga de eixo ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Bx[i], self.By[i], f'Manga de eixo ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Cx, self.Cy, f'({self.Cx:.2f}, {self.Cy:.2f})', fontsize=8, ha='right')

            # Adicionando os ângulos de slip
            plt.text((self.Cx + self.Bx[i]) / 2, (self.Cy + self.By[i]) / 2, f'{self.outer_slip:.2f}°', fontsize=10, ha='center')
            plt.text((self.Ox + self.Ax[i]) / 2, (self.Oy + self.Ay[i]) / 2, f'{self.inner_slip:.2f}°', fontsize=10, ha='center')

            # Achando primeiro ponto de Ackerman

            plt.plot([(self.Ax[i] + self.Ox)/2, (self.Ax[i] + self.Ox)/2 - (self.L3/2 - self.WB)/np.sin(self.theta2_rad[i] + np.pi/2)], [(self.Ay[i] + self.Oy)/2, (self.L3 - self.WB)], 'g', linewidth=1, linestyle = 'dotted')  # Barra do entre-eixos
            plt.plot([(self.Bx[i] + self.Cx)/2, (self.Bx[i] + self.Cx)/2 - (self.L3/2 - self.WB)/np.sin(self.theta4 + np.pi/2)], [(self.By[i] + self.Cy)/2, (self.L3 - self.WB)], 'r', linewidth=1, linestyle = 'dotted')


            if i == 0 or i == (len(self.theta2) - 1):
                plt.show()

            plt.grid()
            plt.axis('equal')
            plt.axis([-300, 1800, -1000, 1000])
            plt.draw()
            plt.pause(0.01)
            plt.clf()

        for i in range(len(self.theta2)-1, -1, -1):

            self.angle = np.radians(self.theta2[i])
            Ax_i, Ay_i, Bx_i, By_i, outer_slip = self.calculate_kinematics()

            plt.plot([self.Ox, self.Ax[i]], [self.Oy, self.Ay[i]], 'r', linewidth=1)
            plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'k', linewidth=2)
            plt.plot([self.Bx[i], self.Cx], [self.By[i], self.Cy], 'r', linewidth=1)
            plt.plot([self.Ox, self.Cx], [self.Oy, self.Cy], 'r', linewidth=1, linestyle = 'dotted')

            # Desenho da barra do entre-eixos
            midpoint_x = (self.Ox + self.Cx) / 2
            rear_x1 = midpoint_x - (self.rear_axle_length / 2) * np.sin(self.alpha)
            rear_x2 = midpoint_x + (self.rear_axle_length / 2) * np.sin(self.alpha)
            

            plt.plot([midpoint_x, midpoint_x], [self.L3, (self.L3 - self.WB)], 'g', linewidth=1)  # Barra do entre-eixos
            plt.plot([rear_x1, rear_x2], [(self.L3 - self.WB), (self.L3 - self.WB)], 'g', linewidth=1)  # Barra do eixo traseiro

            # Adicionando as coordenadas dos pontos
            plt.text(self.Ox, self.Oy, f'({self.Ox:.2f}, {self.Oy:.2f})', fontsize=8, ha='right')
            plt.text(self.Ax[i], self.Ay[i], f'Manga de eixo ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='center')
            plt.text(self.Bx[i], self.By[i], f'Manga de eixo ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='center')

            # Adicionando os ângulos de slip
            plt.text((self.Cx + self.Bx[i]) / 2, (self.Cy + self.By[i]) / 2, f'{self.outer_slip:.2f}°', fontsize=10, ha='center')
            plt.text((self.Ox + self.Ax[i]) / 2, (self.Oy + self.Ay[i]) / 2, f'{self.inner_slip:.2f}°', fontsize=10, ha='center')

            # Achando primeiro ponto de Ackerman

            plt.plot([(self.Ax[i] + self.Ox)/2, (self.Ax[i] + self.Ox)/2 - (self.L3/2 - self.WB)/np.sin(self.theta2_rad[i] + np.pi/2)], [(self.Ay[i] + self.Oy)/2, (self.L3 - self.WB)], 'g', linewidth=1, linestyle = 'dotted')  # Barra do entre-eixos
            plt.plot([(self.Bx[i] + self.Cx)/2, (self.Bx[i] + self.Cx)/2 - (self.L3/2 - self.WB)/np.sin(self.theta4 + np.pi/2)], [(self.By[i] + self.Cy)/2, (self.L3 - self.WB)], 'r', linewidth=1, linestyle = 'dotted')

            plt.grid()
            plt.axis('equal')
            plt.axis([-300, 1800, -1000, 1100])
            plt.draw()
            plt.pause(0.01)
            plt.clf()

            if i == 0:
                print(f'O Slip Angle externo é : {self.outer_slip:.2f}°')

        if self.static_slip_angle is not None:
            print(f"O ângulo de toe é : {self.static_slip_angle:.2f}°")
        else:
            print("Não foi possível determinar um ângulo estático para as rodas dentro do intervalo fornecido.")

    def tire_example():
        """Executa a simulação do mecanismo e exibe os resultados.

        Esta função calcula a cinemática do mecanismo e chama a função
        `plot_mechanism` para animar o movimento. Além disso, exibe o ângulo
        de slip externo e o ângulo de toe estático, se disponível.

        Methods Utilizados
        ------------------
        self.calculate_kinematics() : method
        >>> Calcula e retorna as coordenadas dos pontos do mecanismo para o ângulo
        atual de entrada `self.angle`.

        self.plot_mechanism() : method
        >>> Anima o mecanismo de quatro barras.

        Returns
        -------
        None
        >>> Esta função não retorna valores, mas imprime os ângulos de slip e toe.

        Examples
        --------
        >>> suspension_system.run()
        """

        # Parâmetros de Pacejka
        result = [(0.3336564873588197), (1.6271741344929977), (10), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]

        # Criação de arrays com valores para plotagem e análise
        slip_ratio = np.linspace(-1, 1, 1000)  # Razão de escorregamento variando de -1 a 1
        slip_angles = np.linspace(-9, 9, 1000)  # Ângulos de deslizamento variando de -9 a 9 graus

        # Dados experimentais
        ratio = np.linspace(-1, 1, 19)
        angles = np.array([-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
        tire_lateral_forces_1 = np.array([-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890, 2000, 2090])
        tire_auto_align_moment_1 = np.array([-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09, 50.90, 50.10, 50.81, 48.12, 48.83])

        camber_experimental = [-2060.0, -1950.0, -1840.0, -1700.0, -1540.0, -1350.0, -1130.0, -860.0, -480.0, -30.0, 460.0, 880.0, 1230.0, 1490.0, 1720.0, 1910.0, 2090.0, 2230.0, 2310.0]
        camber_experimental_1 = [-1940.0, -1860.0, -1750.0, -1610.0, -1450.0, -1260.0, -1050.0, -760.0, -400.0, 60.0, 540.0, 970.0, 1290.0, 1550.0, 1780.0, 1980.0, 2150.0, 2280.0, 2370.0]
        camber_experimental_15 = [-1840.0, -1750.0, -1670.0, -1520.0, -1370.0, -1180.0, -960.0, -680.0, -310.0, 130.0, 610.0, 1020.0, 1360.0, 1630.0, 1850.0, 2040.0, 2220.0, 2360.0, 2430.0]

        #Instâncias iniciais
        dynamics_instance = Tire(tire_Fz=1500, tire_Sa=slip_angles, tire_Ls=slip_ratio, tire_friction_coef=1.45, tire_Ca=0)
        mecanism = Tire(B0 = 1393, B1 = 300, B2 = 1400, B3 = 300)

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
        dynamics_instance.plot_graph(predicted_tire_lateral_forces, predicted_tire_auto_align_moment, predicted_tire_longitudinal_forces,
                                    tire_lateral_experimental= tire_lateral_forces_1, tire_auto_align_experimental= tire_auto_align_moment_1, angles=angles, ratio=ratio)

        # Criando uma instância da classe e executando os métodos
        rigidez_deformacao = Tire(track_y=0, tire_k =120)
        rigidez_deformacao.calcular_forca()
        rigidez_deformacao.plotar_deformacao(track_y_values=np.linspace(0, 25, 1000))

        mecanism.plot_mechanism()

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

        Parâmetros
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

        Parâmetros
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

        # Cálculo do ângulo de câmber
        camber_angle = self.calcular_camber(Ax_i, Ay_i, Bx_i, By_i)

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

        Parâmetros
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

        Parâmetros
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

        Parâmetros
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

    def __init__(self, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, rpm, torque, cp):
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

    def CalculateOutputs(self):
        '''
        Método de cálculo de saídas da transmissão

        Return:
        ------------
        - Peso:                             Peso total do veículo.
        - Reação no eixo traseiro:           Força vertical no eixo traseiro.
        - Reação no eixo dianteiro:          Força vertical no eixo dianteiro.
        - Força trativa:                    Força exercida pelo veículo na superfície de contato.
        - Transferência de carga longitudinal: Quantidade de carga transferida longitudinalmente para o eixo traseiro.
        - Carga no eixo traseiro:            Carga total no eixo traseiro, incluindo a reação e a transferência de carga longitudinal.
        - Pico de torque no eixo traseiro:   Torque máximo no eixo traseiro.
        - Carga no pneu:                    Carga suportada por cada pneu do eixo traseiro.
        - Torque no pneu:                   Torque exercido em cada pneu do eixo traseiro.
        - Redução final:                    Relação de redução final da transmissão.
        - Torque necessário do motor:       Torque que o motor deve fornecer para atingir o pico de torque no eixo traseiro.
        - Aceleração primária real:         Aceleração real do veículo com base na força trativa e na massa.
        - Aceleração primária ideal:        Aceleração ideal do veículo.
        - Aceleração real final:            Aceleração real do veículo ajustada para o valor da gravidade.
        - Força trativa ideal:             Força trativa ideal baseada na aceleração ideal e na massa.
        - Torque no pneu ideal:            Torque ideal exercido em cada pneu com base na força trativa ideal.
        - Torque do motor ideal:           Torque ideal que o motor deve fornecer considerando as reduções da transmissão.
        - Transferência de carga ideal:     Transferência de carga longitudinal ideal considerando a força trativa ideal.
        - Transferência de carga real:      Transferência de carga longitudinal real baseada na força trativa real.
        - Carga traseira ideal:            Carga total no eixo traseiro ideal, incluindo a reação e a transferência de carga ideal.
        '''
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
        
        '''
        Método para exibir os resultados do cálculo de saídas da transmissão

        Este método utiliza os valores calculados pelo método `CalculateOutputs` para exibir de forma formatada os resultados relacionados à dinâmica e desempenho do veículo. Os resultados incluem informações sobre peso, forças, torques, acelerações, e cargas envolvidas na operação do veículo.

        Retorna:
        ------------
        - Peso:                             Peso total do veículo em Newtons (N).
        - Reação no eixo traseiro:           Força vertical no eixo traseiro em Newtons (N).
        - Reação no eixo dianteiro:          Força vertical no eixo dianteiro em Newtons (N).
        - Força trativa:                    Força exercida pelo veículo na superfície de contato em Newtons (N).
        - Transferência de carga longitudinal: Quantidade de carga transferida longitudinalmente para o eixo traseiro em Newtons (N).
        - Carga no eixo traseiro:            Carga total no eixo traseiro, incluindo a reação e a transferência de carga longitudinal, em Newtons (N).
        - Pico de torque no eixo traseiro:   Torque máximo no eixo traseiro em Newton-metros (Nm).
        - Carga no pneu:                    Carga suportada por cada pneu do eixo traseiro em Newtons (N).
        - Torque no pneu:                   Torque exercido em cada pneu do eixo traseiro em Newton-metros (Nm).
        - Redução final:                    Relação de redução final da transmissão.
        - Torque necessário do motor:       Torque que o motor deve fornecer para atingir o pico de torque no eixo traseiro em Newton-metros (Nm).
        - Aceleração primária real (g):     Aceleração real do veículo com base na força trativa e na massa, expressa como múltiplo da aceleração da gravidade (g).
        - Aceleração primária ideal:        Aceleração ideal do veículo em metros por segundo ao quadrado (m/s²).
        - Aceleração real final:            Aceleração real do veículo ajustada para o valor da gravidade, em metros por segundo ao quadrado (m/s²).
        - Força trativa ideal:             Força trativa ideal baseada na aceleração ideal e na massa em Newtons (N).
        - Torque no pneu ideal:            Torque ideal exercido em cada pneu com base na força trativa ideal em Newton-metros (Nm).
        - Torque do motor ideal:           Torque ideal que o motor deve fornecer considerando as reduções da transmissão em Newton-metros (Nm).
        - Transferência de carga ideal:     Transferência de carga longitudinal ideal considerando a força trativa ideal em Newtons (N).
        - Transferência de carga real:      Transferência de carga longitudinal real baseada na força trativa real em Newtons (N).
        - Carga traseira ideal:            Carga total no eixo traseiro ideal, incluindo a reação e a transferência de carga ideal, em Newtons (N).
        
        Example:
        ------
        dt_model = Drivetrain(cgx = 853, cgy = 294,  massa = 347, entre_eixos = 1567, coeficiente_atrito = 0.9, raio_pneu = 259, aceleracao_ideal = 1.2, reducao_primaria = 2.12, reducao_unica = 2.76, rpm = 2625, torque= 244.14, cp = 2.22)
        
        dt_model.show_results()        
        
        '''
        
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
        
        '''
        Método de cálculo do desempenho do veículo, o cálculo é potual, a menos que haja variação do rpm.

        Return(Lista):
        ------------
        - forca_trativa:                 Força trativa do veículo, calculada com base no torque do motor, nas reduções da transmissão, no raio do pneu e no rendimento da transmissão.

        - va:                           Velocidade angular das rodas em radianos por segundo, calculada a partir da rotação do motor (RPM) e das reduções da transmissão.

        - velocidade_linear:            Velocidade linear do veículo em km/h, derivada da velocidade angular das rodas e ajustada pelo fator de transmissão motor-roda.

        - fa:                           Força de arrasto do veículo, calculada com base na densidade do ar, na velocidade linear, no coeficiente de arrasto e na área frontal do veículo.

        - rr:                           Resistência de rolamento, calculada considerando a força de resistência básica, a força de resistência dependente da velocidade e o peso do veículo.

        - forca_final:                  Força final disponível para o veículo, resultante da diferença entre a força trativa e as resistências de arrasto e rolamento.
        
                Example:
        ------
        dt_model = Drivetrain(cgx = 853, cgy = 294,  massa = 347, entre_eixos = 1567, coeficiente_atrito = 0.9, raio_pneu = 259, aceleracao_ideal = 1.2, reducao_primaria = 2.12, reducao_unica = 2.76, rpm = 2625, torque= 244.14, cp = 2.22)
        
        performance = dt_model.CarPerformance()
        
        Lidando com variação de rpm:
        ------
        
        dt_model.new_rpm = 2625
        
        performance dt_model.CarPerformance()
        
        - Return: numpy.array()
        
        
        '''
        
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
        

        for rpm in variacao:
            # Cálculo da força trativa (N)
            forca_trativa = ((self.torque * self.reducao_primaria * self.reducao_unica * self.cp) / (self.raio_pneu * 0.001)) * rendimento_transmissao
        
            # Cálculo da velocidade angular (rad/s)
            velocidade_angular = (rpm * 2 * np.pi) / (60 * self.reducao_primaria * self.reducao_unica * self.cp)
            
            # Cálculo da velocidade linear (km/h)
            velocidade_linear = ((velocidade_angular * (self.raio_pneu * 0.001)) * transmissao_motor_roda) * 3.6

            # Cálculo da força de arrasto (N)
            fa = (densidade_ar * velocidade_linear ** 2 * coeficiente_arrasto * area_frontal) / 2

            # Cálculo da resistência de rolamento (N)
            rr = (basic_f + (3.24 * speed_f * ((velocidade_linear / 100 * 0.44704) ** 2.5))) * peso

            # Cálculo da força final (N)
            forca_final = forca_trativa - fa - rr

            # Armazenar os parâmetros calculados em um dicionário
            parametro = {"forca_trativa": forca_trativa,"va": velocidade_angular,"velocidade_linear": velocidade_linear,"fa": fa,"rr": rr,"forca_final": forca_final}

            parametros.append(parametro)
        
        # Retornar os parâmetros calculados
        return parametros, variacao
    
    def printCarPerformance(self):
            
            '''
            Método para exibir o desempenho do veículo

            Este método utiliza os valores calculados  pelo método `CarPerformance` para apresentar o desempenho do veículo. Os resultados incluem informações sobre força trativa, velocidade angular e linear, forças de arrasto, resistência ao rolamento e a força final do veículo. 

            Retorna:
            ------------
            - Força Trativa:                    Força exercida pelo veículo na superfície de contato, em Newtons (N).
            - Velocidade Angular:               Velocidade angular das rodas ou do sistema de transmissão, em radianos por segundo (rad/s).
            - Velocidade Linear:                Velocidade do veículo em quilômetros por hora (Km/h).
            - Força de Arrasto:                 Força que age contra o movimento do veículo devido à resistência do ar, em Newtons (N).
            - Resistência ao Rolamento:         Força que age contra o movimento do veículo devido ao atrito dos pneus com a superfície, em Newtons (N).
            - Força Final:                      Força resultante que o veículo exerce, considerando a força trativa, a força de arrasto e a resistência ao rolamento, em Newtons (N).
        
            Example:
            ------
            dt_model = Drivetrain(cgx = 853, cgy = 294,  massa = 347, entre_eixos = 1567, coeficiente_atrito = 0.9, raio_pneu = 259, aceleracao_ideal = 1.2, reducao_primaria = 2.12, reducao_unica = 2.76, rpm = 2625, torque= 244.14, cp = 2.22)
            
            dt_model.printCarPerformance()       
            '''        
            
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
        
        '''
        Método de dimensionamento dos semieixos

        Parametros:
        ------------
        fsi:                             Fator de segurança ideal aplicado ao cálculo do torque máximo de projeto (valor padrão: 1.25).
        tet:                             Tensão de escoamento do material do semieixo em MPa (valor padrão: 786).
        tec:                             Tensão de cisalhamento do material do semieixo em MPa (valor padrão: 471.6).
        dif:                             Fator de aumento para o torque máximo nos semieixos (valor padrão: 1).

        Return:
        ------------
        - Torque máximo do motor:           Torque máximo disponível do motor, utilizado como base para cálculos subsequentes.

        - Torque máximo nos semieixos:      Torque máximo transmitido aos semieixos considerando as reduções da transmissão e o fator de aumento.

        - Torque máximo de projeto:         Torque máximo que o semieixo deve suportar, ajustado pelo fator de segurança ideal.

        - Diâmetro dos semieixos:           Diâmetro necessário dos semieixos em milímetros, calculado com base no torque máximo de projeto e na tensão de cisalhamento.

        - Fator de segurança ideal:         Fator de segurança ideal utilizado no cálculo do torque máximo de projeto.

        - Fator de segurança obtido:        Fator de segurança real obtido para o diâmetro calculado dos semieixos, comparando com o torque máximo.

        - Fator de segurança para 1 polegada: Fator de segurança obtido para um semieixo com diâmetro de 1 polegada (25.4 mm), usado como referência para comparação.
        
        Example:
        ------
        dt_model = Drivetrain(cgx = 853, cgy = 294,  massa = 347, entre_eixos = 1567, coeficiente_atrito = 0.9, raio_pneu = 259, aceleracao_ideal = 1.2, reducao_primaria = 2.12, reducao_unica = 2.76, rpm = 2625, torque= 244.14, cp = 2.22)
        
        dt_model.HalfShaftsSizing()
        
        '''

        
        # Obtendo o maior torque do motor a partir dos dados experimentais 
        torque_max_motor = self.torque
        
        # Calculando o torque máximo nos semieixos
        torque_max_semieixo = torque_max_motor * self.reducao_primaria * self.reducao_unica * self.cp * dif
        
        # Calculando o torque máximo de projeto
        torque_max_projeto = torque_max_semieixo * fsi
        
        # Calculando o diâmetro dos semieixos (mm)
        diametro_semieixo = (((2 * torque_max_projeto) / (np.pi * tec * 10 ** 6)) ** (1 / 3)) * 2000
        
        # Calculando o fator de segurança obtido
        fator_seguranca_obtido = (np.pi * (((diametro_semieixo / 1000) / 2) ** 3) * tec * (10 ** 6)) / (2 * torque_max_semieixo)
        
        # Calculando o fator de segurança para 1 polegada
        fs1p = (np.pi * ((0.0254 / 2) ** 3) * tec * (10 ** 6)) / (2 * torque_max_semieixo)
        
        # Print dos resultados obtidos
        print("Dimensionamento dos Semieixos:")
        print("Torque máximo do motor:", torque_max_motor, "Nm")
        print("Torque máximo nos semieixos:", torque_max_semieixo, "Nm")
        print("Torque máximo de projeto:", torque_max_projeto, "Nm")
        print("Diâmetro dos semieixos:", diametro_semieixo, "mm")
        print("Fator de segurança ideal:", fsi)
        print("Fator de segurança obtido:", fator_seguranca_obtido)
        print("Fator de segurança para 1 polegada:", fs1p)

class BrakeSystem:
    '''
    Classe que modela o sistema de frenagem de um veículo e calcula diversos parâmetros relacionados à 
    frenagem, incluindo forças, torques, desaceleração, pressões necessárias e forças de atrito.
    '''
    # Inicializa a classe com os parâmetros necessários
    def __init__(self, params):
        self.params = params
        self.RedP = params['RedP']  # Redução do pedal
        self.a = params['a']  # Desaceleração [g]
        self.psi = params['psi']  # Distribuição de carga estática por eixo [adm]
        self.μl = params['μl']  # Coeficiente de atrito do contato pastilha/disco
        self.pi = params['pi']  # Valor de pi
        self.HCG = params['HCG']  # Altura do centro de gravidade [m]
        self.μ = params['μ']  # Coeficiente de atrito do pneu/solo
        self.FzF = params['FzF']  # Força de reação estática na dianteira [N]
        self.FzR = params['FzR']  # Força de reação estática na traseira [N]
        self.Rdp = params['Rdp']  # Raio do pneu [m]
        self.Dcm = params['Dcm']  # Diâmetro do cilindro mestre em metros
        self.Dwc = params['Dwc']  # Diâmetro do cilindro da roda em metros
        self.Npast = params['Npast']  # Número de pastilhas por disco
        self.atrito_coeficiente = params['atrito_coeficiente']  # Coeficiente de atrito da pastilha
        self.red = params['red']  # Raio efetivo do disco [m]
        self.Mt = params['Mt']  # Massa total do veículo [kg]
        self.L = params['L']  # Distância entre eixos [m]
        self.c_rr = params['c_rr']  # Coeficiente de resistência ao rolamento
        self.m_tire = params['m_tire']  # Massa do pneu
        self.m_wheel = params['m_wheel']  # Massa da roda
    
    def calculate_params(self, pedal_force=300):
        """
        Aplica o freio e calcula diversos parâmetros relacionados à frenagem com base na força aplicada ao pedal.

        Returns
        -------
        BF : float
            Fator de freio.
        χ : float
            Razão entre a altura do centro de gravidade (HCG) e a distância entre eixos (L).
        W : float
            Peso do veículo.
        FzF_dyn : float
            Força de reação dinâmica dianteira.
        FzR_dyn : float
            Força de reação dinâmica traseira.
        τF : float
            Torque na dianteira.
        τR : float
            Torque na traseira.
        FnF : float
            Força normal das pastilhas dianteira.
        FnR : float
            Força normal das pastilhas traseira.
        Awc : float
            Área do cilindro de roda.
        Acm : float
            Área do cilindro mestre.
        PF : float
            Pressão necessária na dianteira.
        PR : float
            Pressão necessária na traseira.
        FaCM : float
            Força necessária para acionar o cilindro mestre.
        lF : float
            Distância do eixo ao centro de gravidade em relação à distribuição de carga dianteira.
        lR : float
            Distância do eixo ao centro de gravidade em relação à distribuição de carga traseira.

        Examples
        --------
        >>> resultados = {
                'BF': BF, 'χ': χ, 'W': W, 'FzF_dyn': FzF_dyn, 'FzR_dyn': FzR_dyn,
                'τF': τF, 'τR': τR, 'FnF': FnF, 'FnR': FnR, 'Awc': Awc, 'Acm': Acm,
                'PF': PF, 'PR': PR, 'FaCM': FaCM, 'lF': lF, 'lR': lR
            }
        >>> print("Resultados Calculados:")
            for i, result in resultados.items():
                print(f"{i}: {result}")
                
        Instanciando:        
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.calculate_params()
        >>> Resultados Calculados:
        (BF: 0.9, χ: 0.3333333333333333, W: 2452.5, FzF_dyn: 1929.3, FzR_dyn: 523.1999999999999, τF: 868.185,
        τR: 235.43999999999997, FnF: 1302.2775, FnR: 353.15999999999997, Awc: 0.00080384, Acm: 0.00031400000000000004
        PF: 1620070.536425159, PR: 439341.16242038214, FaCM: 508.7021484375, lF: 0.72, lR: 0.78)

        References
        ----------
        title={Brake Design and Safety},
        author={Limpert, Rudolf},
        year={2011},
        publisher={SAE International},
        edition={3},
        url={https://www.sae.org/publications/books/content/r-397/}
        """
        BF = 2 * self.μl  # Fator de freio
        χ = self.HCG / self.L  # Razão entre HCG e L
        W = self.Mt * 9.81  # Peso do veículo
        FzF_dyn = (1 - self.psi + self.a * χ) * W  # Força de reação dinâmica dianteira
        FzR_dyn = (self.psi - self.a * χ) * W  # Força de reação dinâmica traseira
        τF = FzF_dyn * self.μ * self.Rdp  # Torque na dianteira
        τR = FzR_dyn * self.μ * self.Rdp  # Torque na traseira
        FnF = τF / self.Npast * self.RedP * self.red  # Força normal das pastilhas dianteira
        FnR = τR / self.Npast * self.RedP * self.red  # Força normal das pastilhas traseira
        Awc = (self.pi * (self.Dwc ** 2)) / 4  # Área do cilindro de roda
        Acm = (self.pi * (self.Dcm ** 2)) / 4  # Área do cilindro mestre
        PF = FnF / Awc  # Pressão necessária na dianteira
        PR = FnR / Awc  # Pressão necessária na traseira
        FaCM = PF * Acm  # Força necessária para acionar o cilindro mestre
        lF = self.psi * self.L  # Distância do eixo ao centro de gravidade em relação à distribuição de carga dianteira
        lR = (1 - self.psi) * self.L  # Distância do eixo ao centro de gravidade em relação à distribuição de carga traseira
        
        """# Exemplo:
        resultados = {
                'BF': BF, 'χ': χ, 'W': W, 'FzF_dyn': FzF_dyn, 'FzR_dyn': FzR_dyn,
                'τF': τF, 'τR': τR, 'FnF': FnF, 'FnR': FnR, 'Awc': Awc, 'Acm': Acm,
                'PF': PF, 'PR': PR, 'FaCM': FaCM, 'lF': lF, 'lR': lR
            }

        print("Resultados Calculados de calculate_params:")
        for i, result in resultados.items():
            print(f"{i}: {result}")"""

        return BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR

    def apply_brake(self, pedal_force=300):
        """
        Aplica o freio e calcula os resultados com base na força aplicada ao pedal.

        Returns
        -------
        forca_frenagem : float
            Força de frenagem considerando todos os fatores.

        torque_ajustado : float
            Torque de freio ajustado considerando a resistência ao rolamento.

        forca_f : float
            Força gerada pelo disco de freio ajustado.

        torque_disco_freio : float
            Torque do disco de freio.

        resistencia_rolamento : float
            Resistência ao rolamento.

        torque_resistencia_rolamento : float
            Torque gerado pela resistência ao rolamento.

        Examples
        --------
        >>> resultados = { 
                'resultados': resultados,'forca de frenagem': forca_frenagem, 'torque ajustado': torque_ajustado, 
                'forca gerada pelo disco': forca_f, 'torque do disco de freio': torque_disco_freio, 
                'resistencia ao rolamento': resistencia_rolamento, 'torque de resistencia ao rolamento': torque_resistencia_rolamento
            }
    
        >>> print("Resultados Calculados:")
            for i, result in resultados.items():
                print(f"{i}: {result}")

        Instanciando:
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.apply_brake()
        >>> (Forca de Frenagem: 687.46640625, Torque Ajustado: 190.56374999999997, 
            Forca Gerada Pelo Disco: 635.2125, Torque do Disco de Freio: 201.59999999999997, 
            Resistencia ao Rolamento: 36.7875, Torque de Resistencia ao Rolamento: 11.03625)

        References
        ----------
        title={Brake Design and Safety},
        author={Limpert, Rudolf},
        year={2011},
        publisher={SAE International},
        edition={3},
        url={https://www.sae.org/publications/books/content/r-397/}
        """
        resultados = self.calculate_params(pedal_force)
        BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR = resultados

        # Calcula a pressão do cilindro mestre
        pressao_cilindro_mestre = pedal_force / Acm

        # Calcula a pressurização do fluido
        pressao_fluido = pressao_cilindro_mestre

        # Calcula a transmissão de pressão da linha de freio
        transmissao_pressao = pressao_fluido * Awc

        # Calcula a força de aperto da pinça
        forca_aperto_pinca = transmissao_pressao

        # Calcula a força de atrito da pastilha de freio
        forca_atrito_pastilha = forca_aperto_pinca * self.atrito_coeficiente

        # Calcula o torque do disco de freio
        torque_disco_freio = forca_atrito_pastilha * self.red

        # Calcula a resistência ao rolamento
        resistencia_rolamento = self.c_rr * W

        # Calcula o torque gerado pela resistência ao rolamento
        torque_resistencia_rolamento = resistencia_rolamento * self.Rdp

        # Calcula o torque de freio ajustado considerando a resistência ao rolamento
        torque_ajustado = torque_disco_freio - torque_resistencia_rolamento

        # Calcula a força gerada pelo disco de freio ajustado
        forca_f = torque_ajustado / self.Rdp

        # Calcula a força de frenagem considerando todos os fatores
        forca_frenagem = (FnF + FnR) / self.Npast
        forca_frenagem *= self.atrito_coeficiente
        forca_frenagem *= self.red
        forca_frenagem /= self.Rdp
        forca_frenagem -= resistencia_rolamento

        """# Exemplo
        resultados = { 
                'Forca de Frenagem': forca_frenagem, 'Torque Ajustado': torque_ajustado, 
                'Forca Gerada Pelo Disco': forca_f, 'Torque do Disco de Freio': torque_disco_freio, 
                'Resistencia ao Rolamento': resistencia_rolamento, 'Torque de Resistencia ao Rolamento': torque_resistencia_rolamento
            }

        print("Resultados Calculados de apply_brake:")
        for i, result in resultados.items():
            print(f"{i}: {result}")"""

        return forca_frenagem, torque_ajustado, forca_f, torque_disco_freio, resistencia_rolamento, torque_resistencia_rolamento

    # Calcula a velocidade angular das rodas durante a frenagem
    def calculate_angular_velocity(self, torque_ajustado=190):
        ''' 
        Velocidade angular é uma medida da rapidez com que um objeto gira em torno de um eixo específico. 
        Ela é definida como a taxa de mudança do ângulo em relação ao tempo e é frequentemente medida 
        em radianos por segundo (rad/s). Nessa função, temos como principal objetivo determinar a velocidade 
        angular das rodas durante a frenagem do veículo.

        Returns
        -------
        angular_desaceleration: float
            A desaceleração angular calculada.

        angular_velocity: float
            A velocidade angular final após a iteração completa dos intervalos de tempo.

        angular_velocities: float
            A lista de velocidades angulares ao longo do tempo.

        inertia_wheel: float
            A inércia da roda, que é calculada inicialmente e pode ser útil para outros cálculos ou análises.

        Examples
        --------

        >>> resultados = { 
                'Desaceleração angular': angular_desaceleration,'Velocidade angular': angular_velocity,
                'Primeiro resultado da lista de Velocidades angulares': angular_velocities[0], 
                'Inércia do pneu': inertia_wheel 
            }
            
        >>> print("Resultados Calculados de calculate_angular_velocity:")
            for i, result in resultados.items():
                print(f"{i}: {result}")

        Instanciando:
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.calculate_angular_velocity()
        >>> (Desaceleração angular: 7.916666666666667, Velocidade angular: 50.67340067340082,
            Primeiro resultado da lista de Velocidades angulares: 66.506734006734, Inércia do pneu: 6.0 )

        References
        ----------
        title={Projeto e dimensionamento de um sistema de freios aplicado a um veículo fórmula SAE},
        author={Gustavo Carvalho Martins dos Santos},
        year={2014},
        publisher={Escola Politécnica da Universidade Federal do Rio de Janeiro - UFRJ},
        url={https://www.monografias.poli.ufrj.br/monografias/monopoli10011351.pdf}
        '''
        # Calcula o intervalo de tempo entre cada ponto de 'time_intervals'
        time_step = 0.001

        # Cria uma série de intervalos de tempo igualmente espaçados de 0 até 'tempo' com 100 pontos
        time_intervals = np.arange(0, tempo, time_step)
        
        # Calcula a inércia da roda considerando a massa do pneu e da roda e o raio do pneu
        inertia_wheel = 0.5 * (self.m_tire + self.m_wheel) * (self.Rdp * 2)
        
        # Calcula a desaceleração angular com base no torque ajustado e na inércia da roda
        angular_desaceleration = (torque_ajustado / 4) / inertia_wheel
        
        # Calcula a velocidade angular inicial das rodas com base na velocidade inicial e no raio do pneu
        initial_angular_velocity = initial_speed / self.Rdp
        
        # Define a velocidade angular inicial para o cálculo
        angular_velocity = initial_angular_velocity
        
        # Inicializa uma lista para armazenar as velocidades angulares ao longo do tempo
        angular_velocities = []

        # Itera sobre cada intervalo de tempo
        for i in range(len(time_intervals)):
            # Atualiza a velocidade angular subtraindo a desaceleração angular multiplicada pelo intervalo de tempo
            angular_velocity -= angular_desaceleration * time_step
            
            # Adiciona a velocidade angular atual à lista de velocidades angulares
            angular_velocities.append(angular_velocity)

            if angular_velocities[i] <= 0:
                angular_desaceleration = 0

        """
        # Exemplo
        resultados = { 
                'Desaceleração angular': angular_desaceleration,'Velocidade angular': angular_velocity,
                'Primeiro resultado da lista de Velocidades angulares': angular_velocities[0], 
                'Inércia do pneu': inertia_wheel 
            }

        print("Resultados Calculados de calculate_angular_velocity:")
        for i, result in resultados.items():
            print(f"{i}: {result}")"""

        return angular_desaceleration, angular_velocity, angular_velocities, inertia_wheel
    
    # Imprime os resultados calculados
    def print_resultados(self, pedal_force):
        '''
        Essa função nos mostra os valores de alguns parâmetros que foram caculados anteriormente. Isso serve 
        para monitorarmos os valores e o comportamento da Classe BrakeSystem.
    
        Examples
        --------
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.print_resultados()
        '''

        forca_frenagem, torque_ajustado, forca_f, torque_disco_freio, resistencia_rolamento, torque_resistencia_rolamento = self.apply_brake(pedal_force=pedal_force)
        angular_desaceleration, angular_velocity, angular_velocities, inertia_wheel = self.calculate_angular_velocity(torque_ajustado)

        print("Força de frenagem teórica:", forca_frenagem, "N")
        print("Torque ajustado:", torque_ajustado, "N.m")
        print("Forca Gerada Pelo Disco:", forca_f, "N")
        print("Torque do disco de freio:", torque_disco_freio, "N.m")
        print("Resistência ao rolamento:", resistencia_rolamento, "N")
        print("Torque de resistência ao rolamento:", torque_resistencia_rolamento, "N.m")
        print("Desaceleração Angular:", angular_desaceleration, "rad/s^2")
        print("Inércia da roda:", inertia_wheel, "kg.m^2")
        print("Velocidade Angular:", angular_velocity, "rad/s")

    # Plota o gráfico da força longitudinal em relação à força aplicada ao pedal
    def show_graph(self, pedal_force, time_intervals):
        '''
        Utilizamos a biblioteca matplotlib para plotar um gráfico da força longitudinal em relação à força 
        aplicada ao pedal. A função show_graph recebe dois parâmetros: longitudinal_force (que é uma lista ou 
        array contendo as forças longitudinais) e pedal_forces (que é uma lista ou array contendo as forças 
        aplicadas ao pedal). O gráfico gerado mostra como a força longitudinal varia conforme a força aplicada
        no pedal aumenta. 
        
        Examples
        --------
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.show_graph(forca_longitudinal, pedal_forces)
        '''
        # Definindo os parâmetros
        params = {
        'RedP': 4,  # Redução do pedal
        'a': 0.8,  # Desaceleração [g]
        'psi': 0.48,  # Distribuição de carga estática por eixo [adm]
        'μl': 0.45,  # Coeficiente de atrito do contato pastilha/disco
        'pi': 3.14,  # Valor de pi
        'HCG': 0.5,  # Altura do centro de gravidade [m]
        'μ': 1.5,  # Coeficiente de atrito do pneu/solo
        'FzF': 1471.5,  # Força de reação estática na dianteira [N]
        'FzR': 981.0,  # Força de reação estática na traseira [N]
        'Rdp': 0.30,  # Raio do pneu [m]
        'Dcm': 0.02,  # Diâmetro do cilindro mestre em metros
        'Dwc': 0.032,  # Diâmetro do cilindro da roda em metros
        'Npast': 2,  # Número de pastilhas por disco
        'atrito_coeficiente': 0.35,  # Coeficiente de atrito da pastilha
        'red': 0.75,  # Raio efetivo do disco [m]
        'Mt': 250,  # Massa total do veículo [kg]
        'm_wheel': 10,  # Massa da roda [kg] <-- estimativa
        'm_tire': 10,  # Massa do pneu [kg] <-- estimativa
        'L': 1.5,  # Distância entre eixos [m]]
        'c_rr': 0.015  # Coeficiente de resistência ao rolamento
        }

        # Parâmetros de Pacejka
        result = [(0.3336564873588197), (1.6271741344929977), (10), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]
        torque_ajustado = self.apply_brake(pedal_force = pedal_force)[1]
        angular_velocity = self.calculate_angular_velocity(torque_ajustado)[2]
        
        slip_ratio = []
        for velocity in angular_velocity:
            ratio = Tire.slip_ratio_1(velocidade_angular = velocity, raio_pneu = params['Rdp'])
            slip_ratio.append(ratio)

        #
        longitudinal_forces = []
        for slip in slip_ratio:
            brakes_instance = Tire(tire_Fz=1500, tire_Sa=0, tire_Ls=slip, tire_friction_coef=1.45, tire_Ca=0)
            longitudinal_force = brakes_instance.Tire_forces(result)[2]
            longitudinal_forces.append(longitudinal_force)

        plt.figure(figsize=(10, 6))
        plt.plot(time_intervals, angular_velocity, label='Longitudinal Force', color='blue')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Velocidade Angular (rad/s)')
        plt.title('Velocidade Angular vs. Tempo')
        plt.legend()
        plt.grid(True)
        plt.show()    

        plt.figure(figsize=(10, 6))
        plt.plot(time_intervals, slip_ratio, label='Longitudinal Force', color='blue')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Slip Ratio (Admensional)')
        plt.title('Slip Ratio vs. Tempo')
        plt.legend()
        plt.grid(True)
        plt.show()
        
        plt.figure(figsize=(10, 6))
        plt.plot(time_intervals, longitudinal_forces, label='Longitudinal Force', color='blue')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal vs. Tempo')
        plt.legend()
        plt.grid(True)
        plt.show()

    def brake_system_example():
        params = {
        'RedP': 4,  # Redução do pedal
        'a': 0.8,  # Desaceleração [g]
        'psi': 0.48,  # Distribuição de carga estática por eixo [adm]
        'μl': 0.45,  # Coeficiente de atrito do contato pastilha/disco
        'pi': 3.14,  # Valor de pi
        'HCG': 0.5,  # Altura do centro de gravidade [m]
        'μ': 1.5,  # Coeficiente de atrito do pneu/solo
        'FzF': 1471.5,  # Força de reação estática na dianteira [N]
        'FzR': 981.0,  # Força de reação estática na traseira [N]
        'Rdp': 0.30,  # Raio do pneu [m]
        'Dcm': 0.02,  # Diâmetro do cilindro mestre em metros
        'Dwc': 0.032,  # Diâmetro do cilindro da roda em metros
        'Npast': 2,  # Número de pastilhas por disco
        'atrito_coeficiente': 0.35,  # Coeficiente de atrito da pastilha
        'red': 0.75,  # Raio efetivo do disco [m]
        'Mt': 250,  # Massa total do veículo [kg]
        'm_wheel': 10,  # Massa da roda [kg] <-- estimativa
        'm_tire': 10,  # Massa do pneu [kg] <-- estimativa
        'L': 1.5,  # Distância entre eixos [m]]
        'c_rr': 0.015  # Coeficiente de resistência ao rolamento
        }

        time_step = 0.001
        time_intervals = np.arange(0, tempo, time_step)
        pedal_forces = 1000

        brake = BrakeSystem(params)
        brake.print_resultados(pedal_forces)
        brake.show_graph(pedal_force=pedal_forces, time_intervals=time_intervals)

def dynamics_example(pedal_forces, rpm, torque, slip_angle):
        
        params = {
        'RedP': 4,  # Redução do pedal
        'a': 0.8,  # Desaceleração [g]
        'psi': 0.48,  # Distribuição de carga estática por eixo [adm]
        'μl': 0.45,  # Coeficiente de atrito do contato pastilha/disco
        'pi': 3.14,  # Valor de pi
        'HCG': 0.5,  # Altura do centro de gravidade [m]
        'μ': 1.5,  # Coeficiente de atrito do pneu/solo
        'FzF': 1471.5,  # Força de reação estática na dianteira [N]
        'FzR': 981.0,  # Força de reação estática na traseira [N]
        'Rdp': 0.30,  # Raio do pneu [m]
        'Dcm': 0.02,  # Diâmetro do cilindro mestre em metros
        'Dwc': 0.032,  # Diâmetro do cilindro da roda em metros
        'Npast': 2,  # Número de pastilhas por disco
        'atrito_coeficiente': 0.35,  # Coeficiente de atrito da pastilha
        'red': 0.75,  # Raio efetivo do disco [m]
        'Mt': 250,  # Massa total do veículo [kg]
        'm_wheel': 10,  # Massa da roda [kg] <-- estimativa
        'm_tire': 10,  # Massa do pneu [kg] <-- estimativa
        'L': 1.5,  # Distância entre eixos [m]]
        'c_rr': 0.015  # Coeficiente de resistência ao rolamento
        }

        dt_model = Drivetrain(
        cgx = 853,  # mm
        cgy = 294,  # mm
        massa = 347,  # kg
        entre_eixos = 1567,  # mm
        coeficiente_atrito = 0.9 , # coeficiente de atrito
        raio_pneu = params['Rdp'],  # mm
        aceleracao_ideal = 1.2,  # g
        reducao_primaria = 2.12,  # redução primária
        reducao_unica = 2.76,
        rpm = rpm,
        torque= torque,
        cp = 2.22)

        dt_model.new_rpm = None

        result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]
        # Calculando forças do pneu a partir do input de freio e de direção
        brake = BrakeSystem(params)
        torque_ajustado = brake.apply_brake(pedal_force = pedal_forces)[1]
        angular_velocity = brake.calculate_angular_velocity(torque_ajustado)[1]
        slip_ratio_brake = Tire.slip_ratio_1(angular_velocity, params['Rdp'])

        Tire_instance_b = Tire(tire_Fz=1500, tire_Sa=slip_angle, tire_Ls=slip_ratio_brake, tire_friction_coef=1.45, tire_Ca=0)
        lateral_force, auto_align_moment, brake_longitudinal_forces = Tire_instance_b.Tire_forces(result)
        
        # Calculando forças do pneu a partir do input de drivetrain e de direção
        performance_veiculo, rpm_faixa = dt_model.CarPerformance()
        velocidade_angular_list = [dado['va'] for dado in performance_veiculo]
        for i in velocidade_angular_list:
            velocidade_angular = i
        slip_ratio_acel = Tire.slip_ratio_1(velocidade_angular, params['Rdp'])

        Tire_instance_d = Tire(tire_Fz=1500, tire_Sa=slip_angle, tire_Ls=slip_ratio_acel, tire_friction_coef=1.45, tire_Ca=0)
        lateral_force, auto_align_moment, drivetrain_longitudinal_forces = Tire_instance_d.Tire_forces(result)
        print(f'Resultados de Pneu com slip angle e slip ratio de aceleração e frenagem:\nForça Lateral: {lateral_force} N\nTorque auto-alinhante: {auto_align_moment}N.m\nForça de frenagem: {brake_longitudinal_forces} N\nForça de aceleração: {drivetrain_longitudinal_forces} N')


# Define a velocidade inicial e o tempo para o cálculo
initial_speed = 10  # Velocidade inicial do veículo [m/s]
tempo = 4  # Tempo para o cálculo [s]

dynamics_example(pedal_forces = 1000, rpm = 1000, torque = 100, slip_angle = 9)

BrakeSystem.brake_system_example()
Kinematics.kinematics_example()
Tire.tire_example()
