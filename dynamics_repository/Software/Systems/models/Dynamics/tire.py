import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import io, base64

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

        self.slip_angle_start = slip_angle_start
        self.slip_angle_end = slip_angle_end
        self.angle_step = angle_step
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

    def slip_ratio_1(velocidade_angular, raio_pneu, velocidade_linear):
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

        #velocidade_linear = initial_speed
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

    @staticmethod
    def show_slip_ratio(rpm_values, slip_ratio, velocidade_angular):
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

    @staticmethod
    def show_slip_ratio_base_64(rpm_values, slip_ratio, velocidade_angular):
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


        buffers = []
        imagens_base64 = []
        
        texto = "Valores do Slip Ratio: "

        for dado in slip_ratio:
            texto += f"""{dado}\n"""
            
        plt.figure(figsize=(15, 5))
        
        plt.subplot(1, 2, 1)
        plt.plot(rpm_values, slip_ratio, label = 'Slip Ratio', color = 'blue')
        plt.xlabel('RPM')
        plt.ylabel('Slip Ratio')
        plt.title('Slip Ratio x RPM') 
        plt.legend()
        buff1 = io.BytesIO()
        plt.savefig(buff1, format="png",bbox_inches="tight")
        buff1.seek(0)
        plt.close()    
        img_64 = base64.b64encode(buff1.read()).decode('utf-8')
        imagens_base64.append(img_64)
        
        plt.subplot(1, 2, 2)
        plt.plot(rpm_values, velocidade_angular, label = 'Velocidade Angular', color = 'red')
        plt.xlabel('RPM')
        plt.ylabel('Velocidade Angular (rad/s)')
        plt.title('Velocidade Angular x RPM')
        plt.legend()
        plt.tight_layout()
        buff2 = io.BytesIO()
        plt.savefig(buff2, format="png",bbox_inches="tight")
        buff2.seek(0)
        plt.close()    
        img_64 = base64.b64encode(buff2.read()).decode('utf-8')
        imagens_base64.append(img_64)

        return texto, imagens_base64

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

    def plotar_deformacao_64(self, track_y_values):
        """
        Plota o gráfico de carregamento vs deformação e retorna a imagem em Base64 para uso no Flet.
        """

        force_values = []

        for y in track_y_values:
            self.track_y = y
            force_values.append(self.calcular_forca())

        # Criando o gráfico
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(track_y_values, force_values, color='b')
        ax.set_title('Carregamento Aplicado vs Deformação')
        ax.set_xlabel('Deformação (mm)')
        ax.set_ylabel('Carregamento (N)')
        ax.grid(True)

        # Converter para Base64
        buff = io.BytesIO()
        plt.savefig(buff, format="png", bbox_inches="tight")
        plt.close(fig)
        buff.seek(0)
        img_b64 = base64.b64encode(buff.read()).decode("utf-8")

        return img_b64


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

    def plot_camber_base64(self, predicted_tire_lateral_forces, predicted_tire_lateral_forces_1, predicted_tire_lateral_forces_2, tire_lateral_experimental=None, tire_lateral_experimental_1=None, tire_lateral_experimental_2=None, angles=None, ratio=None):
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

        images_b64 = []

        # --- Gráfico 1: Força lateral ---
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(self.tire_Sa, predicted_tire_lateral_forces, label="Curva com Camber")
        if angles is not None and tire_lateral_experimental is not None:
            ax.scatter(angles, tire_lateral_experimental, color="red", label="Dados Experimentais")
        ax.set_xlabel("Ângulo de Deslizamento Lateral (graus)")
        ax.set_ylabel("Força Lateral do Pneu (N)")
        ax.set_title("Curva com Camber e Dados Experimentais")
        ax.legend()
        ax.grid(True)

        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        images_b64.append(base64.b64encode(buf.read()).decode("utf-8"))

        # --- Gráfico 2: Torque auto-alinhante ---
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(self.tire_Sa, predicted_tire_lateral_forces_1, label="Curva com Camber")
        if angles is not None and tire_lateral_experimental_1 is not None:
            ax.scatter(angles, tire_lateral_experimental_1, color="red", label="Dados Experimentais")
        ax.set_xlabel("Ângulo de Deslizamento Lateral (graus)")
        ax.set_ylabel("Força Lateral do Pneu (N)")
        ax.set_title("Curva com Camber e Dados Experimentais")
        ax.legend()
        ax.grid(True)

        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        images_b64.append(base64.b64encode(buf.read()).decode("utf-8"))

        # --- Gráfico 3: Força longitudinal ---
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(self.tire_Sa, predicted_tire_lateral_forces_2, label="Curva com Camber")
        if angles is not None and tire_lateral_experimental_2 is not None:
            ax.scatter(angles, tire_lateral_experimental_2, color="red", label="Dados Experimentais")
        ax.set_xlabel("Ângulo de Deslizamento Lateral (graus)")
        ax.set_ylabel("Força Lateral do Pneu (N)")
        ax.set_title("Curva com Camber e Dados Experimentais")
        ax.legend()
        ax.grid(True)

        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        images_b64.append(base64.b64encode(buf.read()).decode("utf-8"))

        return images_b64

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

    def plot_graph_base_64(self, predicted_tire_lateral_forces, predicted_tire_auto_align_moment, predicted_tire_longitudinal_forces, tire_lateral_experimental=None, tire_auto_align_experimental=None, angles=None, ratio=None):
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

        buffers = []
        imagens_base64 = []
        

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
        buff1 = io.BytesIO()
        plt.savefig(buff1, format="png",bbox_inches="tight")
        buff1.seek(0)
        img_64 = base64.b64encode(buff1.read()).decode('utf-8')
        imagens_base64.append(img_64)
        plt.close()

        

        # Subplot para torque auto-alinhante
        plt.subplot(1, 3, 2)
        plt.plot(self.tire_Sa, predicted_tire_auto_align_moment, label='Curva Otimizada')
        plt.scatter(angles, tire_auto_align_experimental, color='blue', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Torque auto-alinhante (N.m)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)
        buff2 = io.BytesIO()
        plt.savefig(buff2, format="png",bbox_inches="tight")
        buff2.seek(0)
        plt.close()    
        img_64 = base64.b64encode(buff2.read()).decode('utf-8')
        imagens_base64.append(img_64)

        # Subplot para força longitudinal do pneu
        plt.subplot(1, 3, 3)
        plt.plot(self.tire_Ls, predicted_tire_longitudinal_forces, label='Curva Sem Otimizar')
        plt.xlabel('Taxa de Escorregamento Longitudinal (Admensional)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal - Sem Dados Experimentais')
        plt.legend()
        plt.grid(True)
        buff3 = io.BytesIO()
        plt.savefig(buff3, format="png",bbox_inches="tight")
        buff3.seek(0)
        plt.close()    
        buffers.append(buff3)
        img_64 = base64.b64encode(buff3.read()).decode('utf-8')
        imagens_base64.append(img_64)


        return imagens_base64

        """# Ajusta a disposição dos subplots para evitar sobreposição
        plt.tight_layout(pad=3.0)
        plt.show()"""

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


    def plot_mechanism_base64(self):
        frames_base64 = []

        self.Ax, self.Ay, self.Bx, self.By, self.Px = [], [], [], [], []

        # Animação de ida
        for i in range(len(self.theta2)):
            self.angle = np.radians(self.theta2[i])
            Ax_i, Ay_i, Bx_i, By_i, outer_slip = self.calculate_kinematics()

            self.Ax.append(Ax_i)
            self.Ay.append(Ay_i)
            self.Bx.append(Bx_i)
            self.By.append(By_i)

            fig, ax = plt.subplots(figsize=(6, 4))
            
            # Barras principais
            ax.plot([self.Ox, self.Ax[i]], [self.Oy, self.Ay[i]], 'r', linewidth=1)
            ax.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'k', linewidth=2)
            ax.plot([self.Bx[i], self.Cx], [self.By[i], self.Cy], 'r', linewidth=1)
            ax.plot([self.Ox, self.Cx], [self.Oy, self.Cy], 'r', linewidth=1, linestyle='dotted')

            # Barra do entre-eixos
            midpoint_x = (self.Ox + self.Cx) / 2
            rear_x1 = midpoint_x - (self.rear_axle_length / 2)
            rear_x2 = midpoint_x + (self.rear_axle_length / 2)
            ax.plot([midpoint_x, midpoint_x], [self.L3, (self.L3 - self.WB)], 'g', linewidth=1)
            ax.plot([rear_x1, rear_x2], [(self.L3 - self.WB), (self.L3 - self.WB)], 'g', linewidth=1)

            # Textos originais
            ax.text(self.Ox, self.Oy, f'({self.Ox:.2f}, {self.Oy:.2f})', fontsize=8, ha='right')
            ax.text(self.Ax[i], self.Ay[i], f'Manga de eixo ({self.Ax[i]:.2f}, {self.Ay[i]:.2f})', fontsize=8, ha='center')
            ax.text(self.Bx[i], self.By[i], f'Manga de eixo ({self.Bx[i]:.2f}, {self.By[i]:.2f})', fontsize=8, ha='center')
            ax.text(self.Cx, self.Cy, f'({self.Cx:.2f}, {self.Cy:.2f})', fontsize=8, ha='right')

            # Slip angles
            ax.text((self.Cx + self.Bx[i]) / 2, (self.Cy + self.By[i]) / 2, f'{self.outer_slip:.2f}°', fontsize=10, ha='center')
            ax.text((self.Ox + self.Ax[i]) / 2, (self.Oy + self.Ay[i]) / 2, f'{self.inner_slip:.2f}°', fontsize=10, ha='center')

            # Pontos de Ackerman
            ax.plot(
                [(self.Ax[i] + self.Ox)/2, (self.Ax[i] + self.Ox)/2 - (self.L3/2 - self.WB)/np.sin(self.theta2_rad[i] + np.pi/2)],
                [(self.Ay[i] + self.Oy)/2, (self.L3 - self.WB)],
                'g', linewidth=1, linestyle='dotted'
            )
            ax.plot(
                [(self.Bx[i] + self.Cx)/2, (self.Bx[i] + self.Cx)/2 - (self.L3/2 - self.WB)/np.sin(self.theta4 + np.pi/2)],
                [(self.By[i] + self.Cy)/2, (self.L3 - self.WB)],
                'r', linewidth=1, linestyle='dotted'
            )

            # Configuração de plot
            ax.grid()
            ax.set_aspect('equal')
            ax.set_xlim(-300, 1800)
            ax.set_ylim(-1000, 1000)

            # Salvar para base64
            buff = io.BytesIO()
            plt.savefig(buff, format="png", bbox_inches="tight")
            plt.close(fig)
            buff.seek(0)
            img_b64 = base64.b64encode(buff.read()).decode("utf-8")
            frames_base64.append(img_b64)

        # Aqui você poderia repetir para o "caminho de volta" como no código original

        return frames_base64

        
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

