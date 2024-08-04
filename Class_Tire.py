import numpy as np
import matplotlib.pyplot as plt

class Tire:
    """
    Modelo de Pneu.

    Esta classe representa um modelo de pneu com vários parâmetros para calcular
    forças de pneus, cinemática e plotar gráficos relevantes. O modelo se baseia
    nas equações descritas em "Tire and Vehicle Dynamics" de Hans B. Pacejka, um 
    dos principais recursos sobre a dinâmica de pneus.

    Returns:
    -------
    Tire : objeto
        Uma instância da classe Tire.

    Examples:
    --------
    >>> tire = Tire(tire_Fz=5000, tire_Sa=0.1, tire_Ls=0.05, tire_friction_coef=1.0, tire_Ca=0.02)
    >>> tire.run()
    
    References
    ----------
    Pacejka, H. B. (2002). Tire and Vehicle Dynamics. Elsevier.
    """

    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=0, B0=0, B1=None, 
                B2=0, B3=None, omega=315, slip_angle_start=-10, slip_angle_end=10, angle_step=0.5, WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
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

        # Entradas de rigidez do pneu
        self.track_y = track_y
        self.tire_k = tire_k

    def calculate_kinematics(self):
        """
        Calcular Cinemática.

        Calcula a cinemática do mecanismo de quatro barras para cada ângulo de deslizamento.

        Returns:
        -------
        Nenhum : NoneType
            O método atualiza os atributos da instância com valores calculados.

        Examples:
        --------
        >>> tire = Tire(tire_Fz=5000, tire_Sa=0.1, tire_Ls=0.05, tire_friction_coef=1.0, tire_Ca=0.02)
        >>> tire.calculate_kinematics()

        References
        ----------
        Textbook reference for four-bar linkage mechanisms and calculations.
        """

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
        """
        Calcular Forças do Pneu.

        Calcula as forças do pneu e o momento de auto-alinhamento usando os parâmetros do modelo Pacejka.

        Parâmetros
        ----------
        params : tuple
            Os parâmetros para o modelo de Pacejka.

        Returns:
        -------
        tuple : (float, float, float)
            Força lateral, momento de auto-alinhamento e força longitudinal do pneu.

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
        """
        Calcular Força.

        Calcula a força com base na deformação do pneu.

        Returns:
        -------
        force : float
            A força calculada com base em track_y e tire_k.

        Examples:
        --------
        >>> tire = Tire(track_y=10, tire_k=100)
        >>> tire.calcular_forca()
        1000
        
        References
        ----------
        Custom calculation based on Hooke's Law for spring-like behavior.
        """
        force = self.track_y * self.tire_k
        return force

    def slip_ratio_1(self, velocidade_angular, raio_pneu):
        """
        Calcular Razão de Escorregamento.

        Calcula a razão de escorregamento com base na velocidade angular e no raio do pneu.

        Parâmetros
        ----------
        velocidade_angular : float
            A velocidade angular do pneu.
        raio_pneu : float
            O raio do pneu.

        Returns:
        -------
        slip_ratio : float
            A razão de escorregamento calculada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.slip_ratio_1(100, 0.3)
        0.33333333333333337
        
        References
        ----------
        Basic slip ratio calculations from tire dynamics literature.
        """
        velocidade_linear = initial_speed
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1
        return value

    def plot_graph_slip_ratio(self, value=None):
        """
        Plotar Gráfico de Razão de Escorregamento.

        Plota a razão de escorregamento em relação à força no pedal.

        Parâmetros
        ----------
        value : array-like
            Os valores de razão de escorregamento.

        Returns:
        -------
        Nenhum : NoneType
            O método plota o gráfico e não retorna nada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.plot_graph_slip_ratio([0.1, 0.2, 0.3])
        
        References
        ----------
        Custom plotting function for visualizing slip ratio.
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
        """
        Mostrar Razão de Escorregamento.

        Exibe os valores da razão de escorregamento em relação ao RPM e à velocidade angular.

        Parâmetros
        ----------
        rpm_values : array-like
            Os valores de RPM.
        slip_ratio : array-like
            Os valores da razão de escorregamento.
        velocidade_angular : array-like
            Os valores da velocidade angular.

        Returns:
        -------
        Nenhum : NoneType
            O método plota o gráfico e imprime os valores da razão de escorregamento.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.show_slip_ratio([1000, 2000], [0.1, 0.2], [10, 20])
        
        References
        ----------
        Custom method for displaying slip ratio in graphical form.
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
        """
        Mostrar Força Longitudinal vs. RPM.

        Plota a força longitudinal em relação ao RPM.

        Parâmetros
        ----------
        rpm_values : array-like
            Os valores de RPM.
        tire_longitudinal_force : array-like
            Os valores da força longitudinal.

        Returns:
        -------
        Nenhum : NoneType
            O método plota o gráfico e não retorna nada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.show_longitudinal_rpm([1000, 2000], [1500, 1600])
        
        References
        ----------
        Basic longitudinal force vs. RPM calculation method.
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
        """
        Plotar Deformação.

        Plota o gráfico com base em uma lista de valores fornecidos que representam a variação da pista.

        Parâmetros
        ----------
        track_y_values : array-like
            Os valores que representam a deformação da pista.

        Returns:
        -------
        Nenhum : NoneType
            O método plota o gráfico e não retorna nada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.plotar_deformacao([0, 5, 10])
        
        References
        ----------
        Custom plotting method for track deformation vs. force calculation.
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
        """
        Plotar Camber.

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
        Nenhum : NoneType
            O método plota o gráfico e não retorna nada.

        Examples:
        --------
        >>> tire = Tire(tire_Sa=np.linspace(-0.1, 0.1, 100))
        >>> tire.plot_camber([0, 5, 10], [1, 6, 11], [2, 7, 12])
        
        References
        ----------
        Custom plotting function for visualizing camber effects on tire forces.
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
        """
        Plotar Gráfico de Forças do Pneu.

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
        Nenhum : NoneType
            O método plota o gráfico e não retorna nada.

        Examples:
        --------
        >>> tire = Tire(tire_Sa=np.linspace(-0.1, 0.1, 100), tire_Ls=np.linspace(0, 0.1, 100))
        >>> tire.plot_graph([0, 5, 10], [1, 6, 11], [2, 7, 12])
        
        References
        ----------
        Custom plotting method for visualizing tire forces.
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
        """
        Plotar Mecanismo.

        Plota o movimento do mecanismo com base na cinemática calculada.

        Returns:
        -------
        Nenhum : NoneType
            O método plota o mecanismo e não retorna nada.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.calculate_kinematics()
        >>> tire.plot_mechanism()
        
        References
        ----------
        Custom plotting method for visualizing four-bar linkage mechanism kinematics.
        """
        for i in range(len(self.theta2)):
            plt.figure(1)
            plt.plot([self.Ox[i], self.Ax[i]], [self.Oy[i], self.Ay[i]], 'r', linewidth=1)
            plt.plot([self.Ax[i], self.Bx[i]], [self.Ay[i], self.By[i]], 'k', linewidth=2)
            plt.plot([self.Bx[i], self.Cx[i]], [self.By[i], self.Cy[i]], 'r', linewidth=1)
            plt.plot([self.Ox[i], self.Cx[i]], [self.Oy[i], self.Cy[i]], 'r', linewidth=1, linestyle = 'dotted')

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
            plt.text((self.Ox[i] + self.Ax[i])/2, (self.Oy[i] + self.Ay[i])/2, f'{self.inner_slip[i]:.2f}°', fontsize=10, ha='center')

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
            plt.plot([self.Ox[i], self.Cx[i]], [self.Oy[i], self.Cy[i]], 'r', linewidth=1, linestyle = 'dotted')

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
            plt.text((self.Ox[i] + self.Ax[i])/2, (self.Oy[i] + self.Ay[i])/2, f'{self.inner_slip[i]:.2f}°', fontsize=10, ha='center')

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
        """
        Executar Cálculos e Plotar Mecanismo.

        Executa os cálculos cinemáticos e plota o movimento do mecanismo.

        Returns:
        -------
        Nenhum : NoneType
            O método executa os cálculos, plota o mecanismo e imprime o ângulo de deslizamento estático se determinado.

        Examples:
        --------
        >>> tire = Tire()
        >>> tire.run()
        O ângulo de toe é : 0.00°

        References
        ----------
        Custom method combining kinematic calculations and visualization.
        """
        self.calculate_kinematics()
        self.plot_mechanism()

        if self.static_slip_angle is not None:
            print(f"O ângulo de toe é : {self.static_slip_angle:.2f}°")
        else:
            print("Não foi possível determinar um ângulo estático para as rodas dentro do intervalo fornecido.")
