import math
import matplotlib.pyplot as plt
import numpy as np

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
        
        dt_model.showResults()        
        
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
            variacao = range(self.rpm, self.new_rpm, 25)
            self.rpm = self.new_rpm
            self.new_rpm = 0
            
        else:
            variacao = [self.rpm]
        
        parametros = []
        

        for rpm in variacao:
            # Cálculo da força trativa (N)
            forca_trativa = ((self.torque * self.reducao_primaria * self.reducao_unica * self.cp) / (self.raio_pneu * 0.001)) * rendimento_transmissao
           
            # Cálculo da velocidade angular (rad/s)
            velocidade_angular = (rpm * 2 * math.pi) / (60 * self.reducao_primaria * self.reducao_unica * self.cp)
            
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

        Este método utiliza os valores calculados pelo método `CarPerformance` para apresentar o desempenho do veículo. Os resultados incluem informações sobre força trativa, velocidade angular e linear, forças de arrasto, resistência ao rolamento e a força final do veículo. 

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
        diametro_semieixo = (((2 * torque_max_projeto) / (math.pi * tec * 10 ** 6)) ** (1 / 3)) * 2000
        
        # Calculando o fator de segurança obtido
        fator_seguranca_obtido = (math.pi * (((diametro_semieixo / 1000) / 2) ** 3) * tec * (10 ** 6)) / (2 * torque_max_semieixo)
        
        # Calculando o fator de segurança para 1 polegada
        fs1p = (math.pi * ((0.0254 / 2) ** 3) * tec * (10 ** 6)) / (2 * torque_max_semieixo)
        
        # Print dos resultados obtidos
        print("Dimensionamento dos Semieixos:")
        print("Torque máximo do motor:", torque_max_motor, "Nm")
        print("Torque máximo nos semieixos:", torque_max_semieixo, "Nm")
        print("Torque máximo de projeto:", torque_max_projeto, "Nm")
        print("Diâmetro dos semieixos:", diametro_semieixo, "mm")
        print("Fator de segurança ideal:", fsi)
        print("Fator de segurança obtido:", fator_seguranca_obtido)
        print("Fator de segurança para 1 polegada:", fs1p)


class Tire:

    """Modelo de Pneu.

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

    @staticmethod
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
        
    def plot_graph(self, tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces, tire_lateral_experimental=None, tire_auto_align_experimental=None, angles=None, ratio=None):

        # Definindo um tamanho para a figura
        plt.figure(figsize=(20, 7))

        # Plotagem força lateral
        plt.subplot(1, 3, 1)
        plt.plot(self.tire_Sa, tire_lateral_forces, label='Curva Otimizada')
        plt.scatter(angles, tire_lateral_experimental, color='red', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Força Lateral do Pneu (N)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Plotagem torque auto-alinhante
        plt.subplot(1, 3, 2)
        plt.plot(self.tire_Sa, tire_auto_align_moment, label='Curva Otimizada')
        plt.scatter(angles, tire_auto_align_experimental, color='blue', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Torque auto-alinhante (N.m)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Plotagem força longitudinal
        plt.subplot(1, 3, 3)
        plt.plot(self.tire_Ls, tire_longitudinal_forces, label='Curva Sem Otimizar')
        plt.xlabel('Taxa de Escorregamento Longitudinal (Admensional)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal - Sem Dados Experimentais')
        plt.legend()
        plt.grid(True)

        plt.tight_layout(pad=3.0)  # Aumentando a distância entre os subplots
        plt.show()


def uniteExample():
    dt_model = Drivetrain(
        cgx = 853,  # mm
        cgy = 294,  # mm
        massa = 347,  # kg
        entre_eixos = 1567,  # mm
        coeficiente_atrito = 0.9 , # coeficiente de atrito
        raio_pneu = 259,  # mm
        aceleracao_ideal = 1.2,  # g
        reducao_primaria = 2.12,  # redução primária
        reducao_unica = 2.76,
        rpm = 0,
        torque= 244.14,
        cp = 2.22)
    
    dt_model.showResults()
    
    dt_model.HalfShaftsSizing()
    
    dt_model.new_rpm = 2635
    
    #Exemplo em lista
    dt_model.printCarPerformance()    
            
    #Recebendo os dados da performance do carro. Aqui que se encontra dados da velocidade angular
    performance_veiculo, rpm_faixa = dt_model.CarPerformance()

    #Filtrando a velocidade angular
    velocidade_angular = [dado['va'] for dado in performance_veiculo]
    
    #Printandas as velocidades registradas
    for dado in velocidade_angular:
        print(f'velocidade: {dado}')
        
    #Transformando num array    
    velocidade_angular = np.array(velocidade_angular)

    #Transformando num array  
    rpm_faixa = np.array(rpm_faixa)

    #Calcular o slip ratio
    slip_ratio = Tire.slip_ratio_1(velocidade_angular, 0.259)

    #Plotagem de gráfico do slip ratio e saídas de seus dados no terminal
    Tire.show_slip_ratio(rpm_faixa, slip_ratio, velocidade_angular)

    #Salvando os dados como array para cálculo de força longitudinal
    slip_ratios = np.array(slip_ratio)
    
    #Criando instância da classe Tire
    Slip_model = Tire(tire_Fz=1500, tire_Sa=0, tire_Ls=slip_ratios, tire_friction_coef=1.45, tire_Ca=0)
    
    #Dados experimentais para instância em Tire
    result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]
    
    #Recebendo valores de força lateral, torque auto allinhante e força longitudinal
    tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = Slip_model.Tire_forces(result)
    
    #Plotagem de gráficos
    Slip_model.plot_graph( tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces)

initial_speed = 1

uniteExample()
