import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import pi
from Software.Systems.models.Dynamics.tire import Tire
import io
import base64


class BrakeSystem:
        
    # Define a velocidade inicial e o tempo para o cálculo
    initial_speed = 10  # Velocidade inicial do veículo [m/s]
    tempo = 4  # Tempo para o cálculo [s]
    
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
        self.pi = pi
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
        ----------                
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

        return BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR

    def apply_brake(self, pedal_force=300):
        """
        Aplica o freio e calcula os resultados com base na força aplicada ao pedal.

        Methods Utilizados
        ------------------
        self.calculate_params(pedal_force) : method
        >>> Realiza os cálculos detalhados relacionados à frenagem, incluindo forças dinâmicas, torques 
        e resistências. Retorna os parâmetros necessários para determinar os resultados finais.
        
        Parâmetros
        ----------
        pedal_force : float
        >>> Força aplicada ao pedal do freio, em Newtons.

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
        time_intervals = np.arange(0, self.tempo, time_step)
        
        # Calcula a inércia da roda considerando a massa do pneu e da roda e o raio do pneu
        inertia_wheel = 0.5 * (self.m_tire + self.m_wheel) * (self.Rdp * 2)
        
        # Calcula a desaceleração angular com base no torque ajustado e na inércia da roda
        angular_desaceleration = (torque_ajustado / 4) / inertia_wheel
        
        # Calcula a velocidade angular inicial das rodas com base na velocidade inicial e no raio do pneu
        initial_angular_velocity = self.initial_speed / self.Rdp
        
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
            
            """if angular_velocities[i] <=0:
                angular_desaceleration = 0"""

        return angular_desaceleration, angular_velocity, angular_velocities, inertia_wheel
    
    # Imprime os resultados calculados
    def print_resultados(self, pedal_force):
        '''
        Essa função nos mostra os valores de alguns parâmetros que foram caculados anteriormente. Isso serve 
        para monitorarmos os valores e o comportamento da Classe BrakeSystem.

        Parâmetros
        ----------
        pedal_force : float
            >>> Força aplicada no pedal do freio, em Newtons.

        Methods Utilizados
        ------------------
        BrakeSystem.apply_brake(pedal_force) : method
            >>> Calcula o torque de frenagem ajustado com base na força do pedal de freio.
        BrakeSystem.calculate_angular_velocity(torque) : method
            >>> Determina a velocidade angular da roda após a aplicação do torque de frenagem.
    
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

    def results(self, pedal_force):
        '''
        Essa função nos mostra os valores de alguns parâmetros que foram caculados anteriormente. Isso serve 
        para monitorarmos os valores e o comportamento da Classe BrakeSystem.

        Parâmetros
        ----------
        pedal_force : float
            >>> Força aplicada no pedal do freio, em Newtons.

        Methods Utilizados
        ------------------
        BrakeSystem.apply_brake(pedal_force) : method
            >>> Calcula o torque de frenagem ajustado com base na força do pedal de freio.
        BrakeSystem.calculate_angular_velocity(torque) : method
            >>> Determina a velocidade angular da roda após a aplicação do torque de frenagem.
    
        Examples
        --------
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.print_resultados()
        '''

        forca_frenagem, torque_ajustado, forca_f, torque_disco_freio, resistencia_rolamento, torque_resistencia_rolamento = self.apply_brake(pedal_force=pedal_force)
        angular_desaceleration, angular_velocity, angular_velocities, inertia_wheel = self.calculate_angular_velocity(torque_ajustado)
        text = f"""
        Força de frenagem teórica: {forca_frenagem} N
        Torque ajustado: {torque_ajustado} N.m
        Forca Gerada Pelo Disco: {forca_f} N
        Torque do disco de freio: {torque_disco_freio} N.m
        Resistência ao rolamento: {resistencia_rolamento} N
        Torque de resistência ao rolamento: {torque_resistencia_rolamento} N.m
        Desaceleração Angular: {angular_desaceleration} rad/s^2
        Inércia da roda: {inertia_wheel} kg.m^2
        Velocidade Angular: {angular_velocity} rad/s

        """
        return text


    def graph_2(self, pedal_force, time_intervals, vehicle_speed):

        buffers = []
        imagens_base64 = []
        # Parâmetros de Pacejka
        result = [(0.3336564873588197), (1.6271741344929977), (10), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]
        torque_ajustado = self.apply_brake(pedal_force = pedal_force)[1]
        angular_velocity = self.calculate_angular_velocity(torque_ajustado)[2]
        
        slip_ratio = []
        for velocity in angular_velocity:
            ratio = Tire.slip_ratio_1(velocidade_angular = velocity, raio_pneu = self.Rdp, velocidade_linear = vehicle_speed)
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
        buff1 = io.BytesIO()
        plt.savefig(buff1, format="png",bbox_inches="tight")
        buff1.seek(0)
        img_64 = base64.b64encode(buff1.read()).decode('utf-8')
        imagens_base64.append(img_64)
        plt.close()

        plt.figure(figsize=(10, 6))
        plt.plot(time_intervals, slip_ratio, label='Longitudinal Force', color='blue')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Slip Ratio (Admensional)')
        plt.title('Slip Ratio vs. Tempo')
        plt.legend()
        plt.grid(True)
        buff2 = io.BytesIO()
        plt.savefig(buff2, format="png",bbox_inches="tight")
        buff2.seek(0)
        plt.close()    
        img_64 = base64.b64encode(buff2.read()).decode('utf-8')
        imagens_base64.append(img_64)
        
        plt.figure(figsize=(10, 6))
        plt.plot(time_intervals, longitudinal_forces, label='Longitudinal Force', color='blue')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal vs. Tempo')
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

    # Plota o gráfico da força longitudinal em relação à força aplicada ao pedal
    def show_graph(self, pedal_force, time_intervals, vehicle_speed):
        """Gera gráfico da dinâmica de frenagem.
    
        Utiliza a biblioteca `matplotlib` para criar um gráfico que relaciona a força longitudinal gerada durante 
        a frenagem com a força aplicada no pedal ao longo do tempo. Permite visualizar como a resposta do sistema 
        de freios varia conforme a aplicação de força.
    
        Parâmetros
        ----------
        pedal_force : float
            >>> Força aplicada no pedal do freio, em Newtons.
        time_intervals : numpy.ndarray
            >>> Array de intervalos de tempo para a simulação, em segundos.
    
        Methods Utilizados
        ------------------
        BrakeSystem.apply_brake(pedal_force) : method
            >>> Calcula o torque de frenagem ajustado com base na força do pedal de freio.
        BrakeSystem.calculate_angular_velocity(torque) : method
            >>> Determina a velocidade angular da roda após a aplicação do torque de frenagem.
    
        Returns
        -------
        None
            >>> A função não retorna valores, mas exibe o gráfico gerado.
    
        Example
        --------------
        >>> brake_system = BrakeSystem(params)
        >>> brake_system.show_graph(pedal_force=1000, time_intervals=np.arange(0, 10, 0.01))
        """
       
       
        # Parâmetros de Pacejka
        result = [(0.3336564873588197), (1.6271741344929977), (10), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]
        torque_ajustado = self.apply_brake(pedal_force = pedal_force)[1]
        angular_velocity = self.calculate_angular_velocity(torque_ajustado)[2]
        
        slip_ratio = []
        for velocity in angular_velocity:
            ratio = Tire.slip_ratio_1(velocidade_angular = velocity, raio_pneu = self.Rdp, velocidade_linear = vehicle_speed)
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

    def brake_system_example(self):
        """Exemplo de simulação do sistema de freio.
    
        Esta função demonstra o uso do sistema de freio ao calcular os resultados relacionados 
        à frenagem, como forças aplicadas, torque gerado e outros parâmetros baseados em um 
        conjunto de configurações do veículo e entradas de força no pedal. A função também 
        exibe gráficos relacionados à dinâmica da frenagem.
    
        Parâmetros (internos à função)
        ------------------------------
        params : dict
            >>> Dicionário contendo os parâmetros do sistema de freio e do veículo:
            - 'RedP' : float
                Redução do pedal.
            - 'a' : float
                Desaceleração em termos de gravidade (g).
            - 'psi' : float
                Distribuição de carga estática por eixo (adimensional).
            - 'μl' : float
                Coeficiente de atrito entre pastilha e disco.
            - 'pi' : float
                Valor de pi.
            - 'HCG' : float
                Altura do centro de gravidade do veículo, em metros.
            - 'μ' : float
                Coeficiente de atrito entre pneu e solo.
            - 'FzF' : float
                Força de reação estática no eixo dianteiro, em Newtons.
            - 'FzR' : float
                Força de reação estática no eixo traseiro, em Newtons.
            - 'Rdp' : float
                Raio do pneu, em metros.
            - 'Dcm' : float
                Diâmetro do cilindro mestre, em metros.
            - 'Dwc' : float
                Diâmetro do cilindro da roda, em metros.
            - 'Npast' : int
                Número de pastilhas por disco.
            - 'atrito_coeficiente' : float
                Coeficiente de atrito da pastilha.
            - 'red' : float
                Raio efetivo do disco, em metros.
            - 'Mt' : float
                Massa total do veículo, em kg.
            - 'm_wheel' : float
                Massa da roda, em kg.
            - 'm_tire' : float
                Massa do pneu, em kg.
            - 'L' : float
                Distância entre eixos, em metros.
            - 'c_rr' : float
                Coeficiente de resistência ao rolamento.
    
        time_step : float
            >>> Passo de tempo para os cálculos da simulação, em segundos.
        time_intervals : numpy.ndarray
            >>> Array com os intervalos de tempo para a simulação.
        pedal_forces : float
            >>> Força aplicada no pedal do freio, em Newtons.
    
        Methods Utilizados
        ------------------
        BrakeSystem.print_resultados(pedal_force) : method
            >>> Calcula e exibe os resultados da simulação, como torque e forças geradas pelo sistema de freio.
        BrakeSystem.show_graph(pedal_force, time_intervals) : method
            >>> Gera e exibe gráficos relacionados à frenagem ao longo do tempo, como torque aplicado, 
            velocidade angular e outros parâmetros.
    
        Returns
        -------
        None
            >>> A função não retorna valores, mas exibe resultados no console e gera gráficos explicativos.
    
        Exemplo de Uso
        --------------
        >>> brake_system_example()
        """
        params = {
        'RedP': 4,  # Redução do pedal
        'a': 0.8,  # Desaceleração [g]
        'psi': 0.48,  # Distribuição de carga estática por eixo [adm]
        'μl': 0.45,  # Coeficiente de atrito do contato pastilha/disco
     
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
        time_intervals = np.arange(0, 4, time_step)
        pedal_forces = 1000

        brake = BrakeSystem(params)
        brake.print_resultados(pedal_forces)
        brake.show_graph(pedal_force=pedal_forces, time_intervals=time_intervals)

    def initialSpeedAndTime(self):
    
        print("Initial speed fixed:\t", self.initial_speed)
        print("Time:\t", self.tempo)
        
        change = input("Change time and speed?[Y/N]\n")
        if change.upper() == "Y":
            self.initial_speed = float(input("Insert the vehice's initial speed:\t"))
            self.tempo = float(input("Time of the interval:\t"))
        elif change.upper() != "Y" and change.upper() != "N":
            print("Invalid Input")
            
            return