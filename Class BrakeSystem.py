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
    
    def calculate_params(self, pedal_force):
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
            
        >>> pedal_force = 300
        
        >>> print("Resultados Calculados:")
            for i, result in resultados.items():
                print(f"{i}: {result}")
                
        Instanciando:        
        >>> BrakeSystem = BrakeSystem(params)
        >>> brake_system.calculate_params(300)
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
        
        # Exemplo:
        resultados = {
                'BF': BF, 'χ': χ, 'W': W, 'FzF_dyn': FzF_dyn, 'FzR_dyn': FzR_dyn,
                'τF': τF, 'τR': τR, 'FnF': FnF, 'FnR': FnR, 'Awc': Awc, 'Acm': Acm,
                'PF': PF, 'PR': PR, 'FaCM': FaCM, 'lF': lF, 'lR': lR
            }

        print("Resultados Calculados de calculate_params:")
        for i, result in resultados.items():
            print(f"{i}: {result}")
   
        return BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR

    def apply_brake(self, pedal_force):
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
            
        >>> pedal_force = 300
        
        >>> print("Resultados Calculados:")
            for i, result in resultados.items():
                print(f"{i}: {result}")

        Instanciando:
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.apply_brake(300)
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

        # Exemplo
        resultados = { 
                'Forca de Frenagem': forca_frenagem, 'Torque Ajustado': torque_ajustado, 
                'Forca Gerada Pelo Disco': forca_f, 'Torque do Disco de Freio': torque_disco_freio, 
                'Resistencia ao Rolamento': resistencia_rolamento, 'Torque de Resistencia ao Rolamento': torque_resistencia_rolamento
            }

        print("Resultados Calculados de apply_brake:")
        for i, result in resultados.items():
            print(f"{i}: {result}")

        return forca_frenagem, torque_ajustado, forca_f, torque_disco_freio, resistencia_rolamento, torque_resistencia_rolamento

    # Calcula a velocidade angular das rodas durante a frenagem
    def calculate_angular_velocity(self, torque_ajustado):
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
            
        >>> torque_ajustado = 190
        
        >>> print("Resultados Calculados de calculate_angular_velocity:")
            for i, result in resultados.items():
                print(f"{i}: {result}")

        Instanciando:
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.calculate_angular_velocity(190)
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
        # Cria uma série de intervalos de tempo igualmente espaçados de 0 até 'tempo' com 100 pontos
        time_intervals = np.linspace(0, tempo, 100)
        
        # Calcula a inércia da roda considerando a massa do pneu e da roda e o raio do pneu
        inertia_wheel = 0.5 * (self.m_tire + self.m_wheel) * (self.Rdp * 2)
        
        # Calcula a desaceleração angular com base no torque ajustado e na inércia da roda
        angular_desaceleration = (torque_ajustado / 4) / inertia_wheel
        
        # Calcula a velocidade angular inicial das rodas com base na velocidade inicial e no raio do pneu
        initial_angular_velocity = initial_speed / self.Rdp
        
        # Calcula o intervalo de tempo entre cada ponto de 'time_intervals'
        time_step = time_intervals[1] - time_intervals[0]
        
        # Define a velocidade angular inicial para o cálculo
        angular_velocity = initial_angular_velocity
        
        # Inicializa uma lista para armazenar as velocidades angulares ao longo do tempo
        angular_velocities = []

        # Itera sobre cada intervalo de tempo
        for i in time_intervals:
            # Atualiza a velocidade angular subtraindo a desaceleração angular multiplicada pelo intervalo de tempo
            angular_velocity -= angular_desaceleration * time_step
            # Adiciona a velocidade angular atual à lista de velocidades angulares
            angular_velocities.append(angular_velocity)
        
        # Exemplo
        resultados = { 
                'Desaceleração angular': angular_desaceleration,'Velocidade angular': angular_velocity,
                'Primeiro resultado da lista de Velocidades angulares': angular_velocities[0], 
                'Inércia do pneu': inertia_wheel 
            }

        print("Resultados Calculados de calculate_angular_velocity:")
        for i, result in resultados.items():
            print(f"{i}: {result}")

        return angular_desaceleration, angular_velocity, angular_velocities, inertia_wheel


    # Imprime os resultados calculados
    def print_resultados(self):
        '''
        Essa função nos mostra os valores de alguns parâmetros que foram caculados anteriormente. Isso serve 
        para monitorarmos os valores e o comportamento da Classe BrakeSystem.
    
        Examples
        --------
        >>> BrakeSystem = BrakeSystem(params)
        >>> BrakeSystem.print_resultados()
        '''
        print("Resultados Calculados:")
        print("Força de frenagem teórica:", BrakeSystem.apply_brake(pedal_force=823)[3], 'N')
        print("Inércia da roda:", inertia_wheel, "kg.m^2")
        print("Resistência ao rolamento:", resistencia_rolamento, "N")
        print("Torque do disco de freio:", BrakeSystem.apply_brake(pedal_force=823)[4], "N.m")
        print("Torque de resistência ao rolamento:", torque_resistencia_rolamento, "N.m")
        print("Torque ajustado:", BrakeSystem.apply_brake(pedal_force=823)[2], "N.m")

    # Plota o gráfico da força longitudinal em relação à força aplicada ao pedal
    def show_graph(self, longitudinal_force, pedal_forces):
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
        plt.figure(figsize=(10, 6))
        plt.plot(pedal_forces, longitudinal_force, label='Longitudinal Force', color='blue')
        plt.xlabel('Pedal Force (N)')
        plt.ylabel('Longitudinal Force (N)')
        plt.title('Longitudinal Force vs. Pedal Force')
        plt.legend()
        plt.grid(True)
        plt.show()
