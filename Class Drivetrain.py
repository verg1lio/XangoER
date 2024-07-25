class Drivetrain:
    
    #Construtores com as variáveis mais recorrentes como entrada do código
    def __init__(self, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, matriz_dados):
        self.cgx = cgx
        self.cgy = cgy
        self. massa = massa
        self.entre_eixos =entre_eixos
        self.coeficiente_atrito = coeficiente_atrito
        self.raio_pneu = raio_pneu
        self.aceleracao_ideal = aceleracao_ideal
        self.reducao_primaria = reducao_primaria
        self. reducao_unica = reducao_unica
        self.matriz_dados = matriz_dados
        self.cp = 1
        self.objective_function = 0

    def CalculateOutputs(self):
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

    def show_results(self):
        peso, reacao_traseira, reacao_dianteira, forca_trativa, transferencia_longitudinal, carga_traseira, pico_torque_traseiro, carga_pneu, torque_pneu, reducao_final, torque_necessario_motor, aceleracao_primaria_real, aceleracao_primaria_ideal, aceleraco_real_final, forca_trativa_ideal, torque_pneu_ideal, torque_motor_ideal, transferencia_carga_ideal, transferencia_carga_real, carga_traseira_ideal = Drivetrain.CalculateOutputs(self)
        rpm_values, torque_values, power_values = Drivetrain.CurveTorquePower(self)
        
        # Print dos resultados obtidos
        print(f'''Resultados:
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
        
        print("\nMatriz de RPM, Torque e Potência:")
        print("RPM\t\tTorque (Nm)\tPotência (kW)")
        
        for data in self.matriz_dados:
            rpm = data["rpm"]
            torque = data["trq"]
            potencia = data["ptc"]
            print("{:.2f}\t\t{:.2f}\t\t{:.2f}".format(rpm, torque, potencia))

        # Plotando o gráfico
        plt.figure(figsize=(10, 6))
        plt.plot(rpm_values, torque_values, label='Torque [Nm]', color='blue')
        plt.plot(rpm_values, power_values, label='Power [kW]', color='orange')
        plt.title("Curva de Torque e Potência")
        plt.xlabel('RPM')
        plt.ylabel('Torque [Nm] / Power [kW]')
        plt.legend()
        plt.grid(True)
        plt.show()

    # Definindo a função objetivo para otimização de cp
    def objective_function(cp, self): #Por conta da funcão de otimização, o self é passado como segundo argumento
        
        self.cp = cp
        
        # Calcula os outputs
        outputs = Drivetrain.CalculateOutputs(self)
        torque_necessario_motor = outputs[10]  # Torque necessário no motor
        target_torque_necessario_motor = 80.0  # Alvo para o torque necessário no motor
        
        # Calcular o desempenho do carro com o cp atual
        performance = Drivetrain.CarPerformance(self)
        
        # Extrair a velocidade linear máxima e a menor força final
        max_velocidade_linear = max([p['velocidade_linear'] for p in performance])
        min_forca_final = min([p['forca_final'] for p in performance])
        
        # Penalidades para as novas condições
        penalty = 0
        if max_velocidade_linear < 100:
            penalty += (100 - max_velocidade_linear) ** 2
        if min_forca_final <= 0:
            penalty += abs(min_forca_final) * 1000
        
        
        
        # A função objetivo será a diferença absoluta entre torque_necessario_motor e o valor alvo com penalidade
        return abs(torque_necessario_motor - target_torque_necessario_motor) + penalty

    # Função que realiza a otimização do parâmetro cp
    def optimize_cp(self):
        
        initial_guess = 1.5  # Chute inicial para cp
        result = minimize(Drivetrain.objective_function, initial_guess, args = (self), method='BFGS')
        
        self.cp = result.x[0]
        self.objective_function = result.fun
        
    # Função de desempenho do carro
    def CarPerformance(self):
        peso = self.massa * 9.81
        rendimento_transmissao = 0.9
        transmissao_motor_roda = 0.9
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        basic_f = 0.015
        speed_f = 0.012

        parametros = []

        for dado in self.matriz_dados:
            # Cálculo da força trativa (N)
            forca_trativa = ((dado["trq"] * self.reducao_primaria * self.reducao_unica * self.cp) / (self.raio_pneu * 0.001)) * rendimento_transmissao
        
            # Cálculo da velocidade angular (rad/s)
            velocidade_angular = (dado["rpm"] * 2 * math.pi) / (60 * self.reducao_primaria * self.reducao_unica * self.cp)
            
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
        return parametros 
    
    def print_car_performance(self):
        performance = Drivetrain.CarPerformance(self)
        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]")
        
        for param in performance:
            print(f"{param['forca_trativa']}\t{param['va']}\t{param['velocidade_linear']}")
        
        print("\nForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for param in performance:
            print(f"{param['fa']}\t{param['rr']}\t{param['forca_final']}")

    def HalfShaftsSizing(self, fsi=1.25, tet=786, tec=471.6, dif=1):
        # Obtendo o maior torque do motor a partir dos dados experimentais 
        torque_max_motor = max(data["trq"] for data in self.matriz_dados)
        
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
        
        return torque_max_motor, torque_max_semieixo, torque_max_projeto, diametro_semieixo, fator_seguranca_obtido, fs1p

    def CurveTorquePower(self):
        matriz_dados = self.matriz_dados
        rpm_values = [data["rpm"] for data in matriz_dados]
        torque_values = [data["trq"] for data in matriz_dados]
        power_values = [data["ptc"] for data in matriz_dados]
        return rpm_values, torque_values, power_values

    @classmethod
    def generate_model(cls, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, matriz_dados):
        model = cls(cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, matriz_dados)
        model.optimize_cp()
        return model
    
    def __str__(self):
        return f'Relação Coroa-Pinhão: {self.cp}\nFunção Objetivo: {self.objective_function}'
