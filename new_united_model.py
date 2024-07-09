import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np

class Drivetrain:
    
    #Construtores com as variáveis mais recorrentes como entrada do código
    def __init__(self, *, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, cp = 1):
        self.cgx = cgx
        self.cgy = cgy
        self. massa = massa
        self.entre_eixos =entre_eixos
        self.coeficiente_atrito = coeficiente_atrito
        self.raio_pneu = raio_pneu
        self.aceleracao_ideal = aceleracao_ideal
        self.reducao_primaria = reducao_primaria
        self. reducao_unica = reducao_unica
        self.cp = cp
    # Valores experimentais para a curva de torque e potência
    matriz_dados = [
        {"rpm": 0, "potencia": 0.0, "torque": 245.89},
        {"rpm": 375, "potencia": 9.645583738270778, "torque": 245.61},
        {"rpm": 750, "potencia": 19.252680965147455, "torque": 245.12},
        {"rpm": 1125, "potencia": 28.866061704088473, "torque": 245.01},
        {"rpm": 1500, "potencia": 38.46766085790885, "torque": 244.88},
        {"rpm": 1875, "potencia": 47.988359793900806, "torque": 244.39},
        {"rpm": 2250, "potencia": 57.54597436327078, "torque": 244.22},
        {"rpm": 2625, "potencia": 67.11497779825737, "torque": 244.14},
        {"rpm": 3000, "potencia": 76.65570542895443, "torque": 243.99},
        {"rpm": 3375, "potencia": 86.08922063505362, "torque": 243.57},
        {"rpm": 3750, "potencia": 95.60756325402146, "torque": 243.45},
        {"rpm": 4125, "potencia": 105.02144248491958, "torque": 243.11},
        {"rpm": 4500, "potencia": 114.40861678954424, "torque": 242.77},
        {"rpm": 4875, "potencia": 123.84566647117964, "torque": 242.58},
        {"rpm": 5250, "potencia": 133.2403024463807, "torque": 242.34},
        {"rpm": 5625, "potencia": 142.61019709282843, "torque": 242.09},
        {"rpm": 6000, "potencia": 152.06727546916892, "torque": 242.01},
        {"rpm": 6375, "potencia": 161.53809902815016, "torque": 241.96},
        {"rpm": 6750, "potencia": 170.60206518096516, "torque": 241.34},
        {"rpm": 7125, "potencia": 178.6025469168901, "torque": 239.36},
        {"rpm": 7500, "potencia": 186.73026977211796, "torque": 237.74},
        {"rpm": 7875, "potencia": 194.7307515080429, "torque": 236.12},
        {"rpm": 8250, "potencia": 203.42477588806972, "torque": 235.45},
        {"rpm": 8625, "potencia": 210.73839121146113, "torque": 233.31},
        {"rpm": 9000, "potencia": 217.41265918230565, "torque": 230.67},
        {"rpm": 9375, "potencia": 221.6508880697051, "torque": 225.76},
        {"rpm": 9750, "potencia": 224.86019185656838, "torque": 220.22},
        {"rpm": 10125, "potencia": 227.2738459282842, "torque": 214.34},
        {"rpm": 10500, "potencia": 228.22501256702415, "torque": 207.55},
        {"rpm": 10875, "potencia": 229.2806425938338, "torque": 201.32},
        {"rpm": 11250, "potencia": 226.73660564678286, "torque": 192.45},
        {"rpm": 11625, "potencia": 219.9531616538204, "torque": 180.67},
        {"rpm": 12000, "potencia": 218.2263739946381, "torque": 173.65},
        {"rpm": 12375, "potencia": 209.6627324899464, "torque": 161.78},
        {"rpm": 12750, "potencia": 198.2173152647453, "torque": 148.45},
        {"rpm": 13125, "potencia": 189.38112642426276, "torque": 137.78},
        {"rpm": 13500, "potencia": 179.38170241286863, "torque": 126.88},
        {"rpm": 13875, "potencia": 170.95276369805632, "torque": 117.65},
        {"rpm": 14250, "potencia": 163.32104557640753, "torque": 109.44},
        {"rpm": 14625, "potencia": 152.94618171916892, "torque": 99.86},
        {"rpm": 15000, "potencia": 137.40470006702415, "torque": 87.47}
    ]

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

    def show_results(self): #Alterei o nome para diferenciar da classse
        peso, reacao_traseira, reacao_dianteira, forca_trativa, transferencia_longitudinal, carga_traseira, pico_torque_traseiro, carga_pneu, torque_pneu, reducao_final, torque_necessario_motor, aceleracao_primaria_real, aceleracao_primaria_ideal, aceleraco_real_final, forca_trativa_ideal, torque_pneu_ideal, torque_motor_ideal, transferencia_carga_ideal, transferencia_carga_real, carga_traseira_ideal = Drivetrain.CalculateOutputs(self)
        rpm_values, torque_values, power_values = Drivetrain.CurveTorquePower()
        
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
        
        for data in Drivetrain.matriz_dados:
            rpm = data["rpm"]
            torque = data["torque"]
            potencia = data["potencia"]
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
        
        return result.fun

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

        for dado in Drivetrain.matriz_dados:
            # Cálculo da força trativa (N)
            forca_trativa = ((dado["torque"] * self.reducao_primaria * self.reducao_unica * self.cp) / (self.raio_pneu * 0.001)) * rendimento_transmissao
           
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
        torque_max_motor = max(data["torque"] for data in Drivetrain.matriz_dados)
        
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

    def CurveTorquePower():
        matriz_dados = Drivetrain.matriz_dados
        rpm_values = [data["rpm"] for data in matriz_dados]
        torque_values = [data["torque"] for data in matriz_dados]
        power_values = [data["potencia"] for data in matriz_dados]
        return rpm_values, torque_values, power_values


class Dynamics:
    
    def __init__(self, spring_type=None, spring_k=None, spring_F=None, spring_non_lin_coef=None,
                 tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, damper_type=None, 
                 damper_V=None, damper_F_viscous=None, damper_K_friction=None, damper_F_static=None):
        # Modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # rigidez da mola [N/m]
        self.spring_F = spring_F  # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # coeficiente de ganho não-linear
        # Modelo de pneu
        self.tire_Fz = tire_Fz  # carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # slip angle do pneu [rad]
        self.tire_Ls = tire_Ls  # longitudinal slip do pneu [Admensional]
        self.tire_type = 'Default'
        self.tire_friction_coef = tire_friction_coef # coeficiente de fricção entre o pneu e a pista
        # Modelo de amortecedor
        self.damper_type = damper_type # Coulumb, Integrated, Stribeck
        self.damper_V = damper_V # velocidade relativa amortecedor [m/s]
        self.damper_F_viscous = damper_F_viscous # força viscosa do fluído [N]
        self.damper_F_static = damper_F_static # coeficiente de fricção estática de coulumb [N]
        self.damper_K_friction = damper_K_friction # rigidez de fricção [N/m]

    def Tire(self, params):
        E, Cy, Cx, Cz, c1, c2 = params
        Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))
        D = self.tire_friction_coef * self.tire_Fz
        Bz = Cs / (Cz * D)
        Bx = Cs / (Cx * D)
        By = Cs / (Cy * D)
        tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
        tire_longitudinal_force = D * np.sin(Cx * np.arctan(9 * Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
        tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))


        return tire_lateral_force, (10 + (tire_auto_align_moment/55)), tire_longitudinal_force
    
    #criação de uma função para o slip ratio dentro da classe Dynamics tendo em vista a reutilização com valores de freio e transmissão
    def slip_ratio(velocidade_angular, velocidade_linear, raio_pneu):
        slip_ratio = []
        for i in range(len(velocidade_angular)):
            if velocidade_linear != 0:
                value = (velocidade_angular[i] * raio_pneu /velocidade_linear) - 1
                
                slip_ratio.append(value)
        return slip_ratio
    
    def show_slip_ratio(rpm_values, slip_ratio, velocidade_angular):#,
                        #longitudinal_force):
        
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

#######Instanciando a classe Transmission#######

transmission_model = Drivetrain(
        cgx = 853,  # mm
        cgy = 294,  # mm
        massa = 347,  # kg
        entre_eixos = 1567,  # mm
        coeficiente_atrito = 0.9 , # coeficiente de atrito
        raio_pneu = 259,  # mm
        aceleracao_ideal = 1.2,  # g
        reducao_primaria = 2.12,  # redução primária
        reducao_unica = 2.76,  # redução da marcha única
        )


#Após instanciar a classe, chamar o método optimize_cp() para os demais terem um retorno válido, caso contrário terá valor 1
transmission_model.optimize_cp()
transmission_model.show_results()
transmission_model.HalfShaftsSizing()
transmission_model.print_car_performance()
####### Cálculo de Slip Ratio ############

#Recebendo os dados da performance do carro. Aqui que se encontra dados da velocidade angular
performance_veiculo = transmission_model.CarPerformance()

#Filtrando a velocidade angular
velocidade_angular = [dado["va"] for dado in performance_veiculo]

#Atrinuindo um valor para o raio do pneu(em m)
raio_pneu = 0.259

#Chama a função de slip ratio e salva seus valores numa lista
#Utilizando uma velocidade arbitrária e fixa para o veículo
slip_ratio = Dynamics.slip_ratio(velocidade_angular, 30, raio_pneu)

rpm = [dado["rpm"] for dado in Drivetrain.matriz_dados]

Dynamics.show_slip_ratio(rpm, slip_ratio, velocidade_angular)
