import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
import scipy.optimize as opt


####################


class Transmission:
    
    #Construtores com as variáveis mais recorrentes como entrada do código
    def __init__(self, cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1):
        self.cgx = cgx
        self.cgy = cgy
        self. massa = massa
        self.etex =etex
        self.cfat = cfat
        self.rpneu = rpneu
        self.acpi = acpi
        self.redp = redp
        self. red1 = red1
        
    # Valores experimentais para a curva de torque e potência
    matriz_dados = [
        {"rpm": 0, "ptc": 0.0, "trq": 245.89},
        {"rpm": 375, "ptc": 9.645583738270778, "trq": 245.61},
        {"rpm": 750, "ptc": 19.252680965147455, "trq": 245.12},
        {"rpm": 1125, "ptc": 28.866061704088473, "trq": 245.01},
        {"rpm": 1500, "ptc": 38.46766085790885, "trq": 244.88},
        {"rpm": 1875, "ptc": 47.988359793900806, "trq": 244.39},
        {"rpm": 2250, "ptc": 57.54597436327078, "trq": 244.22},
        {"rpm": 2625, "ptc": 67.11497779825737, "trq": 244.14},
        {"rpm": 3000, "ptc": 76.65570542895443, "trq": 243.99},
        {"rpm": 3375, "ptc": 86.08922063505362, "trq": 243.57},
        {"rpm": 3750, "ptc": 95.60756325402146, "trq": 243.45},
        {"rpm": 4125, "ptc": 105.02144248491958, "trq": 243.11},
        {"rpm": 4500, "ptc": 114.40861678954424, "trq": 242.77},
        {"rpm": 4875, "ptc": 123.84566647117964, "trq": 242.58},
        {"rpm": 5250, "ptc": 133.2403024463807, "trq": 242.34},
        {"rpm": 5625, "ptc": 142.61019709282843, "trq": 242.09},
        {"rpm": 6000, "ptc": 152.06727546916892, "trq": 242.01},
        {"rpm": 6375, "ptc": 161.53809902815016, "trq": 241.96},
        {"rpm": 6750, "ptc": 170.60206518096516, "trq": 241.34},
        {"rpm": 7125, "ptc": 178.6025469168901, "trq": 239.36},
        {"rpm": 7500, "ptc": 186.73026977211796, "trq": 237.74},
        {"rpm": 7875, "ptc": 194.7307515080429, "trq": 236.12},
        {"rpm": 8250, "ptc": 203.42477588806972, "trq": 235.45},
        {"rpm": 8625, "ptc": 210.73839121146113, "trq": 233.31},
        {"rpm": 9000, "ptc": 217.41265918230565, "trq": 230.67},
        {"rpm": 9375, "ptc": 221.6508880697051, "trq": 225.76},
        {"rpm": 9750, "ptc": 224.86019185656838, "trq": 220.22},
        {"rpm": 10125, "ptc": 227.2738459282842, "trq": 214.34},
        {"rpm": 10500, "ptc": 228.22501256702415, "trq": 207.55},
        {"rpm": 10875, "ptc": 229.2806425938338, "trq": 201.32},
        {"rpm": 11250, "ptc": 226.73660564678286, "trq": 192.45},
        {"rpm": 11625, "ptc": 219.9531616538204, "trq": 180.67},
        {"rpm": 12000, "ptc": 218.2263739946381, "trq": 173.65},
        {"rpm": 12375, "ptc": 209.6627324899464, "trq": 161.78},
        {"rpm": 12750, "ptc": 198.2173152647453, "trq": 148.45},
        {"rpm": 13125, "ptc": 189.38112642426276, "trq": 137.78},
        {"rpm": 13500, "ptc": 179.38170241286863, "trq": 126.88},
        {"rpm": 13875, "ptc": 170.95276369805632, "trq": 117.65},
        {"rpm": 14250, "ptc": 163.32104557640753, "trq": 109.44},
        {"rpm": 14625, "ptc": 152.94618171916892, "trq": 99.86},
        {"rpm": 15000, "ptc": 137.40470006702415, "trq": 87.47}
    ]

    def CalculateOutputs(self, optimized_cp):
        peso = self.massa * 9.81
        rnet = (peso * (self.cgx * 0.001)) / (self.etex * 0.001)
        rned = self.massa * 9.81 - rnet
        ftr = (rnet * self.cfat) / (1 - ((self.cgy * 0.001) * self.cfat) / (self.etex * 0.001))
        tdcl = (ftr * self.cgy * 0.001) / (self.etex * 0.001)
        cnet = rnet + tdcl
        ptet = cnet * self.rpneu * 0.001 * self.cfat
        cpneu = cnet / 2
        tpneu = ftr * self.rpneu * 0.001
        redf = self.redp * self.red1
        tpwt = ptet / (self.red1 * self.redp * optimized_cp)
        acpr = (ftr / self.massa) / 9.81
        acfi = self.acpi * 9.81
        acfr = acpr * 9.81
        fti = self.massa * acfi
        tpi = fti * self.rpneu * 0.001
        tpwti = tpi / (self.red1 * self.redp * optimized_cp)
        tci = (fti * self.cgy) / self.etex
        tcr = (ftr * self.cgy) / self.etex
        cteti = rnet + tci
        return peso, rnet, rned, ftr, tdcl, cnet, ptet, cpneu, tpneu, redf, tpwt, acpr, acfi, acfr, fti, tpi, tpwti, tci, tcr, cteti

    def show_results(self, optimized_cp): #Alterei o nome para diferenciar da classse
        peso, rnet, rned, ftr, tdcl, cnet, ptet, cpneu, tpneu, redf, tpwt, acpr, acfi, acfr, fti, tpi, tpwti, tci, tcr, cteti = Transmission.CalculateOutputs(self, optimized_cp)
        rpm_values, torque_values, power_values = Transmission.CurveTorquePower()
        
        # Print dos resultados obtidos
        print(f'''Resultados:
            Peso: {peso}N
            Reação no eixo traseiro: {rnet}N
            Reação no eixo dianteiro: {rned}N
            Força trativa: {ftr}N
            Transferência de carga longitudinal: {tdcl}N
            Carga no eixo traseiro: {cnet}N
            Pico de torque no eixo traseiro: {ptet}Nm
            Carga no pneu: {cpneu}N
            Torque no pneu: {tpneu}Nm
            Redução final: {redf}
            Torque necessário no motor: {tpwt}Nm
            Aceleração primária real (g): {acpr}
            Aceleração final ideal: {acfi}m/s²
            Aceleração final real: {acfr}m/s²
            Força trativa ideal: {fti}N
            Torque no pneu ideal: {tpi}Nm
            Torque no motor ideal: {tpwti}Nm
            Transferência de carga ideal: {tci}N
            Transferência de carga real: {tcr}N
            Carga total no eixo traseiro ideal: {cteti}N\n
            ''')
        
        print("\nMatriz de RPM, Torque e Potência:")
        print("RPM\t\tTorque (Nm)\tPotência (kW)")
        
        for data in Transmission.matriz_dados:
            rpm = data["rpm"]
            trq = data["trq"]
            ptc = data["ptc"]
            print("{:.2f}\t\t{:.2f}\t\t{:.2f}".format(rpm, trq, ptc))

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
        # Calcula os outputs
        outputs = Transmission.CalculateOutputs(self, cp)
        tpwt = outputs[10]  # Torque necessário no motor
        target_tpwt = 80.0  # Alvo para o torque necessário no motor
        
        # Calcular o desempenho do carro com o cp atual
        performance = Transmission.CarPerformance(self, cp)
        
        # Extrair a velocidade linear máxima e a menor força final
        max_vl = max([p['vl'] for p in performance])
        min_ff = min([p['ff'] for p in performance])
        
        # Penalidades para as novas condições
        penalty = 0
        if max_vl < 100:
            penalty += (100 - max_vl) ** 2
        if min_ff <= 0:
            penalty += abs(min_ff) * 1000
        
        
        
        # A função objetivo será a diferença absoluta entre tpwt e o valor alvo com penalidade
        return abs(tpwt - target_tpwt) + penalty

    # Função que realiza a otimização do parâmetro cp
    def optimize_cp(self):
        
        initial_guess = 1.5  # Chute inicial para cp
        result = minimize(Transmission.objective_function, initial_guess, args = (self),
                          #args = (self.cgx, self.cgy, self.massa, self.etex, self.cfat, self.rpneu, self.acpi, self.redp, self.red1),
                          method='BFGS')
        return result.x[0], result.fun

    # Função de desempenho do carro
    def CarPerformance(self, optimized_cp):
        peso = self.massa * 9.81
        rdt = 0.9
        mtrd = 0.9
        cfar = 0.54
        da = 1.162
        af = 1.06
        bscf = 0.015
        spdf = 0.012

        parametros = []

        for dado in Transmission.matriz_dados:
            # Cálculo da força trativa (N)
            ftf = ((dado["trq"] * self.redp * self.red1 * optimized_cp) / (self.rpneu * 0.001)) * rdt

            # Cálculo da velocidade angular (rad/s)
            va = (dado["rpm"] * 2 * math.pi) / (60 * self.redp * self.red1 * optimized_cp)

            # Cálculo da velocidade linear (km/h)
            vl = ((va * (self.rpneu * 0.001)) * mtrd) * 3.6

            # Cálculo da força de arrasto (N)
            fa = (da * vl ** 2 * cfar * af) / 2

            # Cálculo da resistência de rolamento (N)
            rr = (bscf + (3.24 * spdf * ((vl / 100 * 0.44704) ** 2.5))) * peso

            # Cálculo da força final (N)
            ff = ftf - fa - rr

            # Armazenar os parâmetros calculados em um dicionário
            parametro = {"ftf": ftf,"va": va,"vl": vl,"fa": fa,"rr": rr,"ff": ff}

            parametros.append(parametro)

        # Retornar os parâmetros calculados
        return parametros 
    
    def print_car_performance(self, optimized_cp):
        performance = Transmission.CarPerformance(self, optimized_cp)
        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]\tForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for param in performance:
            print(f"{param['ftf']}\t{param['va']}\t{param['vl']}\t{param['fa']}\t{param['rr']}\t{param['ff']}")

    def HalfShaftsSizing(self, optimized_cp, fsi=1.25, tet=786, tec=471.6, dif=1):
        # Obtendo o maior torque do motor a partir dos dados experimentais 
        tmax = max(data["trq"] for data in Transmission.matriz_dados)
        
        # Calculando o torque máximo nos semieixos
        tmsx = tmax * self.redp * self.red1 * optimized_cp * dif
        
        # Calculando o torque máximo de projeto
        tmp = tmsx * fsi
        
        # Calculando o diâmetro dos semieixos (mm)
        dsx = (((2 * tmp) / (math.pi * tec * 10 ** 6)) ** (1 / 3)) * 2000
        
        # Calculando o fator de segurança obtido
        fso = (math.pi * (((dsx / 1000) / 2) ** 3) * tec * (10 ** 6)) / (2 * tmsx)
        
        # Calculando o fator de segurança para 1 polegada
        fs1p = (math.pi * ((0.0254 / 2) ** 3) * tec * (10 ** 6)) / (2 * tmsx)
        
        # Print dos resultados obtidos
        print("Dimensionamento dos Semieixos:")
        print("Torque máximo do motor:", tmax, "Nm")
        print("Torque máximo nos semieixos:", tmsx, "Nm")
        print("Torque máximo de projeto:", tmp, "Nm")
        print("Diâmetro dos semieixos:", dsx, "mm")
        print("Fator de segurança ideal:", fsi)
        print("Fator de segurança obtido:", fso)
        print("Fator de segurança para 1 polegada:", fs1p)
        
        return tmax, tmsx, tmp, dsx, fso, fs1p

    def CurveTorquePower():
        matriz_dados = Transmission.matriz_dados
        rpm_values = [data["rpm"] for data in matriz_dados]
        torque_values = [data["trq"] for data in matriz_dados]
        power_values = [data["ptc"] for data in matriz_dados]
        return rpm_values, torque_values, power_values

    def main():
        
        transmission = Transmission(
        cgx = 853,  # mm
        cgy = 294,  # mm
        massa = 347,  # kg
        etex = 1567,  # mm
        cfat = 0.9 , # coeficiente de atrito
        rpneu = 259,  # mm
        acpi = 1.2,  # g
        redp = 2.12,  # redução primária
        red1 = 2.76,  # redução da marcha única
        )
        # Otimizar o parâmetro cp
        optimized_cp, objective_value = transmission.optimize_cp()
        print(f"Relação coroa-pinhão ideal: {optimized_cp}")
        print(f"Valor da função objetivo: {objective_value}")

        # Usar o parâmetro cp otimizado no restante do código
        transmission.show_results(optimized_cp)
        transmission.print_car_performance(optimized_cp)
        transmission.HalfShaftsSizing(optimized_cp)

#Transmission.main() #Pode ser útil quando quiser o retorno de todas as informações de transmissão

############### Testes ###########
car = Transmission(
        cgx = 853,  # mm
        cgy = 294,  # mm
        massa = 347,  # kg
        etex = 1567,  # mm
        cfat = 0.9 , # coeficiente de atrito
        rpneu = 259,  # mm
        acpi = 1.2,  # g
        redp = 2.12,  # redução primária
        red1 = 2.76,  # redução da marcha única
        )

optimized_cp, objective_value = car.optimize_cp() #É necessário inicializar esse método para ter retorno para os demais
print(optimized_cp)
#car.show_results(optimized_cp)
#car.HalfShaftsSizing(optimized_cp)
#car.print_car_performance(optimized_cp)

#################

class Dynamics:
    
    def __init__(self, 
                 spring_type=None,
                 spring_k=None,
                 spring_F=None,
                 spring_non_lin_coef=None,
                 tire_Fz=None,
                 tire_Sa=None,
                 tire_Ls=None,
                 damper_type=None,
                 damper_V=None,
                 damper_F_viscous=None,
                 damper_K_friction=None,
                 damper_F_static=None):
        
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
        # Modelo de amortecedor
        self.damper_type = damper_type # Coulumb, Integrated, Stribeck
        self.damper_V = damper_V # velocidade relativa amortecedor [m/s]
        self.damper_F_viscous = damper_F_viscous # força viscosa do fluído [N]
        self.damper_F_static = damper_F_static # coeficiente de fricção estática de coulumb [N]
        self.damper_K_friction = damper_K_friction # rigidez de fricção [N/m]
        
    def Tire(self, params):
        
        
        
        E, Cy, Cx, Cz, c1, c2 = params
        Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))
        D = 1.5 * self.tire_Fz
        Bz = Cs / (Cz * D)
        Bx = Cs / (Cx * D)
        By = Cs / (Cy * D)
        tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
        tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))

        return tire_lateral_force, (12 + (tire_auto_align_moment/58))


info = car.CalculateOutputs(optimized_cp)
#Retira as informações que serão entradas para a classe dynamics. É necessário criar o parâmetro otimizado primeiro

print(info) 
#Transforma em um dicioná´rio para facilitar o acesso pelas chaves
info = dict(zip(['peso', 'rnet', 'rned', 'ftr', 'tdcl', 'cnet', 'ptet', 'cpneu', 'tpneu', 'redf', 'tpwt', 'acpr',
                 'acfi', 'acfr', 'fti', 'tpi', 'tpwti', 'tci', 'tcr', 'cteti'], info))
print(info)

#Instancia a classe Dynamics com o parâmetro da carga vertical retirada do dict
dynamics = Dynamics(tire_Fz = info['cpneu'])

print(dynamics.tire_Fz)
