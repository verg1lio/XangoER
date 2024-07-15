import pandas as pd
from scipy.optimize import minimize
import numpy as np
import math
import matplotlib.pyplot as plt

class Drivetrain:
    def __init__(self, matriz_dados):
        self.matriz_dados = matriz_dados
        self.rpm = np.array([item['rpm'] for item in matriz_dados])
        self.ptc = np.array([item['ptc'] for item in matriz_dados])
        self.trq = np.array([item['trq'] for item in matriz_dados])

    def CalculateOutputs(self, cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, optimized_cp):
        peso = massa * 9.81
        rnet = (peso * (cgx * 0.001)) / (etex * 0.001)
        rned = massa * 9.81 - rnet
        ftr = (rnet * cfat) / (1 - ((cgy * 0.001) * cfat) / (etex * 0.001))
        tdcl = (ftr * cgy * 0.001) / (etex * 0.001)
        cnet = rnet + tdcl
        ptet = cnet * rpneu * 0.001 * cfat
        cpneu = cnet / 2
        tpneu = ftr * rpneu * 0.001
        redf = redp * red1
        tpwt = ptet / (red1 * redp * optimized_cp)
        acpr = (ftr / massa) / 9.81
        acfi = acpi * 9.81
        acfr = acpr * 9.81
        fti = massa * acfi
        tpi = fti * rpneu * 0.001
        tpwti = tpi / (red1 * redp * optimized_cp)
        tci = (fti * cgy) / etex
        tcr = (ftr * cgy) / etex
        cteti = rnet + tci
        return peso, rnet, rned, ftr, tdcl, cnet, ptet, cpneu, tpneu, redf, tpwt, acpr, acfi, acfr, fti, tpi, tpwti, tci, tcr, cteti

    def Transmission(self, cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, optimized_cp):
        outputs = self.CalculateOutputs(cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, optimized_cp)
        rpm_values, torque_values, power_values = self.CurveTorquePower(self.matriz_dados)

        # Print dos resultados obtidos
        print("Resultados:")
        labels = ["Peso", "Reação no eixo traseiro", "Reação no eixo dianteiro", "Força trativa",
                  "Transferência de carga longitudinal", "Carga no eixo traseiro", "Pico de torque no eixo traseiro",
                  "Carga no pneu", "Torque no pneu", "Redução final", "Torque necessário no motor",
                  "Aceleração primária real (g)", "Aceleração final ideal", "Aceleração final real",
                  "Força trativa ideal", "Torque no pneu ideal", "Torque no motor ideal", "Transferência de carga ideal",
                  "Transferência de carga real", "Carga total no eixo traseiro ideal"]

        for label, output in zip(labels, outputs):
            print(f"{label}: {output}")

        print("\nMatriz de RPM, Torque e Potência:")
        print("RPM\t\tTorque (Nm)\tPotência (kW)")
        for data in self.matriz_dados:
            rpm = data["rpm"]
            trq = data["trq"]
            ptc = data["ptc"]
            print(f"{rpm:.2f}\t\t{trq:.2f}\t\t{ptc:.2f}")

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

    def objective_function(self, cp, cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1):
        outputs = self.CalculateOutputs(cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, cp)
        tpwt = outputs[10]  # Torque necessário no motor
        target_tpwt = 80.0  # Alvo para o torque necessário no motor

        performance = self.CarPerformance(massa, rpneu, redp, red1, cp)
        max_vl = max([p['vl'] for p in performance])
        min_ff = min([p['ff'] for p in performance])

        penalty = 0
        if max_vl < 100:
            penalty += (100 - max_vl) ** 2
        if min_ff <= 0:
            penalty += abs(min_ff) * 1000

        return abs(tpwt - target_tpwt) + penalty

    def optimize_cp(self, cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1):
        initial_guess = 1.5  # Chute inicial para cp
        result = minimize(self.objective_function, initial_guess, args=(cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1), method='BFGS')
        return result.x[0], result.fun

    def CarPerformance(self, massa, rpneu, redp, red1, optimized_cp):
        peso = massa * 9.81
        rdt = 0.9
        mtrd = 0.9
        cfar = 0.54
        da = 1.162
        af = 1.06
        bscf = 0.015
        spdf = 0.012

        parametros = []

        for dado in self.matriz_dados:
            ftf = ((dado["trq"] * redp * red1 * optimized_cp) / (rpneu * 0.001)) * rdt
            va = (dado["rpm"] * 2 * math.pi) / (60 * redp * red1 * optimized_cp)
            vl = ((va * (rpneu * 0.001)) * mtrd) * 3.6
            fa = (da * vl ** 2 * cfar * af) / 2
            rr = (bscf + (3.24 * spdf * ((vl / 100 * 0.44704) ** 2.5))) * peso
            ff = ftf - fa - rr

            parametro = {
                "ftf": ftf,
                "va": va,
                "vl": vl,
                "fa": fa,
                "rr": rr,
                "ff": ff
            }

            parametros.append(parametro)

        return parametros

    def print_car_performance(self, performance):
        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]\tForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for param in performance:
            print(f"{param['ftf']}\t{param['va']}\t{param['vl']}\t{param['fa']}\t{param['rr']}\t{param['ff']}")

    def HalfShaftsSizing(self, redp, red1, optimized_cp, fsi=1.25, tet=786, tec=471.6, dif=1):
        tmax = max(data["trq"] for data in self.matriz_dados)
        tmsx = tmax * redp * red1 * optimized_cp * dif
        tmp = tmsx * fsi
        dsx = (((2 * tmp) / (math.pi * tec * 10 ** 6)) ** (1 / 3)) * 2000
        fso = (math.pi * (((dsx / 1000) / 2) ** 3) * tec * (10 ** 6)) / (2 * tmsx)
        fs1p = (math.pi * ((0.0254 / 2) ** 3) * tec * (10 ** 6)) / (2 * tmsx)

        print("Dimensionamento dos Semieixos:")
        print("Torque máximo do motor:", tmax, "Nm")
        print("Torque máximo nos semieixos:", tmsx, "Nm")
        print("Torque máximo de projeto:", tmp, "Nm")
        print("Diâmetro dos semieixos:", dsx, "mm")
        print("Fator de segurança ideal:", fsi)
        print("Fator de segurança obtido:", fso)
        print("Fator de segurança para 1 polegada:", fs1p)

        return tmax, tmsx, tmp, dsx, fso, fs1p

    def CurveTorquePower(self, matriz_dados):
        rpm_values = [data["rpm"] for data in matriz_dados]
        torque_values = [data["trq"] for data in matriz_dados]
        power_values = [data["ptc"] for data in matriz_dados]

        return rpm_values, torque_values, power_values

# Carregue os dados da matriz e os parâmetros de uma planilha Excel
arquivo_excel = r'G:\.shortcut-targets-by-id\1aB0nU7DQ_mchppOOiVqn5R-n9pnlVeYw\XANGO E-RACING\7 - DINÂMICA\Técnico\Softwares\Datasheet - Dinâmica.xlsx'

# Planilha com Matriz de Torque e Potência"
df_matriz = pd.read_excel(arquivo_excel, sheet_name='Matriz Torque x Potência')
matriz_dados = df_matriz.dropna().to_dict('records')  # Remove linhas com células vazias

# Parâmetros do sistema
df_parametros = pd.read_excel(arquivo_excel, sheet_name='Transmissão')

# Parâmetros gerais
df_gerais = pd.read_excel(arquivo_excel, sheet_name='Variáveis Globais')

# Extraia os parâmetros 
cgx = df_gerais.at[1, 'Valor']
cgy = df_gerais.at[2, 'Valor']
massa = df_gerais.at[0, 'Valor']
etex = df_gerais.at[8, 'Valor']
cfat = df_parametros.at[0, 'Valor']
rpneu = df_gerais.at[6, 'Valor']
acpi = df_parametros.at[1, 'Valor']
redp = df_parametros.at[2, 'Valor']
red1 = df_parametros.at[3, 'Valor']

# Crie uma instância da classe Drivetrain
drivetrain = Drivetrain(matriz_dados)

# Otimize o valor de cp
optimized_cp, fun_value = drivetrain.optimize_cp(cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1)

print(f"O valor otimizado de cp é: {optimized_cp}")
print(f"O valor da função objetivo é: {fun_value}")

# Calcule a transmissão
drivetrain.Transmission(cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, optimized_cp)

# Calcule o desempenho do carro
performance = drivetrain.CarPerformance(massa, rpneu, redp, red1, optimized_cp)
drivetrain.print_car_performance(performance)

# Dimensione os semi-eixos
tmax, tmsx, tmp, dsx, fso, fs1p = drivetrain.HalfShaftsSizing(redp, red1, optimized_cp)
