import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

class Motor:

    def __init__(self):
        self.rpm_interp = None
        self.torque_continuo_interp = None
        self.torque_pico_interp = None
        self.potencia_continua_interp = None
        self.potencia_pico_interp = None
        self._setup_data()

    def _setup_data(self): 
        rpm = np.array([0, 1000, 2000, 3000, 3500, 4000, 4500, 5000, 5500])
        potencia_continua = np.array([0, 20, 40, 60, 75, 85, 90, 95, 97])
        potencia_pico = np.array([0, 50, 100, 150, 180, 200, 210, 215, 215])
        torque_continuo = np.array([200, 200, 200, 200, 190, 170, 150, 130, 110])
        torque_pico = np.array([500, 500, 500, 500, 480, 450, 400, 370, 350])

        self.rpm_interp = np.linspace(0, 5500, 1000)
        self.torque_continuo_interp = interp1d(rpm, torque_continuo, kind='cubic')(self.rpm_interp)
        self.torque_pico_interp = interp1d(rpm, torque_pico, kind='cubic')(self.rpm_interp)
        self.potencia_continua_interp = interp1d(rpm, potencia_continua, kind='cubic')(self.rpm_interp)
        self.potencia_pico_interp = interp1d(rpm, potencia_pico, kind='cubic')(self.rpm_interp)

    def plot(self):
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.set_xlabel('Motor speed (RPM)', fontsize=12)
        ax1.set_ylabel('Power (kW)', color='tab:green', fontsize=12)
        ax1.plot(self.rpm_interp, self.potencia_continua_interp, label='Continuous power', color='tab:green', linestyle='--', linewidth=2)
        ax1.plot(self.rpm_interp, self.potencia_pico_interp, label='Peak power', color='tab:green', linewidth=2)
        ax1.set_ylim(0, 220)
        ax1.tick_params(axis='y', labelcolor='tab:green')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Torque (Nm)', color='tab:blue', fontsize=12)
        ax2.plot(self.rpm_interp, self.torque_continuo_interp, label='Continuous torque', color='tab:blue', linestyle='--', linewidth=2)
        ax2.plot(self.rpm_interp, self.torque_pico_interp, label='Peak torque', color='tab:blue', linewidth=2)
        ax2.set_ylim(0, 550)
        ax2.tick_params(axis='y', labelcolor='tab:blue')

        fig.suptitle('EMRAX 268MV LC - Torque and Power Curves', fontsize=14)
        ax1.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98))
        ax2.legend(loc="upper right", bbox_to_anchor=(0.98, 0.98))
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def get_torque():
        motor = Motor()
        return motor.torque_continuo_interp[:] * 2

    @staticmethod
    def get_potencia():
        motor = Motor()
        return motor.potencia_continua_interp[:] 
    @staticmethod
    def get_rpm():
        motor = Motor()
        return motor.rpm_interp[:]
    
# Teste de execução

# valores_t = Motor.get_torque()
# valores_p = Motor.get_potencia()
# print(f'torque:{valores_t}')
# print(f'potencia:{valores_p}')
class Drivetrain:

    def __init__(self, cgx, cgy, massa, entre_eixos, raio_pneu, reducao_primaria, reducao_final, cp, tempo_i):
        self.cgx = cgx
        self.cgy = cgy
        self.massa = massa
        self.entre_eixos =entre_eixos
        self.raio_pneu = raio_pneu
        self.reducao_primaria = reducao_primaria
        self.reducao_final = reducao_final
        self.cp = cp
        self.tempo_i = tempo_i
        self.tempo_f = 0

    def Reduções(self, eficiencia_corrente=0.96):
     
        # Dados do motor
        rpm_motor = Motor.get_rpm()
        torque_motor = Motor.get_torque()
        potencia_motor = Motor.get_potencia()


        # Etapa 1: Redução primária (Redutor)
        rpm_pos_redutor = rpm_motor / self.reducao_primaria
        torque_pos_redutor = torque_motor * self.reducao_primaria
        potencia_pos_redutor = potencia_motor  # desprezando perdas aqui

        # Etapa 2: Perdas por corrente de rolos (efeito poligonal e atrito)
        rpm_pos_corrente = rpm_pos_redutor  # corrente não altera rotação, só transmite
        torque_pos_corrente = torque_pos_redutor * eficiencia_corrente
        potencia_pos_corrente = potencia_pos_redutor * eficiencia_corrente
        
        # Etapa 3: Redução final (Diferencial)
        rpm_rodas = rpm_pos_corrente / self.reducao_final
        torque_rodas = torque_pos_corrente * self.reducao_final
        potencia_rodas = potencia_pos_corrente  # perdas desprezadas nesta etapa

        rendimento_transmissao = np.nanmean((potencia_rodas[1:] / potencia_motor[1:]))

        return  [torque_rodas, potencia_rodas, rpm_rodas,rendimento_transmissao]
    
    def CarPerformance(self):
        # Parâmetros do veículo
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        c_r = 0.025
        f_v = 0.012

        # Dados já reduzidos pelo sistema de transmissão
        torque_reduzido, potencia_reduzida, rpm_reduzido, rendimento_transmissao = self.Reduções()

        # Tempo de simulação
        self.tempo_f = len(rpm_reduzido) * 0.01  # tempo final automático baseado na resolução
        variacao_tempo = np.linspace(self.tempo_i, self.tempo_f, len(rpm_reduzido))

        parametros = []
        velocidade_veiculo = 0  # m/s inicial

        for rpm, torque, potencia, tempo in zip(rpm_reduzido, torque_reduzido, potencia_reduzida, variacao_tempo):
            
            # Velocidade angular da roda (rad/s)
            velocidade_angular = (rpm * 2 * np.pi) / 60

            # Velocidade linear na roda (m/s)
            velocidade_linear_roda = velocidade_angular * (self.raio_pneu * 0.001)

            # Força trativa 
            forca_trativa = torque / (self.raio_pneu * 0.001) * rendimento_transmissao

            # Resistência ao rolamento
            resistencia_rolamento = (c_r + 3.24 * f_v * velocidade_linear_roda ** 2.5) * peso

            # Força de arrasto
            forca_arrasto = 0.5 * densidade_ar * coeficiente_arrasto * area_frontal * velocidade_veiculo ** 2

            # Força líquida resultante
            forca_final = np.maximum(forca_trativa - forca_arrasto - resistencia_rolamento, 0)

            # Aceleração (F = ma)
            aceleracao = forca_final / self.massa

            # Integração da velocidade com método de Euler
            if len(parametros) == 0:
                dt = 0
            else:
                dt = tempo - parametros[-1]["tempo"]

            velocidade_veiculo += aceleracao * dt  # m/s
            velocidade_linear_carro = velocidade_veiculo * 3.6  # km/h

            parametro = {
                "ft": forca_trativa,
                "va": velocidade_angular,
                "vlc": velocidade_linear_carro,
                "fa": forca_arrasto,
                "rr": resistencia_rolamento,
                "ff": forca_final,
                "tempo": tempo
            }

            parametros.append(parametro)

        return parametros, rpm_reduzido, variacao_tempo


    def printCarPerformance(self, performance):
        

        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear do Carro [km/h]")
        
        for i in range(0, len(performance), 10):  
            param = performance[i]
            print(f"{param['ft']}\t{param['va']}\t{param['vlc']}")

        print("\nForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        
        for i in range(0, len(performance), 10): 
            param = performance[i]
            print(f"{param['fa']}\t{param['rr']}\t{param['ff']}")


def Instancias():
    # Instância do modelo Drivetrain com os parâmetros fornecidos
    dt_model = Drivetrain(
        cgx=853,             # Centro de gravidade no eixo X [mm]
        cgy=294,             # Centro de gravidade no eixo Y [mm]
        massa=347,           # Massa do veículo [kg]
        entre_eixos=1567,    # Entre-eixos [mm]
        raio_pneu=259,       # Raio do pneu [mm]
        reducao_primaria=2.12,  # Redução primária
        reducao_final=2.76,     # Redução final (diferencial)
        cp=2.22,             # Constante de potência (pode ser usada futuramente)
        tempo_i=0            # Tempo inicial [s]
    )

    # Instância do motor (para plot)
    motor = Motor()
    motor.plot()

    # Cálculo da performance do veículo com os valores reduzidos
    performance_veiculo, variacao_rpm, variacao_tempo = dt_model.CarPerformance()

    # Impressão dos resultados
    dt_model.printCarPerformance(performance_veiculo)

Instancias()
