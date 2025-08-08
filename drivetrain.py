"""
Simulação de desempenho de um veículo Fórmula SAE elétrico.
----------------------------------------------------------
Este código implementa:
1. Modelo do trem de força (motor + transmissão).
2. Cálculo das forças resistivas (arrasto e rolamento).
3. Simulação da aceleração do veículo ao longo do tempo.
4. Modelo simplificado de pneu para slip ratio.
5. Geração de gráficos de desempenho e slip ratio.

Autor: Marco Affonso e Igor Maia
"""

import numpy as np
import matplotlib.pyplot as plt
from pwt1 import Motor                  # Classe Motor definida em outro módulo
from tire2 import Tire

# ============================================================================
# Classe Drivetrain: modelo do sistema de transmissão e veículo
# ============================================================================
class Drivetrain:

    def __init__(self, cgx, cgy, massa, massa_roda, entre_eixos,
                 raio_pneu, reducao_primaria, reducao_final, tempo_i,
                 tire_friction_coef, tire_Fz):

        self.cgx = cgx
        self.cgy = cgy
        self.massa = massa
        self.massa_roda = massa_roda
        self.entre_eixos = entre_eixos
        self.raio_pneu = raio_pneu
        self.reducao_primaria = reducao_primaria
        self.reducao_final = reducao_final
        self.tempo_i = tempo_i
        self.tempo_f = 0

        self.tire_friction_coef = tire_friction_coef
        self.tire_Fz = tire_Fz

        self.motor = Motor(
            rs=0.04585,
            ld=0.00067,
            lq=0.00067,
            jm=0.05769,
            kf=0.1,
            lambda_m=0.13849,
            p=10,
            valor_mu=0.9
        )
        self.tire = Tire(
            tire_Fz=self.tire_Fz,
            tire_Sa=0,
            tire_Ls=1,  # será atualizado depois
            tire_friction_coef=self.tire_friction_coef,
            tire_Ca=0
        )

    def CarPerformance(self):
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        raio = self.raio_pneu * 0.001
        c_r = 0.012
        b = 0.01
        eficiencia_transmissao = 0.95

        if self.tempo_f == 0:
            raise ValueError("Defina o tempo final: dt_model.tempo_f = <valor>")

        dt = 0.001
        variacao_tempo = np.arange(self.tempo_i, self.tempo_f, dt)

        performance = []
        a_linear = 0.0
        v_linear_ms = 0.0

        for tempo in variacao_tempo:
            i_total = self.reducao_primaria * self.reducao_final

            forca_rolamento = c_r * peso
            forca_arrasto = 0.5 * densidade_ar * v_linear_ms**2 * coeficiente_arrasto * area_frontal
            forca_resistiva_total = forca_rolamento + forca_arrasto

            torque_rr = forca_rolamento * raio
            torque_fa = forca_arrasto * raio
            torque_atr_mec = b * (self.motor.wm)
            torque_carga = torque_rr + torque_fa + torque_atr_mec
            torque_carga_motor = torque_carga / (i_total * eficiencia_transmissao)

            self.motor.set_load(torque_carga_motor)
            torque_motor_gerado = self.motor.step(dt)

            torque_trativo_rodas = torque_motor_gerado * i_total * eficiencia_transmissao
            forca_trativa_potencial = torque_trativo_rodas / raio

            # Parâmetros de Pacejka
            params = [0.333, 1.627, 1, 4.396, 931.4, 366.4]
            _, _, forca_longitudinal_maxima = self.tire.Tire_forces(params)
            print(forca_longitudinal_maxima)
            forca_trativa_real = min(forca_trativa_potencial, forca_longitudinal_maxima)
            forca_liquida = forca_trativa_real - forca_resistiva_total
            a_linear = forca_liquida / self.massa
            v_linear_ms += a_linear * dt

            velocidade_angular_roda = self.motor.wm / i_total

            self.motor.store_history(tempo)
            performance.append({
                "tempo": tempo,
                "vlk": v_linear_ms * 3.6,
                "vlm": v_linear_ms,
                "va": velocidade_angular_roda,
                "fa": forca_arrasto,
                "rr": forca_rolamento,
                "ff": forca_trativa_real,
                "fp": forca_trativa_potencial
            })

        return performance, variacao_tempo


    # ------------------------------------------------------------------------
    def printCarPerformance(self, performance):
        """
        Gera gráficos de desempenho do veículo com base no histórico do motor
        e nos dados de simulação armazenados em performance.
        """

        # Histórico do motor já simulado
        motor_history = self.motor.history
        perf_dict = {k: [dic[k] for dic in performance] for k in performance[0]}

        # --- Figura 1: Resumo geral ---
        plt.figure(figsize=(14, 10))

        # Velocidade linear
        plt.subplot(2, 2, 1)
        plt.plot(perf_dict['tempo'], perf_dict['vlk'], color='purple')
        plt.title("Velocidade Linear do Veículo")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade [km/h]")
        plt.grid(True)

        # Forças
        plt.subplot(2, 2, 2)
        plt.plot(perf_dict['tempo'], perf_dict['ff'], label='Força Trativa', color='green')
        plt.plot(perf_dict['tempo'], perf_dict['fa'], label='Força de Arrasto', color='red')
        plt.plot(perf_dict['tempo'], perf_dict['rr'], label='Resist. Rolamento', color='orange')
        plt.title("Forças Atuantes")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Força [N]")
        plt.legend()
        plt.grid(True)

        # RPM do motor
        plt.subplot(2, 2, 3)
        rpm_motor = np.array(motor_history['velocidade']) * 60 / (2 * np.pi)
        plt.plot(motor_history['tempo'], rpm_motor, color='tab:pink')
        plt.title("RPM do Motor")
        plt.xlabel("Tempo [s]")
        plt.ylabel("RPM")
        plt.grid(True)

        # Torques do motor
        plt.subplot(2, 2, 4)
        plt.plot(motor_history['tempo'], motor_history['torque_eletrico'],
                label='Torque Gerado', color='blue')
        plt.plot(motor_history['tempo'], motor_history['torque_carga'],
                label='Torque de Carga', color='cyan', linestyle='--')
        plt.title("Torques no Motor")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Torque [Nm]")
        plt.legend()
        plt.grid(True)

        plt.tight_layout()
        plt.show()

        # Impressão dos parâmetros principais
        print("Tempo [s]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]")
        for i in range(0, len(performance), 150):
            t = performance[i]["tempo"]
            va = performance[i]["va"]
            vlk = performance[i]["vlk"]
            print(f"{t:.2f}\t\t{va:.2f}\t\t{vlk:.2f}")

        print("\nForça de Arrasto [N]\tForça Trativa Potencial [N]\tForça Trativa [N]")
        for i in range(0, len(performance), 150):
            fa = performance[i]["fa"]
            fp = performance[i]["fp"]
            ff = performance[i]["ff"]
            print(f"{fa:.2f}\t\t{fp:.2f}\t\t{ff:.2f}")

        # Obtenção dos vetores principais
        tempo = perf_dict["tempo"]
        velocidade_angular = perf_dict["va"]
        velocidade_linear = perf_dict["vlk"]

        # --- Figura 2: Velocidade Angular e Linear x Tempo ---
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        plt.plot(tempo, velocidade_angular, color='tab:red')
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade Angular [rad/s]")
        plt.title("Velocidade Angular x Tempo")
        plt.grid(True)

        plt.subplot(1, 2, 2)
        plt.plot(tempo, velocidade_linear, color='tab:purple')
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade Linear [km/h]")
        plt.title("Velocidade Linear x Tempo")
        plt.grid(True)

        plt.tight_layout()
        plt.show()

        # --- Figura 3: Velocidade Angular e Linear x RPM ---
        plt.figure(figsize=(12, 6))

        plt.subplot(1, 2, 1)
        plt.plot(rpm_motor[:len(tempo)], velocidade_angular, color='tab:red')
        plt.xlabel("RPM Motor")
        plt.ylabel("Velocidade Angular [rad/s]")
        plt.title("Velocidade Angular x RPM")
        plt.grid(True)

        plt.subplot(1, 2, 2)
        plt.plot(rpm_motor[:len(tempo)], velocidade_linear, color='tab:purple')
        plt.xlabel("RPM Motor")
        plt.ylabel("Velocidade Linear [km/h]")
        plt.title("Velocidade Linear x RPM")
        plt.grid(True)

        plt.tight_layout()
        plt.show()


def Instancias():
    """
    Cria instâncias do veículo e pneu e executa a simulação completa.
    """

    # Modelo do veículo
    dt_model = Drivetrain(
        cgx=853,                    # Centro de gravidade no eixo X [mm]
        cgy=294,                    # Centro de gravidade no eixo Y [mm]
        massa=347,                  # Massa do veículo [kg]
        massa_roda= 6,              # Massa da Roda [kg]    
        entre_eixos=1567,           # Entre-eixos [mm]
        raio_pneu=220,              # Raio do pneu [mm]
        reducao_primaria=2.12,      # Redução primária
        reducao_final=2.76,         # Redução final (diferencial)
        tempo_i=0,                  # Tempo inicial de simulação (s)
        tire_friction_coef=1.45,    # Coeficiente de fricção 
        tire_Fz = 850               # Carga vertical no pneu [N]
    )
    dt_model.tempo_f = 20

    # Simulação de desempenho
    performance_veiculo, variacao_tempo = dt_model.CarPerformance()
    dt_model.printCarPerformance(performance_veiculo)

    # Velocidades para cálculo do slip ratio
    velocidade_angular = np.array([dado["va"] for dado in performance_veiculo])
    velocidade_linear = np.array([dado["vlm"] for dado in performance_veiculo])

    # Slip ratio
    slip_ratio = Tire.SlipRatio(velocidade_angular, dt_model.raio_pneu * 0.001, velocidade_linear)

    # Modelo do pneu
    tire = Tire(
        tire_Fz= 850,               # Carga vertical no pneu [N]
        tire_Sa=0,                  # Ângulo de escorregamento lateral [rad]
        tire_Ls=slip_ratio,         # Slip ratio já calculado
        tire_friction_coef=1.45,    # Coeficiente de fricção
        tire_Ca=0                   # Ângulo de camber
    )

    result = [0.333, 1.627, 1, 4.396, 931.4, 366.4]
    _, _, tire_longitudinal_forces = tire.Tire_forces(result)

    # Gráficos de slip ratio
    tire.printSlipRatio(variacao_tempo, slip_ratio, tire_longitudinal_forces)


# ============================================================================
# Execução
# ============================================================================
if __name__ == "__main__":
    Instancias()
