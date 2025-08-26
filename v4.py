"""
Simulação de desempenho de um veículo Fórmula SAE elétrico -

Autor: Marco Affonso e Igor Maia 
"""

import numpy as np
import matplotlib.pyplot as plt
from pwt import Motor               
from tire import Tire

class Drivetrain:

    def __init__(self, massa, massa_roda,
                 raio_pneu, i, tempo_i,
                 tire_friction_coef, tire_Fz):

        # geometric / mass
        self.massa = massa
        self.massa_roda = massa_roda
        self.raio_pneu = raio_pneu

        # transmissão
        self.i = i
        
        # tempo
        self.tempo_i = tempo_i

        # pneus
        self.tire_friction_coef = tire_friction_coef
        self.tire_Fz = tire_Fz

        # motor 
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

        # instância base do pneu (será atualizada a cada passo com Ls e Fz)
        self.tire = Tire(
            tire_Fz=850,
            tire_Sa=0,
            tire_Ls=1,  # inicial, atualizado em loop
            tire_friction_coef=self.tire_friction_coef,
            tire_Ca=0
        )

    def CarPerformanceDistancia(self, distancia_final=75.0):
        """
        Roda a simulação até o veículo percorrer 'distancia_final' [m].
        Retorna: performance (lista de dicts), vetor de tempo (numpy array), tempo_final
        """

        # --- parâmetros fixos / físicos ---
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06

        # raio em metros
        raio = self.raio_pneu * 0.001
        c_r = 0.012
        b = 0.01
        eficiencia_transmissao = 0.95

        num_roda_motriz = 2
        num_rodas = 4
        Fz_por_roda = self.tire_Fz if self.tire_Fz is not None else (peso / num_rodas)

        pacejka_params = [0.333, 1.627, 1, 4.396, 931.4, 366.4]

        # integração explícita Euler 
        dt = 0.001
        tempo = self.tempo_i

        # estados iniciais
        v_linear_ms = 0.0
        omega_roda = 0.0
        deslocamento = 0.0   # posição [m]

        Jr = 0.5 * self.massa_roda * raio**2
        Jc =  0.5 * self.massa * raio**2
        Jt = Jr + Jc

        eps_v = 1e-3
        performance = []

        while deslocamento < distancia_final:
            

            # forças resistivas
            forca_rolamento = c_r * peso
            forca_arrasto = 0.5 * densidade_ar * coeficiente_arrasto * area_frontal * v_linear_ms**2
            forca_resistiva_total = forca_rolamento + forca_arrasto

            # slip ratio
            denom_v = max(v_linear_ms, eps_v)
            slip_ratio_roda = (omega_roda * raio - v_linear_ms) / denom_v

            # modelo do pneu
            self.tire.tire_Fz = Fz_por_roda
            self.tire.tire_Ls = slip_ratio_roda
            _, _, Fx_roda = self.tire.Tire_forces(pacejka_params)

            # torque de carga no motor
            torque_rr = forca_rolamento * raio
            torque_fa = forca_arrasto * raio
            torque_atr_mec = b * (self.motor.wm) 
            torque_carga = torque_rr + torque_fa + torque_atr_mec
            torque_carga_motor = torque_carga / (self.i * eficiencia_transmissao)

            self.motor.set_load(torque_carga_motor)
            torque_motor_gerado = self.motor.step(dt)

            # torque na roda
            torque_trativo_rodas_total = torque_motor_gerado * self.i * eficiencia_transmissao
            torque_por_roda = torque_trativo_rodas_total / num_roda_motriz

            # EDOs
            forca_trativa_real = num_roda_motriz * Fx_roda
            dvdt = (forca_trativa_real - forca_resistiva_total) / self.massa
            domegadt = (torque_por_roda - raio * Fx_roda) / Jt

            # integração
            v_linear_ms += dvdt * dt
            omega_roda += domegadt * dt
            deslocamento += v_linear_ms * dt
            tempo += dt

            # salvar histórico
            self.motor.store_history(tempo)
            performance.append({
                "tempo": float(tempo),
                "s": float(deslocamento),           # posição acumulada [m]
                "vlk": float(v_linear_ms * 3.6),    # km/h
                "vlm": float(v_linear_ms),          # m/s
                "va": float(omega_roda),            # rad/s
                "fa": float(forca_arrasto),
                "rr": float(forca_rolamento),
                "ff": float(forca_trativa_real),
                "fp": float(torque_por_roda / raio * num_roda_motriz),
                "sr": float(slip_ratio_roda),
                "ft_roda": float(Fx_roda)
            })

        return performance, np.array([p["tempo"] for p in performance]), tempo

    # ------------------------------------------------------------------------
    def printCarPerformance(self, performance):
        """
        Gera gráficos de desempenho do veículo com base no histórico do motor
        e nos dados de simulação armazenados em performance.
        """

        if len(performance) == 0:
            print("Sem dados de performance para plotar.")
            return

        motor_history = self.motor.history
        perf_dict = {k: [dic[k] for dic in performance] for k in performance[0]}

        # --- Figura 1: Resumo geral ---
        plt.figure(figsize=(14, 10))

        # Velocidade linear
        plt.subplot(2, 2, 1)
        plt.plot(perf_dict['tempo'], perf_dict['vlk'])
        plt.title("Velocidade Linear do Veículo")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade [km/h]")
        plt.grid(True)

        # Forças
        plt.subplot(2, 2, 2)
        plt.plot(perf_dict['tempo'], perf_dict['ff'], label='Força Trativa')
        plt.plot(perf_dict['tempo'], perf_dict['fa'], label='Força de Arrasto')
        plt.plot(perf_dict['tempo'], perf_dict['rr'], label='Resist. Rolamento')
        plt.title("Forças Atuantes")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Força [N]")
        plt.legend()
        plt.grid(True)

        # RPM do motor
        plt.subplot(2, 2, 3)
        rpm_motor = np.array(motor_history['velocidade']) * 60 / (2 * np.pi)
        plt.plot(motor_history['tempo'], rpm_motor)
        plt.title("RPM do Motor")
        plt.xlabel("Tempo [s]")
        plt.ylabel("RPM")
        plt.grid(True)

        # Torques do motor
        plt.subplot(2, 2, 4)
        plt.plot(motor_history['tempo'], motor_history['torque_eletrico'], label='Torque Gerado')
        plt.plot(motor_history['tempo'], motor_history['torque_carga'], label='Torque de Carga', linestyle='--')
        plt.title("Torques no Motor")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Torque [Nm]")
        plt.legend()
        plt.grid(True)

        plt.tight_layout()
        plt.show()

         # Impressões resumidas
        print("Tempo [s]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]")
        for i in range(0, len(performance), max(1, len(performance)//10)):
            t = performance[i]["tempo"]
            va = performance[i]["va"]
            vlk = performance[i]["vlk"]
            print(f"{t:.2f}\t\t{va:.2f}\t\t{vlk:.2f}")

        print("\nForça de Arrasto [N]\tForça Trativa Potencial [N]\tForça Trativa [N]")
        for i in range(0, len(performance), max(1, len(performance)//10)):
            fa = performance[i]["fa"]
            fp = performance[i]["fp"]
            ff = performance[i]["ff"]
            print(f"{fa:.2f}\t\t{fp:.2f}\t\t{ff:.2f}")

        # Plots secundários: velocidade angular e linear x tempo
        tempo = perf_dict["tempo"]
        velocidade_angular = perf_dict["va"]
        velocidade_linear = perf_dict["vlk"]

        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.plot(tempo, velocidade_angular)
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade Angular [rad/s]")
        plt.title("Velocidade Angular x Tempo")
        plt.grid(True)

        plt.subplot(1, 2, 2)
        plt.plot(tempo, velocidade_linear)
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade Linear [km/h]")
        plt.title("Velocidade Linear x Tempo")
        plt.grid(True)

        plt.tight_layout()
        plt.show()


# ---------------------------------------------------------------------------
def Instancias():
    """
    Cria instâncias do veículo e pneu e executa a simulação completa até 75 m.
    """

    # Modelo do veículo - parâmetros originais mantidos
    dt_model = Drivetrain(
     
        massa=347,                  # Massa do veículo [kg]
        massa_roda=6,               # Massa da Roda [kg]
        raio_pneu=220,              # Raio do pneu [mm]
        i=5,                     # Redução total
        tempo_i=0,                  # Tempo inicial de simulação (s)
        tire_friction_coef=1.45,    # Coeficiente de fricção 
        tire_Fz=850                 # Carga vertical no pneu [N]
    )

    # Simulação de desempenho até 75 m
    performance_veiculo, variacao_tempo, tempo_final = dt_model.CarPerformanceDistancia(distancia_final=75)
    dt_model.printCarPerformance(performance_veiculo)

    # Opcional: calcular slip ratio final e forças para plot separado (vetorizado)
    velocidade_angular = np.array([dado["va"] for dado in performance_veiculo])
    velocidade_linear = np.array([dado["vlm"] for dado in performance_veiculo])
    slip_ratio = Tire.SlipRatio(velocidade_angular, dt_model.raio_pneu * 0.001, velocidade_linear)

    # Modelo do pneu para o vetor slip_ratio (apenas para gráficos)
    tire = Tire(
        tire_Fz=850,
        tire_Sa=0,
        tire_Ls=slip_ratio,
        tire_friction_coef=1.45,
        tire_Ca=0
    )
    result = [0.333, 1.627, 1, 4.396, 931.4, 366.4]
    _, _, tire_longitudinal_forces = tire.Tire_forces(result)

    # Gráficos de slip ratio (usa função existente da classe Tire)
    tire.printSlipRatio(variacao_tempo, slip_ratio, tire_longitudinal_forces)

    print(f"Tempo final para percorrer 75 m: {tempo_final:.3f} s")

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    Instancias()
