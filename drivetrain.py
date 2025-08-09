"""
Simulação de desempenho de um veículo Fórmula SAE elétrico -
Autor: Marco Affonso e Igor Maia
"""

import numpy as np
import matplotlib.pyplot as plt
from pwt import Motor               
from tire2 import Tire

class Drivetrain:

    def __init__(self, cgx, cgy, massa, massa_roda, entre_eixos,
                 raio_pneu, reducao_primaria, reducao_final, tempo_i,
                 tire_friction_coef, tire_Fz):

        # geometric / mass
        self.cgx = cgx
        self.cgy = cgy
        self.massa = massa
        self.massa_roda = massa_roda
        self.entre_eixos = entre_eixos
        # raio em mm na entrada — guardo valor original e valor em m quando necessário
        self.raio_pneu = raio_pneu

        # transmissão
        self.reducao_primaria = reducao_primaria
        self.reducao_final = reducao_final

        # tempo
        self.tempo_i = tempo_i
        self.tempo_f = 0

        # pneus
        self.tire_friction_coef = tire_friction_coef
        self.tire_Fz = tire_Fz

        # motor (mantive os mesmos parâmetros do seu código original)
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
            tire_Fz=self.tire_Fz,
            tire_Sa=0,
            tire_Ls=1,  # inicial, atualizado em loop
            tire_friction_coef=self.tire_friction_coef,
            tire_Ca=0
        )

    def CarPerformance(self):
        """
        Roda a simulação usando EDOs acopladas para v(t) e omega_w(t).
        Retorna: performance (lista de dicts) e vetor de tempo (numpy array)
        """

        # --- parâmetros fixos / físicos ---
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        # converto raio mm -> m
        raio = self.raio_pneu * 0.001
        c_r = 0.012
        b = 0.01                     # coeficiente de atrito/atrito mecânico no rotor (mantido)
        eficiencia_transmissao = 0.95

        # número de rodas motrizes (RWD)
        num_roda_motriz = 2
        num_rodas = 4

        # carga por roda 
        Fz_por_roda = self.tire_Fz if self.tire_Fz is not None else (peso / num_rodas)

        # parâmetros Pacejka placeholder
        pacejka_params = [0.333, 1.627, 1, 4.396, 931.4, 366.4]

        # verificações de tempo
        if self.tempo_f == 0:
            raise ValueError("Defina o tempo final: dt_model.tempo_f = <valor>")

        # integração explícita Euler 
        dt = 0.001
        variacao_tempo = np.arange(self.tempo_i, self.tempo_f + dt/2, dt)

        # --- estados iniciais ---
        v_linear_ms = 0.0      # velocidade linear do CM (m/s)  -> estado 1
        omega_roda = 0.0       # velocidade angular da roda (rad/s) -> estado 2

        # inércia efetiva de cada roda 
        Iw = 0.5 * self.massa_roda * raio**2

        # numerics
        eps_v = 1e-3           # evita divisão por zero em slip ratio
        performance = []

        for tempo in variacao_tempo:
            # --- cálculo da redução total entre motor e roda ---
            i_total = self.reducao_primaria * self.reducao_final

            # --- forças resistivas ---
            forca_rolamento = c_r * peso
            forca_arrasto = 0.5 * densidade_ar * coeficiente_arrasto * area_frontal * v_linear_ms**2
            forca_resistiva_total = forca_rolamento + forca_arrasto

            # --- slip ratio dinâmico por roda ---
            denom_v = max(v_linear_ms, eps_v)
            slip_ratio_roda = (omega_roda * raio - v_linear_ms) / denom_v

            # --- consulta ao modelo de pneu (força longitudinal por roda) ---
            self.tire.tire_Fz = Fz_por_roda
            self.tire.tire_Ls = slip_ratio_roda
            # suponho que Tire_forces(params) retorna (Fy, Mz, Fx)
            _, _, Fx_roda = self.tire.Tire_forces(pacejka_params)

            Fx_roda = np.clip(Fx_roda, -self.tire_friction_coef*Fz_por_roda, self.tire_friction_coef*Fz_por_roda)

            # --- torque de carga no motor ---
            # torque por forças resistivas refletidas ao eixo do motor:
            torque_rr = forca_rolamento * raio
            torque_fa = forca_arrasto * raio
            torque_atr_mec = b * (self.motor.wm) 
            torque_carga = torque_rr + torque_fa + torque_atr_mec

            # torque_carga_motor é o torque no motor que representa a carga externa
            torque_carga_motor = torque_carga / (i_total * eficiencia_transmissao)

            self.motor.set_load(torque_carga_motor)
            torque_motor_gerado = self.motor.step(dt)   # Nm no eixo do motor

            # torque disponível na superfície da roda 
            torque_trativo_rodas_total = torque_motor_gerado * i_total * eficiencia_transmissao
            # torque por roda motriz
            torque_por_roda = torque_trativo_rodas_total / num_roda_motriz

            # --- EDOs: cálculo das derivadas ---
            # translacional:
            forca_trativa_real = num_roda_motriz * Fx_roda
            dvdt = (forca_trativa_real - forca_resistiva_total) / self.massa

            # rotacional (por roda): 
            domegadt = (torque_por_roda - raio * Fx_roda) / Iw

            # integração explícita (Euler)
            v_linear_ms += dvdt * dt
            omega_roda += domegadt * dt
    
            # armazenar histórico para saída
            self.motor.store_history(tempo)   # manter histórico do motor (sua função)
            performance.append({
                "tempo": float(tempo),
                "vlk": float(v_linear_ms * 3.6),     # km/h
                "vlm": float(v_linear_ms),           # m/s
                "va": float(omega_roda),             # rad/s
                "fa": float(forca_arrasto),
                "rr": float(forca_rolamento),
                "ff": float(forca_trativa_real),
                "fp": float(torque_por_roda / raio * num_roda_motriz),  # força trativa potencial aproximada
                "sr": float(slip_ratio_roda),
                "ft_roda": float(Fx_roda)
            })

        return performance, variacao_tempo

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
    Cria instâncias do veículo e pneu e executa a simulação completa.
    """

    # Modelo do veículo - parâmetros originais mantidos
    dt_model = Drivetrain(
        cgx=853,                    # Centro de gravidade no eixo X [mm]
        cgy=294,                    # Centro de gravidade no eixo Y [mm]
        massa=347,                  # Massa do veículo [kg]
        massa_roda=6,               # Massa da Roda [kg]
        entre_eixos=1567,           # Entre-eixos [mm]
        raio_pneu=220,              # Raio do pneu [mm]
        reducao_primaria=2.12,      # Redução primária
        reducao_final=2.76,         # Redução final (diferencial)
        tempo_i=0,                  # Tempo inicial de simulação (s)
        tire_friction_coef=1.45,    # Coeficiente de fricção 
        tire_Fz=850                 # Carga vertical no pneu [N]
    )
    dt_model.tempo_f = 20

    # Simulação de desempenho (agora com EDOs acopladas)
    performance_veiculo, variacao_tempo = dt_model.CarPerformance()
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

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    Instancias()
