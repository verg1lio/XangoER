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


# ============================================================================
# Classe Drivetrain: modelo do sistema de transmissão e veículo
# ============================================================================
class Drivetrain:

    def __init__(self, cgx, cgy, massa, massa_roda, entre_eixos,
                 raio_pneu, reducao_primaria, reducao_final, tempo_i):
        
        #Inicializa os parâmetros do veículo

        self.cgx = cgx                             # cgx, cgy: coordenadas do centro de gravidade                   
        self.cgy = cgy                       
        self.massa = massa                         # massa: massa total do veículo (kg)
        self.massa_roda = massa_roda               # massa_roda: massa de cada roda (kg)
        self.entre_eixos = entre_eixos             # entre_eixos: distância entre eixos (mm)
        self.raio_pneu = raio_pneu                 # raio_pneu: raio dinâmico do pneu (mm)
        self.reducao_primaria = reducao_primaria   # reducao_primaria: relação de redução da corrente
        self.reducao_final = reducao_final         # reducao_final: relação de redução do diferencial
        self.tempo_i = tempo_i                     # tempo_i: tempo inicial de simulação (s)
        self.tempo_f = 0

       
        # Inicializa motor com parâmetros fixos
        
        self.motor = Motor(
            rs=0.04585,                            # Resistência do estator
            ld=0.00067,                            # Indutância d
            lq=0.00067,                            # Indutância q
            jm=0.05769,                            # Inércia do rotor
            kf=0.1,                                # Coeficiente de atrito
            lambda_m=0.13849,                      # Fluxo do ímã
            p=10,                                  # Número de pares de polos
            valor_mu=0.9
        )

  
    def CarPerformance(self):
        """
        Simula o desempenho do veículo ao longo do tempo.
        Retorna:
        - performance: lista de dicionários com resultados
        - variacao_tempo: vetor de tempo
        """

        # Parâmetros físicos do veículo
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        raio = self.raio_pneu * 0.001  # mm → m
        c_r = 0.012                    # Coef. de rolamento
        b = 0.01                       # Atrito mecânico
        eficiencia_transmissao = 0.95

        if self.tempo_f == 0:
            raise ValueError("Defina o tempo final: dt_model.tempo_f = <valor>")

        # Discretização temporal
        dt = 0.001
        variacao_tempo = np.arange(self.tempo_i, self.tempo_f, dt)

        # Condições iniciais
        performance = []
        a_linear =0
        v_linear_ms = 0.0

     
        for tempo in variacao_tempo:
          

            # Relação de Transmissão total
            i_total = self.reducao_primaria * self.reducao_final

            # Forças resistivas
           
            forca_rolamento = c_r * peso
            forca_arrasto = 0.5 * densidade_ar * v_linear_ms**2 * coeficiente_arrasto * area_frontal
            forca_inercia = self.massa * a_linear 
            forca_resistiva_total = forca_rolamento + forca_arrasto + forca_inercia

            # Torque resistivo (convertido para motor)
            torque_i = forca_inercia * raio
            torque_rr = forca_rolamento * raio
            torque_fa = forca_arrasto * raio
            torque_atr_mec = b * (self.motor.wm)

            torque_carga = torque_rr + torque_fa + torque_atr_mec + torque_i
          
            torque_carga_motor = torque_carga / (i_total * eficiencia_transmissao)

            # Simulação do motor
            self.motor.set_load(torque_carga_motor)
            torque_motor_gerado = self.motor.step(dt)

            # Força trativa nas rodas
            torque_trativo_rodas = torque_motor_gerado * i_total * eficiencia_transmissao
            forca_trativa = torque_trativo_rodas / raio
            
            # Dinâmica do veículo
            forca_liquida = forca_trativa - forca_resistiva_total
            a_linear = np.abs(forca_liquida) / self.massa
            v_linear_ms += a_linear * dt
            
            # Velocidade angular da roda
            velocidade_angular_roda = self.motor.wm / i_total

            # Armazena histórico
            self.motor.store_history(tempo)
            performance.append({
                "tempo": tempo,
                "vlk": v_linear_ms * 3.6,   # km/h
                "vlm": v_linear_ms,         # m/s
                "va": velocidade_angular_roda,
                "fa": forca_arrasto,
                "rr": forca_rolamento,
                "ff": forca_trativa
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
        for i in range(0, len(performance), 10):
            t = performance[i]["tempo"]
            va = performance[i]["va"]
            vlk = performance[i]["vlk"]
            print(f"{t:.2f}\t\t{va:.2f}\t\t{vlk:.2f}")

        print("\nForça de Arrasto [N]\tResistência ao Rolamento [N]\tForça Trativa [N]")
        for i in range(0, len(performance), 10):
            fa = performance[i]["fa"]
            rr = performance[i]["rr"]
            ff = performance[i]["ff"]
            print(f"{fa:.2f}\t\t{rr:.2f}\t\t{ff:.2f}")

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

# ============================================================================
# Classe Tire: modelo simplificado de pneu
# ============================================================================
class Tire:
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None,
                 tire_friction_coef=None, tire_Ca=0,
                 B0=0, B1=0, B2=1, B3=1, omega=315,
                 slip_angle_start=-9, slip_angle_end=9, angle_step=0.5,
                 WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
        """
        Inicializa parâmetros do pneu para simulação de forças e slip ratio.
        """
        self.tire_Fz = tire_Fz
        self.tire_Sa = tire_Sa
        self.tire_Ls = tire_Ls
        self.tire_type = 'Default'
        self.tire_friction_coef = tire_friction_coef
        self.tire_Ca = tire_Ca

        # Constantes adicionais e variáveis geométricas
        self.L0 = B0
        self.L1 = B1
        self.L2 = B2
        self.L3 = B3
        self.alpha = np.radians(omega)
        self.theta2 = np.arange(slip_angle_start + 90, slip_angle_end + 91, angle_step)
        self.theta2_rad = np.radians(self.theta2)
        self.angle = self.theta2_rad[0]
        self.AC = []
        self.beta = []
        self.psi = []
        self.lamda = []
        self.theta3 = []
        self.theta4 = []
        self.Ox, self.Oy = 0, 0
        self.Ax, self.Ay = [], []
        self.Bx, self.By = [], []
        self.Cx, self.Cy = B0, 0
        self.w = []
        self.om2, self.om4 = [], []
        self.alpha_dot = []
        self.outer_slip = []
        self.inner_slip = []
        self.static_slip_angle = None
        self.WB = WB
        self.rear_axle_length = rear_axle_length
        self.rear_axle_center = (B0 / 2, 0)
        self.track_y = track_y
        self.tire_k = tire_k

  
    def Tire_forces(self, params):
        """
        Calcula forças do pneu (laterais, longitudinais e momento de alinhamento)
        com base em um modelo simplificado de Pacejka.
        """
        E, Cy, Cx, Cz, c1, c2 = params
        Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))
        D = self.tire_friction_coef * self.tire_Fz

        Bz = Cs / (Cz * D) if (Cz * D) != 0 else 0
        Bx = Cs / (Cx * D) if (Cx * D) != 0 else 0
        By = Cs / (Cy * D) if (Cy * D) != 0 else 0

        tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
        tire_longitudinal_force = D * np.sin(Cx * np.arctan(9 * Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
        tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))
        camber_thrust = D * np.sin(Cy * np.arctan(By * self.tire_Ca))

        return (tire_lateral_force + 0.5 * camber_thrust,
                (10 + (tire_auto_align_moment / 55)),
                tire_longitudinal_force)


    @staticmethod
    def SlipRatio(velocidade_angular, raio_pneu, velocidade_linear):
        """
        Calcula o slip ratio: diferença relativa entre
        velocidade da roda e velocidade do veículo.
        """
        return (velocidade_angular * raio_pneu / (velocidade_linear + 1e-6)) - 1


    def printSlipRatio(self, tempo, slip_ratio, tire_longitudinal_forces):
        """
        Plota gráficos do slip ratio e imprime alguns valores no console.
        """
        passo = 10
        slip_ratio_constante = max(slip_ratio) if len(slip_ratio) > 0 else 0

        # Mantém slip ratio constante no início (0-0.05s)
        i = 0
        while i < len(tempo) and tempo[i] <= 0.05:
            slip_ratio[i] = slip_ratio_constante
            i += 1

        print("Valores do Slip Ratio:")
        for i in range(0, len(slip_ratio), passo):
            print(f"{slip_ratio[i]:.2f}")

        # Gráficos
        plt.figure(figsize=(14, 6))

        # Slip ratio vs tempo
        plt.subplot(1, 2, 1)
        plt.plot(tempo, slip_ratio, label='Slip Ratio', color='blue')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Slip Ratio')
        plt.title('Slip Ratio x Tempo')
        plt.grid(True)
        plt.legend()

        # Força longitudinal vs slip ratio
        plt.subplot(1, 2, 2)
        plt.plot(slip_ratio, tire_longitudinal_forces, label='Força Longitudinal', color='green')
        plt.xlabel('Slip Ratio')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal x Slip Ratio')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()


# ============================================================================
# Função principal de execução
# ============================================================================
def Instancias():
    """
    Cria instâncias do veículo e pneu e executa a simulação completa.
    """
    # Modelo do veículo
    dt_model = Drivetrain(
        cgx=853, cgy=294, massa=340, massa_roda=6,
        entre_eixos=1567, raio_pneu=220,
        reducao_primaria=2.12, reducao_final=2.76, tempo_i=0
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
    slip_model = Tire(
        tire_Fz=851, tire_Sa=0, tire_Ls=slip_ratio,
        tire_friction_coef=1.45, tire_Ca=0
    )
    result = [0.333, 1.627, 1, 4.396, 931.4, 366.4]
    _, _, tire_longitudinal_forces = slip_model.Tire_forces(result)

    # Gráficos de slip ratio
    slip_model.printSlipRatio(variacao_tempo, slip_ratio, tire_longitudinal_forces)


# ============================================================================
# Execução
# ============================================================================
if __name__ == "__main__":
    Instancias()
