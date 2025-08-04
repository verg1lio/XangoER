import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control as clt
from pwt import Motor

class Drivetrain:

    def __init__(self, cgx, cgy, massa, massa_roda, entre_eixos, raio_pneu, reducao_primaria, reducao_final, cp, tempo_i):
        self.cgx = cgx
        self.cgy = cgy
        self.massa = massa
        self.massa_roda = massa_roda
        self.entre_eixos =entre_eixos
        self.raio_pneu = raio_pneu
        self.reducao_primaria = reducao_primaria
        self.reducao_final = reducao_final
        self.cp = cp
        self.tempo_i = tempo_i
        self.tempo_f = 0

        # Criar motor corretamente com os parâmetros
        self.motor = Motor(
            rs=0.04585,      # Resistência do estator
            ld=0.00067,      # Indutância d
            lq=0.00067,      # Indutância q
            jm=0.05769,      # Inércia do rotor
            kf=0.1,          # Coeficiente de atrito
            lambda_m=0.13849, # Fluxo do ímã
            p=10,            # Nº pares de polos
            valor_mu=0.9
        )

    # <<< OBSERVAÇÃO: Este método está definido mas não é utilizado no fluxo atual do código.
    def Reduções(self, eficiencia_corrente=0.96):

        self.motor.speed_ref = 600
        self.motor.simulate()

        rpm_motor = np.array(self.motor.velocidade) * 60 / (2 * np.pi)
        torque_motor = np.array(self.motor.torque_mecanico)
        potencia_motor = torque_motor * self.motor.wm

        rpm_pos_corrente1 = rpm_motor / self.reducao_primaria
        torque_pos_corrente1 = torque_motor * self.reducao_primaria
        potencia_pos_corrente1 = potencia_motor

        rpm_pos_corrente = rpm_pos_corrente1
        torque_pos_corrente = torque_pos_corrente1 * eficiencia_corrente
        potencia_pos_corrente = potencia_pos_corrente1 * eficiencia_corrente

        rpm_dif = rpm_pos_corrente / self.reducao_final
        torque_dif = torque_pos_corrente * self.reducao_final
        potencia_dif = potencia_pos_corrente

        rendimento_transmissao = (potencia_dif / (potencia_motor + 1e-9))
        torque_ef = torque_dif * rendimento_transmissao

        return [torque_ef, potencia_dif, rpm_dif]
        
    def CarPerformance(self):
        # Parâmetros do veículo
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        raio = self.raio_pneu * 0.001
        c_r = 0.012  # Coeficiente de rolamento
        b = 0.01
        eficiencia_transmissao = 0.95

        if self.tempo_f == 0:
            raise ValueError("Defina o tempo final: dt_model.tempo_f = <valor>")

        dt = 0.001
        variacao_tempo = np.arange(self.tempo_i, self.tempo_f, dt)
        
        performance = []
        v_linear_ms = 0.0
        v_angular = 0.0

        for tempo in variacao_tempo:
            # Inércias
            J_roda = 0.5 * self.massa_roda * raio**2        # inércia de uma roda         
            J_carro_equivalente = self.massa * raio**2      # Inércia equivalente da massa do veículo transferida ao eixo
           
            J_total = ((4 * J_roda) + J_carro_equivalente)  # Inércia total sentida pelo motor/trenó de força

            # 1. Calcular forças resistivas com base na velocidade atual
            forca_rr = c_r * peso
            forca_arrasto = 0.5 * densidade_ar * v_linear_ms**2 * coeficiente_arrasto * area_frontal

            forca_resistiva_total = forca_rr + forca_arrasto
            
            # 2. Calcular torque de carga no eixo do motor
            torque_rr = forca_rr * raio
            torque_fa = forca_arrasto * raio
            torque_atr_mec = b * v_angular

            torque_carga = torque_rr + torque_fa + torque_atr_mec
            # Relação de transmissão (correntes e diferencial)
            i_total = self.reducao_primaria * self.reducao_final
            # Torque de carga transferido para o motor
            torque_carga_motor = torque_carga / (i_total * eficiencia_transmissao)
            
            self.motor.set_load(torque_carga_motor)
            torque_motor_gerado = self.motor.step(dt)
            
            torque_trativo_rodas = torque_motor_gerado * i_total * eficiencia_transmissao
            forca_trativa = torque_trativo_rodas / raio
            
            forca_liquida = forca_trativa - forca_resistiva_total
            a_linear = forca_liquida / self.massa
            v_linear_ms += a_linear * dt
            
            if v_linear_ms < 0:
                v_linear_ms = 0

            # <<< CORREÇÃO: Adicionando as chaves "vlm" e "va" que eram necessárias depois.
            # A velocidade angular (va) aqui é a da RODA.
            velocidade_angular_roda = self.motor.wm / i_total

            self.motor.store_history(tempo)
            performance.append({
                "tempo": tempo,
                "vlk": v_linear_ms * 3.6,
                "vlm": v_linear_ms, 
                "va": velocidade_angular_roda,
                "fa": forca_arrasto,
                "rr": forca_rr,
                "ff": forca_trativa
            })

        # <<< CORREÇÃO: Retornando os 2 valores que serão utilizados.
        return performance, variacao_tempo

    def printCarPerformance(self, performance):
        motor_history = self.motor.history
        
        # Converte a lista de dicionários para um dicionário de listas para facilitar a plotagem
        perf_dict = {k: [dic[k] for dic in performance] for k in performance[0]}

        plt.figure(figsize=(14, 10))
        
        plt.subplot(2, 2, 1)
        plt.plot(perf_dict['tempo'], perf_dict['vlk'], color='purple')
        plt.title("Velocidade Linear do Veículo")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade [km/h]")
        plt.grid(True)
        
        plt.subplot(2, 2, 2)
        plt.plot(perf_dict['tempo'], perf_dict['ff'], label='Força Trativa', color='green')
        plt.plot(perf_dict['tempo'], perf_dict['fa'], label='Força de Arrasto', color='red')
        plt.plot(perf_dict['tempo'], perf_dict['rr'], label='Resist. Rolamento', color='orange')
        plt.title("Forças Atuantes")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Força [N]")
        plt.legend()
        plt.grid(True)

        plt.subplot(2, 2, 3)
        rpm_motor = np.array(motor_history['velocidade']) * 60 / (2 * np.pi)
        plt.plot(motor_history['tempo'], rpm_motor, color='tab:pink')
        plt.title("RPM do Motor")
        plt.xlabel("Tempo [s]")
        plt.ylabel("RPM")
        plt.grid(True)

        plt.subplot(2, 2, 4)
        plt.plot(motor_history['tempo'], motor_history['torque_eletrico'], label='Torque Gerado', color='blue')
        plt.plot(motor_history['tempo'], motor_history['torque_carga'], label='Torque de Carga', color='cyan', linestyle='--')
        plt.title("Torques no Motor")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Torque [Nm]")
        plt.legend()
        plt.grid(True)

        plt.tight_layout()
        plt.show()

# A classe Tire não precisa de alterações
class Tire:
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=0, B0=0, B1=0, B2=1, B3=1, omega=315, slip_angle_start=-9, slip_angle_end=9, angle_step=0.5, WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
        self.tire_Fz = tire_Fz; self.tire_Sa = tire_Sa; self.tire_Ls = tire_Ls; self.tire_type = 'Default'; self.tire_friction_coef = tire_friction_coef; self.tire_Ca = tire_Ca; self.L0 = B0; self.L1 = B1; self.L2 = B2; self.L3 = B3; self.alpha = np.radians(omega); self.theta2 = np.arange(slip_angle_start + 90, slip_angle_end + 91, angle_step); self.theta2_rad = np.radians(self.theta2); self.angle = self.theta2_rad[0]; self.AC = []; self.beta = []; self.psi = []; self.lamda = []; self.theta3 = []; self.theta4 = []; self.Ox, self.Oy = 0, 0; self.Ax, self.Ay = [], []; self.Bx, self.By = [], []; self.Cx, self.Cy = B0, 0; self.w = []; self.om2, self.om4 = [], []; self.alpha_dot = []; self.outer_slip = []; self.inner_slip = []; self.static_slip_angle = None; self.WB = WB; self.rear_axle_length = rear_axle_length; self.rear_axle_center = (B0 / 2, 0); self.track_y = track_y; self.tire_k = tire_k
    def Tire_forces(self, params):
        E, Cy, Cx, Cz, c1, c2 = params; Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2)); D = self.tire_friction_coef * self.tire_Fz; Bz = Cs / (Cz * D) if (Cz * D) != 0 else 0; Bx = Cs / (Cx * D) if (Cx * D) != 0 else 0; By = Cs / (Cy * D) if (Cy * D) != 0 else 0; tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa)))); tire_longitudinal_force = D * np.sin(Cx * np.arctan(9 * Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls)))); tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa)))); camber_thrust = D * np.sin(Cy * np.arctan(By * self.tire_Ca)); return tire_lateral_force + 0.5 * camber_thrust, (10 + (tire_auto_align_moment / 55)), tire_longitudinal_force
    @staticmethod
    def SlipRatio(velocidade_angular, raio_pneu, velocidade_linear):
        # Adicionar proteção contra divisão por zero
        return (velocidade_angular * raio_pneu / (velocidade_linear + 1e-6)) - 1
    def printSlipRatio(self, tempo, slip_ratio, tire_longitudinal_forces):
        passo = 10; slip_ratio_constante = max(slip_ratio) if len(slip_ratio) > 0 else 0; i = 0
        while i < len(tempo) and tempo[i] <= 0.05: slip_ratio[i] = slip_ratio_constante; i += 1
        print("Valores do Slip Ratio:");
        for i in range(0, len(slip_ratio), passo): print(f"{slip_ratio[i]:.2f}")
        plt.figure(figsize=(14, 6)); plt.subplot(1, 2, 1); plt.plot(tempo, slip_ratio, label='Slip Ratio', color='blue'); plt.xlabel('Tempo [s]'); plt.ylabel('Slip Ratio'); plt.title('Slip Ratio x Tempo'); plt.grid(True); plt.legend(); plt.subplot(1, 2, 2); plt.plot(slip_ratio, tire_longitudinal_forces, label='Força Longitudinal', color='green'); plt.xlabel('Slip Ratio'); plt.ylabel('Força Longitudinal (N)'); plt.title('Força Longitudinal x Slip Ratio'); plt.grid(True); plt.legend(); plt.tight_layout(); plt.show()
        
def Instancias():
    dt_model = Drivetrain(
        cgx=853, cgy=294, massa=347, massa_roda= 6, entre_eixos=1567, 
        raio_pneu=220, reducao_primaria=2.12, reducao_final=2.76, 
        cp=2.22, tempo_i=0
    )

    dt_model.tempo_f = 20

    # <<< CORREÇÃO: A chamada agora espera e recebe 2 valores.
    performance_veiculo, variacao_tempo = dt_model.CarPerformance()

    # <<< CORREÇÃO: Passando a lista de dicionários diretamente.
    dt_model.printCarPerformance(performance_veiculo)

    # Extração das velocidades dos dados de performance
    # O código agora encontra as chaves "va" e "vlm"
    velocidade_angular = np.array([dado["va"] for dado in performance_veiculo])
    velocidade_linear = np.array([dado["vlm"] for dado in performance_veiculo])
 
    # Cálculo do slip ratio
    slip_ratio = Tire.SlipRatio(velocidade_angular, dt_model.raio_pneu * 0.001, velocidade_linear)
        
    slip_model = Tire(
        tire_Fz= 851, tire_Sa=0, tire_Ls=slip_ratio, 
        tire_friction_coef=1.45, tire_Ca=0
    )

    result = [0.333, 1.627, 1, 4.396, 931.4, 366.4]

    tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = slip_model.Tire_forces(result)

    slip_model.printSlipRatio(variacao_tempo, slip_ratio, tire_longitudinal_forces)
    
Instancias()
