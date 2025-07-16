import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
import control as ctl

class Motor:

    def __init__(self, rs, rr, ls, lr, mrs, jm, kf):
        # Constants
        self.pi23 = 2 * np.pi / 3
        self.rq23 = np.sqrt(2 / 3)
        self.rq3 = np.sqrt(3)

        # Machine parameters
        self.rs = rs  # Stator resistance (ohms)
        self.rr = rr  # Rotor resistance (ohms)
        self.ls = ls  # Stator inductance (henries)
        self.lr = lr  # Rotor inductance (henries)
        self.msr = mrs  # Mutual inductance between stator and rotor (henries)
        self.lso = 0.1 * self.ls  # Stator leakage inductance (henries)
        
        self.jm = jm  # Moment of inertia (kg*m^2)
        self.kf = kf  # Friction coefficient (N*m*s)
        self.cte_tempo_mec = self.jm / self.kf  # Mechanical time constant (s)
        self.idt = 1 / (self.ls * self.lr - self.msr * self.msr)  # Inverse of the determinant
        self.p = 2  # Number of pole pairs
        self.amsr = self.p * self.idt * self.msr  # Constant for torque calculation
        self.J_eff = 0         # inércia efetiva (transmissão + rodas) refletida ao motor
        self.prev_wm = 0.0     # velocidade angular anterior (rad/s)

        # Simulation parameters
        self.h = 1.e-5  # Time step (s)
        self.tmax = 1  # Maximum simulation time (s)
        self.hp = self.tmax / 2000  # Plotting time step (s)
        if self.hp < self.h:
            self.hp = self.h

        # Initial conditions
        self.reset_initial_conditions()
        
        # Storage for output variables
        self.tempo = []  # Time (s)
        self.corrented = []  # Direct-axis current (A)
        self.correnteq = []  # Quadrature-axis current (A)
        self.corrente1 = []  # Phase 1 current (A)
        self.corrente2 = []  # Phase 2 current (A)
        self.corrente3 = []  # Phase 3 current (A)
        self.tensao1 = []  # Phase 1 voltage (V)
        self.tensao2 = []  # Phase 2 voltage (V)
        self.tensao3 = []  # Phase 3 voltage (V)
        self.tensaosd = []  # Direct-axis voltage (V)
        self.tensaosq = []  # Quadrature-axis voltage (V)
        self.fluxord = []  # Direct-axis rotor flux (Wb)
        self.fluxorq = []  # Quadrature-axis rotor flux (Wb)
        self.fluxos1 = []  # Phase 1 stator flux (Wb)
        self.fluxos2 = []  # Phase 2 stator flux (Wb)
        self.fluxos3 = []  # Phase 3 stator flux (Wb)
        self.fluxosd = []  # Direct-axis stator flux (Wb)
        self.fluxosq = []  # Quadrature-axis stator flux (Wb)
        self.fluxos = []   # Zero-sequence stator flux (Wb)
        self.conjugado = []  # Electromagnetic torque (N*m)
        self.velocidade = []  # Mechanical speed (rad/s)
        self.frequencia = []  # Electrical frequency (rad/s)
        self.conjcarga = []  # Load torque (N*m)
        self.correnteo = []  # Zero-sequence current (A)
        self.torque_mecanico = []  # Mechanical torque (N*m)

    def reset_initial_conditions(self):
        # Initialize conditions
        self.cl = 0  # Load torque (N*m)
        self.wm = 0.0  # Mechanical speed (rad/s)
        self.prev_wm = 0.0    # inicializa velocidade anterior
        self.t = 0  # Time (s)
        self.tp = 0  # Plotting time (s)
        self.j = 0  # Plotting index
        self.ce = 0  # Electromagnetic torque (N*m)
        self.ws = 377  # Synchronous speed (rad/s)
        self.Vsm = 220 * np.sqrt(2)  # Peak stator voltage (V)
        self.Vs = self.Vsm  # Stator voltage (V)
        self.tete = 0  # Electrical angle (rad)
        self.fsd = 0  # Direct-axis stator flux (Wb)
        self.fsq = 0  # Quadrature-axis stator flux (Wb)
        self.frd = 0  # Direct-axis rotor flux (Wb)
        self.frq = 0  # Quadrature-axis rotor flux (Wb)
        self.isd = 0  # Direct-axis stator current (A)
        self.isq = 0  # Quadrature-axis stator current (A)
        self.ird = 0  # Direct-axis rotor current (A)
        self.irq = 0  # Quadrature-axis rotor current (A)
        self.iso = 0  # Zero-sequence stator current (A)
        self.rg = 0  # Rotor angle (rad)

    def source_voltage(self,):
        self.tete += self.h * self.ws
        if self.tete >= 2 * np.pi:
            self.tete -= 2 * np.pi
            
        vs1 = self.Vs * np.cos(self.tete)
        vs2 = self.Vs * np.cos(self.tete - self.pi23)
        vs3 = self.Vs * np.cos(self.tete + self.pi23)
        return vs1, vs2, vs3

    def load_torque(self,):
            if self.t >= self.tmax / 2:
                self.cl = 10        
    
    def direct_voltage_and_quadrature(self, vs1, vs2, vs3):
        vsd = self.rq23 * (vs1 - vs2 / 2 - vs3 / 2)
        vsq = self.rq23 * (vs2 * self.rq3 / 2 - vs3 * self.rq3 / 2)
        vso = (vs1 + vs2 + vs3) / self.rq3
        return vsd, vsq, vso

    def calculate_derivatives(self, vsd, vsq, vso):
        dervfsd = vsd - self.rs * self.isd
        dervfsq = vsq - self.rs * self.isq
        dervfrd = -self.rr * self.ird - self.frq * self.wm
        dervfrq = -self.rr * self.irq + self.frd * self.wm
        deriso = (vso - self.rs * self.iso) / self.lso
        return dervfsd, dervfsq, dervfrd, dervfrq, deriso
    
    def update_fluxes_and_currents(self, dervfsd, dervfsq, dervfrd, dervfrq, deriso):
        self.fsd += dervfsd * self.h
        self.fsq += dervfsq * self.h
        self.frd += dervfrd * self.h
        self.frq += dervfrq * self.h
        self.iso += deriso * self.h
        fso = self.lso * self.iso
        return fso
    
    def calculate_electromagnetic_torque(self,):
        self.ce = self.amsr * (self.fsq * self.frd - self.fsd * self.frq)
            
        self.isd = self.idt * (self.lr * self.fsd - self.msr * self.frd)
        self.isq = self.idt * (self.lr * self.fsq - self.msr * self.frq)
            
        self.ird = self.idt * (-self.msr * self.fsd + self.ls * self.frd)
        self.irq = self.idt * (-self.msr * self.fsq + self.ls * self.frq)
        return

    def currents_and_fluxes_phases(self, fso):
        is1 = self.rq23 * self.isd + self.iso / self.rq3
        is2 = self.rq23 * (-self.isd / 2 + self.rq3 * self.isq / 2) + self.iso / self.rq3
        is3 = self.rq23 * (-self.isd / 2 - self.rq3 * self.isq / 2) + self.iso / self.rq3
            
        fs1 = self.rq23 * self.fsd + fso / self.rq3
        fs2 = self.rq23 * (-self.fsd / 2 + self.rq3 * self.fsq / 2) + fso / self.rq3
        fs3 = self.rq23 * (-self.fsd / 2 - self.rq3 * self.fsq / 2) + fso / self.rq3
        return is1, is2, is3, fs1, fs2, fs3

    def mechanical_speed(self,):
        wm = self.wm + (self.ce - self.cl - self.wm * self.kf) * self.h / self.jm
        return wm
    
    def mechanical_torque(self,):
        cm = self.ce - self.cl
        return cm

    def outputs(self, is1, is2, is3, fs1, fs2, fs3, fso, cm, vso, vsd, vsq):
        self.tempo.append(self.t)
        self.corrented.append(self.isd)
        self.correnteq.append(self.isq)
        self.corrente1.append(is1)
        self.corrente2.append(is2)
        self.corrente3.append(is3)
        self.tensao1.append(vso)
        self.tensao2.append(vsd)
        self.tensao3.append(vsq)
        self.fluxord.append(self.frd)
        self.fluxorq.append(self.frq)
        self.fluxosd.append(self.fsd)
        self.fluxosq.append(self.fsq)
        self.fluxos1.append(fs1)
        self.fluxos2.append(fs2)
        self.fluxos3.append(fs3)
        self.fluxos.append(fso)
        self.conjugado.append(self.ce)
        self.velocidade.append(self.wm)
        self.correnteo.append(self.iso)
        self.frequencia.append(self.ws)
        self.torque_mecanico.append(cm)
        self.conjcarga.append(self.cl) 

    def simulate(self):
        while self.t < self.tmax:
            self.t += self.h
            vs1, vs2, vs3 = self.source_voltage()
            self.load_torque()
            vsd, vsq, vso = self.direct_voltage_and_quadrature(vs1, vs2, vs3)
            dervfsd, dervfsq, dervfrd, dervfrq, deriso = self.calculate_derivatives(vsd, vsq, vso)
            fso = self.update_fluxes_and_currents(dervfsd, dervfsq, dervfrd, dervfrq, deriso)
            self.calculate_electromagnetic_torque()
            is1, is2, is3, fs1, fs2, fs3 = self.currents_and_fluxes_phases(fso)
            self.wm = self.mechanical_speed()
            self.prev_wm = self.wm             # guarda para próximo passo
            cm = self.mechanical_torque()
            if self.t >= self.tp:
                self.tp += self.hp
                self.outputs(is1, is2, is3, fs1, fs2, fs3, fso, cm, vso, vsd, vsq)

    """def plot_motor(self):
        # Plotting currents
        plt.figure(1)
        plt.plot(self.tempo, self.corrente1, label='Current 1 (A)')
        plt.plot(self.tempo, self.corrente2, label='Current 2 (A)')
        plt.plot(self.tempo, self.corrente3, label='Current 3 (A)')
        plt.title('Currents (A)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Current (A)')
        plt.show()

        # Plotting voltages
        plt.figure(2)
        plt.plot(self.tempo, self.tensao1, label='Voltage 1 (V)')
        plt.plot(self.tempo, self.tensao2, label='Voltage 2 (V)')
        plt.plot(self.tempo, self.tensao3, label='Voltage 3 (V)')
        plt.title('Voltages (V)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Voltage (V)')
        plt.show()

        # Plotting fluxes
        plt.figure(3)
        plt.plot(self.tempo, self.fluxos1, label='Flux 1 (Wb)')
        plt.plot(self.tempo, self.fluxos2, label='Flux 2 (Wb)')
        plt.plot(self.tempo, self.fluxos3, label='Flux 3 (Wb)')
        plt.title('Fluxes (Wb)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Flux (Wb)')
        plt.show()

        # Plotting zero-sequence current
        plt.figure(5)
        plt.plot(self.tempo, self.correnteo, label='Current o (A)')
        plt.title('Current o (A)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Current (A)')
        plt.show()

        # Plotting multiple graphs
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 2, 1)
        plt.plot(self.tempo, self.conjcarga, label='Load Torque (N*m)')
        plt.title('Load Torque (N*m)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Torque (N*m)')

        plt.subplot(2, 2, 2)
        plt.plot(self.tempo, self.velocidade, label='Speed (rad/s)')
        plt.title('Speed (rad/s)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Speed (rad/s)')

        plt.subplot(2, 2, 3)
        plt.plot(self.tempo, self.conjugado, label='Electromagnetic Torque (N*m)')
        plt.title('Electromagnetic Torque (N*m)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Torque (N*m)')

        plt.subplot(2, 2, 4)
        plt.plot(self.tempo, self.torque_mecanico, label='Mechanical Torque (N*m)')
        plt.title('Mechanical Torque (N*m)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Torque (N*m)')

        plt.tight_layout()
        plt.show()

    def transfer_function(self):
        num = [self.msr* self.p]
        den = [self.jm, 2 * self.kf, self.rs + self.ls]
        
        return num, den

    def plot_bode(self):
        num, den = self.transfer_function()
        system = signal.TransferFunction(num, den)
        w, mag, phase = signal.bode(system)

        plt.figure(figsize=(10, 6))
        plt.subplot(2, 1, 1)
        plt.semilogx(w, mag)
        plt.title('Diagrama de Bode - Motor')
        plt.ylabel('Magnitude (dB)')
        plt.grid(which="both", axis="both")

        plt.subplot(2, 1, 2)
        plt.semilogx(w, phase)
        plt.xlabel('Frequência (rad/s)')
        plt.ylabel('Fase (graus)')
        plt.grid(which="both", axis="both")

        plt.tight_layout()
        plt.show()

    def plot_nyquist(self):
        num, den = self.transfer_function()
        motor_system = ctl.TransferFunction(num, den)

        # Frequências para o Diagrama de Nyquist
        w_start = 1e-2
        w_stop = 1e3
        num_points = 1000
        frequencies = np.logspace(np.log10(w_start), np.log10(w_stop), num_points)

        plt.figure()
        ctl.nyquist_plot(motor_system, omega=frequencies)
        plt.title("Diagrama de Nyquist - Motor Trifásico")
        plt.grid(True)
        plt.show()

    def state_space_representation(self):
        # Coeficientes da função de transferência
        num, den = self.transfer_function()

        # Sistema de segunda ordem: Numerador e denominador
        # Exemplo: num = [b0], den = [a2, a1, a0]
        a2 = den[0]
        a1 = den[1]
        a0 = den[2]
        b0 = num[0]

        # Matrizes A, B, C, D no espaço de estados
        A = np.array([[0, 1],
                      [-a0/a2, -a1/a2]])

        B = np.array([[0],
                      [b0/a2]])

        C = np.array([[1, 0]])

        D = np.array([[0]])

        return A, B, C, D

    def print_state_space(self):
        A, B, C, D = self.state_space_representation()
        print("Matriz A:")
        print(A)
        print("\nMatriz B:")
        print(B)
        print("\nMatriz C:")
        print(C)
        print("\nMatriz D:")
        print(D)

    def step_response(self):
        num, den = self.transfer_function()
        system = signal.TransferFunction(num, den)
        t, response = signal.step(system)

        plt.figure(figsize=(10, 6))
        plt.plot(t, response, label='Resposta ao Degrau Unitário')
        plt.title('Resposta ao Degrau Unitário - Sistema de Segunda Ordem')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Amplitude')
        plt.grid(True)
        plt.legend()
        plt.show()

    def example():
        
        motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        motor.simulate()
        motor.plot_motor()
        motor.plot_bode()
        motor.plot_nyquist() 
        motor.print_state_space()
        motor.step_response()"""

# exemplos
#Motor.example()

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

    def Reduções(self, eficiencia_corrente=0.96):

        motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)

        motor.simulate()
        
        # Dados do motor
        rpm_motor = np.array(motor.velocidade) * 60 / (2 * np.pi)
        torque_eletromag = np.array(motor.conjugado)
        potencia_motor = 11

        # Etapa 1: Redução primária (Redutor)
        rpm_pos_redutor = rpm_motor / self.reducao_primaria
        torque_pos_redutor = torque_eletromag * self.reducao_primaria
        potencia_pos_redutor = potencia_motor  # desprezando perdas aqui

        # Etapa 2: Perdas por corrente de rolos (efeito poligonal e atrito)
        rpm_pos_corrente = rpm_pos_redutor  # corrente não altera rotação, só transmite
        torque_pos_corrente = torque_pos_redutor * eficiencia_corrente
        potencia_pos_corrente = potencia_pos_redutor * eficiencia_corrente
        
        # Etapa 3: Redução final (Diferencial)
        rpm_dif = rpm_pos_corrente / self.reducao_final
        torque_dif = torque_pos_corrente * self.reducao_final
        potencia_dif = potencia_pos_corrente  # perdas desprezadas nesta etapa

        rendimento_transmissao = (potencia_dif / potencia_motor)

        torque_ef = torque_dif * rendimento_transmissao

        return  [torque_ef, potencia_dif, rpm_dif]
    
    def CarPerformance(self):

        # Parâmetros do veículo
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        r = self.raio_pneu * 0.001      # raio da roda (m)
        c_r = 0.012                     # coeficiente de rolamento
        b = 0.01                        # coeficiente de atrito (transmissão)
        torque = 0                      # torque mecânico (será calculado depois que o torque de carga for obtido)

        # inicializa as velocidades
        v_angular = 0                   # velocidade angular inicial (rad/s)
        v_linear_ms = 0                 # velocidade linear inicial (m/s)

        # Sinais do motor já reduzidos pelo sistema de transmissão
        torque_reduzido_raw, _, rpm_reduzido_raw = self.Reduções()

        # Verifica se tempo final foi definido
        if self.tempo_f == 0:
            raise ValueError("Defina o tempo final manualmente: dt_model.tempo_f = <valor_em_segundos>")

        # Gerar vetor de tempo com resolução de 0.01s (100 Hz)
        variacao_tempo = np.arange(self.tempo_i, self.tempo_f, 0.01)
        n_amostras = len(variacao_tempo)

        # Expandir sinais do motor para cobrir todo o tempo
        def expandir(vetor):
            if len(vetor) >= n_amostras:
                return vetor[:n_amostras]
            else:
                ultima = vetor[-1]
                repeticoes = n_amostras - len(vetor)
                return np.concatenate((vetor, np.full(repeticoes, ultima)))

        # vetores expandidos
        torque_reduzido = expandir(torque_reduzido_raw)
        rpm_reduzido = expandir(rpm_reduzido_raw)

        parametros = []
        
        for _, torque_e, tempo in zip(rpm_reduzido, torque_reduzido, variacao_tempo):

            # passo de tempo
            if not parametros:
                dt = 0
            else:
                dt = tempo - parametros[-1]["tempo"]

            # Inércias
            J_roda = 0.5 * self.massa_roda * r**2           # inércia de uma roda         
            J_carro_equivalente = self.massa * r**2         # Inércia equivalente da massa do veículo transferida ao eixo
            J_total = (4 * J_roda) + J_carro_equivalente    # Inércia total sentida pelo motor/trenó de força

            # Torque de carga (rolamento, arrasto e atrito)
            torque_rr = (c_r * peso) * r
            torque_fa = ((densidade_ar * v_linear_ms**2 * coeficiente_arrasto * area_frontal) / 2) * r
            torque_atr_mec = b * v_angular
            
            # Torque de carga e resultante
            torque_carga = torque_rr + torque_fa + torque_atr_mec
            torque = torque_e - torque_carga

            # Forças
            forca_trativa = torque / r
            forca_arrasto = torque_fa / r 
            forca_rr = torque_rr / r

            forca_resistiva = forca_arrasto + forca_rr
            forca_efetiva = forca_trativa - forca_resistiva

            # Acelerações
            a_angular = torque / J_total  # usando a inércia total (roda + carro)
            a_linear = forca_efetiva / self.massa

            # Velocidades
            v_angular += a_angular * dt
            v_linear_ms += a_linear * dt
            v_linear_kmh = v_linear_ms * 3.6

            parametros.append({
                "ft": forca_trativa,
                "va": v_angular,
                "vlm": v_linear_ms,
                "vlk": v_linear_kmh,
                "fa": forca_arrasto,
                "rr": forca_rr,
                "ff": forca_efetiva,
                "tempo": tempo
            })

        return parametros, rpm_reduzido, variacao_tempo

    def printCarPerformance(self, performance):
        
        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear do Carro [km/h]")
        for i in range(0, len(performance), 10):  
            param = performance[i]
            print(f"{param['ft']:.2f}\t\t\t{param['va']:.2f}\t\t\t\t{param['vlk']:.2f}")

        print("\nForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for i in range(0, len(performance), 10): 
            param = performance[i]
            print(f"{param['fa']:.2f}\t\t\t{param['rr']:.2f}\t\t\t{param['ff']:.2f}")

        # Corrigir tamanho do rpm_motor para bater com o tempo
        tempo = [param["tempo"] for param in performance]
        n_amostras = len(tempo)

        motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        motor.simulate()

        rpm_motor_raw = np.array(motor.velocidade) * 60 / (2 * np.pi)
        if len(rpm_motor_raw) >= n_amostras:
            rpm_motor = rpm_motor_raw[:n_amostras]
        else:
            rpm_motor = np.concatenate([
                rpm_motor_raw,
                np.full(n_amostras - len(rpm_motor_raw), rpm_motor_raw[-1])
            ])

        velocidade_angular = [param["va"] for param in performance]
        velocidade_linear = [param["vlk"] for param in performance]

        # --- Figura 1: RPM x Tempo ---
        plt.figure(figsize=(8, 4))
        plt.plot(tempo, rpm_motor, color='tab:pink')
        plt.xlabel("Tempo [s]")
        plt.ylabel("RPM do Motor")
        plt.title("RPM x Tempo")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

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

        # Plot dos gráficos
        plt.figure(figsize=(12, 6))

        # Plot Velocidade Angular (rad/s) x RPM
        plt.subplot(1, 2, 1)
        plt.plot(rpm_motor, velocidade_angular, color='tab:red')
        plt.xlabel("RPM Motor")
        plt.ylabel("Velocidade Angular [rad/s]")
        plt.title("Velocidade Angular x RPM")
        plt.grid(True)

        # Plot Velocidade Linear (km/h) x RPM
        plt.subplot(1, 2, 2)
        plt.plot(rpm_motor, velocidade_linear, color='tab:purple')
        plt.xlabel("RPM Motor")
        plt.ylabel("Velocidade Linear [km/h]")
        plt.title("Velocidade Linear x RPM")
        plt.grid(True)

        plt.tight_layout()
        plt.show()

class Tire:
    
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=0, B0=0, B1=0, B2=1, B3=1, omega=315, slip_angle_start=-9, slip_angle_end=9, angle_step=0.5, WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
        
        # Inicializando parâmetros para o modelo de pneu
        self.tire_Fz = tire_Fz  # Carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # Ângulo de deslizamento lateral do pneu [rad]
        self.tire_Ls = tire_Ls  # Escorregamento longitudinal do pneu [Adimensional]
        self.tire_type = 'Default'  # Tipo de pneu
        self.tire_friction_coef = tire_friction_coef  # Coeficiente de fricção entre pneu e pista
        self.tire_Ca = tire_Ca  # Ângulo de camber do pneu
            
        # Comprimentos das barras do mecanismo de quatro barras
        self.L0 = B0  # Comprimento da barra de direção
        self.L1 = B1  # Comprimento do braço de direção
        self.L2 = B2  # Comprimento da bitola
        self.L3 = B3  # Comprimento do braço de direção

        # Ângulo de orientação da barra longitudinal em relação ao eixo horizontal
        self.alpha = np.radians(omega)

        # Array de ângulos de deslizamento lateral do pneu em graus
        self.theta2 = np.arange(slip_angle_start + 90, slip_angle_end + 91, angle_step)

        # Convertendo os ângulos de deslizamento lateral do pneu para radianos
        self.theta2_rad = np.radians(self.theta2)
        self.angle = self.theta2_rad[0]

        # Inicialização das listas de resultados para armazenar dados ao longo do cálculo
        self.AC = []  # Lista para armazenar AC
        self.beta = []  # Lista para armazenar beta
        self.psi = []  # Lista para armazenar psi
        self.lamda = []  # Lista para armazenar lambda
        self.theta3 = []  # Lista para armazenar theta3
        self.theta4 = []  # Lista para armazenar theta4
        self.Ox, self.Oy = 0, 0  # Lista para armazenar coordenadas de O
        self.Ax, self.Ay = [], []  # Listas para armazenar coordenadas de A
        self.Bx, self.By = [], []  # Listas para armazenar coordenadas de B
        self.Cx, self.Cy = B0, 0 # Listas para armazenar coordenadas de C
        self.w = []  # Lista para armazenar w
        self.om2, self.om4 = [], []  # Listas para armazenar om2 e om4
        self.alpha_dot = []  # Lista para armazenar alpha_dot
        self.outer_slip = []  # Lista para armazenar ângulos de deslizamento lateral externo
        self.inner_slip = []  # Lista para armazenar ângulos de deslizamento lateral interno
        self.static_slip_angle = None  # Variável para armazenar o ângulo de deslizamento estático (inicialmente não definido)

        # Coordenadas da barra longitudinal (entre-eixos)
        self.WB = WB  # Entre-eixos (wheelbase) fixo
        self.rear_axle_length = rear_axle_length  # Comprimento do eixo traseiro fixo
        self.rear_axle_center = (B0 / 2, 0)

        # Entradas de rigidez do pneu
        self.track_y = track_y
        self.tire_k = tire_k   
        
    def Tire_forces(self, params):
        
            # Desembalando os parâmetros de Pacejka
            E, Cy, Cx, Cz, c1, c2 = params

            # Calculando parâmetros intermediários
            Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))
            D = self.tire_friction_coef * self.tire_Fz
            Bz = Cs / (Cz * D)
            Bx = Cs / (Cx * D)
            By = Cs / (Cy * D)

            # Calculando a força lateral do pneu
            tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))

            # Calculando a força longitudinal do pneu
            tire_longitudinal_force = D * np.sin(Cx * np.arctan(9 * Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))

            # Calculando o momento de auto-alinhamento do pneu
            tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))

            # Calculando a força de camber
            camber_thrust = D * np.sin(Cy * np.arctan(By * self.tire_Ca))

            # Retornando as forças calculadas e o momento de auto-alinhamento
            return tire_lateral_force + 0.5 * camber_thrust, (10 + (tire_auto_align_moment / 55)), tire_longitudinal_force
    
    @staticmethod
    def SlipRatio(velocidade_angular, raio_pneu, velocidade_linear):
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1
        return value

    def printSlipRatio(self, tempo, slip_ratio, tire_longitudinal_forces):
        passo = 10  # Intervalo para impressão
        slip_ratio_constante = max(slip_ratio)  # <<< Defina aqui o valor constante desejado
        i = 0

        # Loop while para corrigir os valores de slip_ratio nos primeiros instantes
        while i < len(tempo) and tempo[i] <= 0.05:
            slip_ratio[i] = slip_ratio_constante
            i += 1

        print("Valores do Slip Ratio:")
        for i in range(0, len(slip_ratio), passo):
            print(f"{slip_ratio[i]:.2f}")

        # Figura com dois subgráficos
        plt.figure(figsize=(14, 6))

        # Gráfico 1: Slip Ratio x Tempo
        plt.subplot(1, 2, 1)
        plt.plot(tempo, slip_ratio, label='Slip Ratio', color='blue')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Slip Ratio')
        plt.title('Slip Ratio x Tempo')
        plt.grid(True)
        plt.legend()

        # Gráfico 2: Força Longitudinal x Slip Ratio
        plt.subplot(1, 2, 2)
        plt.plot(slip_ratio, tire_longitudinal_forces, label='Força Longitudinal', color='green')
        plt.xlabel('Slip Ratio')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal x Slip Ratio')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()
        
def Instancias():
    # Instância do modelo Drivetrain com os parâmetros fornecidos
    dt_model = Drivetrain(
        cgx=853,                # Centro de gravidade no eixo X [mm]
        cgy=294,                # Centro de gravidade no eixo Y [mm]
        massa=347,              # Massa do veículo [kg]
        massa_roda= 6,          # Massa da Roda [kg]    
        entre_eixos=1567,       # Entre-eixos [mm]
        raio_pneu=220,          # Raio do pneu [mm]
        reducao_primaria=2.12,  # Redução primária
        reducao_final=2.76,     # Redução final (diferencial)
        cp=2.22,                # Constante de potência (pode ser usada futuramente)
        tempo_i=0               # Tempo inicial [s]
    )

    dt_model.tempo_f = 10  # segundos

    # Cálculo da performance do veículo com os valores reduzidos
    performance_veiculo, variacao_rpm, variacao_tempo = dt_model.CarPerformance()

    # Impressão dos resultados
    dt_model.printCarPerformance(performance_veiculo)

    # Extração das velocidadesdos dados de performance
    velocidade_angular = np.array([dado["va"] for dado in performance_veiculo])
    velocidade_linear = np.array([dado["vlm"] for dado in performance_veiculo])
  
    # Cálculo do slip ratio
    slip_ratio = Tire.SlipRatio(velocidade_angular, dt_model.raio_pneu * 0.001, velocidade_linear)
        
    # Instância da classe Tire com os slip ratios calculados
    slip_model = Tire(
        tire_Fz= 850,               # Carga vertical no pneu [N]
        tire_Sa=0,                 # Ângulo de escorregamento lateral [rad]
        tire_Ls=slip_ratio,       # Slip ratio já calculado
        tire_friction_coef=1.45,   # Coeficiente de fricção
        tire_Ca=0                  # Ângulo de camber
    )

    # Parâmetros experimentais fornecidos
    result = [
        0.3336564873588197,
        1.6271741344929977,
        1,
        4.3961693695846655,
        931.4055775279057,
        366.4936818126405,
    ]

    # Cálculo das forças e do momento de auto-alinhamento
    tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = slip_model.Tire_forces(result)

    # Plotagem dos gráficos de slip ratio
    slip_model.printSlipRatio(variacao_tempo, slip_ratio, tire_longitudinal_forces)
    
Instancias()
