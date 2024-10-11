import math
import matplotlib.pyplot as plt
import numpy as np


class Torque:
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
            cm = self.mechanical_torque()
            if self.t >= self.tp:
                self.tp += self.hp
                self.outputs (is1, is2, is3, fs1, fs2, fs3, fso, cm, vso, vsd, vsq)




    def example():
        motor = Torque(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        motor.simulate()
        val_torque = np.array(motor.torque_mecanico)
        return (val_torque [:167] )
        
Torque.example()


class Drivetrain:

    def __init__(self, cgx, cgy, massa, entre_eixos, coeficiente_atrito, raio_pneu, aceleracao_ideal, reducao_primaria, reducao_unica, rpm, cp, tempo_i):
        self.cgx = cgx
        self.cgy = cgy
        self. massa = massa
        self.entre_eixos =entre_eixos
        self.coeficiente_atrito = coeficiente_atrito
        self.raio_pneu = raio_pneu
        self.aceleracao_ideal = aceleracao_ideal
        self.reducao_primaria = reducao_primaria
        self.reducao_unica = reducao_unica
        self.rpm = rpm
        self.cp = cp
        self.new_rpm = 0
        self.tempo_i = tempo_i
        self.tempo_f = 0


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

    def showResults(self):
        
        peso, reacao_traseira, reacao_dianteira, forca_trativa, transferencia_longitudinal, carga_traseira, pico_torque_traseiro, carga_pneu, torque_pneu, reducao_final, torque_necessario_motor, aceleracao_primaria_real, aceleracao_primaria_ideal, aceleraco_real_final, forca_trativa_ideal, torque_pneu_ideal, torque_motor_ideal, transferencia_carga_ideal, transferencia_carga_real, carga_traseira_ideal = Drivetrain.CalculateOutputs(self)
                
        # Print dos resultados obtidos
        print(f'''
            Resultados:
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
        
        
    def CarPerformance(self):
        # Definindo parâmetros do veículo
        peso = self.massa * 9.81  # Força peso
        rendimento_transmissao = 0.9
        transmissao_motor_roda = 0.9
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        basic_f = 0.015
        speed_f = 0.012
        c_r = 0.015  # Coeficiente de resistência ao rolamento 
        Iw = (self.massa * self.raio_pneu) / 2  # Momento de inércia da roda

        # Listas para armazenar os dados
        variacao_rpm = []
        variacao_tempo = []
        variacao_torque = []

        
        if self.new_rpm:
            variacao_rpm = range(self.rpm, self.new_rpm, 30)
            variacao_tempo = np.linspace(self.tempo_i, self.tempo_f, len(variacao_rpm))
            variacao_torque = Torque.example()
            self.tempo_i = self.tempo_f
            self.rpm = self.new_rpm
            self.new_rpm = 0
            self.tempo_f = 0
           
        else:
            variacao_rpm = [self.rpm]
            variacao_tempo = [self.tempo_i]
            

        parametros = []

        for rpm, tempo, torque in zip(variacao_rpm, variacao_tempo, variacao_torque):

            # Inicializa a velocidade angular (rad/s) a partir da RPM
            velocidade_angular = (rpm * 2 * math.pi) / (60 * self.reducao_primaria * self.reducao_unica * self.cp)

            # Cálculo da força trativa (N)
            forca_trativa = ((torque * self.reducao_primaria * self.reducao_unica * self.cp) / (self.raio_pneu * 0.001)) * rendimento_transmissao

            # Cálculo da resistência ao rolamento (N)
            F_r = (basic_f + (3.24 * speed_f * (((velocidade_angular) * (self.raio_pneu * 0.001)) ** 2.5))) * peso

            # Cálculo do momento My (Nm)
            My = c_r * peso * self.raio_pneu * 0.001  

            # Cálculo da aceleração angular (rad/s^2)
            aceleracao_angular = (torque - My) / Iw

            # Atualizando a velocidade angular (rad/s)
            velocidade_angular += aceleracao_angular * tempo  # Atualiza a velocidade angular com base no tempo

            # Cálculo da velocidade linear (km/h)
            velocidade_linear = ((velocidade_angular * (self.raio_pneu * 0.001)) * transmissao_motor_roda) * 3.6

            # Cálculo da força de arrasto (N)
            fa = (densidade_ar * velocidade_linear ** 2 * coeficiente_arrasto * area_frontal) / 2

            # Cálculo da força final (N)
            forca_final = forca_trativa - fa - F_r

            # Armazenar os parâmetros calculados em um dicionário
            parametro = {
                "forca_trativa": forca_trativa,
                "va": velocidade_angular,
                "velocidade_linear": velocidade_linear,
                "fa": fa,
                "rr": F_r,
                "My": My,
                "forca_final": forca_final
            }

            parametros.append(parametro)

        return parametros, variacao_rpm, variacao_tempo

    def printCarPerformance(self):
                  
        rpm = self.rpm
        new_rpm = self.new_rpm
        
        performance, variacao_rpm, variacao_tempo = self.CarPerformance()
        
        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]")
        
        for param in performance:
            print(f"{param['forca_trativa']}\t{param['va']}\t{param['velocidade_linear']}")
        
        print("\nForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for param in performance:
            print(f"{param['fa']}\t{param['rr']}\t{param['forca_final']}")

        self.rpm = rpm
        self.new_rpm = new_rpm


    def HalfShaftsSizing(self, fsi=1.25, tet=786, tec=471.6, dif=1):
        
        # Obtendo o maior torque do motor a partir dos dados experimentais 
        torque_max_motor = max(Torque.example())
        
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


class Tire:
    
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=0, B0=0, B1=0, 
                B2=1, B3=1, omega=315, slip_angle_start=-9, slip_angle_end=9, angle_step=0.5, WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
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
    def slip_ratio_1(velocidade_angular, raio_pneu):
       
        velocidade_linear = initial_speed
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1

        return value

    @staticmethod
    def show_slip_ratio(rpm_values, slip_ratio, velocidade_angular):
        
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
        
    def plot_graph(self, tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces, tire_lateral_experimental=None, tire_auto_align_experimental=None, angles=None, ratio=None):

        # Definindo um tamanho para a figura
        plt.figure(figsize=(20, 7))

        # Plotagem força lateral
        plt.subplot(1, 3, 1)
        plt.plot(self.tire_Sa, tire_lateral_forces, label='Curva Otimizada')
        plt.scatter(angles, tire_lateral_experimental, color='red', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Força Lateral do Pneu (N)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Plotagem torque auto-alinhante
        plt.subplot(1, 3, 2)
        plt.plot(self.tire_Sa, tire_auto_align_moment, label='Curva Otimizada')
        plt.scatter(angles, tire_auto_align_experimental, color='blue', label='Dados Experimentais')
        plt.xlabel('Ângulo de Deslizamento Lateral (graus)')
        plt.ylabel('Torque auto-alinhante (N.m)')
        plt.title('Curva Otimizada com os Dados Experimentais')
        plt.legend()
        plt.grid(True)

        # Plotagem força longitudinal
        plt.subplot(1, 3, 3)
        plt.plot(self.tire_Ls, tire_longitudinal_forces, label='Curva Sem Otimizar')
        plt.xlabel('Taxa de Escorregamento Longitudinal (Admensional)')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal - Sem Dados Experimentais')
        plt.legend()
        plt.grid(True)

        plt.tight_layout(pad=3.0)  # Aumentando a distância entre os subplots
        plt.show()


def uniteExample():
    dt_model = Drivetrain(
        cgx = 853,  # mm
        cgy = 294,  # mm
        massa = 347,  # kg
        entre_eixos = 1567,  # mm
        coeficiente_atrito = 0.9 , # coeficiente de atrito
        raio_pneu = 259,  # mm
        aceleracao_ideal = 1.2,  # g
        reducao_primaria = 2.12,  # redução primária
        reducao_unica = 2.76,
        rpm = 0,
        cp = 2.22,
        tempo_i=0)
    
    dt_model.showResults()
    
    dt_model.HalfShaftsSizing()
    
    dt_model.new_rpm = 5000
    dt_model.tempo_f = 5

    #Exemplo em lista
    dt_model.printCarPerformance() 

    #Recebendo os dados da performance do carro. Aqui que se encontra dados da velocidade angular
    performance_veiculo, variacao_rpm, variacao_tempo = dt_model.CarPerformance()

    #Filtrando a velocidade angular
    velocidade_angular = [dado['va'] for dado in performance_veiculo]
       
    #Transformando num array    
    velocidade_angular = np.array(velocidade_angular)

    #Transformando num array  
    variacao_rpm = np.array(variacao_rpm)

    #Calcular o slip ratio
    slip_ratio = Tire.slip_ratio_1(velocidade_angular, 0.259)

    #Plotagem de gráfico do slip ratio e saídas de seus dados no terminal
    Tire.show_slip_ratio(variacao_rpm, slip_ratio, velocidade_angular)

    #Salvando os dados como array para cálculo de força longitudinal
    slip_ratios = np.array(slip_ratio)
    
    #Criando instância da classe Tire
    Slip_model = Tire(tire_Fz=1500, tire_Sa=0, tire_Ls=slip_ratios, tire_friction_coef=1.45, tire_Ca=0)
    
    #Dados experimentais para instância em Tire
    result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]
    
    #Recebendo valores de força lateral, torque auto allinhante e força longitudinal
    tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = Slip_model.Tire_forces(result)
    
    #Plotagem de gráficos
    Slip_model.plot_graph( tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces)

initial_speed = 1

uniteExample()

Tire()
