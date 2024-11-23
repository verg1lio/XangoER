import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control as ctl



__all__ = ["Motor", "ModulacaoEscalar", "Peso"]


class Motor:
    def __init__(self, rs, rr, ls, lr, mrs, jm, kf, q1, q2, q3, valor_mu):
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
        self.f_onda_p = 50 # Frequência da onda portadora do sinal PWM
        self.f_onda_m = 5 # Frequência da onda modulante do sinal PWM
        self.q1 = q1 # Chave de comutação do inversor
        self.q2 = q2 # Chave de comutação do inversor
        self.q3 = q3 # Chave de comutação do inversor
        self.jm = jm  # Moment of inertia (kg*m^2)
        self.kf = kf  # Friction coefficient (N*m*s)
        self.cte_tempo_mec = self.jm / self.kf  # Mechanical time constant (s)
        self.idt = 1 / (self.ls * self.lr - self.msr * self.msr)  # Inverse of the determinant
        self.p = 2  # Number of pole pairs
        self.amsr = self.p * self.idt * self.msr  # Constant for torque calculation
        self.valor_mu = valor_mu  # Escalar


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
        self.e_star = self.Vs

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
                self.outputs(is1, is2, is3, fs1, fs2, fs3, fso, cm, vso, vsd, vsq)

    def plot_motor(self):
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
        
     
    def chaves(self):
            
        """
        Determina a configuração das seis chaves que 
        compõem o do inversor de frequência (q1, q2, q3, q4, q5, q6). 
        Sendo q4, q5 e q6 os complementares de q1, q2 e q3, respectivamente.
        
        Returns
        -------
        str
            A configuração de chaves do inversor.
        """
        
        if self.q1 == 1: chave_1 = True # Retorno chave 1 fechada
        else: chave_1 = False # Retorno chave 1 aberta

        if self.q2 == 1: chave_2 = True # Retorno chave 2 fechada
        else: chave_2 = False # Retorno chave 2 aberta

        if self.q3 == 1: chave_3 = True # Retorno chave 3 fechada
        else: chave_3 = False # Retorno chave 3 aberta

        q1_bar = 1 - self.q1
        q2_bar = 1 - self.q2
        q3_bar = 1 - self.q3
            
        if q1_bar == 1: chave_4 = True # Retorno chave 1_bar fechada
        else: chave_4 = False # Retorno chave 1_bar aberta

        if q2_bar == 1: chave_5 = True # Retorno chave 2_bar fechada
        else: chave_5 = False # Retorno chave 2_bar aberta

        if q3_bar == 1: chave_6 = True # Retorno chave 3_bar fechada
        else: chave_6 = False # Retorno chave 3_bar aberta    

        return print(f'A configuração de chaves do inversor é: C1={chave_1}, C2={chave_2}, C3={chave_3}, C4={chave_4}, C5={chave_5}, C6={chave_6}')
    
    
    def controle_pwm(self):
        
        """
        Essa função gera o sinal PWM para controle do inversor com modulação escalar conforme o capítulo 7.4 do livro 
        "Sistemas de Acionamento Estático de Máquina Elétrica" de Cursino Brandão Jacobina.
        
        Parameters
        Self.Vs e Self.t_pwm (Passo de tempo)
        None
        
        Returns
        Plots do gráfico do sinal PWM, onda triangular portadora, correntes, 
        tensão de entrada e tensão modulada
    
        
        Notes
        -----
        A tensão modulada é calculada com base na equação (7.4) do livro "Sistemas de Acionamento Estático de Máquina Elétrica" 
        de Cursino Brandão Jacobina.
        A onda triangular portadora é gerada com base na equação (7.5) do mesmo livro.
        O sinal PWM puro é gerado com base na equação (7.6) do mesmo livro.
        
        """
        # Define o passo de tempo do controlador
        self.t_pwm = np.linspace(0, 2, 1000)  

        # Tensões de entrada
        self.v1 = self.Vs * np.sin(2 * np.pi * self.t_pwm)
        self.v2 = self.Vs * np.sin(2 * np.pi * self.t_pwm + 2 * np.pi / 3)
        self.v3 = self.Vs * np.sin(2 * np.pi * self.t_pwm + 4 * np.pi / 3)

        # Correntes
        self.i1 = self.Vs * np.sin(2 * np.pi * self.t_pwm - np.pi / 2)
        self.i2 = self.Vs * np.sin(2 * np.pi * self.t_pwm + 2 * np.pi / 3 - np.pi / 2)
        self.i3 = self.Vs * np.sin(2 * np.pi * self.t_pwm + 4 * np.pi / 3 - np.pi / 2)
       
        # Calculo max e min das tensoes        
        self.vN0max_star = (self.Vs / 2) - np.maximum.reduce([self.v1, self.v2, self.v3])
        self.vN0mim_star = -self.Vs / 2 - np.minimum.reduce([self.v1, self.v2, self.v3])
        
        # Tensão homopolar
        self.vN0_star =  self.valor_mu * self.vN0max_star + (1 - self.valor_mu) * self.vN0mim_star  
        
        # Tensões moduladas
        self.v10 = self.v1 + self.vN0_star
        self.v20 = self.v2 + self.vN0_star
        self.v30 = self.v3 + self.vN0_star
        
        # Geração da onda portadora triangular:
        periodo_port = 1 / self.f_onda_p  
        onda_port = 1 - 2 * np.abs((self.t_pwm % periodo_port) * self.f_onda_p - 0.5)  

        #Sinal pwm
        PWM_signal = np.where(self.v10 >= onda_port, self.Vs, 0)

         #Correntes
        plt.figure(figsize=(8, 4))
        plt.plot(self.t_pwm, self.i1, label='Current PWM 1 (A)')
        plt.plot(self.t_pwm, self.i2, label='Current PWM 2 (A)')
        plt.plot(self.t_pwm, self.i3, label='Current PWM 3 (A)')
        plt.title('Currents PWM (A)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Current PWM (A)')
        
        
        # Gráficos dos resultados
        plt.figure(figsize=(10, 8))

        # Tensão modulada
        plt.subplot(4, 1, 1)
        plt.plot(self.t_pwm, self.v10, label='V10', color='black')
        plt.title('Tensão modulada')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')
        plt.legend()


        # Onda triangular portadora
        plt.subplot(4, 1, 2)
        plt.plot(self.t_pwm, onda_port, label='Onda portadora')
        plt.title('Sinal da onda portadora')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')

        
        # Gráfico do sinal PWM puro
        plt.subplot(4, 1, 3)
        plt.step(self.t_pwm, PWM_signal, label='PWM', color='red')
        plt.title('Sinal PWM para controle do torque')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')


        # Onda de referência
        plt.subplot(4, 1, 4)
        plt.step(self.t_pwm, self.v1, label='Tensão 1', color='green')
        plt.title('Tensão de referência')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')
        
        plt.tight_layout()
        plt.show()   

    
    def exec_pwm(self):
        self.mu_values = np.linspace(0, 10e-1, 1)  # Variação de mu
        
        if self.valor_mu != 1:
            self.controle_pwm()  
        else:
            for mu in self.mu_values:
                
                self.valor_mu = mu  
                self.controle_pwm()  

    def example():
        motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01, q1=1, q2=1, q3=0, valor_mu=1) # Varia o valor de mu entre 0 e 1
        motor.simulate()
        # motor.plot_motor()
        # motor.plot_bode()
        # motor.plot_nyquist() 
        # motor.print_state_space()
        # motor.step_response()
        motor.chaves()
        motor.exec_pwm()
       

        

        

class Peso:
    def __init__(self, peso_bateria, peso_inversor, peso_motor, peso_chicote):
        self.peso_bateria = peso_bateria
        self.peso_inversor = peso_inversor
        self.peso_motor = peso_motor
        self.peso_chicote = peso_chicote

    def peso_total(self):
        """Peso total dos componentes do sistema de powertrain

        Calcula o peso total dos componentes do sistema de powertrain, somando o peso da bateria, do inversor e do motor.

        Returns
        -------
        peso_total : float
            Somatorio dos pesos dos componentes do sistema

        Examples
        --------
        >>> peso_pwt = peso(10, 10, 65, 5)
            90
        """
        peso_total = self.peso_bateria + self.peso_inversor + self.peso_motor + self.peso_chicote
        return peso_total

    def example():
        peso_pwt = Peso(10, 10, 65, 5)
        total = peso_pwt.peso_total()
        print(f"O peso total é {total} Kg")



# exemplos
Motor.example()
#Peso.example()