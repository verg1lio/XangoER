import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control as clt

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

    def source_voltage(self,):
        """Tensão da fonte

        Calcula a tensão das três fases em função do ângulo elétrico do estator.

        Returns
        -------
        vs1 : float
            Tensão da fase 1
        vs2 : float
            Tensão da fase 2
        vs3 : float
            Tensão da fase 3

        Examples
        --------
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> v1, v2, v3 = motor.source_voltage()
        >>> v1, v2, v3
        (311.1247727163462, -154.5465853680912, -156.57818734825486)
        """
        
        # Atualiza o ângulo elétrico do estator
        self.tete += self.h * self.ws
        if self.tete >= 2 * np.pi:
            self.tete -= 2 * np.pi
            
        # Calcula as tensões de cada fase 
        vs1 = self.Vs * np.cos(self.tete)
        vs2 = self.Vs * np.cos(self.tete - self.pi23)
        vs3 = self.Vs * np.cos(self.tete + self.pi23)

        return vs1, vs2, vs3

    def load_torque(self,):
        """Torque de carga

        Ajusta o torque de carga (cl) com base no tempo atual (t).
        Se o tempo for maior ou igual à metade do tempo máximo (tmax), o torque de carga é definido como 10.

        Examples
        --------
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> motor.t = 5
        >>> motor.tmax = 10
        >>> motor.load_torque()
        >>> motor.cl
        10
        """
        if self.t >= self.tmax / 2:
            self.cl = 10        
    
    def direct_voltage_and_quadrature(self, vs1, vs2, vs3):
        """Tensão direta e em quadratura

        Converte as tensões das três fases em coordenadas de eixo direto (d) e de quadratura (q), além de calcular a tensão de sequência zero (vso).

        Parameters
        ----------
        vs1 : float
            Tensão da fase 1.
        vs2 : float
            Tensão da fase 2.
        vs3 : float
            Tensão da fase 3.

        Returns
        -------
        vsd : float
            Tensão na coordenada d (direta).
        vsq : float
            Tensão na coordenada q (quadratura).
        vso : float
            Tensão de sequência zero.

        Examples
        --------
        >>> v1, v2, v3 = 311.1247727163462, -154.5465853680912, -156.57818734825486
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> vsd1, vsq1, vs01 = motor.direct_voltage_and_quadrature(v1, v2, v3)
        >>> vsd1, vsq1, vs01
        (381.0484697472187, 1.4365595368457367, 6.563712636189232e-14)
        """

        vsd = self.rq23 * (vs1 - vs2 / 2 - vs3 / 2)
        vsq = self.rq23 * (vs2 * self.rq3 / 2 - vs3 * self.rq3 / 2)
        vso = (vs1 + vs2 + vs3) / self.rq3

        return vsd, vsq, vso

    def calculate_derivatives(self, vsd, vsq, vso):
        """Cálculo de derivadas

        Calcula as derivadas dos fluxos e correntes com base nas tensões das fases e parâmetros do sistema.

        Parameters
        ----------
        vsd : float
            Tensão na coordenada d (direta).
        vsq : float
            Tensão na coordenada q (quadratura).
        vso : float
            Tensão de sequência zero.

        Returns
        -------
        dervfsd : float
            Derivada do fluxo do estator na coordenada d.
        dervfsq : float
            Derivada do fluxo do estator na coordenada q.
        dervfrd : float
            Derivada do fluxo do rotor na coordenada d.
        dervfrq : float
            Derivada do fluxo do rotor na coordenada q.
        deriso : float
            Derivada da corrente iso.

        Examples
        --------
        >>> vsd, vsq, vso = 381.04, 1.43, 6.56e-14
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> defsd, defsq, defrd, defrq, deiso = motor.calculate_derivatives(vsd, vsq, vso)
        >>> defsd, defsq, defrd, defrq, deiso
        (381.04, 1.43, -0.0, 0.0, 6.98e-12)
        """

        dervfsd = vsd - self.rs * self.isd
        dervfsq = vsq - self.rs * self.isq
        dervfrd = -self.rr * self.ird - self.frq * self.wm
        dervfrq = -self.rr * self.irq + self.frd * self.wm
        deriso = (vso - self.rs * self.iso) / self.lso

        return dervfsd, dervfsq, dervfrd, dervfrq, deriso
    
    def update_fluxes_and_currents(self, dervfsd, dervfsq, dervfrd, dervfrq, deriso):
        """Atualiza fluxos e correntes

        Atualiza os valores de fluxo e corrente no sistema, com base nas suas respectivas derivadas.

        Parameters
        ----------
        dervfsd : float
            Derivada do fluxo do estator na coordenada d.
        dervfsq : float
            Derivada do fluxo do estator na coordenada q.
        dervfrd : float
            Derivada do fluxo do rotor na coordenada d.
        dervfrq : float
            Derivada do fluxo do rotor na coordenada q.
        deriso : float
            Derivada da corrente iso.

        Returns
        -------
        fso : float
            Fluxo associado à corrente 'iso'.

        Examples
        --------
        >>> dervfsd, dervfsq, dervfrd, dervfrq, deriso = 381.048, 1.437, 0.0, 0.0, 6.983e-12
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> fso1 = motor.update_fluxes_and_currents(dervfsd, dervfsq, dervfrd, dervfrq, deriso)
        >>> fso1
        (6.563712636189233e-19)
        """

        self.fsd += dervfsd * self.h
        self.fsq += dervfsq * self.h
        self.frd += dervfrd * self.h
        self.frq += dervfrq * self.h
        self.iso += deriso * self.h
        fso = self.lso * self.iso

        return fso
    
    def calculate_electromagnetic_torque(self,):
        """Cálculo do torque eletromagnético

        Calcula o torque eletromagnético com base nos fluxos do estator e rotor.

        Returns
        -------
        ce : float
            Torque eletromagnético calculado.

        Examples
        --------
        >>> torque = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> ce = torque.calculate_electromagnetic_torque()
        >>> ce
        (0.0)
        """

        # Calcula o torque eletromagnético utilizando a fórmula baseada nos fluxos
        self.ce = self.amsr * (self.fsq * self.frd - self.fsd * self.frq)
            
        # Atualiza as correntes do estator nas coordenadas d e q    
        self.isd = self.idt * (self.lr * self.fsd - self.msr * self.frd)
        self.isq = self.idt * (self.lr * self.fsq - self.msr * self.frq)

        # Atualiza as correntes do rotor nas coordenadas d e q    
        self.ird = self.idt * (-self.msr * self.fsd + self.ls * self.frd)
        self.irq = self.idt * (-self.msr * self.fsq + self.ls * self.frq)
        
        return self.ce

    def currents_and_fluxes_phases(self, fso):
        """Fases de correntes e fluxos

        Calcula as fases das correntes e dos fluxos com base nas correntes e fluxos do sistema.

        Parameters
        ----------
        fso : float
            Fluxo associado à corrente 'iso'.

        Returns
        -------
        is1 : float
            Corrente fase 1.
        is2 : float
            Corrente fase 2.
        is3 : float
            Corrente fase 3.
        fs1 : float
            Fluxo fase 1.
        fs2 : float
            Fluxo fase 2.
        fs3 : float
            Fluxo fase 3.

        Examples
        --------
        >>> calculo_fluxo_corrente = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> is1, is2, is3, fs1, fs2, fs3 = calculo_fluxo_corrente.currents_and_fluxes_phases(fso=1)
        >>> is1, is2, is3, fs1, fs2, fs3
        (0.0, 0.0, 0.0, 0.5773502691896258, 0.5773502691896258, 0.5773502691896258)
        """

        # Calcula as correntes para as fases 1, 2 e 3
        is1 = self.rq23 * self.isd + self.iso / self.rq3
        is2 = self.rq23 * (-self.isd / 2 + self.rq3 * self.isq / 2) + self.iso / self.rq3
        is3 = self.rq23 * (-self.isd / 2 - self.rq3 * self.isq / 2) + self.iso / self.rq3

        # Calcula os fluxos para as fases 1, 2 e 3    
        fs1 = self.rq23 * self.fsd + fso / self.rq3
        fs2 = self.rq23 * (-self.fsd / 2 + self.rq3 * self.fsq / 2) + fso / self.rq3
        fs3 = self.rq23 * (-self.fsd / 2 - self.rq3 * self.fsq / 2) + fso / self.rq3
        
        return is1, is2, is3, fs1, fs2, fs3

    def mechanical_speed(self,):
        """Velocidade mecânica

        Calcula a velocidade mecânica do motor com base no torque eletromagnético, no torque de carga e nas características do motor.

        Returns
        -------
        wm : float
            Velocidade mecânica do motor.

        Examples
        --------
        >>> velocidade_mecanica = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> wm = velocidade_mecanica.mechanical_speed()
        >>> wm
        (0.0)
        """
        
        # Calcula a nova velocidade mecânica usando a equação do movimento
        wm = self.wm + (self.ce - self.cl - self.wm * self.kf) * self.h / self.jm
        return wm
    
    def mechanical_torque(self,):
        """Torque mecânico

        Calcula o torque mecânico do motor, que é a diferença entre o torque eletromagnético e o torque de carga.

        Returns
        -------
        cm : float
            Torque mecânico do motor.

        Examples
        --------
        >>> torque_mecanica = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> cm = torque_mecanica.mechanical_torque()
        >>> cm
        (0)
        """
        
        # Calcula o torque mecânico como a diferença entre o torque eletromagnético e o torque de carga
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

        # Plota as correntes das fases
        plt.figure(1)
        plt.plot(self.tempo, self.corrente1, label='Current 1 (A)')
        plt.plot(self.tempo, self.corrente2, label='Current 2 (A)')
        plt.plot(self.tempo, self.corrente3, label='Current 3 (A)')
        plt.title('Currents (A)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Current (A)')
        plt.show()

        # Plota as tensões das fases
        plt.figure(2)
        plt.plot(self.tempo, self.tensao1, label='Voltage 1 (V)')
        plt.plot(self.tempo, self.tensao2, label='Voltage 2 (V)')
        plt.plot(self.tempo, self.tensao3, label='Voltage 3 (V)')
        plt.title('Voltages (V)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Voltage (V)')
        plt.show()

        # Plota os fluxos das fases
        plt.figure(3)
        plt.plot(self.tempo, self.fluxos1, label='Flux 1 (Wb)')
        plt.plot(self.tempo, self.fluxos2, label='Flux 2 (Wb)')
        plt.plot(self.tempo, self.fluxos3, label='Flux 3 (Wb)')
        plt.title('Fluxes (Wb)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Flux (Wb)')
        plt.show()

        # Plotando a corrente homopolar
        plt.figure(5)
        plt.plot(self.tempo, self.correnteo, label='Current o (A)')
        plt.title('Current o (A)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Current (A)')
        plt.show()

        # Plota múltiplos gráficos em uma única figura
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


    def example():
        motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01, q1=1, q2=1, q3=0, valor_mu=1) # Varia o valor de mu entre 0 e 1
        motor.simulate()
        motor.plot_motor()

Motor.example()
