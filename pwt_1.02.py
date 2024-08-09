import numpy as np
import matplotlib.pyplot as plt

__all__ = ["Motor", "ModulacaoEscalar", "Peso"]


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

def example_motor():
    motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
    motor.simulate()
    motor.plot_motor()


# NAO ESTA COMPLETO !!!!!!!!!!!!!!!!!
class ModulacaoEscalar:
    def __init__(self, e_star, mu_values, num_cycles=1):
        self.e_star = e_star
        self.mu_values = mu_values
        self.num_cycles = num_cycles
        self.t = None
        self.vs10_star = None
        self.vs20_star = None
        self.vs30_star = None
        self.i1 = None
        self.i2 = None
        self.i3 = None
        self.thd_values = []

    def t_entrada(self):
        self.vs1_star = np.sin(2 * np.pi * self.t)
        self.vs2_star = np.sin(2 * np.pi * self.t - 2 * np.pi / 3)
        self.vs3_star = np.sin(2 * np.pi * self.t + 2 * np.pi / 3)

    def calculate_vN0max_vN0mim(self):
        self.vN0max_star = self.e_star / 2 - np.maximum.reduce([self.vs1_star, self.vs2_star, self.vs3_star])
        self.vN0mim_star = -self.e_star / 2 - np.minimum.reduce([self.vs1_star, self.vs2_star, self.vs3_star])
        
    def calculates_vN0(self, mu):
        self.vN0_star = mu * self.vN0max_star + (1 - mu) * self.vN0mim_star

    def tm(self):
        self.vs10_star = self.vs1_star + self.vN0_star
        self.vs20_star = self.vs2_star + self.vN0_star
        self.vs30_star = self.vs3_star + self.vN0_star

    def c(self):
        self.i1 = np.sin(2 * np.pi * self.t)
        self.i2 = np.sin(2 * np.pi * self.t - 2 * np.pi / 3)
        self.i3 = np.sin(2 * np.pi * self.t + 2 * np.pi / 3)

    def calculate_thd(self):# Placeholder para cálculo real do THD. Aqui estou assumindo um valor fictício para exemplo.
        thd = np.random.random() * 10  # Exemplo: THD fictício entre 0 e 10%
        return thd

    def simulate_mod(self, mu):
        self.t = np.linspace(0, self.num_cycles, 1000)  # Vetor de tempo
        self.t_entrada()  # Calcular tensões de entrada
        self.calculate_vN0max_vN0mim()  # Calcular vN0max* e vN0mim*
        self.calculates_vN0(mu)  # Calcular vN0*
        self.tm()  # Calcular tensões moduladas
        self.c()  # Calcular correntes

    def plot_mod(self, mu):
        plt.figure(figsize=(10, 6))
        
        plt.subplot(2, 1, 1)
        plt.plot(self.t, self.vs10_star, label=r'$v_{s10}^*$', color='black')
        plt.title(f'Tensão Modulada $v_{{s10}}^*$ (µ = {mu})')
        plt.ylabel('Tensão')
        plt.legend()
        
        plt.subplot(2, 1, 2)
        plt.plot(self.t, self.i1, label=r'$i_1$', color='blue')
        plt.plot(self.t, self.i2, label=r'$i_2$', color='red')
        plt.plot(self.t, self.i3, label=r'$i_3$', color='green')
        plt.title('Correntes')
        plt.ylabel('Corrente')
        plt.legend()
        
        plt.tight_layout()
        plt.show()

    def plot_thd_vs_mu(self):
        plt.figure(figsize=(8, 6))
        plt.plot(self.mu_values, self.thd_values, label='THD', color='black', marker='o')
        plt.title('THD em função de µ')
        plt.xlabel('µ')
        plt.ylabel('THD (%)')
        plt.legend()
        plt.grid(True)
        plt.show()

    def run(self):
        for mu in self.mu_values:
            self.simulate_mod(mu)
            thd = self.calculate_thd()
            self.thd_values.append(thd)
            self.plot_mod(mu)
        
        self.plot_thd_vs_mu()

def exemplo_modulacao_escalar():
    e_star = 1.0  # Tensão DC do inversor
    mu_values = [0.5, 1.0]  # Valores de µ
    simulation = ModulacaoEscalar(e_star, mu_values)
    simulation.run()


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

def example_peso():
    peso_pwt = Peso(10, 10, 65, 5)
    total = peso_pwt.peso_total()
    print(f"O peso total é {total} Kg")



#exemplo
example_motor()
example_peso()
exemplo_modulacao_escalar()
