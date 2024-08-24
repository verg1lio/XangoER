import numpy as np
# from scipy.integrate import odeint
import matplotlib.pyplot as plt
#import cmath  # Para operações com números complexos
#from matplotlib.widgets import Slider
#from scipy import signal
class ControleFluxoEstatoricoQuadratura:
 
    def __init__(self, regu_l=0, time_end=10, num_points=100):
        # Initialize constants with more realistic values
        self.regu_l = regu_l
        self.war = 50.0  # Velocidade angular relativa (rad/s)
        self.wr = 314.0  # Velocidade angular do rotor (rad/s)
        self.lm = 0.15  # Indutância magnetizante (H)
        self.rr = 0.5   # Resistência do rotor (Ohm)
        self.flr = 1.0  # Fluxo de referência do rotor (Wb)
        self.fkp = 10.0  # Fator de ganho proporcional
        self.fki = 5.0  # Fator de ganho integral
        self.rs = 0.1   # Resistência do estator (Ohm)
        self.te = 0.001  # Período de amostragem (s)
        self.eki_d = 0.0
        self.eki_q = 0.0
        self.anga = 0.0
        self.pid = 2 * np.pi
        self.fsa = 0.0
        self.fsb = 0.0
        self.dei = [0.0, 0.0]
        self.flsr = 1.0  # Fluxo do estator de referência (Wb)
        self.esdl = 0.1
        self.fr_d = 1.0  # Componente de Fluxo no eixo d (Wb)
        self.fr_q = 1.0  # Componente de Fluxo no eixo q (Wb)

        # Simulation parameters
        self.time = np.linspace(0, time_end, num_points)
        self.h = self.time[1] - self.time[0]  # Time step

        # Storage for output variables
        self.vsd = []
        self.vsq = []
        self.vsa_r = []
        self.vsb_r = []

    def compute(self):
        for t in self.time:
            # Example current values (replace these with actual simulation data)
            isa = np.sin(t)
            #isb = np.cos(t)

            if self.regu_l == 0:
                # P controller strategy with quadrature field
            
                
                self.war = self.lm * self.rr * isa / (self.flr * self.rr)
                self.anga += (self.wr + self.war) * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid
                
                cga = np.cos(self.anga)
                sga = np.sin(self.anga)

                fsd = self.fsa * cga + self.fsb * sga
                fsq = self.fsb * cga - self.fsa * sga

                self.dei[0] = self.flsr - fsd
                self.dei[1] = -fsq

                vsd = self.fkp * self.dei[0]
                vsq = self.fkp * self.dei[1]

                vsa_r = vsd * cga - vsq * sga
                vsb_r = vsq * cga + vsd * sga

                self.vsa_r.append(vsa_r)
                self.vsb_r.append(vsb_r)

            elif self.regu_l == 1:
                # PI controller strategy with quadrature field

                
                self.anga += (self.wr + self.war) * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid

                cga = np.cos(self.anga)
                sga = np.sin(self.anga)

                fsd = self.fsa * cga + self.fsb * sga
                fsq = self.fsb * cga - self.fsa * sga
                frd = self.fr_d * cga + self.fr_q * sga 
                frq = self.fr_q * cga - self.fr_q * sga  

                self.dei[0] = self.flsr - fsd
                self.dei[1] = -fsq

                vsd = self.eki_d + (self.fkp + self.fki) * self.dei[0] # SAIDA PROPORCIONAL 
                self.eki_d += self.fki * self.dei[0]# ACUMULO DO ERRO
                vsd -= self.esdl * frd - (self.war + self.wr) * fsd 

                vsq = self.eki_q + (self.fkp + self.fki) * self.dei[1] # SAIDA INTEGRAL
                self.eki_q += self.fki * self.dei[1] # ACUMULO DO ERRRO
                vsq -= self.esdl * frq + (self.war + self.wr) * fsd

                # Saturation
                if self.eki_d >= 180.:
                    self.eki_d = 180.
                if self.eki_d <= -180.:
                    self.eki_d = -180.
                if self.eki_q >= 180.:
                    self.eki_q = 180.
                if self.eki_q <= -180.:
                    self.eki_q = -180.

                vsa_r = vsd * cga - vsq * sga
                vsb_r = vsq * cga + vsd * sga

                self.vsa_r.append(vsa_r)
                self.vsb_r.append(vsb_r)

            # Saturation
            if len(self.vsa_r) > 0:
                 self.vsa_r[-1] = np.clip(self.vsa_r[-1], -311., 311.)
                 self.vsb_r[-1] = np.clip(self.vsb_r[-1], -311., 311.)
             
                
    def plot(self):
        """
        Plots the reference voltage vsa_r and vsb_r over time.
        """
        plt.figure(figsize=(12, 6))

        plt.plot(self.time, self.vsa_r, label='vsa_r')
        plt.plot(self.time, self.vsb_r, label='vsb_r')

        # Set title, labels, and legend
        if self.regu_l == 1: 
            plt.axhline(311.0, color='red', linestyle='--', label='Limite superior (+311 V)')
            plt.axhline(-311.0, color='red', linestyle='--', label='Limite inferior (-311 V)')
        
        plt.title('Tensões de Referência vsa_r e vsb_r ao Longo do Tempo')
        plt.xlabel('Tempo')
        plt.ylabel('Tensão de Referência')
        plt.grid(True)
        plt.legend()

        # Show the plot
        plt.show()  # plt.savefig('ControleQuadratura.png')


    def example_ControleFluxoEstatoricoQuadratura():
        control = ControleFluxoEstatoricoQuadratura(regu_l=1)  # Set regu_l=0 for P control or 1 for PI control
        control.compute()
        control.plot()




ControleFluxoEstatoricoQuadratura.example_ControleFluxoEstatoricoQuadratura()

