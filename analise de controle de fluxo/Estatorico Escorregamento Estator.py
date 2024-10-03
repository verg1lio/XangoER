import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#import cmath  # Para operações com números complexos
#from matplotlib.widgets import Slider
# from scipy import signal

class ControleFluxoEstatoricoEscorregamentoEstator:
    def __init__(self, regu_l=1, time_end=0.3, num_points=1000000):
        # Initialize constants
        self.regu_l = regu_l
        self.war = 0.0
        self.wr = 0.0
        self.flr = 1.0
        self.fkp = 1.0
        self.fki = 0.5
        self.rs = 0.1
        self.te = 1.0
        self.eki_d = 0.0
        self.eki_q = 0.0
        self.anga = 0.0
        self.pid = 2 * np.pi
        self.flsr = 1.0
        self.esdl = 0.1
        self.fsa = 0.0
        self.fsb = 0.0
        self.dei = [0.0, 0.0]

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
            # Example current values (you should replace these with actual simulation data)
            isa = np.sin(t)
            isb = np.cos(t)

            if self.regu_l == 0:
                # P controller strategy
                #geração do angulo de referencia 
                self.anga += (self.wr + self.war) * self.te 
                if self.anga >= self.pid:
                    self.anga -= self.pid
                
                cga = np.cos(self.anga)
                sga = np.sin(self.anga)
                
                #geração referencia em referencial estatorico 
                fsa_i = self.flsr * cga
                fsb_i = self.flsr * sga

                #erros de fluxo
                self.dei[0] = fsa_i - self.fsa
                self.dei[1] = fsb_i - self.fsb

                #reguladores - saída de tensões de referência
                #eixo d - controle de fluxo
                vsd = self.fkp * self.dei[0] + self.rs * isa  
                #eixo q - controle do conjugado
                vsq = self.fkp * self.dei[1] + self.rs * isb  #duvida (ao inves de isb, deveria ser isa?)

                #* *#
                self.vsa_r.append(vsd) 
                self.vsb_r.append(vsq)
#-------------------------------------------------------------------------------------------#
            elif self.regu_l == 1:
                # Estratégia de controle usando reguladores PI (eixos dq)
                # Geração do angulo de referencia 
                self.anga += (self.wr + self.war) * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid

                # Definindo variáveis sen e cos
                cga = np.cos(self.anga)
                sga = np.sin(self.anga)

                # Geração referencia em referencial estatórico
                fsa_i = self.flsr * cga
                fsb_i = self.flsr * sga

                #Calculo dos erros do fluxo e do torque 
                self.dei[0] = fsa_i - self.fsa #(flxuo)
                self.dei[1] = fsb_i - self.fsb #(torque)

                # Reguladores - saídas de tensões de referencia 
                # Eixo d - controle de fluxo 
                vsd = self.eki_d + (self.fkp + self.fki) * self.dei[0]
                self.eki_d += self.fki * self.dei[0]
                vsd -= self.esdl * self.fsa #fra?

                # Eixo q - controle de torque 
                vsq = self.eki_q + (self.fkp + self.fki) * self.dei[1]
                self.eki_q += self.fki * self.dei[1]
                vsq -= self.esdl * self.fsb  #frb?

                # Limitação da saída regulador PI parte integral 
                if self.eki_d >= 180.:
                    self.eki_d = 180.
                if self.eki_d <= -180.:
                    self.eki_d = -180.
                if self.eki_q >= 180.:
                    self.eki_q = 180.
                if self.eki_q <= -180.:
                    self.eki_q = -180.

                #Transformação de referencia campo/estator
                self.vsa_r.append(vsd)
                self.vsb_r.append(vsq)

            # Saturation
            if len(self.vsa_r) > 0:
                self.vsa_r[-1] = np.clip(self.vsa_r[-1], -311., 311.)
                self.vsb_r[-1] = np.clip(self.vsb_r[-1], -311., 311.)

    def plot(self):
        plt.figure(figsize=(12, 6))

        plt.plot(self.time, self.vsa_r, label='vsa_r')
        plt.plot(self.time, self.vsb_r, label='vsb_r')
        plt.title('Tensões de Referência vsa_r e vsb_r ao Longo do Tempo')
        plt.xlabel('Tempo')
        plt.ylabel('Tensão de Referência')
        plt.grid(True)
        plt.legend()
        plt.show()

def example_ControleFluxoEstatoricoEscorregamentoEstator():
    control = ControleFluxoEstatoricoEscorregamentoEstator(regu_l=1)  # Set regu_l=0 for P control or 1 for PI control
    control.compute()
    control.plot()

if __name__ == "__main__":
    example_ControleFluxoEstatoricoEscorregamentoEstator()
