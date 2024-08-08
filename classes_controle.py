import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import cmath  # Para operações com números complexos
from matplotlib.widgets import Slider
from scipy import signal

__all__ = ["bateria", "motor_gaiola", "inversor", "controle_fluxo", "controle_inversor", "peso"]




# patricio caio
class ControleFluxoRotoricoEscorregamentoEstator:
    def __init__(self, regu_l=0, time_end=10, num_points=100):
        # Initialize constants
        self.regu_l = regu_l
        self.wbr = 0.0
        self.flrr = 1.0
        self.fra = 1.0
        self.frb = 0.0
        self.tr = 1.0
        self.lm = 1.0
        self.lr = 1.0
        self.rs = 0.1
        self.dkr = 1.0
        self.dka = 0.5
        self.hinv = 1.0
        self.te = 1.0
        self.eki_d = 0.0
        self.eki_q = 0.0
        self.fkp = 1.0
        self.fki = 0.5
        self.erf = 0.0
        self.erw = 0.0
        self.wr = 0.0
        self.angb = 0.0
        self.cg_b = 1.0
        self.sg_b = 0.0
        self.cgb = 1.0
        self.sgb = 0.0
        self.deico = [0.0, 0.0]
        
        # Simulation parameters
        self.time = np.linspace(0, time_end, num_points)
        self.h = self.time[1] - self.time[0]  # Time step

        # Storage for output variables
        self.vsd = []
        self.vsq = []
        self.vsa_r = []
        self.vsb_r = []
        self.isd_i = []
        self.isq_i = []

    def compute(self):
        for t in self.time:
            # Example current values (you should replace these with actual simulation data)
            isa = np.sin(t)
            isb = np.cos(t)
            flrr = self.flrr
            angb = self.angb

            if self.regu_l == 0:
                # Case 0
                self.wbi = self.wbr + self.wr
                self.angb = self.angb + self.wbi * self.h
                if self.angb >= 2 * np.pi:
                    self.angb -= 2 * np.pi

                fsd_i = flrr * np.cos(self.angb)
                fsq_i = flrr * np.sin(self.angb)

                vsd = -self.wbi * fsq_i
                vsq = self.wbi * fsd_i
                vsa_r = vsd
                vsb_r = vsq

            elif self.regu_l == 1:
                # Case 1
                frd_i = flrr * np.cos(self.angb)
                frq_i = flrr * np.sin(self.angb)

                self.deico[2] = frd_i - self.fra
                self.deico[3] = frq_i - self.frb

                isd_i = self.erf + (self.fkp + self.fki) * self.deico[2]
                self.erf += self.fki * self.deico[2]
                if self.erf >= 15.:
                    self.erf = 15.
                if self.erf <= -15.:
                    self.erf = -15.
                isd_i += (self.wr * frd_i * self.tr / self.lm)

                isq_i = self.erw + (self.fkp + self.fki) * self.deico[3]
                self.erw += self.fki * self.deico[3]
                if self.erw >= 15.:
                    self.erw = 15.
                if self.erw <= -15.:
                    self.erw = -15.
                isq_i -= (self.wr * frq_i * self.tr / self.lm)

                vsa_r = vsd
                vsb_r = vsq

            elif self.regu_l == 2:
                # Case 2 (Predictive regulation)
                isd_i = self.hinv * frd_i - self.hinv * (1. - self.te / self.tr) * self.fra + self.hinv * self.wr * self.te * frq_i
                isq_i = self.hinv * frq_i - self.hinv * self.wr * self.te * self.fra - self.hinv * (1. - self.te / self.tr) * frq_i

                if isd_i >= 8.:
                    isd_i = 8.
                if isd_i <= -8.:
                    isd_i = -8.
                if isq_i >= 8.:
                    isq_i = 8.
                if isq_i <= -8.:
                    isq_i = -8.

                # Store results
                self.isd_i.append(isd_i)
                self.isq_i.append(isq_i)

            # Storing voltage results
            self.vsd.append(vsd)
            self.vsq.append(vsq)
            self.vsa_r.append(vsa_r)
            self.vsb_r.append(vsb_r)

            # Update angle for next iteration
            self.angb += (self.wr + self.wbr) * self.h
            if self.angb >= 2 * np.pi:
                self.angb -= 2 * np.pi

    def plot(self):
        plt.figure(figsize=(12, 8))
        
        plt.subplot(3, 1, 1)
        plt.plot(self.time, self.vsd, label='vsd')
        plt.title('vsd over Time')
        plt.grid(True)
        plt.legend()
        
        plt.subplot(3, 1, 2)
        plt.plot(self.time, self.vsq, label='vsq')
        plt.title('vsq over Time')
        plt.grid(True)
        plt.legend()
        
        plt.subplot(3, 1, 3)
        plt.plot(self.time, self.vsa_r, label='vsa_r')
        plt.plot(self.time, self.vsb_r, label='vsb_r')
        plt.title('vsa_r and vsb_r over Time')
        plt.grid(True)
        plt.legend()
        
        plt.tight_layout()
        plt.show()

def example_ControleFluxoRotoricoEscorregamentoEstator():
    control = ControleFluxoRotoricoEscorregamentoEstator(regu_l=0)  # Set regu_l=1 or 2 to switch to other strategies
    control.compute()
    control.plot()

if __name__ == "__main__":
    example_ControleFluxoRotoricoEscorregamentoEstator()


# jm igor
class ControleFluxoEstatoricoEscorregamentoEstator:
    def __init__(self, regu_l=0, time_end=10, num_points=100):
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
                self.anga += (self.wr + self.war) * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid
                
                cga = np.cos(self.anga)
                sga = np.sin(self.anga)

                fsa_i = self.flsr * cga
                fsb_i = self.flsr * sga

                self.dei[0] = fsa_i - self.fsa
                self.dei[1] = fsb_i - self.fsb

                vsd = self.fkp * self.dei[0] + self.rs * isa
                vsq = self.fkp * self.dei[1] + self.rs * isb

                self.vsa_r.append(vsd)
                self.vsb_r.append(vsq)

            elif self.regu_l == 1:
                # PI controller strategy
                self.anga += (self.wr + self.war) * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid
                
                cga = np.cos(self.anga)
                sga = np.sin(self.anga)

                fsa_i = self.flsr * cga
                fsb_i = self.flsr * sga

                self.dei[0] = fsa_i - self.fsa
                self.dei[1] = fsb_i - self.fsb

                vsd = self.eki_d + (self.fkp + self.fki) * self.dei[0]
                self.eki_d += self.fki * self.dei[0]
                vsd -= self.esdl * self.fsa

                vsq = self.eki_q + (self.fkp + self.fki) * self.dei[1]
                self.eki_q += self.fki * self.dei[1]
                vsq -= self.esdl * self.fsb

                # Saturation
                if self.eki_d >= 180.:
                    self.eki_d = 180.
                if self.eki_d <= -180.:
                    self.eki_d = -180.
                if self.eki_q >= 180.:
                    self.eki_q = 180.
                if self.eki_q <= -180.:
                    self.eki_q = -180.

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


# iago alan
class ControleFluxoEstatoricoQuadratura:
    def __init__(self, regu_l=0, time_end=10, num_points=100):
        # Initialize constants
        self.regu_l = regu_l
        self.war = 0.0
        self.wr = 0.0
        self.lm = 1.0
        self.rr = 1.0
        self.flr = 1.0
        self.fkp = 1.0
        self.fki = 0.5
        self.rs = 0.1
        self.te = 1.0
        self.eki_d = 0.0
        self.eki_q = 0.0
        self.anga = 0.0
        self.pid = 2 * np.pi
        self.fsa = 0.0
        self.fsb = 0.0
        self.dei = [0.0, 0.0]
        self.flsr = 1.0
        self.esdl = 0.1
        self.fr_d = 0.0
        self.fr_q = 0.0

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
            isb = np.cos(t)

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
                frq = self.fr_q * cga - self.fr_d * sga

                self.dei[0] = self.flsr - fsd
                self.dei[1] = -fsq

                vsd = self.eki_d + (self.fkp + self.fki) * self.dei[0]
                self.eki_d += self.fki * self.dei[0]
                vsd -= self.esdl * frd - (self.war + self.wr) * fsq

                vsq = self.eki_q + (self.fkp + self.fki) * self.dei[1]
                self.eki_q += self.fki * self.dei[1]
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
        plt.figure(figsize=(12, 6))

        plt.plot(self.time, self.vsa_r, label='vsa_r')
        plt.plot(self.time, self.vsb_r, label='vsb_r')
        plt.title('Tensões de Referência vsa_r e vsb_r ao Longo do Tempo')
        plt.xlabel('Tempo')
        plt.ylabel('Tensão de Referência')
        plt.grid(True)
        plt.legend()
        plt.show()

def example_ControleFluxoEstatoricoQuadratura():
    control = ControleFluxoEstatoricoQuadratura(regu_l=1)  # Set regu_l=0 for P control or 1 for PI control
    control.compute()
    control.plot()

if __name__ == "__main__":
    example_ControleFluxoEstatoricoQuadratura()

    
# gabriel nikolas
class ControleFluxoRotoricoEscorregamentoCampo:
    def __init__(self, regu_l=0, fcor_r=0, time_end=10, num_points=100):
        # Initialize constants
        self.regu_l = regu_l
        self.fcor_r = fcor_r
        self.wbr = 0.0
        self.wr = 0.0
        self.lm = 1.0
        self.rr = 1.0
        self.flr = 1.0
        self.fkp = 1.0
        self.fki = 0.5
        self.rs = 0.1
        self.te = 1.0
        self.eki_d = 0.0
        self.eki_q = 0.0
        self.anga = 0.0
        self.pid = 2 * np.pi
        self.fsa = 0.0
        self.fsb = 0.0
        self.dei = [0.0, 0.0, 0.0, 0.0]
        self.flr_r = 1.0
        self.esdl = 0.1
        self.fr_d = 0.0
        self.fr_q = 0.0
        self.sigma = 0.01
        self.ls = 1.0

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
            isb = np.cos(t)

            if self.regu_l == 0:
                # P controller strategy with rotor slip
                self.anga += (self.wr + self.wbr) * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid

                cgb_r = np.cos(self.anga)
                sgb_r = np.sin(self.anga)

                frd = self.fr_d * cgb_r + self.fr_q * sgb_r
                frq = self.fr_q * cgb_r - self.fr_d * sgb_r

                self.dei[2] = self.flr_r - frd
                self.dei[3] = -frq

                isd_i = self.fkp * self.dei[2]
                isd_i -= (self.lm * self.wbr * frq) / (self.rr * self.lm)
                isq_i = self.fkp * self.dei[3]
                isq_i += (self.lm * self.wbr * frd) / (self.rr * self.lm)

                # Saturation
                isd_i = np.clip(isd_i, -55.0, 55.0)
                isq_i = np.clip(isq_i, -55.0, 55.0)

            elif self.regu_l == 1:
                # PI controller strategy with rotor slip
                self.anga += (self.wr + self.wbr) * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid

                cgb_r = np.cos(self.anga)
                sgb_r = np.sin(self.anga)

                frd = self.fr_d * cgb_r + self.fr_q * sgb_r
                frq = self.fr_q * cgb_r - self.fr_d * sgb_r

                self.dei[2] = self.flr_r - frd
                self.dei[3] = -frq

                isd_i = self.eki_d + (self.fkp + self.fki) * self.dei[2]
                self.eki_d += self.fki * self.dei[2]
                isd_i -= (self.lm * self.wbr * self.flr * sgb_r) / (self.rr * self.lm)

                isq_i = self.eki_q + (self.fkp + self.fki) * self.dei[3]
                self.eki_q += self.fki * self.dei[3]
                isq_i += (self.lm * self.wbr * self.flr * cgb_r) / (self.rr * self.lm)

                # Saturation
                self.eki_d = np.clip(self.eki_d, -55.0, 55.0)
                self.eki_q = np.clip(self.eki_q, -55.0, 55.0)

            if self.fcor_r == 0:
                # Stator current source
                self.dei[0] = isa - self.eki_d
                self.dei[1] = isb - self.eki_q

                vsd = self.eki_d + (self.fkp + self.fki) * self.dei[0] - (self.wr + self.wbr) * self.lm * self.flr * sgb_r / self.ls
                vsq = self.eki_q + (self.fkp + self.fki) * self.dei[1] + (self.wr + self.wbr) * self.lm * self.flr * cgb_r / self.ls

                self.eki_d += self.fki * self.dei[0]
                self.eki_q += self.fki * self.dei[1]

                vsa_r = vsd * cgb_r - vsq * sgb_r
                vsb_r = vsq * cgb_r + vsd * sgb_r

            else:
                # Synchronous current source
                isd_i = isa * cgb_r + isb * sgb_r
                isq_i = -isa * sgb_r + isb * cgb_r
                self.dei[0] = isd_i - isd_i
                self.dei[1] = isq_i - isq_i

                vsd = self.eki_d + (self.fkp + self.fki) * self.dei[0]
                self.eki_d += self.fki * self.dei[0]
                vsd -= (self.wbr + self.wr) * self.sigma * self.ls * isq_i - (self.lm * self.flr * self.rr) / (self.ls * self.ls)

                vsq = self.eki_q + (self.fkp + self.fki) * self.dei[1]
                self.eki_q += self.fki * self.dei[1]
                vsq += (self.wbr + self.wr) * self.sigma * self.ls * isd_i + (self.wr * self.lm * self.flr) / self.ls

                vsa_r = vsd * cgb_r - vsq * sgb_r
                vsb_r = vsq * cgb_r + vsd * sgb_r

            # Saturation
            vsa_r = np.clip(vsa_r, -311.0, 311.0)
            vsb_r = np.clip(vsb_r, -311.0, 311.0)

            self.vsd.append(vsd)
            self.vsq.append(vsq)
            self.vsa_r.append(vsa_r)
            self.vsb_r.append(vsb_r)

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

def example_ControleFluxoRotoricoEscorregamentoCampo():
    control = ControleFluxoRotoricoEscorregamentoCampo(regu_l=1, fcor_r=0)  # Set regu_l=0 for P control or 1 for PI control
    control.compute()
    control.plot()

if __name__ == "__main__":
    example_ControleFluxoRotoricoEscorregamentoCampo()


# samuel leo
class ControleFluxoRotoricoEscorregamentoRotor:
    def __init__(self, regu_l=0, fcor_r=0, time_end=10, num_points=100):
        # Initialize constants
        self.regu_l = regu_l
        self.fcor_r = fcor_r
        self.wbr = 0.0
        self.wr = 0.0
        self.lm = 1.0
        self.rr = 1.0
        self.flr = 1.0
        self.fkp = 1.0
        self.fki = 0.5
        self.rs = 0.1
        self.te = 1.0
        self.eki_d = 0.0
        self.eki_q = 0.0
        self.anga = 0.0
        self.pid = 2 * np.pi
        self.fsa = 0.0
        self.fsb = 0.0
        self.dei = [0.0, 0.0, 0.0, 0.0]
        self.flr_r = 1.0
        self.esdl = 0.1
        self.fr_d = 0.0
        self.fr_q = 0.0
        self.sigma = 0.01
        self.ls = 1.0
        self.ctt = 1.0

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
            isb = np.cos(t)

            if self.regu_l == 0:
                # P controller strategy
                self.anga += self.wbr * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid

                cgb_r = np.cos(self.anga)
                sgb_r = np.sin(self.anga)
                ct_t = self.ctt
                st_t = np.sin(ct_t)

                frd = self.flr * ct_t + self.flr * st_t
                frq = self.flr * ct_t - self.flr * st_t

                frd_i = self.flr_r * cgb_r
                frq_i = self.flr_r * sgb_r

                self.dei[2] = frd_i - frd
                self.dei[3] = frq_i - frq

                isd_i = self.fkp * self.dei[2]
                isq_i = self.fkp * self.dei[3]

            elif self.regu_l == 1:
                # PI controller strategy
                self.anga += self.wbr * self.te
                if self.anga >= self.pid:
                    self.anga -= self.pid

                cgb_r = np.cos(self.anga)
                sgb_r = np.sin(self.anga)
                ct_t = self.ctt
                st_t = np.sin(ct_t)

                frd = self.flr * ct_t + self.flr * st_t
                frq = self.flr * ct_t - self.flr * st_t

                frd_i = self.flr_r * cgb_r
                frq_i = self.flr_r * sgb_r

                self.dei[2] = frd_i - frd
                self.dei[3] = frq_i - frq

                isd_i = self.eki_d + (self.fkp + self.fki) * self.dei[2]
                self.eki_d += self.fki * self.dei[2]
                isd_i = np.clip(isd_i, -15.0, 15.0)

                isq_i = self.eki_q + (self.fkp + self.fki) * self.dei[3]
                self.eki_q += self.fki * self.dei[3]
                isq_i = np.clip(isq_i, -15.0, 15.0)

            # Coordinate transformation
            isa_i = isd_i * np.cos(ct_t) - isq_i * np.sin(ct_t)
            isb_i = isq_i * np.cos(ct_t) + isd_i * np.sin(ct_t)

            if self.fcor_r == 0:
                # Stator current source
                self.dei[0] = isa_i - isa
                self.dei[1] = isb_i - isb

                vsd = self.eki_d + (self.fkp + self.fki) * self.dei[0] - (self.wr + self.wbr) * self.lm * self.flr * sgb_r / self.rr
                vsq = self.eki_q + (self.fkp + self.fki) * self.dei[1] + (self.wr + self.wbr) * self.lm * self.flr * cgb_r / self.rr

                self.eki_d += self.fki * self.dei[0]
                self.eki_q += self.fki * self.dei[1]

                vsa_r = vsd
                vsb_r = vsq

            else:
                # Synchronous current source
                isd_i = isa_i * cgb_r + isb_i * sgb_r
                isq_i = -isa_i * sgb_r + isb_i * cgb_r

                self.dei[0] = isd_i - isd_i
                self.dei[1] = isq_i - isq_i

                vsd = self.eki_d + (self.fkp + self.fki) * self.dei[0]
                self.eki_d += self.fki * self.dei[0]
                vsd -= (self.wbr + self.wr) * self.sigma * self.ls * isq_i - (self.lm * self.flr * self.rr) / (self.rr * self.rr)

                vsq = self.eki_q + (self.fkp + self.fki) * self.dei[1]
                self.eki_q += self.fki * self.dei[1]
                vsq += (self.wbr + self.wr) * self.sigma * self.ls * isd_i + (self.wr * self.lm * self.flr) / self.rr

                vsa_r = vsd * cgb_r - vsq * sgb_r
                vsb_r = vsq * cgb_r + vsd * sgb_r

            # Saturation
            vsa_r = np.clip(vsa_r, -311.0, 311.0)
            vsb_r = np.clip(vsb_r, -311.0, 311.0)

            self.vsd.append(vsd)
            self.vsq.append(vsq)
            self.vsa_r.append(vsa_r)
            self.vsb_r.append(vsb_r)

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

def example_ControleFluxoRotoricoEscorregamentoRotor():
    control = ControleFluxoRotoricoEscorregamentoRotor(regu_l=1, fcor_r=0)  # Set regu_l=0 for P control or 1 for PI control
    control.compute()
    control.plot()

if __name__ == "__main__":
    example_ControleFluxoRotoricoEscorregamentoRotor()



class ModulacaoPWM:
    def __init__(self, te=0.01, rq3=1.0, ec=1.0, pid=2*np.pi, num_points=100):
        self.te = te
        self.rq3 = rq3
        self.ec = ec
        self.pid = pid
        self.num_points = num_points

        # Initialize variables
        self.vsar = 1.0
        self.vsbr = 1.0
        self.vsc = 1.0
        self.cont = 0.0
        self.kor = 0
        self.kli_m = 0
        self.se = 0
        self.kst2 = 0
        self.tk = 0.0
        self.tkl = 0.0
        self.to = 0.0

        # Storage for output variables
        self.vf = np.zeros((3, num_points))
        self.vonf = np.zeros(num_points)
        self.vslf = np.zeros(num_points)
        self.vs2f = np.zeros(num_points)
        self.vs3f = np.zeros(num_points)

        # Time vector
        self.time = np.linspace(0, te*(num_points-1), num_points)
        self.index = 0

    def calcular_pwm(self):
        for t in self.time:
            # Update PWM parameters
            self.cont += 1
            self.kor = (self.kor + 1) % 2
            tck = 2.0 * self.te / self.rq3

            # Determine sector
            vsa_r = self.vsar / self.rq3
            vsb_r = self.vsbr / self.rq3
            aux = [self.rq3 * vsa_r + vsb_r, -self.rq3 * vsa_r + vsb_r, -2.0 * vsb_r]
            
            if aux[0] >= 0 and aux[1] <= 0 and aux[2] <= 0:
                self.se = 1
            elif aux[0] >= 0 and aux[1] >= 0 and aux[2] <= 0:
                self.se = 2
            elif aux[0] <= 0 and aux[1] >= 0 and aux[2] <= 0:
                self.se = 3
            elif aux[0] <= 0 and aux[1] >= 0 and aux[2] >= 0:
                self.se = 4
            elif aux[0] <= 0 and aux[1] <= 0 and aux[2] >= 0:
                self.se = 5
            elif aux[0] >= 0 and aux[1] <= 0 and aux[2] >= 0:
                self.se = 6

            # Calculate application times
            self.tk = tck * (self.vsar * np.sin(self.se * np.pi / 3.0) - self.vsbr * np.cos(self.se * np.pi / 3.0)) / self.ec
            self.tkl = -tck * (self.vsar * np.sin((self.se - 1) * np.pi / 3.0) - self.vsbr * np.cos((self.se - 1) * np.pi / 3.0)) / self.ec
            self.tk = max(self.tk, 0.0)
            self.tkl = max(self.tkl, 0.0)
            self.to = self.te - self.tk - self.tkl
            self.to = min(max(self.to, 0.0), self.te)

            # Determine instantaneous voltage vectors
            if self.kli_m == 0:
                if self.cont <= self.to / 2.0:
                    self.kst2 = 7
                else:
                    if self.kor == 0:
                        if self.cont < (self.to / 2.0 + self.tk):
                            self.kst2 = self.se
                        elif self.cont < (self.to / 2.0 + self.tk + self.tkl):
                            self.kst2 = self.se + 1
                            if self.kst2 > 6:
                                self.kst2 = 1
                        else:
                            self.kst2 = 7
                    else:
                        if self.cont < (self.to / 2.0 + self.tkl):
                            self.kst2 = self.se + 1
                            if self.kst2 > 6:
                                self.kst2 = 1
                        elif self.cont < (self.to / 2.0 + self.tk + self.tkl):
                            self.kst2 = self.se
                        else:
                            self.kst2 = 7

            # Calculate output variables
            gtf = np.zeros(3)
            if self.kst2 == 1:
                gtf = [1, 0, 0]
            elif self.kst2 == 2:
                gtf = [1, 1, 0]
            elif self.kst2 == 3:
                gtf = [0, 1, 0]
            elif self.kst2 == 4:
                gtf = [0, 1, 1]
            elif self.kst2 == 5:
                gtf = [0, 0, 1]
            elif self.kst2 == 6:
                gtf = [1, 0, 1]
            elif self.kst2 == 7:
                gtf = [0, 0, 0]

            vf = np.array([gtf[0] * (self.ec / 2) + (gtf[0] - 1) * (self.ec / 2),
                           gtf[1] * (self.ec / 2) + (gtf[1] - 1) * (self.ec / 2),
                           gtf[2] * (self.ec / 2) + (gtf[2] - 1) * (self.ec / 2)])
            vonf = (1.0 / 3.0) * np.sum(vf)
            self.vslf[self.index] = vf[0] - vonf
            self.vs2f[self.index] = vf[1] - vonf
            self.vs3f[self.index] = vf[2] - vonf

            self.index += 1

    def plot_results(self):
        plt.figure(figsize=(12, 6))

        plt.plot(self.time, self.vslf, label='vslf')
        plt.plot(self.time, self.vs2f, label='vs2f')
        plt.plot(self.time, self.vs3f, label='vs3f')
        plt.title('Resultados da Modulação PWM')
        plt.xlabel('Tempo')
        plt.ylabel('Tensão')
        plt.grid(True)
        plt.legend()
        plt.show()

def example_ModulacaoPWM():
    pwm = ModulacaoPWM(te=0.01, rq3=1.0, ec=1.0, pid=2*np.pi, num_points=100)
    pwm.calcular_pwm()
    pwm.plot_results()

if __name__ == "__main__":
    example_ModulacaoPWM()

