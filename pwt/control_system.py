import numpy as np
import matplotlib.pyplot as plt

class PWT:
    def __init__(self, k, amp_fluxo, r_s, v_s, i_s, r_r, v_r, i_r, ref_torque, ref_fluxo_rotor, ref_fluxo_estator, t):
        self.k = k  # Parâmetro da máquina
        self.amp_fluxo = amp_fluxo  # Amplitude do fluxo
        self.r_s = r_s  # Resistência do estator
        self.v_s = v_s  # Tensão do estator
        self.i_s = i_s  # Corrente do estator
        self.r_r = r_r  # Resistência do rotor
        self.v_r = v_r  # Tensão do rotor
        self.i_r = i_r  # Corrente do rotor
        self.ref_torque = ref_torque  # Referência de torque
        self.ref_fluxo_rotor = ref_fluxo_rotor  # Referência de fluxo do rotor
        self.ref_fluxo_estator = ref_fluxo_estator  # Referência de fluxo do estator
        self.t = t  # Vetor de tempo
        self.dt = t[1] - t[0]  # Intervalo de tempo
        self.w1 = 0  # Frequência do estator
        self.wr = 0  # Frequência do rotor
        self.ce = 0  # Torque eletromagnético
        self.psi_r = np.zeros_like(t)  # Inicialização do fluxo do rotor
        self.psi_s = np.zeros_like(t)  # Inicialização do fluxo do estator
        self.erro_torque = np.zeros_like(t)  # Inicialização do erro de torque
        self.erro_fluxo_rotor = np.zeros_like(t)  # Inicialização do erro de fluxo do rotor
        self.erro_fluxo_estator = np.zeros_like(t)  # Inicialização do erro de fluxo do estator

    def calcular_fluxo_rotorico(self):
        # Cálculo do fluxo do rotor (integração da tensão do rotor)
        self.psi_r = np.cumsum((self.v_r - self.i_r * self.r_r)*self.dt)

    def calcular_fluxo_estatorico(self):
        # Cálculo do fluxo do estator (integração da tensão do estator)
        self.psi_s = np.cumsum((self.v_s - self.i_s * self.r_s)*self.dt)

    def calcular_torque(self):
        # Cálculo do conjugado eletromagnético
        self.ce = self.k * (self.psi_r ** 2) * (self.w1 - self.wr)

    def controle(self):
        self.calcular_fluxo_rotorico()
        self.calcular_fluxo_estatorico()
        self.calcular_torque()
        # Diferença entre o torque de referência e o torque calculado
        self.erro_torque = self.ref_torque - self.ce
        # Diferença entre o fluxo de referência e o fluxo rotorico calculado
        self.erro_fluxo_rotor = self.ref_fluxo_rotor - self.psi_r
        # Diferença entre o fluxo de referência e o fluxo estatorico calculado
        self.erro_fluxo_estator = self.ref_fluxo_estator - self.psi_s
        # Seleção de vetor de tensão baseado nos erros (lógica simplificada)
        v_s_corrigido_rotor = self.v_r + self.erro_torque * 0.1 + self.erro_fluxo_rotor * 0.1 #perguntar
        v_s_corrigido_estator = self.v_s + self.erro_torque * 0.1 + self.erro_fluxo_estator * 0.1 #perguntar
        return v_s_corrigido_rotor, v_s_corrigido_estator, self.erro_torque, self.erro_fluxo_rotor, self.erro_fluxo_estator

# Plotar resultados
def customPlot(subplot, x, y, label, xlabel, ylabel):
    plt.subplot(5, 1, subplot)
    plt.plot(x, y, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)

def generatePlots(t, pwt, erro_torque, erro_fluxo_rotor, saveFig):
    plt.figure()

    # Gráfico do erro do torque
    customPlot(1, t, erro_torque, 'Erro de Torque', 'Tempo (s)', 'Erro de Torque')
    # Gráfico do erro do fluxo rotorico
    customPlot(2, t, erro_fluxo_rotor, 'Erro de Fluxo Rotorico', 'Tempo (s)', 'Erro de Fluxo Rotorico')
    # Gráfico do conjugado eletromagnético
    customPlot(3, t, pwt.ce, 'Conjugado Eletromagnético', 'Tempo (s)', 'Torque')
    # Gráfico do fluxo do rotor
    customPlot(4, t, pwt.psi_r, 'Fluxo do Rotor', 'Tempo (s)', 'Fluxo do Rotor')
    # Gráfico do fluxo do estator
    customPlot(5, t, pwt.psi_s, 'Fluxo do Estator', 'Tempo (s)', 'Fluxo do Estator')

    plt.tight_layout()

    if (saveFig): 
        plt.savefig('pwt.png')
    else:
        plt.show()

def example(saveToFile = False):
    # Simulação
    t = np.linspace(0, 1, 1000)
    v_s = np.sin(2 * np.pi * 50 * t)  # Tensão do estator (exemplo)
    i_s = np.sin(2 * np.pi * 50 * t)  # Corrente do estator (exemplo)
    v_r = np.sin(2 * np.pi * 50 * t)  # Tensão do rotor (exemplo)
    i_r = np.sin(2 * np.pi * 50 * t)  # Corrente do rotor (exemplo)
    ref_torque = 1  # Referência de torque
    ref_fluxo_rotor = 1  # Referência de fluxo do rotor
    ref_fluxo_estator = 1  # Referência de fluxo do estator

    # Instanciar o objeto PWT com os parâmetros necessários
    pwt = PWT(k=1, amp_fluxo=1, r_s=0.1, v_s=v_s, i_s=i_s, r_r=0.1, v_r=v_r, i_r=i_r, ref_torque=ref_torque, ref_fluxo_rotor=ref_fluxo_rotor, ref_fluxo_estator=ref_fluxo_estator, t=t)

    # Controle
    v_s_corrigido_rotor, v_s_corrigido_estator, erro_torque, erro_fluxo_rotor, erro_fluxo_estator = pwt.controle()
    generatePlots(t, pwt, erro_torque, erro_fluxo_rotor, saveToFile)

example(saveToFile = True)