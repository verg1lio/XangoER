import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class Inversor:
    def __init__(self, V_m, f, V_dc, i_dc, phi, theta, f_s, m, peso_inv=None):
        self.V_m = V_m # Tensão de entrada RMS (V)
        self.f = f # Frequência de entrada (Hz)
        self.V_dc = V_dc # Tensão CC (V)
        self.i_dc = i_dc # Corrente CC (A)
        self.phi = phi # Ângulo de fase da tensão de entrada (rad)
        self.theta = theta # Ângulo de fase da corrente de saída (rad)
        self.f_s = f_s # Frequência de chaveamento (Hz)
        self.m = m # Índice de modulação
        self.peso_inv = peso_inv # peso da inversor (kg)

        # Cálculo do período de chaveamento
        self.T_s = 1 / f_s 

        # Cálculo do ângulo de modulação
        self.theta_m = np.arcsin(m)

        # Cálculo da frequência angular de entrada
        self.omega = 2 * np.pi * f

        # Cálculo da tensão de saída fundamental
        self.V_o1 = (2 * V_dc / np.pi) * m

    def gerar_tensoes_saida(self, t):

        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
        Description.

        Examples
        --------
        =>>> example
        """        

        u_a = self._gerar_funcao_comutacao(t, self.theta_m)
        u_b = self._gerar_funcao_comutacao(t, self.theta_m + 2 * np.pi / 3)
        u_c = self._gerar_funcao_comutacao(t, self.theta_m + 4 * np.pi / 3)

        v_sw_a = self.V_dc * u_a
        v_sw_b = self.V_dc * u_b
        v_sw_c = self.V_dc * u_c

        v_a = v_sw_a * np.sin(self.omega * t - self.phi)
        v_b = v_sw_b * np.sin(self.omega * t - self.phi - 2 * np.pi / 3)
        v_c = v_sw_c * np.sin(self.omega * t - self.phi + 2 * np.pi / 3)

        return v_a, v_b, v_c

    def gerar_correntes_saida(self, t):

        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
        Description.

        Examples
        --------
        =>>> example
        """

        # Utilizando a mesma forma das tensões, por simplicidade
        i_a = self.i_dc * np.sin(self.omega * t - self.theta + np.pi/2)
        i_b = self.i_dc * np.sin(self.omega * t - self.theta - 2 * np.pi / 3 + np.pi/2)
        i_c = self.i_dc * np.sin(self.omega * t - self.theta + 2 * np.pi / 3 + np.pi/2)
        return i_a, i_b, i_c

    def calcular_potencias(self, v_a, v_b, v_c, i_a, i_b, i_c):

        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
        Description.

        Examples
        --------
        =>>> example
        """

        V_ef = np.sqrt(np.mean(v_a**2 + v_b**2 + v_c**2))
        I_ef = np.sqrt(np.mean(i_a**2 + i_b**2 + i_c**2))
        P = 3 * V_ef * I_ef * np.cos(0)
        Q = np.mean(v_a * i_a * np.sin(self.phi)) + np.mean(v_b * i_b * np.sin(self.phi)) + np.mean(v_c * i_c * np.sin(self.phi))
        S = np.sqrt(P**2 + Q**2)
        return P, Q, S

    def _gerar_funcao_comutacao(self, t, theta_m):

        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
        Description.

        Examples
        --------
        =>>> example
        """

        k = np.floor((t + self.T_s / 4) / self.T_s)
        u = (t < (k * self.T_s + theta_m)).astype(int)
        return u

def plotar_inversor(t, v_a, v_b, v_c, i_a, i_b, i_c, Pa, Q, S, lambda_m):
    plt.figure(figsize=(10, 12))
    plt.suptitle('Parâmetros de Saída do Inversor')

    inversorCustomPlot(t, [v_a, v_b, v_c], ['Tensão Fase A', 'Tensão Fase B', 'Tensão Fase A'], ylabel='Tensão (V)', subplot=1)
    inversorCustomPlot(t, [i_a, i_b, i_c], ['Corrente Fase A', 'Corrente Fase B', 'Corrente Fase A'], ylabel='Corrente (A)', subplot=2)
    inversorCustomPlot(
        t, 
        [np.full_like(t, Pa), np.full_like(t, Q), np.full_like(t, S)], 
        ['Potência Ativa (W)', 'Potência Reativa (VAR)', 'Potência Aparente (VA)'], 
        ylabel='Potência', 
        subplot=3
    )
    inversorCustomPlot(t, [lambda_m], ['Fluxo de Magnetização ($lambda_m$)'],ylabel='Fluxo de Magnetização (Wb)', subplot=4)

    plt.xlabel('Tempo (s)')
    plt.grid(True)
    plt.show()
    # plt.savefig("inversor.png")

def deriv(y, t, V_d, V_q, omega):
    R_s = 0.435 # resistência do stator
    I_d, I_q, lambda_d, lambda_q = y
    
    d_lambda_d_dt = V_d - R_s * I_d + omega * lambda_q
    d_lambda_q_dt = V_q - R_s * I_q - omega * lambda_d
    
    return [I_d, I_q, d_lambda_d_dt, d_lambda_q_dt]

def inversorCustomPlot(x, ys, titles, ylabel, subplot):
    plt.subplot(4, 1, subplot)
    for i in range(len(ys)):
        plt.plot(x, ys[i], label=titles[i])

    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)

def example_inversor():
    y0 = [0, 0, 0, 0] # I_d0, I_q0, lambda_d0, lambda_q0
    omega = 2 * np.pi * 60 # velocidade síncrona(rad/s)
    
    # Tensão aplicada (exemplo)
    V_d = 10
    V_q = 10

    t_sim = 0.2
    t = np.linspace(0, t_sim, 1000)

    # Resolver as equações diferenciais
    sol = odeint(deriv, y0, t, args=(V_d, V_q, omega))
    
    lambda_d = sol[:, 2]
    lambda_q = sol[:, 3]
    lambda_m = np.sqrt(lambda_d**2 + lambda_q**2)

    my_inversor = Inversor(220, 60, 220, 118, 0, 0, 1, 0.01)

    v_a, v_b, v_c = my_inversor.gerar_tensoes_saida(t) # as tensões aplicadas não deviam ser calculadas a partir destas?
    i_a, i_b, i_c = my_inversor.gerar_correntes_saida(t)
    Pa, Q, S = my_inversor.calcular_potencias(v_a, v_b, v_c, i_a, i_b, i_c)

    plotar_inversor(t, v_a, v_b, v_c, i_a, i_b, i_c, Pa, Q, S, lambda_m)

def plotar_potencia_vs_frequencia():
    f_s_range = np.linspace(4, 76, num=10000)
    P_list = []
    
    for f_s in f_s_range:
        my_inversor = Inversor(220, 60, 220, 118, 0, 0, f_s, 0.01)
        
        t_sim = 0.2
        t = np.linspace(0, t_sim, 1000)
        
        v_a, v_b, v_c = my_inversor.gerar_tensoes_saida(t)
        i_a, i_b, i_c = my_inversor.gerar_correntes_saida(t)
        P, Q, S = my_inversor.calcular_potencias(v_a, v_b, v_c, i_a, i_b, i_c)
        
        P_list.append(P)
    
    plt.figure()
    plt.plot(f_s_range, P_list)
    #plt.xscale('log')
    plt.xlabel('Frequência de Chaveamento (Hz)')
    plt.ylabel('Potência Ativa (W)')
    plt.title('Potência Ativa vs Frequência de Chaveamento')
    plt.grid(True)
    plt.show()

plotar_potencia_vs_frequencia()

example_inversor()