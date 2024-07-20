import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import cmath  # Para operações com números complexos
from matplotlib.widgets import Slider
from scipy import signal

__all__ = ["bateria", "motor_gaiola", "inversor", "controle_fluxo", "controle_inversor", "peso"]

class bateria:
    def __init__(self, bat_type, bat_Idis, bat_Ich, bat_t, bat_Es, bat_K, bat_Q, bat_A, bat_B, bat_Iast, bat_R, time=None):
        self.bat_type = 'Modified Shepherd' if bat_type == 'S' else 'Rint model' if bat_type == 'R' else bat_type
        self.bat_Idis = bat_Idis # Current discharge
        self.bat_Ich = bat_Ich # Current charge
        self.bat_t = bat_t # Sampling time
        self.bat_Es = bat_Es # Open-circuit voltage
        self.bat_K = bat_K # Change in polarization resistance factor (mΩ A/h)
        self.bat_Q = bat_Q # Battery nominal capacity (Ah)
        self.bat_A = bat_A # Voltage factor
        self.bat_B = bat_B # Capacity factor
        self.bat_Iast = bat_Iast # Filtered current
        self.bat_R = bat_R # Internal resistance

        if time is None:
            self.time = self.bat_Q / self.bat_Idis  # A*h/A = h
        else:
            self.time = time

    def battery(self):
        """Bateria do sistema de powertrain

        Função que determinará os valores de tensão de carga e tensão de descarga da bateria de acordo com pontos do tempo.

        Returns
        -------
        time_points : np.ndarray
            Array do NumPy contendo 100 valores, entre '0' e 'self.time'.
        discharge_voltages : list
            Valores das tensões de descarga
        charge_voltages : list
            Valores das tensões de carga
        
        Examples
        --------
        >>> my_battery = bateria('S', 10, 5, 1, 12, 0.1, 100, 1, 0.01, 1, 0.05)
        >>> time, discharge_voltage, charge_voltage = my_battery.battery()
            print(time, discharge_voltage, charge_voltage)
        """
        time_points = np.linspace(0, self.time, 100)  # 100 time points
        discharge_voltages = []
        charge_voltages = []

        for t in time_points:
            if self.bat_type == 'Modified Shepherd':
                V_dis = self.bat_Es - self.bat_R * self.bat_Idis - self.bat_K * (self.bat_Q / (self.bat_Q - (self.bat_Idis * t))) * ((self.bat_Idis * t) + self.bat_Iast) + self.bat_A * np.exp( -self.bat_B * self.bat_Idis * t)
                V_ch = self.bat_Es - self.bat_R * self.bat_Ich - self.bat_K * (self.bat_Q / ((self.bat_Ich * t) - 0.1 * self.bat_Q)) * self.bat_Iast - self.bat_K * (((self.bat_Q) / (self.bat_Q - (self.bat_Ich * t))) * (self.bat_Ich * t)) + self.bat_A * np.exp(-self.bat_B * self.bat_Ich * t)
           
            elif self.bat_type == 'Rint model':
                V_dis = self.bat_Es - self.bat_R * self.bat_Idis
                V_ch = self.bat_Es - self.bat_R * self.bat_Ich

            discharge_voltages.append(V_dis)
            charge_voltages.append(V_ch)

        return time_points, discharge_voltages, charge_voltages

def plot_bateria(time, discharge_voltage, charge_voltage):
    plt.plot(time, discharge_voltage, label='Discharge Voltage')
    plt.plot(time, charge_voltage, label='Charge Voltage')
    plt.xlabel('Time')
    plt.ylabel('Voltage')
    plt.title('Battery Voltage over Time')
    plt.legend()
    plt.grid(True)
    plt.show()

def example_bateria():

    my_battery = bateria('S', 10, 5, 1, 12, 0.1, 100, 1, 0.01, 1, 0.05)
    time, discharge_voltage, charge_voltage = my_battery.battery()

    plot_bateria(time, discharge_voltage, charge_voltage)



class inversor:
    def __init__(self, V_m, f, V_dc, i_dc, phi, theta, f_s, m):
        self.V_m = V_m # Tensão de entrada RMS (V)
        self.f = f # Frequência de entrada (Hz)
        self.V_dc = V_dc # Tensão CC (V)
        self.i_dc = i_dc # Corrente CC (A)
        self.phi = phi # Ângulo de fase da tensão de entrada (rad)
        self.theta = theta # Ângulo de fase da corrente de saída (rad)
        self.f_s = f_s # Frequência de chaveamento (Hz)
        self.m = m # Índice de modulação
        
        self.T_s = 1 / f_s # Cálculo do período de chaveamento
        self.theta_m = np.arcsin(m) # Cálculo do ângulo de modulação
        self.omega = 2 * np.pi * f # Cálculo da frequência angular de entrada
        self.V_o1 = (2 * V_dc / np.pi) * m # Cálculo da tensão de saída fundamental

    def gerar_tensoes_saida(self, t):
        """Gera as tensõs de saída do inversor.

        Returns
        -------
        v_a : np.ndarray
            Tensão da Fase A do inversor
        v_b : np.ndarray
            Tensão da Fase B do inversor
        v_c : np.ndarray
            Tensão da Fase C do inversor

        Examples
        --------
        >>> t_sim = 0.2
        >>> t = np.linspace(0, t_sim, 1000)
        >>> tensoes = Inversor(220, 60, 220, 118, 0, 0, 1, 0.01)
        >>> v_a, v_b, v_c = tensoes.gerar_tensoes_saida(t)
            [...]
        """
        u_a = self.gerar_funcao_comutacao(t, self.theta_m)
        u_b = self.gerar_funcao_comutacao(t, self.theta_m + 2 * np.pi / 3)
        u_c = self.gerar_funcao_comutacao(t, self.theta_m + 4 * np.pi / 3)

        v_sw_a = self.V_dc * u_a
        v_sw_b = self.V_dc * u_b
        v_sw_c = self.V_dc * u_c

        v_a = v_sw_a * np.sin(self.omega * t - self.phi)
        v_b = v_sw_b * np.sin(self.omega * t - self.phi - 2 * np.pi / 3)
        v_c = v_sw_c * np.sin(self.omega * t - self.phi + 2 * np.pi / 3)
        return v_a, v_b, v_c

    def gerar_correntes_saida(self, t):
        """Gera as correntes de saída do inversor.
        
        Returns
        -------
        i_a : np.ndarray
            Corrente da Fase A do inversor
        i_b : np.ndarray
            Corrente da Fase B do inversor
        i_c : np.ndarray
            Corrente da Fase C do inversor.
      
        Examples
        --------
        >>> t_sim = 0.2
        >>> t = np.linspace(0, t_sim, 1000)
        >>> correntes = Inversor(220, 60, 220, 118, 0, 0, 1, 0.01)
            i_A, i_B, i_B = correntes.gerar_correntes_saida(t)
            [...]
        """
        # Utilizando a mesma forma das tensões, por simplicidade
        i_a = self.i_dc * np.sin(self.omega * t - self.theta + np.pi/2)
        i_b = self.i_dc * np.sin(self.omega * t - self.theta - 2 * np.pi / 3 + np.pi/2)
        i_c = self.i_dc * np.sin(self.omega * t - self.theta + 2 * np.pi / 3 + np.pi/2)
        return i_a, i_b, i_c

    def calcular_potencias(self, v_a, v_b, v_c, i_a, i_b, i_c):
        """Calculo da potencia do inversor.

        Returns
        -------
        P : float
            Potencia ativa.
        Q : float
            Potencia reativa.
        S : float
            Potencia aparente. 

        Examples
        --------
        >>> v_a, v_b, v_c, i_a, i_b, i_c = (5, 6, 9, 18, 10, 14)
        >>> potencias = Inversor(220, 60, 220, 118, 0, 0, 1, 0.01)
        >>> p_P, p_Q, p_S = potencias.calcular_potencias(v_a, v_b, v_c, i_a, i_b, i_c)
            890.1460554313545, 0.0, 890.1460554313545
        """
        V_ef = np.sqrt(np.mean(v_a**2 + v_b**2 + v_c**2))
        I_ef = np.sqrt(np.mean(i_a**2 + i_b**2 + i_c**2))
        P = 3 * V_ef * I_ef * np.cos(0)
        Q = np.mean(v_a * i_a * np.sin(self.phi)) + np.mean(v_b * i_b * np.sin(self.phi)) + np.mean(v_c * i_c * np.sin(self.phi))
        S = np.sqrt(P**2 + Q**2)
        return P, Q, S

    def gerar_funcao_comutacao(self, t, theta_m):
        """Função de comutação.

        Gera a função responsável por ligar e desligar a chave eletrônica de potencia.

        Returns
        -------
        u : np.ndarray
            Defini a corrente como contínua ou alternada.

        Examples
        --------
        >>> t_sim = 0.2
        >>> t = np.linspace(0, t_sim, 1000)
        >>> theta_m = np.pi * 45
        >>> comutacao = Inversor(220, 60, 220, 118, 0, 0, 1, 0.01)
        >>> exemple_u = comutacao._gerar_funcao_comutacao(t, theta_m)
            [1, 1, ..., 1]
        """
        k = np.floor((t + self.T_s / 4) / self.T_s)
        u = (t < (k * self.T_s + theta_m)).astype(int)
        return u

def deriv(y, t, V_d, V_q, omega):
        """Calculo da derivada dos fluxos magneticos.

        Calcula as derivadas das componentes do fluxo magnético no referencial dq e as componentes da corrente no eixo d e q em relação ao tempo.

        Returns
        -------
        I_d : float
            Componente da corrente no eixo d.
        I_q : float
            Componente da corrente no eixo q.
        d_lambda_d_dt : float
            Derivada da componente do fluxo magnético no eixo d.
        d_lambda_q_dt : float
            Derivada da componente do fluxo magnético no eixo q.

        Examples
        --------
        >>> t_sim = 0.2
        >>> t = np.linspace(0, t_sim, 1000)
        >>> y = [0, 0, 0, 0] 
        >>> V_d, V_q, omega = (10, 10, 2 * np.pi * 60)
        >>> derivadas = inversor(220, 60, 220, 118, 0, 0, 1, 0.01)
        >>> exemple_id, exemple_iq, exemple_lambda_d, exemple_lambda_q = derivadas.deriv(y, t, V_d, V_q, omega)
            [0, 0, 10.0, 10.0]
        """
        R_s = 0.435 # resistência do stator
        I_d, I_q, lambda_d, lambda_q = y
        
        d_lambda_d_dt = V_d - R_s * I_d + omega * lambda_q
        d_lambda_q_dt = V_q - R_s * I_q - omega * lambda_d
        return [I_d, I_q, d_lambda_d_dt, d_lambda_q_dt]

def plotar_inversor(t, v_a, v_b, v_c, i_a, i_b, i_c, Pa, Q, S, lambda_m):
    plt.figure(figsize=(10, 12))
    plt.suptitle('Parâmetros de Saída do Inversor')

    plt.subplot(4, 1, 1)
    plt.plot(t, v_a, label='Tensão Fase A')
    plt.plot(t, v_b, label='Tensão Fase B')
    plt.plot(t, v_c, label='Tensão Fase C')
    plt.ylabel('Tensão (V)')
    plt.legend()
    plt.grid(True)

    plt.subplot(4, 1, 2)
    plt.plot(t, i_a, label='Corrente Fase A')
    plt.plot(t, i_b, label='Corrente Fase B')
    plt.plot(t, i_c, label='Corrente Fase C')
    plt.ylabel('Corrente (A)')
    plt.legend()
    plt.grid(True)

    plt.subplot(4, 1, 3)
    plt.plot(t, np.full_like(t, Pa), label='Potência Ativa (W)')
    plt.plot(t, np.full_like(t, Q), label='Potência Reativa (VAR)')
    plt.plot(t, np.full_like(t, S), label='Potência Aparente (VA)')
    plt.ylabel('Potência')
    plt.legend()
    plt.grid(True)

    plt.subplot(4, 1, 4)
    plt.plot(t, lambda_m, label='Fluxo de Magnetização ($lambda_m$)')
    plt.ylabel('Fluxo de Magnetização (Wb)')
    plt.legend()
    plt.grid(True)

    plt.xlabel('Tempo (s)')
    plt.grid(True)
    plt.show()

def plotar_potencia_vs_frequencia():
    f_s_range = np.linspace(4, 76, num=10000)
    P_list = []
    
    for f_s in f_s_range:
        my_inversor = inversor(220, 60, 220, 118, 0, 0, f_s, 0.01)
        
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

        my_inversor = inversor(220, 60, 220, 118, 0, 0, 1, 0.01)

        v_a, v_b, v_c = my_inversor.gerar_tensoes_saida(t)
        i_a, i_b, i_c = my_inversor.gerar_correntes_saida(t)
        Pa, Q, S = my_inversor.calcular_potencias(v_a, v_b, v_c, i_a, i_b, i_c)

        plotar_inversor(t, v_a, v_b, v_c, i_a, i_b, i_c, Pa, Q, S, lambda_m)



class motor_gaiola:
    def __init__(self, frequencia, P, R1, X1, R2, X2, Xm, K, V_m):
        self.frequencia = frequencia # frequência em Hz
        self.P = P  # Número de polos
        self.R1 = R1 # Resistência do estator
        self.X1 = X1 # Reatância do estator
        self.R2 = R2 # Resistência do rotor
        self.X2 = X2 # Reatância do rotor
        self.Xm = Xm # Reatância magnética
        self.K = K  # Constante de proporcionalidade para o torque
        self.w_s = 2 * np.pi * self.frequencia / (self.P / 2)  # Velocidade síncrona
        self.V_m = V_m  # Tensão de entrada RMS (V)

    def calcular_impedancia(self, s):
        """Calculo da impedância do motor.

        Representa a oposição que um circuito elétrico oferece à passagem da corrente elétrica alternada.

        Returns
        -------
        Z1 + Z2_prime : np.ndarray
            Somatorio da impedância do rotor e a do estator
        Z1 : np.ndarray
            Valor de impedância do estator

        Examples
        --------
        >>> s = 5
        >>> motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)
        >>> soma_impedancia, impedancia_estator = motor.calcular_impedancia
            (0.14098964096979372+0.8908940265114892j), (0.135+0.768j)
        """
        j = complex(0, 1)
        Z1 = self.R1 + j * self.X1
        Z2 = (self.R2 / s) + j * self.X2
        Zm = j * self.Xm
        Z2_prime = Z2 * Zm / (Z2 + Zm)
        return Z1 + Z2_prime, Z1

    def calcular_corrente(self, V_fase, s):
        """Cálculo da corrente de fase no sistema.
        
        Calcula a corrente de fase para um motor de indução de gaiola de esquilo, dado o valor da tensão de fase e do escorregamento.

        Returns
        -------
        I_fase : complex
            Corrente da fase

        Examples
        --------
        >>> V_fase, s = (10, 5)
        >>> motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)
        >>> corrente = motor.calcular_corrente(V_fase, s)
            (1.73297440237383-10.950425382691304j)
        """
        Z = self.calcular_impedancia(s)[0]
        I_fase = V_fase / Z
        return I_fase

    def calcular_tensao_induzida(self, V_fase, s):
        """Calcula a tensão induzida do motor.

        Calcula a tensão induzida (E2) no motor de indução de gaiola de esquilo, subtraindo a queda de tensão na impedância do estator da tensão de fase aplicada.

        Returns
        -------
        E2 : complex
            Valor total da tensão induzida.

        Examples
        --------
        >>> V_fase, s = (10, 5)
        >>> motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)
        >>> exemple_E2 = motor.calcular_tensao_induzida(V_fase, s)
            (1.3561217617726111+0.1473830856402245j)
        """
        E2 = V_fase - self.calcular_corrente(V_fase, s) * self.calcular_impedancia(s)[1]
        return E2

    def calcular_corrente_de_partida(self, V_fase, s):
        """Calculo da corrente de partida do motor.

        Calcula a corrente necessária para iniciar a operação do motor, considerando a impedância magnetizante.

        Returns
        -------
        I2 : complex
            Valor da corrente necessária para partida.

        Examples
        --------
        >>> V_fase, s = (10, 5)
        >>> motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)
        >>> corrente = motor.calcular_corrente_de_partida(V_fase, s)
            (1.7234443829657302-10.951461103602337j)
        """
        Im = self.calcular_tensao_induzida(V_fase, s) / self.Xm
        I2 = self.calcular_corrente(V_fase, s) - Im
        return I2

    def calcular_torque(self, V_fase, s):
        """Calcuclo do torque do motor

        Calcuclo do torque do motor elétrico aplicando a constante de proporcionalidade

        Returns
        -------
        self.K * torque : float 
            Torque com a constante de proporcionalidade 

        Examples
        --------
        >>> V_fase, s = (10, 5)
        >>> motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)
        >>> calculo_torque = motor.calcular_torque(V_fase, s)
            (0.011149713124255235)
        """
        I2 = self.calcular_corrente_de_partida(V_fase, s) # Corrente do rotor
        P_r = 3 * abs(I2)**2 * (self.R2 / s)# Potência no rotor e torque
        torque = P_r / self.w_s
        return self.K * torque 

    def encontrar_maior_torque(self, V_fase, escorregamentos):
        """Encontrar o maior torque gerado.
        
        Calcula o torque para diferentes valores de escorregamento e retorna o escorregamento correspondente ao maior torque gerado.

        Returns
        -------
        max_s : float
            Escorregamento do torque máximo.
        max_torque : float
            Torque máximo.

        Examples
        --------
        >>> V_fase = 220
        >>> escorregamentos = [0.01, 0.05, 0.1, 0.15]
        >>> motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)
        >>> max_s_exemple, max_torque_exemple = motor.encontrar_maior_torque(V_fase, escorregamentos)
            (0.05, 325.6747000287233)
        """
        torques = [self.calcular_torque(V_fase, s) for s in escorregamentos]
        max_torque = max(torques)
        max_s = escorregamentos[torques.index(max_torque)]
        return max_s, max_torque

    def calcular_velocidade_angular(self, escorregamentos):
        """Calculo da valocidade angular do rotor.

        Função que calcula a velocidade de giro do rotor.

        Returns
        -------
        escorregamentos :  np.ndarray
            Retorna o valor da velocidade de giro do rotor 

        Examples
        --------
        >>> escorregamentos = [0.01, 0.05, 0.1, 0.15]
        >>> motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)
        >>> velocidade_angular_ex = motor.calcular_velocidade_angular(escorregamentos)
            [1782. 1710. 1620.]
        """
        return (self.w_s * (1 - escorregamentos)) * (30 / np.pi)

    def calcular_fator_potencia(self, s):
        """Calcula o fator de potencia.

        calcula o fator de potência e no retorno deve-se subtrair o fator de potencia por 1, para se ter o valor esperado 

        Returns
        -------
        fator_potencia : float
            Fator de potencia menos 1 

        Examples
        --------
        >>> s = (5) 
        >>> motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)
        >>> calculo_fator_potencia = motor.calcular_fator_potencia(s)
            0.8436889515099687
        """
        Z, _ = self.calcular_impedancia(s)
        fator_potencia = np.cos(np.angle(Z))
        return 1-fator_potencia

    def simular_motor(self, t_final, num_steps):
        
        t = np.linspace(0, t_final, num_steps)
        s = np.linspace(0.0001, 1, num_steps)  # Slip varies from 0 to 1 during simulation

        corrente = [self.calcular_corrente(220, si) for si in s]
        torque = [self.calcular_torque(220, si) for si in s]
        fator_potencia = [self.calcular_fator_potencia(si) for si in s]

        return t, torque, corrente, fator_potencia

def plotar_motor(t, torque, corrente, fator_potencia):
    plt.figure(figsize=(10, 10))

    plt.subplot(3, 1, 1)
    plt.plot(t, torque, label='Torque')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Torque (Nm)')
    plt.title('Torque em função do Tempo')
    plt.grid(True)
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(t, corrente, label='Corrente')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Corrente (A)')
    plt.title('Corrente em função do Tempo')
    plt.grid(True)
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(t, fator_potencia, label='Fator de Potência')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Fator de Potência')
    plt.title('Fator de Potência em função do Tempo')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()

def plotar_motor_dinamico(t, torque, corrente, fator_potencia):
    plt.figure(figsize=(10, 10))
    
    for i in range(len(t)):
        plt.subplot(3, 1, 1)
        plt.plot(t[:i], torque[:i], label='Torque')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Torque (Nm)')
        plt.title('Torque em função do Tempo')
        plt.grid(True)
        plt.legend()

        plt.subplot(3, 1, 2)
        plt.plot(t[:i], corrente[:i], label='Corrente')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Corrente (A)')
        plt.title('Corrente em função do Tempo')
        plt.grid(True)
        plt.legend()

        plt.subplot(3, 1, 3)
        plt.plot(t[:i], fator_potencia[:i], label='Fator de Potência')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Fator de Potência')
        plt.title('Fator de Potência em função do Tempo')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.draw()
        plt.pause(0.01)
        plt.clf()

def example_motor():

    motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)

    t_final = 3
    num_steps = 1000
    t, torque, corrente, fator_potencia = motor.simular_motor(t_final, num_steps)

    plotar_motor(t, torque, corrente, fator_potencia)
    plotar_motor_dinamico(t, torque, corrente, fator_potencia)


    
class controle_fluxo:
    def __init__(self, k, amp_fluxo, r_s, v_s, i_s, r_r, v_r, i_r, ref_torque, ref_fluxo_rotor, ref_fluxo_estator, t, w1, wr):
        self.k = k
        self.amp_fluxo = amp_fluxo
        self.r_s = r_s
        self.v_s = v_s
        self.i_s = i_s
        self.r_r = r_r
        self.v_r = v_r
        self.i_r = i_r
        self.ref_torque = ref_torque
        self.ref_fluxo_rotor = ref_fluxo_rotor
        self.ref_fluxo_estator = ref_fluxo_estator
        self.t = t
        self.w1 = w1
        self.wr = wr
        self.dt = t[1] - t[0]
        self.psi_r = np.zeros_like(t)
        self.psi_s = np.zeros_like(t)
        self.erro_torque = np.zeros_like(t)
        self.erro_fluxo_rotor = np.zeros_like(t)
        self.erro_fluxo_estator = np.zeros_like(t)

    def calcular_fluxo_rotorico(self):
        """Calcula o fluxo do rotor por integração da tensão do rotor.

        Utiliza a tensão do rotor e a corrente do rotor ajustada pela resistência do rotor para calcular
        o fluxo do rotor ao longo do tempo.

        Returns
        -------
        psi_r : array_like
            Fluxo do rotor calculado ao longo do tempo.
        """
        self.psi_r = np.cumsum((self.v_r - self.i_r * self.r_r) * self.dt)
        return self.psi_r

    def calcular_fluxo_estatorico(self):
        """Calcula o fluxo do estator por integração da tensão do estator.

        Utiliza a tensão do estator e a corrente do estator ajustada pela resistência do estator para calcular
        o fluxo do estator ao longo do tempo.

        Returns
        -------
        psi_s : array_like
            Fluxo do estator calculado ao longo do tempo.
        """
        self.psi_s = np.cumsum((self.v_s - self.i_s * self.r_s) * self.dt)
        return self.psi_s

    def calcular_torque(self):
        """Calcula o torque eletromagnético com base no fluxo do rotor e nas velocidades.

        Utiliza o fluxo do rotor e as velocidades do estator e do rotor para calcular o torque
        eletromagnético ao longo do tempo.

        Returns
        -------
        ce : array_like
            Torque eletromagnético calculado ao longo do tempo.
        """
        self.ce = self.k * (self.psi_r ** 2) * (self.w1 - self.wr)
        return self.ce

    def controle(self):
        """Realiza o controle ajustando o fluxo e o torque de acordo com os valores de referência.

        Calcula o fluxo do rotor, o fluxo do estator e o torque, e ajusta as tensões baseadas nos erros
        entre os valores calculados e os valores de referência.

        Returns
        -------
        v_s_corrigido_rotor : array_like
            Tensão corrigida do rotor.
        v_s_corrigido_estator : array_like
            Tensão corrigida do estator.
        erro_torque : array_like
            Erro de torque calculado ao longo do tempo.
        erro_fluxo_rotor : array_like
            Erro de fluxo do rotor calculado ao longo do tempo.
        erro_fluxo_estator : array_like
            Erro de fluxo do estator calculado ao longo do tempo.
        """
        self.calcular_fluxo_rotorico()
        self.calcular_fluxo_estatorico()
        self.calcular_torque()
        self.erro_torque = self.ref_torque - self.ce
        self.erro_fluxo_rotor = self.ref_fluxo_rotor - self.psi_r
        self.erro_fluxo_estator = self.ref_fluxo_estator - self.psi_s
        v_s_corrigido_rotor = self.v_r + self.erro_torque * 0.1 + self.erro_fluxo_rotor * 0.1
        v_s_corrigido_estator = self.v_s + self.erro_torque * 0.1 + self.erro_fluxo_estator * 0.1
        return v_s_corrigido_rotor, v_s_corrigido_estator, self.erro_torque, self.erro_fluxo_rotor, self.erro_fluxo_estator

def customPlot(subplot, x, y, label, xlabel, ylabel):
    plt.subplot(5, 1, subplot)
    plt.plot(x, y, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)

def generatePlots(t, controle, erro_torque, erro_fluxo_rotor, saveFig):
    plt.figure(figsize=(10, 20))  # Aumentar o tamanho da figura

    customPlot(1, t, erro_torque, 'Erro de Torque', 'Tempo (s)', 'Erro de Torque') # Gráfico do erro do torque
    customPlot(2, t, erro_fluxo_rotor, 'Erro de Fluxo Rotorico', 'Tempo (s)', 'Erro de Fluxo Rotorico')  # Gráfico do erro do fluxo rotorico
    customPlot(3, t, controle.ce, 'Conjugado Eletromagnético', 'Tempo (s)', 'Torque') # Gráfico do conjugado eletromagnético   
    customPlot(4, t, controle.psi_r, 'Fluxo do Rotor', 'Tempo (s)', 'Fluxo do Rotor') # Gráfico do fluxo do rotor
    customPlot(5, t, controle.psi_s, 'Fluxo do Estator', 'Tempo (s)', 'Fluxo do Estator') # Gráfico do fluxo do estator

    plt.tight_layout(pad=3.0)  # Aumentar o espaço entre os subplots

    if saveFig:
        plt.savefig('controle_fluxo.png')
    else:
        plt.show()

def example_controle_fluxo(saveToFile=False):
    # Simulação
    t = np.linspace(0, 1, 1000)
    v_s = np.sin(2 * np.pi * 50 * t) * np.exp(t)  # Tensão do estator (exemplo)
    i_s = np.sin(2 * np.pi * 50 * t) * np.exp(t)  # Corrente do estator (exemplo)
    v_r = np.sin(2 * np.pi * 50 * t) * np.exp(t)  # Tensão do rotor (exemplo)
    i_r = np.sin(2 * np.pi * 50 * t) * np.exp(t)  # Corrente do rotor (exemplo)
    ref_torque = 1  # Referência de torque
    ref_fluxo_rotor = 1  # Referência de fluxo do rotor
    ref_fluxo_estator = 1  # Referência de fluxo do estator

    # Instanciar o objeto controle_fluxo com os parâmetros necessários
    controle = controle_fluxo(k=1, amp_fluxo=1, r_s=0.1, v_s=v_s, i_s=i_s, r_r=0.1, v_r=v_r, i_r=i_r, ref_torque=ref_torque, ref_fluxo_rotor=ref_fluxo_rotor, ref_fluxo_estator=ref_fluxo_estator, t=t, w1=50, wr=40)

    # Controle
    v_s_corrigido_rotor, v_s_corrigido_estator, erro_torque, erro_fluxo_rotor, erro_fluxo_estator = controle.controle()
    generatePlots(t, controle, erro_torque, erro_fluxo_rotor, saveFig=saveToFile)



class controle_inversor:
    def __init__(self, q1, q2, q3, v_dc, t_s, V, DDP_da_fonte):
        self.q1 = q1  # Interruptor Q1
        self.q2 = q2  # Interruptor Q2
        self.q3 = q3  # Interruptor Q3
        self.v_dc = v_dc
        self.t_s = t_s
        self.V = V  # Tensão de Alimentação
        self.DDP_da_fonte = DDP_da_fonte  # Diferença de Potencial da Fonte

        self.f = 50  # frequency in Hz
        self.T = 1 / self.f  # period
        self.t = np.linspace(0, 2 * self.T, 1000)  # time vector
        self.f_start = 1  # initial frequency in Hz
        self.f_end = 1000  # final frequency in Hz
        self.w_start = 2 * np.pi * self.f_start  # initial frequency in Rad
        self.w_end = 2 * np.pi * self.f_end  # final frequency in Rad
        self.num_points = 400  # Number of frequency points for smooth plot

    def chaves(self):
        """Determinar se as chaves de comutação.

        Função para determinar se as chaves de comutação estão abertas (False) ou fechadas (True). 
        Referido no livro Sistemas de Acionamento Estatico de Maquina Eletrica, Cursino Brandão Jacobina, Pag 89.

        Returns
        -------
        chave_1 : bool
            Retorna True se a chave 1 estiver fechada, False se estiver aberta.
        chave_2 : bool
            Retorna True se a chave 2 estiver fechada, False se estiver aberta.
        chave_3 : bool
            Retorna True se a chave 3 estiver fechada, False se estiver aberta.
    
        Examples
        --------
        >>> self.q1 = 1
        >>> self.q2 = 0
        >>> self.q3 = 1
        >>> chaves()
            (True, False, True)
        """
        chave_1 = self.q1 == 1
        chave_2 = self.q2 == 1
        chave_3 = self.q3 == 1
        return chave_3, chave_2, chave_1

    def tensao_nos_terminais(self):
        """Calcula as tensões trifásicas nos terminais de cargas.

        A função retorna as tensões nos terminais de uma carga trifásica com base no estado das chaves de comutação 
        e na diferença de potencial fornecida pela fonte, conforme descrito nas equações (7.1), (7.2) e (7.3) do livro 
        "Sistemas de Acionamento Estático de Máquina Elétrica" de Cursino Brandão Jacobina.

        Returns
        -------
        terminal_1 : float
            Tensão no terminal 1.
        terminal_2 : float
            Tensão no terminal 2.
        terminal_3 : float
            Tensão no terminal 3.

        Examples
        --------
        >>> self.q1 = 1
        >>> self.q2 = 0
        >>> self.q3 = 1
        >>> self.V = 220
        >>> self.DDP_da_fonte = 10
        >>> tensao_nos_terminais()
            (120.0, -100.0, 120.0)
        """
        terminal_1 = (2 * self.q1 - 1) * self.V / 2 + self.DDP_da_fonte  # Jacobina equação (7.1)
        terminal_2 = (2 * self.q2 - 1) * self.V / 2 + self.DDP_da_fonte  # Jacobina equação (7.2)
        terminal_3 = (2 * self.q3 - 1) * self.V / 2 + self.DDP_da_fonte  # Jacobina equação (7.3)
        return terminal_1, terminal_2, terminal_3

    def componente_direta(self, terminal_1, terminal_2, terminal_3):
        """Retorna o valor da componente direta.

        Calcula a componente direta das tensões trifásicas nos terminais usando a fórmula descrita na equação (7.10) 
        do livro "Sistemas de Acionamento Estático de Máquina Elétrica" de Cursino Brandão Jacobina.

        Returns
        -------
        componente_direta : float
            Valor da componente direta das tensões nos terminais.

        Examples
        --------
        >>> terminal_1 = 120.0
        >>> terminal_2 = -100.0
        >>> terminal_3 = 120.0
        >>> componente_direta(terminal_1, terminal_2, terminal_3)
            146.097563556976
        """
        componente_direta = np.sqrt(2 / 3) * (terminal_1 - (terminal_2 / 2) - (terminal_3 / 2))  # Jacobina equação (7.10)
        return componente_direta

    def componente_quadratura(self, terminal_2, terminal_3):
        """Retorna o valor da componente de quadratura.

        Calcula a componente de quadratura das tensões trifásicas nos terminais usando a fórmula descrita na equação (7.11) 
        do livro "Sistemas de Acionamento Estático de Máquina Elétrica" de Cursino Brandão Jacobina.

        Returns
        -------
        componente_quadratura : float
            Valor da componente de quadratura das tensões nos terminais.

        Examples
        --------
        >>> terminal_2 = -100.0
        >>> terminal_3 = 120.0
        >>> componente_quadratura(terminal_2, terminal_3)
            -170.78251276599332
        """
        componente_quadratura = np.sqrt(2 / 3) * (np.sqrt(3) / 2 * terminal_2 - (np.sqrt(3) / 2) * terminal_3)  # Jacobina equação (7.11)
        return componente_quadratura

    def vetores(self):
        """Retorna os vetores de tensão não nulos para a modulação vetorial.

        Calcula os vetores de tensão que são utilizados na modulação vetorial, conforme descrito nas equações (7.14) a (7.19) 
        do livro "Sistemas de Acionamento Estático de Máquina Elétrica" de Cursino Brandão Jacobina.

        Returns
        -------
        vetor_1 : complex
            Vetor de tensão 1.
        vetor_2 : complex
            Vetor de tensão 2.
        vetor_3 : complex
            Vetor de tensão 3.
        vetor_4 : complex
            Vetor de tensão 4.
        vetor_5 : complex
            Vetor de tensão 5.
        vetor_6 : complex
            Vetor de tensão 6.

        Examples
        --------
        >>> self.V = 220
        >>> vetores()
            (179.61044104991788); 
            (89.80522052495894+155.56349186104046j); 
            (-89.80522052495894+155.56349186104046j);
            (-179.61044104991788);
            (89.80522052495894-155.56349186104046j);
            (-89.80522052495894-155.56349186104046j);
        """
        vetor_1 = np.sqrt(2 / 3) * self.V  # Jacobina equação (7.14)
        vetor_2 = (self.V / np.sqrt(6)) + (1j * self.V / np.sqrt(2))  # Jacobina equação (7.15)
        vetor_3 = (-self.V / np.sqrt(6)) + (1j * self.V / np.sqrt(2))  # Jacobina equação (7.16)
        vetor_4 = - np.sqrt(2 / 3) * self.V  # Jacobina equação (7.17)
        vetor_5 = (self.V / np.sqrt(6)) - (1j * self.V / np.sqrt(2))  # Jacobina equação (7.18)
        vetor_6 = (-self.V / np.sqrt(6)) - (1j * self.V / np.sqrt(2))  # Jacobina equação (7.19)
        
        return vetor_1, vetor_2, vetor_3, vetor_4, vetor_5, vetor_6

    def transform_clarke(self, v_a, v_b, v_c):
        """Transformação de Clarke.

        Converte as tensões trifásicas (v_a, v_b, v_c) para as componentes alfa (v_alpha) e beta (v_beta) 
        usando a Transformação de Clarke.
        
        Returns
        -------
        v_alpha : float
            Componente alfa da transformação de Clarke.
        v_beta : float
            Componente beta da transformação de Clarke.

        Examples
        --------
        >>> v_a = 220
        >>> v_b = 110
        >>> v_c = -110
        >>> transform_clarke(v_a, v_b, v_c)
            (146.66666666666666, 190.5255888325765)
        """
        v_alpha = (2 / 3) * (v_a - 0.5 * (v_b + v_c))
        v_beta = (2 / 3) * (np.sqrt(3) / 2 * (v_b - v_c))
        return v_alpha, v_beta

    def transform_park(self, v_alpha, v_beta):
        """Transformação de Park.

        Converte as componentes alfa (v_alpha) e beta (v_beta) para as componentes direta (D) e de quadratura (Q) 
        usando a Transformação de Park.

        Returns
        -------
        D : float
            Componente direta da transformação de Park.
        Q : float
            Componente de quadratura da transformação de Park.

        Examples
        --------
        >>> v_alpha = 146.66666666666666
        >>> v_beta = 190.5255888325765
        >>> self.f = 50
        >>> self.t = 0.01
        >>> transform_park(v_alpha, v_beta)
            (140.3174570760944, 196.66287296729442)
        """
        theta = 2 * np.pi * self.f * self.t
        D = v_alpha * np.cos(theta) + v_beta * np.sin(theta)
        Q = -v_alpha * np.sin(theta) + v_beta * np.cos(theta)
        return D, Q

    def find_sector(self, v_alpha, v_beta):
        """Determina em qual setor o vetor de referência está localizado.

        Returns
        -------
        sector : int
            O número do setor onde o vetor de referência está localizado.
            Os setores são definidos como:
            1: 0° ≤ setor < 60°
            2: 60° ≤ setor < 120°
            3: 120° ≤ setor < 180°
            4: 180° ≤ setor < 240°
            5: 240° ≤ setor < 300°
            6: 300° ≤ setor < 360°

        Examples
        --------
        >>> find_sector(1, 1)
            1
        >>> find_sector(-1, 1)
            2
        >>> find_sector(-1, -1)
            4
        >>> find_sector(1, -1)
            6
        """
        sector = np.arctan2(v_beta, v_alpha) * 180 / np.pi
        if sector < 0:
            sector += 360
        if 0 <= sector < 60:
            return 1
        elif 60 <= sector < 120:
            return 2
        elif 120 <= sector < 180:
            return 3
        elif 180 <= sector < 240:
            return 4
        elif 240 <= sector < 300:
            return 5
        elif 300 <= sector < 360:
            return 6

    def calculate_times(self, v_alpha, v_beta):
        """Calcula os tempos de comutação T1, T2 e T0.

        A função calcula os tempos de comutação T1, T2 e T0 com base nas componentes alfa (v_alpha) e beta (v_beta) 
        da tensão, no tempo de comutação total (t_s) e na tensão do barramento (v_dc).

        Returns
        -------
        t1 : float
            Tempo de comutação T1.
        t2 : float
            Tempo de comutação T2.
        t0 : float
            Tempo de comutação T0.

        Examples
        --------
        >>> self.t_s = 0.001
        >>> self.v_dc = 400
        >>> v_alpha = 100
        >>> v_beta = 50
        >>> calculate_times(v_alpha, v_beta)
            (0.000125, 0.00025, 0.000625)
        """
        t1 = v_beta * self.t_s / self.v_dc
        t2 = v_alpha * self.t_s / self.v_dc
        t0 = self.t_s - t1 - t2
        return t1, t2, t0

    def svm(self):
        """Modulação por vetor espacial (SVM).

        Executa a modulação por vetor espacial (SVM) para tensões trifásicas de referência.

        O procedimento inclui a transformação de Clarke das tensões trifásicas, a identificação do setor de referência, 
        e o cálculo dos tempos de comutação T1, T2 e T0.

        Returns
        -------
        v_alpha : float
            Componente alfa da transformação de Clarke.
        v_beta : float
            Componente beta da transformação de Clarke.
        sector : int
            Setor onde o vetor de referência está localizado.
        t1 : float
            Tempo de comutação T1.
        t2 : float
            Tempo de comutação T2.
        t0 : float
            Tempo de comutação T0.

        Examples
        --------
        >>> self.v_a = 220
        >>> self.v_b = 110
        >>> self.v_c = -110
        >>> self.t_s = 0.001
        >>> self.v_dc = 400
        >>> svm()
            (146.66666666666666, 190.5255888325765, 1, 0.0004763147220814413, 0.0003666666666666667, 0.00015601861125189196)
        """
        v_alpha, v_beta = self.transform_clarke(self.v_a, self.v_b, self.v_c)
        sector = self.find_sector(v_alpha, v_beta)
        t1, t2, t0 = self.calculate_times(v_alpha, v_beta)
        return v_alpha, v_beta, sector, t1, t2, t0

    def first_order_system(self, omega_n, zeta):
        """Define uma função de transferência de primeira ordem.

        Esta função cria uma função de transferência para um sistema de primeira ordem com frequência natural `omega_n` 
        e fator de amortecimento `zeta`. A função de transferência resultante é dada por: H(s) = omega_n^2 / (s^2 + 2 * zeta * omega_n * s + omega_n^2)

        Parameters
        ----------
        omega_n : float
            Frequência natural do sistema.
        zeta : float
            Fator de amortecimento do sistema.

        Returns
        -------
        transfer_function : scipy.signal.TransferFunction
            A função de transferência de primeira ordem definida pelos parâmetros fornecidos.
            
        """
        num = [omega_n**2]
        den = [1, 2 * zeta * omega_n, omega_n**2]
        return signal.TransferFunction(num, den)

    def second_order_system(self, omega_n, zeta):
        """Define uma função de transferência de segunda ordem.

        Esta função cria um sistema de segunda ordem com a frequência natural `omega_n` e o fator de amortecimento `zeta`.

        Returns
        -------
        system : signal.TransferFunction
            Função de transferência do sistema de segunda ordem.

        Examples
        --------
        >>> omega_n = 5.0
        >>> zeta = 0.7
        >>> system = second_order_system(omega_n, zeta)

        """
        num = [omega_n**2]
        den = [1, 2 * zeta * omega_n, omega_n**2]
        return signal.TransferFunction(num, den)

    def frequency_sweep(self):
        """Gera uma varredura de frequência logarítmica.

        Esta função cria uma série de pontos de frequência em uma escala logarítmica, começando de `w_start` até `w_end` com `num_points` pontos.

        Returns
        -------
        w : np.ndarray
            Array de frequências em uma escala logarítmica.

        Examples
        --------
        >>> sweeper = YourClassName(w_start=1, w_end=1000, num_points=100)
        >>> frequencies = sweeper.frequency_sweep()
            [ 1.   1.04712855   1.0964782   ...   953.03136355    1000. ]
        """
        self.w = np.logspace(np.log10(self.w_start), np.log10(self.w_end), self.num_points)
        return self.w

    def plot_bode(self):
        """
        Plota o diagrama de Bode para sistemas de primeira e segunda ordem.
        """
        # Define sistemas de funções de transferência (exemplo)
        system1 = self.first_order_system(10, 0.5)  # Exemplo de sistema de primeira ordem
        system2 = self.second_order_system(20, 0.7)  # Exemplo de sistema de segunda ordem

        # Gera dados de Bode plot para cada sistema
        w = self.frequency_sweep()
        w, mag1, phase1 = signal.bode(system1, w)
        w, mag2, phase2 = signal.bode(system2, w)

        plt.figure(figsize=(10, 6))

        plt.subplot(2, 1, 1)
        plt.semilogx(w, mag1, label='System 1')
        plt.semilogx(w, mag2, label='System 2')
        plt.title('Bode Plot (Magnitude)')
        plt.ylabel('Magnitude (dB)')
        plt.grid(True, which="both", ls="--")
        plt.legend()

        plt.subplot(2, 1, 2)
        plt.semilogx(w, phase1, label='System 1')
        plt.semilogx(w, phase2, label='System 2')
        plt.xlabel('Frequency (rad/s)')
        plt.ylabel('Phase (degrees)')
        plt.grid(True, which="both", ls="--")
        plt.legend()

        plt.tight_layout()
        plt.show()

    def plot_tensoes(self):
        """Plota as tensões de fase e as componentes Alpha-Beta e D-Q.
        """
        plt.figure(figsize=(10, 8))

        plt.subplot(3, 1, 1)
        plt.plot(self.t, self.v_a, label='Va')
        plt.plot(self.t, self.v_b, label='Vb')
        plt.plot(self.t, self.v_c, label='Vc')
        plt.title('Tensões de Fase')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensão (V)')
        plt.legend()

        plt.subplot(3, 1, 2)
        plt.plot(self.t, self.v_alpha, label='Alpha')
        plt.plot(self.t, self.v_beta, label='Beta')
        plt.title('Componentes Alpha-Beta')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensão (V)')
        plt.legend()

        plt.subplot(3, 1, 3)
        plt.plot(self.t, self.D, label='D')
        plt.plot(self.t, self.Q, label='Q')
        plt.title('Componentes D-Q')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensão (V)')
        plt.legend()

        plt.tight_layout()
        plt.show()

    def plot_signals(self, T_pwm):
        # Calculando a frequência de chaveamento
        f_chaveamento = 1 / T_pwm  # Frequência de chaveamento (Hz)

        # Tempo de simulação
        t_final = 0.1  # Tempo final da simulação (segundos)
        t = np.arange(0, t_final, T_pwm)  # Vetor de tempo

        # Geração da onda senoidal de referência
        V_out_peak = self.v_dc / 2  # Amplitude de pico da tensão de saída CA
        omega = 2 * np.pi * self.f  # Frequência angular
        V_ref = V_out_peak * np.sin(omega * t)  # Onda senoidal de referência

        # Geração do sinal PWM
        V_carrier = V_out_peak * np.sign(np.sin(2 * np.pi * f_chaveamento * t))

        # Sinal PWM
        PWM_signal = np.where(V_ref >= V_carrier, self.v_dc, 0)

        # Plotando os sinais
        self.ax1.clear()
        self.ax1.plot(t, V_ref, label='Onda Senoidal de Referência')
        self.ax1.set_title('Onda Senoidal de Referência')
        self.ax1.set_xlabel('Tempo (s)')
        self.ax1.set_ylabel('Tensão (V)')
        self.ax1.legend()

        self.ax2.clear()
        self.ax2.plot(t, V_carrier, label='Onda Portadora')
        self.ax2.set_title('Onda Portadora')
        self.ax2.set_xlabel('Tempo (s)')
        self.ax2.set_ylabel('Tensão (V)')
        self.ax2.legend()

        self.ax3.clear()
        self.ax3.plot(t, PWM_signal, label='Sinal PWM')
        self.ax3.set_title('Sinal PWM')
        self.ax3.set_xlabel('Tempo (s)')
        self.ax3.set_ylabel('Tensão (V)')
        self.ax3.legend()

        self.fig.canvas.draw_idle()

    def configure_slider_plot(self, T_pwm_initial):
        # Configuração inicial do gráfico
        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, 1, figsize=(12, 8))
        plt.subplots_adjust(left=0.1, bottom=0.25)

        # Plot inicial
        self.plot_signals(T_pwm_initial)

        # Eixo para o slider de controle de T_pwm
        ax_slider = plt.axes([0.1, 0.1, 0.8, 0.03], facecolor='lightgoldenrodyellow')
        T_pwm_slider = Slider(ax_slider, 'T_pwm (s)', 0.00001, 0.1, valinit=T_pwm_initial)

        # Função de atualização do slider
        def update(val):
            T_pwm = T_pwm_slider.val
            self.plot_signals(T_pwm)

        T_pwm_slider.on_changed(update)

        plt.show()

def example_controle_inversor():
    v_a = 220 * np.sin(2 * np.pi * 50 * np.linspace(0, 2 * (1/50), 1000))
    v_b = 220 * np.sin(2 * np.pi * 50 * np.linspace(0, 2 * (1/50), 1000) - 2*np.pi/3)
    v_c = 220 * np.sin(2 * np.pi * 50 * np.linspace(0, 2 * (1/50), 1000) + 2*np.pi/3)

    controle = controle_inversor(q1=1, q2=0, q3=1, v_dc=400, t_s=0.001, V=220, DDP_da_fonte=0)
    controle.v_a = v_a
    controle.v_b = v_b
    controle.v_c = v_c

    # Calcular transformações
    controle.v_alpha, controle.v_beta = controle.transform_clarke(controle.v_a, controle.v_b, controle.v_c)
    controle.D, controle.Q = controle.transform_park(controle.v_alpha, controle.v_beta)

    # Plotar gráficos
    controle.plot_bode()
    controle.plot_tensoes()
    controle.configure_slider_plot(0.0001)



class peso:
    def __init__(self, peso_bateria, peso_inversor, peso_motor):
        self.peso_bateria = peso_bateria
        self.peso_inversor = peso_inversor
        self.peso_motor = peso_motor

    def peso_total(self):
        """Peso total dos componentes do sistema de powertrain

        Calcula o peso total dos componentes do sistema de powertrain, somando o peso da bateria, do inversor e do motor.

        Returns
        -------
        peso_total : float
            Somatorio dos pesos dos componentes do sistema

        Examples
        --------
        >>> peso_pwt = peso(10, 10, 65)
            85
        """
        peso_total = self.peso_bateria + self.peso_inversor + self.peso_motor
        return peso_total

def example_peso():

    peso_pwt = peso(10, 10, 65)
    total = peso_pwt.peso_total()
    
    print(f"O peso total é {total} Kg")



# testes
# example_bateria()
# example_inversor()
# plotar_potencia_vs_frequencia()
# example_motor()
# example_controle_fluxo(saveToFile=False)
# example_controle_inversor()
# example_peso()
