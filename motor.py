import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import cmath  # Para operações com números complexos

__all__ = ["bateria", "motor_gaiola", "inversor", "controle_fluxo", "controle_inversor"]

class bateria:
    def __init__(self, bat_type, bat_Idis, bat_Ich, bat_t, bat_Es, bat_K, bat_Q, bat_A, bat_B, bat_Iast, bat_R, peso_bat=None,time=None):
     
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
        self.peso_bat = peso_bat # peso da bateria (kg)

        if time is None:
            self.time = self.bat_Q / self.bat_Idis  # A*h/A = h
        else:
            self.time = time

    def battery(self):

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

    plt.subplot(4, 1, 1)
    plt.plot(t, v_a, label='Tensão Fase A')
    plt.plot(t, v_b, label='Tensão Fase B')
    plt.plot(t, v_c, label='Tensão Fase C')
    plt.title('Parâmetros de Saída do Inversor')
    plt.ylabel('Tensão (V)')
    plt.legend()

    plt.subplot(4, 1, 2)
    plt.plot(t, i_a, label='Corrente Fase A')
    plt.plot(t, i_b, label='Corrente Fase B')
    plt.plot(t, i_c, label='Corrente Fase C')
    plt.ylabel('Corrente (A)')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.plot(t, np.full_like(t, Pa), label='Potência Ativa (W)')
    plt.plot(t, np.full_like(t, Q), label='Potência Reativa (VAR)')
    plt.plot(t, np.full_like(t, S), label='Potência Aparente (VA)')
    plt.ylabel('Potência')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.plot(t, lambda_m, label='Fluxo de Magnetização (lambda)')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Fluxo de Magnetização (Wb)')
    plt.legend()

    plt.grid(True)
    plt.show()

def example_inversor():

    t_sim = 0.2
    t = np.linspace(0, t_sim, 1000)

    my_inversor = inversor(220, 60, 220, 118, 0, 0, 10000, 0.01)

    v_a, v_b, v_c = my_inversor.gerar_tensoes_saida(t)
    i_a, i_b, i_c = my_inversor.gerar_correntes_saida(t)
    Pa, Q, S = my_inversor.calcular_potencias(v_a, v_b, v_c, i_a, i_b, i_c)

    plotar_inversor(t, v_a, v_b, v_c, i_a, i_b, i_c, Pa, Q, S, np.zeros_like(t))  # Passamos np.zeros_like(t) para lambda_m como exemplo


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
        j = complex(0, 1)
        Z1 = self.R1 + j * self.X1
        Z2 = (self.R2 / s) + j * self.X2
        Zm = j * self.Xm
        Z2_prime = Z2 * Zm / (Z2 + Zm)
        return Z1 + Z2_prime, Z1

    def calcular_corrente(self, V_fase, s):
        Z = self.calcular_impedancia(s)[0]
        I_fase = V_fase / Z
        return I_fase

    def calcular_tensao_induzida(self, V_fase, s):
        E2 = V_fase - self.calcular_corrente(V_fase, s) * self.calcular_impedancia(s)[1]
        return E2

    def calcular_corrente_de_partida(self, V_fase, s):
        Im = self.calcular_tensao_induzida(V_fase, s) / self.Xm
        I2 = self.calcular_corrente(V_fase, s) - Im
        return I2

    def calcular_torque(self, V_fase, s):
        # Corrente do rotor
        I2 = self.calcular_corrente_de_partida(V_fase, s)
        # Potência no rotor e torque
        P_r = 3 * abs(I2)**2 * (self.R2 / s)
        torque = P_r / self.w_s

        return self.K * torque  # Aplica a constante de proporcionalidade

    def encontrar_maior_torque(self, V_fase, escorregamentos):
        torques = [self.calcular_torque(V_fase, s) for s in escorregamentos]
        max_torque = max(torques)
        max_s = escorregamentos[torques.index(max_torque)]
        return max_s, max_torque

    def calcular_velocidade_angular(self, escorregamentos):
        return (self.w_s * (1 - escorregamentos)) * (30 / np.pi)

    def calcular_fator_potencia(self, s):
        """Calculates power factor based on slip."""
        Z, _ = self.calcular_impedancia(s)
        fator_potencia = np.cos(np.angle(Z))
        return 1-fator_potencia

    def simular_motor(self, t_final, num_steps):
        """Simulates motor behavior with a linear slip."""
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

def example_motor():

    motor = motor_gaiola(60, 4, 0.135, 0.768, 0.03, 0.123, 142.3, 0.95, 220)

    t_final = 10
    num_steps = 1000
    t, torque, corrente, fator_potencia = motor.simular_motor(t_final, num_steps)

    plotar_motor(t, torque, corrente, fator_potencia)


    
class controle_fluxo:
    def __init__(self,):
        return
        
class controle_inversor:
    def __init__(self,):
        return

def peso_total():

    battery = bateria(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0)
    inverter = inversor(0, 0, 0, 0, 0, 0, 0, 0, 50)
    motor = motor_gaiola(0, 0, 0, 0, 20)
    total_weight = battery.peso_bat + inverter.peso_inv + motor.peso_mot

    print(f"O peso total é {total_weight} Kg")  # Output: O peso total é 120 Kg



# testes
#example_bateria()
example_inversor()
example_motor()
# peso_total()

# mapamental

# refatoração
# [  ] completar as descriçoes: gabriel, samuel e alan
# [  ] comcertar o fluxo: igor, nikolas e caio
# [ ] juntar os 2 codigos de motor: fabricio e leo
# [  ] concertar o peso: muca

# concluir os codigos 
# [  ] fluxo: caio, alan,
# [  ] pwm: mari, jm e iago
