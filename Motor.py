import numpy as np
import matplotlib.pyplot as plt

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
        s = np.linspace(0, 1, num_steps)  # Slip varies from 0 to 1 during simulation

        corrente = [self.calcular_corrente_de_partida(220, si) for si in s]
        torque = [self.calcular_torque(220, si) for si in s]
        fator_potencia = [self.calcular_fator_potencia(si) for si in s]

        return t, torque, corrente, fator_potencia

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
    motor = motor_gaiola(60, 4, 0.135, 0.768, 0.116, 0.123, 142.3, 0.95, 220)

    t_final = 3
    num_steps = 100
    t, torque, corrente, fator_potencia = motor.simular_motor(t_final, num_steps)

    plotar_motor_dinamico(t, torque, corrente, fator_potencia)

example_motor()
