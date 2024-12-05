import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control as clt

from class_motor import Motor

class Controle:
    def __init__(self, motor: Motor):
        """
        Inicializa a classe Controle com uma instância de Motor.
        
        Parameters
        ----------
        motor : Motor
            Instância da classe Motor.
        """
        self.motor = motor
        self.msr = motor.msr
        self.p = motor.p
        self.jm = motor.jm
        self.kf = motor.kf
        self.rs = motor.rs
        self.ls = motor.ls
        self.q1 = motor.q1
        self.q2 = motor.q2
        self.q3 = motor.q3
        self.valor_mu = motor.valor_mu
        self.Vs = motor.Vs
        self.f_onda_p = motor.f_onda_p

    def transfer_function(self):
        """Função de transferência

        Calcula a função de transferência do sistema, que é a relação entre a entrada e a saída em termos de numerador e denominador.

        Returns
        -------
        num : list
            Coeficientes do numerador da função de transferência.
        den : list
            Coeficientes do denominador da função de transferência.

        Examples
        --------
        >>> x = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> num, dem = x.transfer_function()
        >>> num, dem
        ([0.182], [0.04, 0.02, 0.484])
        """
        
        # Calcula os coeficientes do numerador e denominador da função de transferência
        num = [self.msr* self.p]  # Numerador
        den = [self.jm, 2 * self.kf, self.rs + self.ls]  # Denominador
        
        return num, den

    def plot_bode(self):
        """Gera o diagrama de Bode do sistema.

        Esta função calcula a função de transferência do motor e plota o diagrama de Bode, que inclui a magnitude e a fase em função da frequência.

        Examples
        --------
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> motor.plot_bode()
        """

        num, den = self.transfer_function()
        system = signal.TransferFunction(num, den)
        w, mag, phase = signal.bode(system)

        # Plota o diagrama de Bode
        plt.figure(figsize=(10, 6))

        # Subplot para a magnitude
        plt.subplot(2, 1, 1)
        plt.semilogx(w, mag)
        plt.title('Diagrama de Bode - Motor')
        plt.ylabel('Magnitude (dB)')
        plt.grid(which="both", axis="both")

        # Subplot para a fase
        plt.subplot(2, 1, 2)
        plt.semilogx(w, phase)
        plt.xlabel('Frequência (rad/s)')
        plt.ylabel('Fase (graus)')
        plt.grid(which="both", axis="both")

        plt.tight_layout()
        plt.show()

    def plot_nyquist(self):
        """Gera o diagrama de Nyquist do sistema.

        Esta função calcula a função de transferência do motor e plota o diagrama de Nyquist, que representa a resposta em frequência do sistema no domínio complexo.

        Examples
        --------
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> motor.plot_nyquist()
        """

        num, den = self.transfer_function()
        motor_system = clt.TransferFunction(num, den)

        # Define as frequências para o Diagrama de Nyquist
        w_start = 1e-2
        w_stop = 1e3
        num_points = 1000
        frequencies = np.logspace(np.log10(w_start), np.log10(w_stop), num_points)

        # Plota o diagrama de Nyquist
        plt.figure()
        clt.nyquist_plot(motor_system, omega=frequencies)
        plt.title("Diagrama de Nyquist - Motor Trifásico")
        plt.grid(True)
        plt.show()

    def state_space_representation(self):
        """Representação do espaço de estado

        Calcula e retorna a representação no espaço de estados do motor elétrico.

        Returns:
        A: numpy.ndarray
            Matriz A do sistema no espaço de estado.
        B: numpy.ndarray
            Matriz B do sistema no espaço de estado.
        C: numpy.ndarray
            Matriz C do sistema no espaço de estado.
        D: numpy.ndarray
            Matriz D do sistema no espaço de estado.

        Examples:
        --------
        >>> representacao_espaco_estado = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> A, B, C, D = representacao_espaco_estado.state_space_representation()

        A = [[  0.    1. ]
            [-12.1  -0.5]]
        B = [[0.  ]
            [4.55]]
        C = [[1 0]]
        D = [[0]]
        """
    
        # Coeficientes da função de transferência
        num, den = self.transfer_function()

        # Sistema de segunda ordem: Numerador e denominador
        # Exemplo: num = [b0], den = [a2, a1, a0]
        a2 = den[0]
        a1 = den[1]
        a0 = den[2]
        b0 = num[0]

        # Matrizes A, B, C, D no espaço de estados
        A = np.array([[0, 1],
                      [-a0/a2, -a1/a2]])

        B = np.array([[0],
                      [b0/a2]])

        C = np.array([[1, 0]])

        D = np.array([[0]])

        return A, B, C, D

    def print_state_space(self):
        A, B, C, D = self.state_space_representation()
        print("Matriz A:")
        print(A)
        print("\nMatriz B:")
        print(B)
        print("\nMatriz C:")
        print(C)
        print("\nMatriz D:")
        print(D)

    def step_response(self):
        """Resposta ao Degrau

        Calcula e plota a resposta ao degrau unitário do sistema representado pela função de transferência do motor.

        Examples:
        --------
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> motor.step_response()

        """
        
        num, den = self.transfer_function()
        system = signal.TransferFunction(num, den)
        t, response = signal.step(system)

        # Plota a resposta degrau
        plt.figure(figsize=(10, 6))
        plt.plot(t, response, label='Resposta ao Degrau Unitário')
        plt.title('Resposta ao Degrau Unitário - Sistema de Segunda Ordem')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Amplitude')
        plt.grid(True)
        plt.legend()
        plt.show()
        
    def chaves(self):
        """ A configuração de chaves do inversor.


        Determina a configuração das seis chaves que 
        compõem o do inversor de frequência (q1, q2, q3, q4, q5, q6). 
        Sendo q4, q5 e q6 os complementares de q1, q2 e q3, respectivamente.
        
        Returns:
            Chave 1: Chave 1 Fechada (True) ou Aberta (False)
            Chave 2: Chave 2 Fechada (True) ou Aberta (False)
            Chave 3: Chave 3 Fechada (True) ou Aberta (False)
            Chave 4: Chave 4 Fechada (True) ou Aberta (False)
            Chave 5: Chave 5 Fechada (True) ou Aberta (False)
            Chave 6: Chave 6 Fechada (True) ou Aberta (False)
                   
        Example
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01, q1=1, q2=1, q3=0, valor_mu=1) 
        >>> motor.chaves()
       
        
        """
        
        if self.q1 == 1: chave_1 = True # Retorno chave 1 fechada
        else: chave_1 = False # Retorno chave 1 aberta

        if self.q2 == 1: chave_2 = True # Retorno chave 2 fechada
        else: chave_2 = False # Retorno chave 2 aberta

        if self.q3 == 1: chave_3 = True # Retorno chave 3 fechada
        else: chave_3 = False # Retorno chave 3 aberta

        q1_bar = 1 - self.q1
        q2_bar = 1 - self.q2
        q3_bar = 1 - self.q3
            
        if q1_bar == 1: chave_4 = True # Retorno chave 1_bar fechada
        else: chave_4 = False # Retorno chave 1_bar aberta

        if q2_bar == 1: chave_5 = True # Retorno chave 2_bar fechada
        else: chave_5 = False # Retorno chave 2_bar aberta

        if q3_bar == 1: chave_6 = True # Retorno chave 3_bar fechada
        else: chave_6 = False # Retorno chave 3_bar aberta    

        return print(f'A configuração de chaves do inversor é: C1={chave_1}, C2={chave_2}, C3={chave_3}, C4={chave_4}, C5={chave_5}, C6={chave_6}')
            
    def controle_pwm(self):
        """Sinal PWM para controle do 

        Implementa o controle PWM (Modulação por Largura de Pulso) para o sistema. A função realiza o cálculo das tensões moduladas, das correntes e do sinal PWM, além de gerar gráficos para análise visual dos resultados.

        Returns
        Plots do gráfico do sinal PWM, onda triangular portadora, correntes, tensão de entrada e tensão modulada
        None

        Examples
        --------
        >>> motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        >>> motor.controle_pwm()
        # O código exibirá gráficos de tensões moduladas, sinais PWM e correntes.

        Notes
        -----
        A função considera parâmetros do sistema previamente definidos, como tensão de alimentação, frequência da onda portadora, e fator de modulação.

        """

        # Define o passo de tempo do controlador
        self.t_pwm = np.linspace(0, 2, 1000)

        # Tensões de entrada
        self.v1 = self.Vs * np.sin(2 * np.pi * self.t_pwm)
        self.v2 = self.Vs * np.sin(2 * np.pi * self.t_pwm + 2 * np.pi / 3)
        self.v3 = self.Vs * np.sin(2 * np.pi * self.t_pwm + 4 * np.pi / 3)

        # Correntes
        self.i1 = self.Vs * np.sin(2 * np.pi * self.t_pwm - np.pi / 2)
        self.i2 = self.Vs * np.sin(2 * np.pi * self.t_pwm + 2 * np.pi / 3 - np.pi / 2)
        self.i3 = self.Vs * np.sin(2 * np.pi * self.t_pwm + 4 * np.pi / 3 - np.pi / 2)

        # Cálculo máximo e mínimo das tensões
        self.vN0max_star = (self.Vs / 2) - np.maximum.reduce([self.v1, self.v2, self.v3])
        self.vN0mim_star = -self.Vs / 2 - np.minimum.reduce([self.v1, self.v2, self.v3])

        # Tensão homopolar
        self.vN0_star = self.valor_mu * self.vN0max_star + (1 - self.valor_mu) * self.vN0mim_star

        # Tensões moduladas
        self.v10 = self.v1 + self.vN0_star
        self.v20 = self.v2 + self.vN0_star
        self.v30 = self.v3 + self.vN0_star

        # Geração da onda portadora triangular:
        periodo_port = 1 / self.f_onda_p
        onda_port = 1 - 2 * np.abs((self.t_pwm % periodo_port) * self.f_onda_p - 0.5)

        # Sinal PWM
        PWM_signal = np.where(self.v10 >= onda_port, self.Vs, 0)

        # Correntes
        plt.figure(figsize=(8, 4))
        plt.plot(self.t_pwm, self.i1, label='Current PWM 1 (A)')
        plt.plot(self.t_pwm, self.i2, label='Current PWM 2 (A)')
        plt.plot(self.t_pwm, self.i3, label='Current PWM 3 (A)')
        plt.title('Currents PWM (A)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Current PWM (A)')

        # Gráficos dos resultados
        plt.figure(figsize=(10, 8))

        # Tensão modulada
        plt.subplot(4, 1, 1)
        plt.plot(self.t_pwm, self.v10, label='V10', color='black')
        plt.title('Tensão modulada')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')
        plt.legend()

        # Onda triangular portadora
        plt.subplot(4, 1, 2)
        plt.plot(self.t_pwm, onda_port, label='Onda portadora')
        plt.title('Sinal da onda portadora')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')

        # Gráfico do sinal PWM puro
        plt.subplot(4, 1, 3)
        plt.step(self.t_pwm, PWM_signal, label='PWM', color='red')
        plt.title('Sinal PWM para controle do torque')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')

        # Onda de referência
        plt.subplot(4, 1, 4)
        plt.step(self.t_pwm, self.v1, label='Tensão 1', color='green')
        plt.title('Tensão de referência')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')

        plt.tight_layout()
        plt.show()

    def exec_pwm(self):
        self.mu_values = np.linspace(0, 10e-1, 1)  # Variação de mu
        
        if self.valor_mu != 1:
            self.controle_pwm()  
        else:
            for mu in self.mu_values:
                
                self.valor_mu = mu  
                self.controle_pwm()  
