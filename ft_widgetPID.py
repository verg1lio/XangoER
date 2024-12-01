import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from scipy import signal
import control as ctl



__all__ = ["Motor", "ModulacaoEscalar", "Peso"]
# Escolha de quais listas de gráficos exibir listas em (def setup_plots(self, i))
escolha = [1]  # Altere para [0], [1], ou [0, 1], conforme desejado
velocidade_desejada = 150

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
        self.tmax = 10  # Maximum simulation time (s)
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
        self.erros = []
        self.graficos = []  # Lista para armazenar funções de plotagem


    def reset_initial_conditions(self):
        # Initialize conditions
        self.cl = 0  # Load torque (N*m)
        self.wm = 0.0  # Mechanical speed (rad/s)
        self.t = 0  # Time (s)
        self.tp = 0  # Plotting time (s)
        self.j = 0  # Plotting index
        self.ce = 0  # Electromagnetic torque (N*m)
        
        self.freq = 200 # frequencia  
        self.V = 72 # tensão rms

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
        
        self.erro_acumulado = 0
        self.erro_anterior = 0
        self.wm_anterior = 0
        # Parametros PID
        self.P_cr = 0.4125 #periodo critico
        self.K_cr = 1.e-3/0.5 # ganho crítico
        self.K_p = 5000*self.K_cr # Ganho proporcional PID
        self.K_i = 0.5*self.P_cr # Ganho integrativo
        self.K_d = 0.125*self.P_cr # Ganho derivativo

    def source_voltage(self,):
        self.Vsm = self.V * np.sqrt(2)  # Peak stator voltage (V)
        self.Vs = self.Vsm  # Stator voltage (V)
        self.ws = 2*np.pi*self.freq/self.p # Synchronous speed 6000 RPM (628,32 rad/s)
        self.tete += self.h * self.ws
        if self.tete >= 2 * np.pi:
            self.tete -= 2 * np.pi
            
        vs1 = self.Vs * np.cos(self.tete)
        vs2 = self.Vs * np.cos(self.tete - self.pi23)
        vs3 = self.Vs * np.cos(self.tete + self.pi23)
        return vs1, vs2, vs3

    def load_torque(self,):
        if self.t>=3:
            self.cl = 1
    
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

    def controlePID(self,):
        erro = velocidade_desejada - self.wm 
       
        P = self.K_p * erro # Termo proporcional
        self.erro_acumulado += erro * self.h  # Termo integrativo
        I = self.K_i * self.erro_acumulado
        D = self.K_d*(erro-self.erro_anterior)/self.h # Termo derivativo

         # Atualizar o erro anterior
        self.erro_anterior = erro

        self.V += P + I + D
        self.V = max(24, min(72, self.V))
  
        return erro
    
    def outputs(self, is1, is2, is3, fs1, fs2, fs3, fso, cm, vso, vsd, vsq, erro):
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

        self.erros.append(erro)


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
            erro = self.controlePID()
            if self.t >= self.tp:
                self.tp += self.hp
                self.outputs(is1, is2, is3, fs1, fs2, fs3, fso, cm, vso, vsd, vsq, erro)

    def graph_current(self, ax):
        # Plotting currents
        self.fig.clear()
        ax = self.fig.add_subplot(1, 1, 1)
        ax.plot(self.tempo, self.corrente1, label='Current 1 (A)')
        ax.plot(self.tempo, self.corrente2, label='Current 2 (A)')
        ax.plot(self.tempo, self.corrente3, label='Current 3 (A)')
        ax.set_title('Currents (A)')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Current (A)')
        ax.legend()
        self.recreate_buttons_and_atualization_canvas()

    def graph_voltages(self, ax):
        # Plotting voltages
        self.fig.clear()
        ax = self.fig.add_subplot(1, 1, 1)
        ax.plot(self.tempo, self.tensao1, label='Voltage 1 (V)')
        ax.plot(self.tempo, self.tensao2, label='Voltage 2 (V)')
        ax.plot(self.tempo, self.tensao3, label='Voltage 3 (V)')
        ax.set_title('Voltages (V)')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Voltage (V)')
        ax.legend()
        self.recreate_buttons_and_atualization_canvas()

    def graph_fluxes(self, ax):
        # Plotting fluxes
        self.fig.clear()
        ax = self.fig.add_subplot(1, 1, 1)
        ax.plot(self.tempo, self.fluxos1, label='Flux 1 (Wb)')
        ax.plot(self.tempo, self.fluxos2, label='Flux 2 (Wb)')
        ax.plot(self.tempo, self.fluxos3, label='Flux 3 (Wb)')
        ax.set_title('Fluxes (Wb)')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Flux (Wb)')
        ax.legend()
        self.recreate_buttons_and_atualization_canvas()

    def graph_zero_current(self, ax):
        # Plotting zero-sequence current
        self.fig.clear()
        ax = self.fig.add_subplot(1, 1, 1)
        ax.plot(self.tempo, self.correnteo, label='Current o (A)')
        ax.set_title('Current o (A)')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Current (A)')
        ax.legend()
        self.recreate_buttons_and_atualization_canvas()

    def graph_erros(self, ax):
        self.fig.clear()
        ax = self.fig.add_subplot(1,1,1)
        ax.plot(self.tempo, self.erros, label='Erros')
        ax.set_title('Erro de rastreamento')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Erro')
        ax.legend()
        self.recreate_buttons_and_atualization_canvas()

    def graph_multiples(self, ax):
        """Plota múltiplos gráficos de torque e velocidade."""
        self.fig.clear()
        axs = self.fig.subplots(2, 2)
        self.fig.subplots_adjust(bottom=0.2)

        axs[0, 0].plot(self.tempo, self.conjcarga, label='Load Torque (N*m)')
        axs[0, 0].set_title('Load Torque (N*m)')
        axs[0, 0].legend()
        axs[0, 0].set_xlabel('Time (s)')
        axs[0, 0].set_ylabel('Torque (N*m)')

        axs[0, 1].plot(self.tempo, self.velocidade, label='Speed (rad/s)')
        axs[0, 1].set_title('Speed (rad/s)')
        axs[0, 1].legend()
        axs[0, 1].set_xlabel('Time (s)')
        axs[0, 1].set_ylabel('Speed (rad/s)')

        axs[1, 0].plot(self.tempo, self.conjugado, label='Electromagnetic Torque (N*m)')
        axs[1, 0].set_title('Electromagnetic Torque (N*m)')
        axs[1, 0].legend()
        axs[1, 0].set_xlabel('Time (s)')
        axs[1, 0].set_ylabel('Torque (N*m)')

        axs[1, 1].plot(self.tempo, self.torque_mecanico, label='Mechanical Torque (N*m)')
        axs[1, 1].set_title('Mechanical Torque (N*m)')
        axs[1, 1].legend()
        axs[1, 1].set_xlabel('Time (s)')
        axs[1, 1].set_ylabel('Torque (N*m)')

        self.recreate_buttons_and_atualization_canvas()

    def setup_plots(self, i):
        if i == 0:
            self.graficos = [
                self.graph_current,
                self.graph_voltages,
                self.graph_fluxes,
                self.graph_zero_current,
            ]
        elif i==1:
            self.graficos = [
                self.graph_erros,
                self.graph_multiples
            ]
        elif i == 2:
            self.graficos = [
                self.plot_bode,
                self.plot_nyquist,
                self.step_response,
                self.plot_fourier_transform
            ]

    def recreate_buttons_and_atualization_canvas(self):
        """Recria os botões após um fig.clear()."""
        self.fig.subplots_adjust(bottom=0.2)

        ax_prev = plt.axes([0.3, 0.05, 0.1, 0.075])
        self.btn_prev = Button(ax_prev, 'Prev')
        self.btn_prev.on_clicked(self.prev_plot)

        ax_next = plt.axes([0.6, 0.05, 0.1, 0.075])
        self.btn_next = Button(ax_next, 'Next')
        self.btn_next.on_clicked(self.next_plot)

        #Atualizar canvas
        self.fig.canvas.draw_idle()

    def update_plot(self):
        """Atualiza o gráfico com base no índice atual."""
        ax = self.fig.add_subplot(1, 1, 1)  # Cria ou seleciona o eixo para o gráfico
        self.graficos[self.current_index](ax)  # Chama a função do gráfico desejado
        plt.draw()

    def next_plot(self, event):
        # Avança para o próximo gráfico
        self.current_index = (self.current_index + 1) % len(self.graficos)
        self.update_plot()

    def prev_plot(self, event):
        # Volta para o gráfico anterior
        self.current_index = (self.current_index - 1) % len(self.graficos)
        self.update_plot()

    def plot_motor(self, i):
        self.current_index = 0
        """Inicializa os gráficos e exibe."""
        self.fig = plt.figure(figsize=(12, 8) if i == 1 else (10, 6) if i == 2 else (8, 6))
        plt.subplots_adjust(bottom=0.2)  # Espaço para os botões
        self.setup_plots(i)
        self.update_plot()
        plt.show()

    def transfer_function(self):
        num = [self.msr* self.p]
        den = [self.jm, 2 * self.kf, self.rs + self.ls]
        return num, den

    def plot_bode(self, ax):
        num, den = self.transfer_function()
        system = signal.TransferFunction(num, den)
        w, mag, phase = signal.bode(system)

        # Limpar a figura
        self.fig.clear()

        # Criar subplots: 2 linhas e 1 coluna
        axs = self.fig.subplots(2, 1)

        # Gráfico de magnitude
        axs[0].semilogx(w, mag)
        axs[0].set_title('Diagrama de Bode - Motor')
        axs[0].set_ylabel('Magnitude (dB)')
        axs[0].grid(which="both", axis="both")

        # Gráfico de fase
        axs[1].semilogx(w, phase)
        axs[1].set_xlabel('Frequência (rad/s)')
        axs[1].set_ylabel('Fase (graus)')
        axs[1].grid(which="both", axis="both")

        # Recriar botões e atualizar canvas
        self.recreate_buttons_and_atualization_canvas()
        
    def plot_nyquist(self, ax):
        num, den = self.transfer_function()
        motor_system = ctl.TransferFunction(num, den)

        # Frequências para o Diagrama de Nyquist
        w_start = 1e-2
        w_stop = 1e3
        num_points = 1000
        frequencies = np.logspace(np.log10(w_start), np.log10(w_stop), num_points)

        self.fig.clear()
        ax = self.fig.add_subplot(1, 1, 1)
        ctl.nyquist_plot(motor_system, omega=frequencies)
        ax.set_title("Diagrama de Nyquist - Motor Trifásico")
        ax.grid(True)
        self.recreate_buttons_and_atualization_canvas()

    def state_space_representation(self):
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

    def step_response(self, ax):
        num, den = self.transfer_function()
        system = signal.TransferFunction(num, den)
        t, response = signal.step(system)

        ax = self.fig.add_subplot(1, 1, 1)
        ax.plot(t, response, label='Resposta ao Degrau Unitário')
        ax.set_title('Resposta ao Degrau Unitário - Sistema de Segunda Ordem')
        ax.set_xlabel('Tempo (s)')
        ax.set_ylabel('Amplitude')
        ax.grid(True)
        ax.legend()
        self.recreate_buttons_and_atualization_canvas()

    def plot_fourier_transform(self, ax):
        # Converte a lista de velocidades e tempos em arrays numpy
        velocidade = np.array(self.velocidade)
        tempo = np.array(self.tempo)

        # Calcula a transformada de Fourier da velocidade
        velocidade_fft = np.fft.fft(velocidade)
        frequencias = np.fft.fftfreq(len(velocidade), d=self.h)  # Frequências correspondentes

        # Apenas valores positivos de frequência e amplitude
        indices_positivos = frequencias > 0
        frequencias = frequencias[indices_positivos]
        amplitudes = np.abs(velocidade_fft[indices_positivos])

        # Plot do espectro de frequência
        self.fig.clear()
        ax = self.fig.add_subplot(1, 1, 1)
        ax.plot(frequencias, amplitudes, label="Espectro de Frequência")
        ax.set_title("Transformada de Fourier da Velocidade")
        ax.set_xlabel("Frequência (Hz)")
        ax.set_ylabel("Amplitude")
        ax.legend()
        ax.grid()
         # Limita os eixos até 1000
        ax.set_xlim(0, 10000)  # Limite do eixo X (frequência)
        ax.set_ylim(0, 10000)  # Limite do eixo Y (amplitude)
        self.recreate_buttons_and_atualization_canvas()

    def example():
        motor = Motor(0.39, 1.41, 0.094, 0.094, 0.091, 0.04, 0.01)
        motor.simulate()

        # Loop pelas escolhas selecionadas
        for i in escolha:
            motor.plot_motor(i)


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

    def example():
        peso_pwt = Peso(10, 10, 65, 5)
        total = peso_pwt.peso_total()
        print(f"O peso total é {total} Kg")



# exemplos
Motor.example()
Peso.example()
