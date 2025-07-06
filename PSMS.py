import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control as clt

class Motor:
    """
    Classe que modela o comportamento de um motor síncrono de Ímã Permanente (PMSM).
    """
    
    def __init__(self, rs, ld, lq, jm, kf, lambda_m, p, q1, q2, q3, valor_mu):
        # Constants
        self.pi23 = 2 * np.pi / 3
        self.rq23 = np.sqrt(2 / 3)
        self.rq3 = np.sqrt(3)

        # Machine parameters
        self.rs = rs      # Stator resistance (ohms)
        self.ld = ld      # Indutância direta (d-axis) (H)
        self.lq = lq      # Indutância em quadratura (q-axis) (H)
        self.jm = jm      # Momento de inércia (kg*m^2)
        self.kf = kf      # Coeficiente de atrito (N*m*s)
        self.lambda_m = lambda_m  # Fluxo magnético do ímã permanente (Wb)
        self.p = p        # Número de pares de polos
        self.f_onda_p = 50 # Frequência da onda portadora PWM
        self.f_onda_m = 5  # Frequência da onda modulante PWM
        self.q1 = q1      # Chave de comutação do inversor
        self.q2 = q2      # Chave de comutação do inversor
        self.q3 = q3      # Chave de comutação do inversor
        self.valor_mu = valor_mu  # Escalar

        # Simulation parameters
        self.h = 1.e-5    # Time step (s)
        self.tmax = 1     # Maximum simulation time (s)
        self.hp = self.tmax / 2000  # Plotting time step (s)
        if self.hp < self.h:
            self.hp = self.h

        # Initial conditions
        self.reset_initial_conditions()
        
        # Storage for output variables
        self.tempo = []        # Time (s)
        self.corrented = []    # Corrente direta (d-axis) (A)
        self.correnteq = []    # Corrente em quadratura (q-axis) (A)
        self.corrente1 = []    # Corrente fase 1 (A)
        self.corrente2 = []    # Corrente fase 2 (A)
        self.corrente3 = []    # Corrente fase 3 (A)
        self.tensao1 = []      # Tensão fase 1 (V)
        self.tensao2 = []      # Tensão fase 2 (V)
        self.tensao3 = []      # Tensão fase 3 (V)
        self.tensaosd = []     # Tensão direta (d-axis) (V)
        self.tensaosq = []     # Tensão em quadratura (q-axis) (V)
        self.fluxosd = []      # Fluxo direto (d-axis) (Wb)
        self.fluxosq = []      # Fluxo em quadratura (q-axis) (Wb)
        self.fluxos1 = []      # Fluxo fase 1 (Wb)
        self.fluxos2 = []      # Fluxo fase 2 (Wb)
        self.fluxos3 = []      # Fluxo fase 3 (Wb)
        self.conjugado = []    # Torque eletromagnético (N*m)
        self.velocidade = []   # Velocidade mecânica (rad/s)
        self.frequencia = []   # Frequência elétrica (rad/s)
        self.conjcarga = []    # Torque de carga (N*m)
        self.torque_mecanico = [] # Torque mecânico (N*m)
        self.temperatura = []  # Temperatura (K)

    def reset_initial_conditions(self):
        # Initialize conditions
        self.cl = 0       # Torque de carga (N*m)
        self.wm = 0.0     # Velocidade mecânica (rad/s)
        self.t = 0        # Tempo (s)
        self.tp = 0       # Plotting time (s)
        self.ce = 0       # Torque eletromagnético (N*m)
        self.ws = 377*10     # Velocidade síncrona (rad/s)
        self.Vsm = 340 / np.sqrt(3)   # Tensão de pico (V)
        self.Vs = self.Vsm # Tensão do estator (V)
        self.tete = 0     # Ângulo elétrico (rad)
        self.isd = 0      # Corrente d-axis (A)
        self.isq = 0      # Corrente q-axis (A)
        self.iso = 0      # Corrente de sequência zero (A)
        self.temp = 25    # Temperatura (C°)
        self.m = 22   # Massa do estator (Kg)
        self.C = 0.385    # Capacidade térmica específica (J/(kg·K))

    def source_voltage(self):
        # Usar velocidade elétrica real (sincronismo)
        we = self.p * self.wm
        self.tete += self.h * we
        
        # Garantir ângulo dentro de [0, 2π]
        self.tete %= 2 * np.pi
        
        # Gerar tensões alinhadas ao quadrature-axis (controle vetorial)
        vs1 = self.Vs * np.sin(self.tete)
        vs2 = self.Vs * np.sin(self.tete - self.pi23)
        vs3 = self.Vs * np.sin(self.tete + self.pi23)
        
        return vs1, vs2, vs3

    def load_torque(self):
        """Aplica torque de carga após metade do tempo"""
        if self.t >= self.tmax / 2:
            self.cl = 40        

    def direct_voltage_and_quadrature(self, vs1, vs2, vs3):
        """Transformação ABC para DQ0"""
        vsd = self.rq23 * (vs1 - vs2 / 2 - vs3 / 2)
        vsq = self.rq23 * (vs2 * self.rq3 / 2 - vs3 * self.rq3 / 2)
        vso = (vs1 + vs2 + vs3) / self.rq3

        return vsd, vsq, vso

    def calculate_derivatives(self, vsd, vsq, vso):
        we = self.p * self.wm  # Velocidade elétrica
        
        # Adicionar pequena constante para evitar divisão por zero
        epsilon = 1e-6
        ld_eff = max(self.ld, epsilon)
        lq_eff = max(self.lq, epsilon)
        
        # Calcular derivadas com limitação
        dervisd = (vsd - self.rs * self.isd + we * lq_eff * self.isq) / ld_eff
        dervisq = (vsq - self.rs * self.isq - we * (ld_eff * self.isd + self.lambda_m)) / lq_eff
        derviso = (vso - self.rs * self.iso) / max(0.1 * (ld_eff + lq_eff)/2, epsilon)
        
        return dervisd, dervisq, derviso

    def update_currents(self, dervisd, dervisq, derviso):
        # Limitar taxas de variação antes de atualizar
        dervisd = np.clip(dervisd, -1e6, 1e6)
        dervisq = np.clip(dervisq, -1e6, 1e6)
        
        self.isd += dervisd * self.h
        self.isq += dervisq * self.h
        self.iso += derviso * self.h
        
        # Aplicar limites físicos baseados na folha de dados
        max_current = 220 * np.sqrt(2)  # 220A RMS -> 311A pico
        current_magnitude = np.sqrt(self.isd**2 + self.isq**2)
        
        if current_magnitude > max_current:
            scaling = max_current / current_magnitude
            self.isd *= scaling
            self.isq *= scaling
        
        # Atualiza fluxos
        self.flux_d = self.ld * self.isd + self.lambda_m
        self.flux_q = self.lq * self.isq
        flux_o = 0.1 * (self.ld + self.lq)/2 * self.iso  # Fluxo de sequência zero
        
        return flux_o

    def calculate_electromagnetic_torque(self):
        # Verificar valores válidos antes do cálculo
        if np.isfinite(self.isd) and np.isfinite(self.isq):
            self.ce = 1.5 * self.p * (self.lambda_m * self.isq)
            # Termo adicional apenas se Ld != Lq
            if abs(self.ld - self.lq) > 1e-6:
                self.ce += 1.5 * self.p * (self.ld - self.lq) * self.isd * self.isq
        else:
            self.ce = 0
        return self.ce

    def currents_and_fluxes_phases(self, flux_o):
        """Transformação DQ0 para ABC"""
        # Correntes
        is1 = self.rq23 * self.isd + self.iso / self.rq3
        is2 = self.rq23 * (-self.isd / 2 + self.rq3 * self.isq / 2) + self.iso / self.rq3
        is3 = self.rq23 * (-self.isd / 2 - self.rq3 * self.isq / 2) + self.iso / self.rq3

        # Fluxos
        fs1 = self.rq23 * self.flux_d + flux_o / self.rq3
        fs2 = self.rq23 * (-self.flux_d / 2 + self.rq3 * self.flux_q / 2) + flux_o / self.rq3
        fs3 = self.rq23 * (-self.flux_d / 2 - self.rq3 * self.flux_q / 2) + flux_o / self.rq3
        
        return is1, is2, is3, fs1, fs2, fs3

    def mechanical_speed(self):
        """Dinâmica mecânica"""
        wm = self.wm + (self.ce - self.cl - self.wm * self.kf) * self.h / self.jm
        return wm
    
    def mechanical_torque(self):
        """Torque mecânico resultante"""
        cm = self.ce - self.cl
        return cm

    def calcular_temperatura(self, h):
        try:
            corrente_eficaz = np.sqrt(self.isd**2 + self.isq**2) / np.sqrt(2)
            potencia_perdida = 3 * self.rs * (corrente_eficaz**2)
            dT = (potencia_perdida * h) / (self.m * self.C)
            self.temp += dT
            return np.clip(self.temp, 20, 120)  # Limitar entre 20°C e 120°C
        except:
            return self.temp  # Manter último valor válido

    def outputs(self, is1, is2, is3, fs1, fs2, fs3, flux_o, cm, vso, vsd, vsq):
        self.tempo.append(self.t)
        self.corrented.append(self.isd)
        self.correnteq.append(self.isq)
        self.corrente1.append(is1)
        self.corrente2.append(is2)
        self.corrente3.append(is3)
        self.tensao1.append(vso)
        self.tensao2.append(vsd)
        self.tensao3.append(vsq)
        self.fluxosd.append(self.flux_d)
        self.fluxosq.append(self.flux_q)
        self.fluxos1.append(fs1)
        self.fluxos2.append(fs2)
        self.fluxos3.append(fs3)
        self.conjugado.append(self.ce)
        self.velocidade.append(self.wm)
        self.frequencia.append(self.ws)
        self.torque_mecanico.append(cm)
        self.conjcarga.append(self.cl) 
        self.temperatura.append(self.temp)

    def simulate(self):
        while self.t < self.tmax:
            self.t += self.h
            vs1, vs2, vs3 = self.source_voltage()
            self.load_torque()
            vsd, vsq, vso = self.direct_voltage_and_quadrature(vs1, vs2, vs3)
            dervisd, dervisq, derviso = self.calculate_derivatives(vsd, vsq, vso)
            flux_o = self.update_currents(dervisd, dervisq, derviso)
            self.calculate_electromagnetic_torque()
            is1, is2, is3, fs1, fs2, fs3 = self.currents_and_fluxes_phases(flux_o)
            self.wm = self.mechanical_speed()
            cm = self.mechanical_torque()
            self.temp = self.calcular_temperatura(self.h)
            
            if self.t >= self.tp:
                self.tp += self.hp
                self.outputs(is1, is2, is3, fs1, fs2, fs3, flux_o, cm, vso, vsd, vsq)

    def plot_motor(self):

        # Plota as correntes das fases
        '''plt.figure(1)
        plt.plot(self.tempo, self.corrente1, label='Current 1 (A)')
        plt.plot(self.tempo, self.corrente2, label='Current 2 (A)')
        plt.plot(self.tempo, self.corrente3, label='Current 3 (A)')
        plt.title('Currents (A)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Current (A)')
        plt.show()

        # Plota as tensões das fases
        plt.figure(2)
        plt.plot(self.tempo, self.tensao1, label='Voltage 1 (V)')
        plt.plot(self.tempo, self.tensao2, label='Voltage 2 (V)')
        plt.plot(self.tempo, self.tensao3, label='Voltage 3 (V)')
        plt.title('Voltages (V)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Voltage (V)')
        plt.show()

        # Plota os fluxos das fases
        plt.figure(3)
        plt.plot(self.tempo, self.fluxos1, label='Flux 1 (Wb)')
        plt.plot(self.tempo, self.fluxos2, label='Flux 2 (Wb)')
        plt.plot(self.tempo, self.fluxos3, label='Flux 3 (Wb)')
        plt.title('Fluxes (Wb)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Flux (Wb)')
        plt.show()

        # Plotando a corrente homopolar
        plt.figure(5)
        plt.plot(self.tempo, self.corrented, label='Current d (A)')
        plt.plot(self.tempo, self.correnteq, label='Current q (A)')
        plt.title('Current o (A)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Current (A)')
        plt.show()

       # Plota a temperatura do motor
        plt.figure(6)
        plt.plot(self.tempo, self.temperatura, label='Temperature C°')
        plt.title('Stator Temperature')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (C°)')
        plt.show()'''

        # Plota múltiplos gráficos em uma única figura
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 2, 1)
        plt.plot(self.tempo, self.conjcarga, label='Load Torque (N*m)')
        plt.title('Load Torque (N*m)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Torque (N*m)')

        plt.subplot(2, 2, 2)
        plt.plot(self.tempo, self.velocidade, label='Speed (rad/s)')
        plt.title('Speed (rad/s)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Speed (rad/s)')

        plt.subplot(2, 2, 3)
        plt.plot(self.tempo, self.conjugado, label='Electromagnetic Torque (N*m)')
        plt.title('Electromagnetic Torque (N*m)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Torque (N*m)')

        plt.subplot(2, 2, 4)
        plt.plot(self.tempo, self.torque_mecanico, label='Mechanical Torque (N*m)')
        plt.title('Mechanical Torque (N*m)')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Torque (N*m)')

        plt.tight_layout()
        plt.show()

    @staticmethod

    def example():
        # Parâmetros típicos de um PMSM
        motor = Motor(
        rs=0.04585,          # Resistência do estator (45.85 mΩ)
        ld=0.00067,          # Indutância do eixo direto (670 µH)
        lq=0.00067,          # Indutância do eixo em quadratura (670 µH)
        jm=0.05769,          # Inércia do rotor (0.05769 kg·m²)
        kf=0.05,             # Coeficiente de atrito (estimado)
        lambda_m=0.13849,    # Fluxo magnético do ímã (0.13849 Wb)
        p=10,                # Número de pares de polos
        q1=1, q2=1, q3=0,    # Configuração das chaves do inversor
        valor_mu=1           # Escalar para controle PWM
        )
        motor.simulate()
        motor.plot_motor()

# Classes Controle e Peso permanecem inalteradas
# ...

# Executar exemplo
if __name__ == "__main__":
    Motor.example()
