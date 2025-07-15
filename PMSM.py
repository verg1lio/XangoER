import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control as clt

class PIController:
    """Controlador PI para FOC"""
    def __init__(self, kp, ki, limit):
        self.kp = kp
        self.ki = ki
        self.limit = limit
        self.integral = 0
        self.prev_error = 0
        
    def update(self, error, dt):
        self.integral += error * dt
        derivative = (error - self.prev_error) / dt if dt > 0 else 0
        self.prev_error = error
        
        # Ação PID
        output = self.kp * error + self.ki * self.integral
        
        # Anti-windup
        if output > self.limit:
            output = self.limit
            self.integral -= error * dt  # Anti-windup
        elif output < -self.limit:
            output = -self.limit
            self.integral -= error * dt  # Anti-windup
            
        return output

class Motor:
    """Classe que modela um motor síncrono de Ímã Permanente (PMSM) com FOC"""
    
    def __init__(self, rs, ld, lq, jm, kf, lambda_m, p, valor_mu):
        # Constantes
        self.pi23 = 2 * np.pi / 3
        self.rq23 = np.sqrt(2 / 3)
        self.rq3 = np.sqrt(3)
        
        # Parâmetros do motor EMRAX 268 HV+42%
        self.rs = rs          # Resistência do estator (45.85 mΩ)
        self.ld = ld          # Indutância do eixo direto (670 µH)
        self.lq = lq          # Indutância do eixo em quadratura (670 µH)
        self.jm = jm          # Inércia do rotor (0.05769 kg·m²)
        self.kf = kf          # Coeficiente de atrito
        self.lambda_m = lambda_m  # Fluxo magnético do ímã (0.13849 Wb)
        self.p = p            # Número de pares de polos (10)
        
        # Controladores FOC
        self.id_controller = PIController(kp=0.5, ki=100, limit=1000)
        self.iq_controller = PIController(kp=0.5, ki=100, limit=1000)
        self.speed_controller = PIController(kp=1.0, ki=0, limit=500)
        
        # Referências
        self.id_ref = 0        # Referência de corrente d (campo)
        self.iq_ref = 0        # Referência de corrente q (torque)
        self.speed_ref = 300   # Velocidade de referência (rad/s mecânico)
        
        # Parâmetros de simulação
        self.h = 1e-5          # Passo de tempo (s)
        self.tmax = 2.0        # Tempo máximo de simulação (s)
        self.hp = self.tmax / 2000  # Passo para plotagem
        
        # Inicialização
        self.reset_initial_conditions()
        
        # Armazenamento de dados
        self.initialize_storage()

    def initialize_storage(self):
        self.tempo = []
        self.corrented = []
        self.correnteq = []
        self.corrente1 = []
        self.corrente2 = []
        self.corrente3 = []
        self.tensao1 = []
        self.tensao2 = []
        self.tensao3 = []
        self.tensaosd = []
        self.tensaosq = []
        self.fluxosd = []
        self.fluxosq = []
        self.conjugado = []
        self.velocidade = []
        self.frequencia = []
        self.conjcarga = []
        self.torque_mecanico = []
        self.temperatura = []
        self.vd_control = []
        self.vq_control = []
        self.speed_error = []

    def reset_initial_conditions(self):
        self.cl = 0           # Torque de carga
        self.wm = 0.0         # Velocidade mecânica
        self.t = 0            # Tempo
        self.tp = 0           # Tempo para plotagem
        self.ce = 0           # Torque eletromagnético
        self.isd = 0          # Corrente d-axis
        self.isq = 0          # Corrente q-axis
        self.iso = 0          # Corrente de sequência zero
        self.theta_e = 0      # Ângulo elétrico
        self.theta_m = 0      # Ângulo mecânico
        self.temp = 25        # Temperatura inicial
        self.m = 22           # Massa do motor (kg)
        self.C = 0.385        # Capacidade térmica
        self.Vdc = 830        # Tensão DC do barramento
        self.Vs = self.Vdc / np.sqrt(3)  # Tensão de fase máxima

    def field_oriented_control(self):
        """Implementação do Field-Oriented Control (FOC)"""
        # 1. Controle de velocidade
        speed_error = self.speed_ref - self.wm
        self.speed_error.append(speed_error)
        torque_ref = self.speed_controller.update(speed_error, self.h)
        
        # 2. Referências de corrente (estratégia id=0)
        self.iq_ref = torque_ref / (1.5 * self.p * self.lambda_m)
        self.id_ref = 0  # Maximiza torque por ampère
        
        # Limitar correntes de referência
        max_current = 220 * np.sqrt(2)  # 220A RMS -> 311A pico
        if abs(self.iq_ref) > max_current:
            self.iq_ref = np.sign(self.iq_ref) * max_current
        
        # 3. Controle de corrente
        error_d = self.id_ref - self.isd
        error_q = self.iq_ref - self.isq
        
        # Termos de desacoplamento
        we = self.p * self.wm  # Velocidade elétrica
        decoupling_d = -we * self.lq * self.isq
        decoupling_q = we * (self.ld * self.isd + self.lambda_m)
        
        # Ações de controle
        vd = self.id_controller.update(error_d, self.h) + decoupling_d
        vq = self.iq_controller.update(error_q, self.h) + decoupling_q
        
        # Armazenar para análise
        self.vd_control.append(vd)
        self.vq_control.append(vq)
        
        return vd, vq

    def inverse_park_transform(self, vd, vq):
        """Transformação inversa de Park (dq -> abc)"""
        # Transformação dq -> αβ
        valpha = vd * np.cos(self.theta_e) - vq * np.sin(self.theta_e)
        vbeta = vd * np.sin(self.theta_e) + vq * np.cos(self.theta_e)
        
        # Transformação αβ -> abc
        v0 = 0  # Sem componente de sequência zero
        vs1 = valpha
        vs2 = -0.5 * valpha + (np.sqrt(3)/2 * vbeta)
        vs3 = -0.5 * valpha - (np.sqrt(3)/2 * vbeta)
        
        return vs1, vs2, vs3, v0

    def calculate_derivatives(self, vsd, vsq, vso):
        """Modelo elétrico do PMSM"""
        we = self.p * self.wm  # Velocidade elétrica
        
        # Equações diferenciais
        dervisd = (vsd - self.rs * self.isd + we * self.lq * self.isq) / self.ld
        dervisq = (vsq - self.rs * self.isq - we * (self.ld * self.isd + self.lambda_m)) / self.lq
        derviso = (vso - self.rs * self.iso) / (0.1 * (self.ld + self.lq)/2)
        
        return dervisd, dervisq, derviso

    def update_currents(self, dervisd, dervisq, derviso):
        """Atualiza correntes e fluxos"""
        # Limitar derivadas para estabilidade numérica
        dervisd = np.clip(dervisd, -1e6, 1e6)
        dervisq = np.clip(dervisq, -1e6, 1e6)
        
        self.isd += dervisd * self.h
        self.isq += dervisq * self.h
        self.iso += derviso * self.h
        
        # Limitar correntes fisicamente
        max_current = 220 * np.sqrt(2)  # 220A RMS -> 311A pico
        current_magnitude = np.sqrt(self.isd**2 + self.isq**2)
        
        if current_magnitude > max_current:
            scaling = max_current / current_magnitude
            self.isd *= scaling
            self.isq *= scaling
        
        # Atualizar fluxos
        self.flux_d = self.ld * self.isd + self.lambda_m
        self.flux_q = self.lq * self.isq
        flux_o = 0.1 * (self.ld + self.lq)/2 * self.iso
        
        return flux_o

    def calculate_electromagnetic_torque(self):
        """Calcula torque eletromagnético"""
        self.ce = 1.5 * self.p * (self.lambda_m * self.isq)
        return self.ce

    def mechanical_dynamics(self):
        """Modelo mecânico"""
        # Atualizar velocidade
        dwm = (self.ce - self.cl - self.kf * self.wm) / self.jm
        self.wm += dwm * self.h
        
        # Atualizar posição
        self.theta_m += self.wm * self.h
        self.theta_e = self.p * self.theta_m  # Atualizar ângulo elétrico
        
        # Torque mecânico
        cm = self.ce - self.cl
        return cm

    def thermal_model(self):
        """Modelo térmico simplificado"""
        try:
            # Corrente RMS por fase
            i_rms = np.sqrt(self.isd**2 + self.isq**2) / np.sqrt(2)
            
            # Perdas no cobre (3 fases)
            copper_losses = 3 * self.rs * i_rms**2
            
            # Variação de temperatura
            dT = (copper_losses * self.h) / (self.m * self.C)
            self.temp += dT
            
            # Limitar temperatura entre 20°C e 120°C
            self.temp = np.clip(self.temp, 20, 520)
        except:
            # Manter último valor válido em caso de erro
            pass
        
        return self.temp

    def load_torque(self):
        """Perfil de torque de carga"""
        if self.t < 0.5:
            self.cl = 0
        elif self.t < 1.0:
            self.cl = 100  # 100 Nm
        else:
            self.cl = 250  # 250 Nm (nominal)

    def abc_currents(self, flux_o):
        """Transformação dq0 -> abc"""
        # Correntes
        is1 = self.rq23 * (self.isd * np.cos(self.theta_e) - self.isq * np.sin(self.theta_e))
        is2 = self.rq23 * (self.isd * np.cos(self.theta_e - self.pi23) - self.isq * np.sin(self.theta_e - self.pi23))
        is3 = self.rq23 * (self.isd * np.cos(self.theta_e + self.pi23) - self.isq * np.sin(self.theta_e + self.pi23))
        
        # Fluxos (apenas para registro)
        fs1 = self.rq23 * self.flux_d
        fs2 = self.rq23 * self.flux_d
        fs3 = self.rq23 * self.flux_d
        
        return is1, is2, is3, fs1, fs2, fs3

    def outputs(self, is1, is2, is3, fs1, fs2, fs3, flux_o, cm, vso, vsd, vsq):
        """Armazena dados para plotagem"""
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
        self.conjugado.append(self.ce)
        self.velocidade.append(self.wm)
        self.torque_mecanico.append(cm)
        self.conjcarga.append(self.cl)
        self.temperatura.append(self.temp)

    def simulate(self):
        """Loop principal de simulação"""
        while self.t < self.tmax:
            self.t += self.h
            
            # 1. Atualizar torque de carga
            self.load_torque()
            
            # 2. Executar controle FOC
            vd, vq = self.field_oriented_control()
            vs1, vs2, vs3, vso = self.inverse_park_transform(vd, vq)
            
            # 3. Modelo elétrico
            dervisd, dervisq, derviso = self.calculate_derivatives(vd, vq, vso)
            flux_o = self.update_currents(dervisd, dervisq, derviso)
            self.calculate_electromagnetic_torque()
            
            # 4. Modelo mecânico
            cm = self.mechanical_dynamics()
            
            # 5. Modelo térmico
            self.thermal_model()
            
            # 6. Transformar para ABC
            is1, is2, is3, fs1, fs2, fs3 = self.abc_currents(flux_o)
            
            # 7. Armazenar dados para plotagem
            if self.t >= self.tp:
                self.tp += self.hp
                self.outputs(is1, is2, is3, fs1, fs2, fs3, flux_o, cm, vso, vd, vq)

    def plot_results(self):
        """Plotagem abrangente dos resultados"""
        plt.figure(figsize=(15, 20))
        
        # Velocidade e torque
        plt.subplot(5, 2, 1)
        plt.plot(self.tempo, self.velocidade, 'b-', linewidth=2)
        plt.title('Velocidade Mecânica')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Velocidade (rad/s)')
        plt.grid(True)
        
        plt.subplot(5, 2, 2)
        plt.plot(self.tempo, self.conjugado, 'r-', label='Torque Elétrico')
        plt.plot(self.tempo, self.conjcarga, 'g--', label='Torque de Carga')
        plt.title('Torque do Motor')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Torque (Nm)')
        plt.legend()
        plt.grid(True)
        
        # Correntes
        plt.subplot(5, 2, 3)
        plt.plot(self.tempo, self.corrented, 'b-', label='Id')
        plt.plot(self.tempo, self.correnteq, 'r-', label='Iq')
        plt.title('Correntes dq')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Corrente (A)')
        plt.legend()
        plt.grid(True)
        
        plt.subplot(5, 2, 4)
        plt.plot(self.tempo, self.corrente1, 'b-', label='Fase 1')
        plt.plot(self.tempo, self.corrente2, 'r-', label='Fase 2')
        plt.plot(self.tempo, self.corrente3, 'g-', label='Fase 3')
        plt.title('Correntes de Fase')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Corrente (A)')
        plt.legend()
        plt.grid(True)
        
        # Tensões
        '''plt.subplot(5, 2, 5)
        plt.plot(self.tempo, self.tensaosd, 'b-', label='Vd')
        plt.plot(self.tempo, self.tensaosq, 'r-', label='Vq')
        plt.title('Tensões de Controle dq')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensão (V)')
        plt.legend()
        plt.grid(True)'''
        
        plt.subplot(5, 2, 6)
        plt.plot(self.tempo, self.tensao1, 'b-', label='Fase 1')
        plt.plot(self.tempo, self.tensao2, 'r-', label='Fase 2')
        plt.plot(self.tempo, self.tensao3, 'g-', label='Fase 3')
        plt.title('Tensões de Fase')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensão (V)')
        plt.legend()
        plt.grid(True)
        
        # Fluxos
        plt.subplot(5, 2, 7)
        plt.plot(self.tempo, self.fluxosd, 'b-', label='Fluxo d')
        plt.plot(self.tempo, self.fluxosq, 'r-', label='Fluxo q')
        plt.title('Fluxos Magnéticos')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Fluxo (Wb)')
        plt.legend()
        plt.grid(True)
        
        # Temperatura
        plt.subplot(5, 2, 8)
        plt.plot(self.tempo, self.temperatura, 'm-')
        plt.title('Temperatura do Motor')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Temperatura (°C)')
        plt.grid(True)
        
        # Sinal de controle
        '''plt.subplot(5, 2, 9)
        plt.plot(self.tempo, self.vd_control, 'b-', label='Vd control')
        plt.plot(self.tempo, self.vq_control, 'r-', label='Vq control')
        plt.title('Sinais de Controle FOC')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensão (V)')
        plt.legend()
        plt.grid(True)'''
        
        # Erro de velocidade
        '''plt.subplot(5, 2, 10)
        plt.plot(self.tempo, self.speed_error, 'k-')
        plt.title('Erro de Velocidade')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Erro (rad/s)')
        plt.grid(True)'''
        
        plt.tight_layout()
        plt.show()

# Configuração e execução da simulação
def simulate_emrax268():
    motor = Motor(
        rs=0.04585,       # 45.85 mΩ
        ld=0.00067,       # 670 µH
        lq=0.00067,       # 670 µH
        jm=0.05769,       # 0.05769 kg·m²
        kf=0.1,           # Coeficiente de atrito estimado
        lambda_m=0.13849, # Fluxo do ímã permanente
        p=10,             # 10 pares de polos
        valor_mu=0.9      # Fator de modulação
    )
    
    # Configurar referências
    motor.speed_ref = 600  # Velocidade desejada (rad/s mecânico) ~ 2865 RPM
    
    # Executar simulação
    motor.simulate()
    motor.plot_results()

if __name__ == "__main__":
    simulate_emrax268()