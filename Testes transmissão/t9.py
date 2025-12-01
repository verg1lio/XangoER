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
            self.cl = 0 # 100 Nm
        else:
            self.cl = 0 # 250 Nm (nominal)

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

class Drivetrain:

    def __init__(self, cgx, cgy, massa, massa_roda, entre_eixos, raio_pneu, reducao_primaria, reducao_final, cp, tempo_i):
        self.cgx = cgx
        self.cgy = cgy
        self.massa = massa
        self.massa_roda = massa_roda
        self.entre_eixos =entre_eixos
        self.raio_pneu = raio_pneu
        self.reducao_primaria = reducao_primaria
        self.reducao_final = reducao_final
        self.cp = cp
        self.tempo_i = tempo_i
        self.tempo_f = 0

    def Reduções(self, eficiencia_corrente=0.96):

        # Criar motor corretamente com os parâmetros
        motor = Motor(
            rs=0.04585,       # Resistência do estator
            ld=0.00067,       # Indutância d
            lq=0.00067,       # Indutância q
            jm=0.05769,       # Inércia do rotor
            kf=0.1,           # Coeficiente de atrito
            lambda_m=0.13849, # Fluxo do ímã
            p=10,             # Nº pares de polos
            valor_mu=0.9
        )

        motor.speed_ref = 600  # Pode ajustar a referência
        motor.simulate()

        # Dados do motor
        rpm_motor = np.array(motor.velocidade) * 60 / (2 * np.pi)
        torque_motor = np.array(motor.torque_mecanico)   # <<<< extração correta do torque
        potencia_motor = torque_motor * motor.wm         # potência instantânea

        # Etapa 1: Redução corrente
        rpm_pos_corrente1 = rpm_motor / self.reducao_primaria
        torque_pos_corrente1 = torque_motor * self.reducao_primaria
        potencia_pos_corrente1 = potencia_motor

        # Etapa 2: perdas na corrente
        rpm_pos_corrente = rpm_pos_corrente1
        torque_pos_corrente = torque_pos_corrente1 * eficiencia_corrente
        potencia_pos_corrente = potencia_pos_corrente1 * eficiencia_corrente

        # Etapa 3: Redução final (diferencial)
        rpm_dif = rpm_pos_corrente / self.reducao_final
        torque_dif = torque_pos_corrente * self.reducao_final
        potencia_dif = potencia_pos_corrente

        rendimento_transmissao = (potencia_dif / (potencia_motor + 1e-9))  # evitar div/0
        torque_ef = torque_dif * rendimento_transmissao

        return [torque_ef, potencia_dif, rpm_dif]
        
    def CarPerformance(self):

        # Parâmetros do veículo
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        raio = self.raio_pneu * 0.001               
        c_r = 0.012                         # coeficiente de rolamento
        b = 0.1                            # coeficiente e atrito (transmissão)

        # inicializa as velocidades
        v_angular = 0                       # velocidade angular inicial (rad/s)
        v_linear_ms = 0                     # velocidade linear inicial (m/s)
        

        # Verifica se tempo final foi definido
        if self.tempo_f == 0:
            raise ValueError("Defina o tempo final manualmente: dt_model.tempo_f = <valor_em_segundos>")

        # Vetor de tempo (resolução de 0.01s)
        variacao_tempo = np.arange(self.tempo_i, self.tempo_f, 0.01)
        n_amostras = len(variacao_tempo)

        # Sinais do motor reduzidos pelo sistema de transmissão
        torque_raw, _, rpm_raw = self.Reduções()

        # Função para expandir vetor até n_amostras
        def expandir(vetor):
            return np.resize(vetor, n_amostras) if len(vetor) < n_amostras else vetor[:n_amostras]

        # Vetores expandidos
        torque_reduzido = expandir(torque_raw)
        rpm_reduzido = expandir(rpm_raw)

        parametros = []
        
        for torque, tempo in zip(torque_reduzido, variacao_tempo):

            # passo de tempo
            if not parametros:
                dt = 0
            else:
                dt = tempo - parametros[-1]["tempo"]

            # Inércias
            J_roda = 0.5 * self.massa_roda * raio**2        # inércia de uma roda         
            J_carro_equivalente = self.massa * raio**2      # Inércia equivalente da massa do veículo transferida ao eixo
           
            J_total = ((4 * J_roda) + J_carro_equivalente)  # Inércia total sentida pelo motor/trenó de força

            # Torque de carga (rolamento, arrasto e atrito)
            torque_rr = (c_r * peso) * raio
            torque_fa = ((densidade_ar * v_linear_ms**2 * coeficiente_arrasto * area_frontal) / 2) * raio
            torque_atr_mec = b * v_angular

            torque_carga = torque_rr + torque_fa + torque_atr_mec
            torque_efetivo = torque - torque_carga

            # Forças
            forca_trativa = torque_efetivo / raio
            forca_arrasto = torque_fa / raio 
            forca_rr = torque_rr / raio

            forca_resistiva = torque_carga / raio
            forca_efetiva = forca_trativa - forca_resistiva

            # Acelerações
            a_angular = torque_efetivo / J_total                     
            a_linear = forca_efetiva / self.massa

            # Velocidades
            v_angular += a_angular * dt
            v_linear_ms += a_linear * dt
            v_linear_kmh = v_linear_ms * 3.6

            parametros.append({
                "t": torque,
                "va": v_angular,
                "vlm": v_linear_ms,
                "vlk": v_linear_kmh,
                "fa": forca_arrasto,
                "rr": forca_rr,
                "ff": forca_trativa,
                "tempo": tempo
            })

        return parametros, rpm_reduzido, variacao_tempo

    def printCarPerformance(self, performance):
        
        # Impressão dos parâmetros principais
        print("Torque [N/m]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]")
        for i in range(0, len(performance), 10):
            t = performance[i]["t"]
            va = performance[i]["va"]
            vlk = performance[i]["vlk"]
            print(f"{t:.2f}\t\t\t{va:.2f}\t\t\t\t{vlk:.2f}")

        print("\nForça de Arrasto [N]\tResistência ao Rolamento [N]\tForça Trativa [N]")
        for i in range(0, len(performance), 10):
            fa = performance[i]["fa"]
            rr = performance[i]["rr"]
            ff = performance[i]["ff"]
            print(f"{fa:.2f}\t\t\t{rr:.2f}\t\t\t{ff:.2f}")

        # Obtenção do vetor de tempo
        tempo = [p["tempo"] for p in performance]
        n_amostras = len(tempo)

        # Simulação do motor
        motor = Motor(
        rs=0.04585,       # Resistência do estator
        ld=0.00067,       # Indutância d
        lq=0.00067,       # Indutância q
        jm=0.05769,       # Inércia do rotor
        kf=0.1,           # Coeficiente de atrito
        lambda_m=0.13849, # Fluxo do ímã
        p=10,             # Nº pares de polos
        valor_mu=0.9
    )
        motor.simulate()

        # Conversão de rad/s para RPM
        rpm_motor_raw = np.array(motor.velocidade) * 60 / (2 * np.pi)

        # Ajuste do vetor de RPM do motor ao tamanho de amostras
        rpm_motor = np.resize(rpm_motor_raw, n_amostras)

        # Extração de dados adicionais do desempenho
        velocidade_angular = [p["va"] for p in performance]
        velocidade_linear = [p["vlk"] for p in performance]


        # --- Figura 1: RPM x Tempo ---
        plt.figure(figsize=(8, 4))
        plt.plot(tempo, rpm_motor, color='tab:pink')
        plt.xlabel("Tempo [s]")
        plt.ylabel("RPM do Motor")
        plt.title("RPM x Tempo")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        # --- Figura 2: Velocidade Angular e Linear x Tempo ---
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        plt.plot(tempo, velocidade_angular, color='tab:red')
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade Angular [rad/s]")
        plt.title("Velocidade Angular x Tempo")
        plt.grid(True)

        plt.subplot(1, 2, 2)
        plt.plot(tempo, velocidade_linear, color='tab:purple')
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade Linear [km/h]")
        plt.title("Velocidade Linear x Tempo")
        plt.grid(True)

        plt.tight_layout()
        plt.show()

        # Plot dos gráficos
        plt.figure(figsize=(12, 6))

        # Plot Velocidade Angular (rad/s) x RPM
        plt.subplot(1, 2, 1)
        plt.plot(rpm_motor, velocidade_angular, color='tab:red')
        plt.xlabel("RPM Motor")
        plt.ylabel("Velocidade Angular [rad/s]")
        plt.title("Velocidade Angular x RPM")
        plt.grid(True)

        # Plot Velocidade Linear (km/h) x RPM
        plt.subplot(1, 2, 2)
        plt.plot(rpm_motor, velocidade_linear, color='tab:purple')
        plt.xlabel("RPM Motor")
        plt.ylabel("Velocidade Linear [km/h]")
        plt.title("Velocidade Linear x RPM")
        plt.grid(True)

        plt.tight_layout()
        plt.show()

class Tire:
    
    def __init__(self, tire_Fz=None, tire_Sa=None, tire_Ls=None, tire_friction_coef=None, tire_Ca=0, B0=0, B1=0, B2=1, B3=1, omega=315, slip_angle_start=-9, slip_angle_end=9, angle_step=0.5, WB=1500, rear_axle_length=1600, track_y=0, tire_k=0):
        
        # Inicializando parâmetros para o modelo de pneu
        self.tire_Fz = tire_Fz  # Carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # Ângulo de deslizamento lateral do pneu [rad]
        self.tire_Ls = tire_Ls  # Escorregamento longitudinal do pneu [Adimensional]
        self.tire_type = 'Default'  # Tipo de pneu
        self.tire_friction_coef = tire_friction_coef  # Coeficiente de fricção entre pneu e pista
        self.tire_Ca = tire_Ca  # Ângulo de camber do pneu
            
        # Comprimentos das barras do mecanismo de quatro barras
        self.L0 = B0  # Comprimento da barra de direção
        self.L1 = B1  # Comprimento do braço de direção
        self.L2 = B2  # Comprimento da bitola
        self.L3 = B3  # Comprimento do braço de direção

        # Ângulo de orientação da barra longitudinal em relação ao eixo horizontal
        self.alpha = np.radians(omega)

        # Array de ângulos de deslizamento lateral do pneu em graus
        self.theta2 = np.arange(slip_angle_start + 90, slip_angle_end + 91, angle_step)

        # Convertendo os ângulos de deslizamento lateral do pneu para radianos
        self.theta2_rad = np.radians(self.theta2)
        self.angle = self.theta2_rad[0]

        # Inicialização das listas de resultados para armazenar dados ao longo do cálculo
        self.AC = []  # Lista para armazenar AC
        self.beta = []  # Lista para armazenar beta
        self.psi = []  # Lista para armazenar psi
        self.lamda = []  # Lista para armazenar lambda
        self.theta3 = []  # Lista para armazenar theta3
        self.theta4 = []  # Lista para armazenar theta4
        self.Ox, self.Oy = 0, 0  # Lista para armazenar coordenadas de O
        self.Ax, self.Ay = [], []  # Listas para armazenar coordenadas de A
        self.Bx, self.By = [], []  # Listas para armazenar coordenadas de B
        self.Cx, self.Cy = B0, 0 # Listas para armazenar coordenadas de C
        self.w = []  # Lista para armazenar w
        self.om2, self.om4 = [], []  # Listas para armazenar om2 e om4
        self.alpha_dot = []  # Lista para armazenar alpha_dot
        self.outer_slip = []  # Lista para armazenar ângulos de deslizamento lateral externo
        self.inner_slip = []  # Lista para armazenar ângulos de deslizamento lateral interno
        self.static_slip_angle = None  # Variável para armazenar o ângulo de deslizamento estático (inicialmente não definido)

        # Coordenadas da barra longitudinal (entre-eixos)
        self.WB = WB  # Entre-eixos (wheelbase) fixo
        self.rear_axle_length = rear_axle_length  # Comprimento do eixo traseiro fixo
        self.rear_axle_center = (B0 / 2, 0)

        # Entradas de rigidez do pneu
        self.track_y = track_y
        self.tire_k = tire_k   
        
    def Tire_forces(self, params):
        
            # Desembalando os parâmetros de Pacejka
            E, Cy, Cx, Cz, c1, c2 = params

            # Calculando parâmetros intermediários
            Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))
            D = self.tire_friction_coef * self.tire_Fz
            Bz = Cs / (Cz * D)
            Bx = Cs / (Cx * D)
            By = Cs / (Cy * D)

            # Calculando a força lateral do pneu
            tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))

            # Calculando a força longitudinal do pneu
            tire_longitudinal_force = D * np.sin(Cx * np.arctan(9 * Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))

            # Calculando o momento de auto-alinhamento do pneu
            tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))

            # Calculando a força de camber
            camber_thrust = D * np.sin(Cy * np.arctan(By * self.tire_Ca))

            # Retornando as forças calculadas e o momento de auto-alinhamento
            return tire_lateral_force + 0.5 * camber_thrust, (10 + (tire_auto_align_moment / 55)), tire_longitudinal_force
    
    @staticmethod
    def SlipRatio(velocidade_angular, raio_pneu, velocidade_linear):
        value = (velocidade_angular * raio_pneu / velocidade_linear) - 1
        return value

    def printSlipRatio(self, tempo, slip_ratio, tire_longitudinal_forces):
        passo = 10  # Intervalo para impressão
        slip_ratio_constante = max(slip_ratio)  # <<< Defina aqui o valor constante desejado
        i = 0

        # Loop while para corrigir os valores de slip_ratio nos primeiros instantes
        while i < len(tempo) and tempo[i] <= 0.05:
            slip_ratio[i] = slip_ratio_constante
            i += 1

        print("Valores do Slip Ratio:")
        for i in range(0, len(slip_ratio), passo):
            print(f"{slip_ratio[i]:.2f}")

        # Figura com dois subgráficos
        plt.figure(figsize=(14, 6))

        # Gráfico 1: Slip Ratio x Tempo
        plt.subplot(1, 2, 1)
        plt.plot(tempo, slip_ratio, label='Slip Ratio', color='blue')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Slip Ratio')
        plt.title('Slip Ratio x Tempo')
        plt.grid(True)
        plt.legend()

        # Gráfico 2: Força Longitudinal x Slip Ratio
        plt.subplot(1, 2, 2)
        plt.plot(slip_ratio, tire_longitudinal_forces, label='Força Longitudinal', color='green')
        plt.xlabel('Slip Ratio')
        plt.ylabel('Força Longitudinal (N)')
        plt.title('Força Longitudinal x Slip Ratio')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()
        
def Instancias():

    # Instância do modelo Drivetrain com os parâmetros fornecidos
    dt_model = Drivetrain(
        cgx=853,                # Centro de gravidade no eixo X [mm]
        cgy=294,                # Centro de gravidade no eixo Y [mm]
        massa=347,              # Massa do veículo [kg]
        massa_roda= 6,          # Massa da Roda [kg]    
        entre_eixos=1567,       # Entre-eixos [mm]
        raio_pneu=220,          # Raio do pneu [mm]
        reducao_primaria=2.12,  # Redução primária
        reducao_final=2.76,     # Redução final (diferencial)
        cp=2.22,                # Constante de potência (pode ser usada futuramente)
        tempo_i=0               # Tempo inicial [s]
    )

    # tempo (s)
    dt_model.tempo_f = 20

    # Cálculo da performance do veículo com os valores reduzidos
    performance_veiculo, variacao_rpm, variacao_tempo = dt_model.CarPerformance()

    # Impressão dos resultados
    dt_model.printCarPerformance(performance_veiculo)

    # Extração das velocidadesdos dados de performance
    velocidade_angular = np.array([dado["va"] for dado in performance_veiculo])
    velocidade_linear = np.array([dado["vlm"] for dado in performance_veiculo])
  
    # Cálculo do slip ratio
    slip_ratio = Tire.SlipRatio(velocidade_angular, dt_model.raio_pneu * 0.001, velocidade_linear)
        
    # Instância da classe Tire com os slip ratios calculados
    slip_model = Tire(
        tire_Fz= 850,               # Carga vertical no pneu [N]
        tire_Sa=0,                 # Ângulo de escorregamento lateral [rad]
        tire_Ls=slip_ratio,       # Slip ratio já calculado
        tire_friction_coef=1.45,   # Coeficiente de fricção
        tire_Ca=0                  # Ângulo de camber
    )

    # Parâmetros experimentais fornecidos
    result = [
        0.3336564873588197,
        1.6271741344929977,
        1,
        4.3961693695846655,
        931.4055775279057,
        366.4936818126405,
    ]

    # Cálculo das forças e do momento de auto-alinhamento
    tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = slip_model.Tire_forces(result)

    # Plotagem dos gráficos de slip ratio
    slip_model.printSlipRatio(variacao_tempo, slip_ratio, tire_longitudinal_forces)
    
Instancias()
