import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# 1. BIBLIOTECA DE FÍSICA
# Funções puras que podem ser reutilizadas em qualquer parte do código
# ==============================================================================

def calc_reynolds(rho, velocidade, diametro_hidraulico, viscosidade_dinamica):
    """Calcula o Número de Reynolds (Adimensional)."""
    # Evita divisão por zero
    if viscosidade_dinamica <= 0: return 0.0
    return (rho * velocidade * diametro_hidraulico) / viscosidade_dinamica

def calc_nusselt_interno(Re, Pr):
        """
        Calcula Nu automaticamente baseado no Reynolds.
        - Se Re < 2300: Regime Laminar (Nu Constante = 3.66)
        - Se Re >= 2300: Regime Transição/Turbulento (Correlação de Gnielinski)
        
        Nota: Gnielinski é escolhida pois Dittus-Boelter tem erros grandes
        na faixa de transição (2300 < Re < 10000), comum em variações de vazão.
        """
        # Proteção contra Re negativo ou zero
        if Re <= 1.0: 
            return 3.66

        if Re < 2300.0:
            # --- REGIME LAMINAR ---
            # Para tubos circulares com temperatura superficial constante (T_s = cte)
            # Valor teórico exato: 3.66
            return 3.66 
        
        else:
            # --- REGIME TRANSIÇÃO E TURBULENTO ---
            # Correlação de Gnielinski (Válida para 2300 < Re < 5x10^6)
            # 1. Fator de atrito de Darcy (f) para tubos lisos (Fórmula de Petukhov)
            f = (0.79 * np.log(Re) - 1.64)**(-2)
            
            # 2. Equação de Gnielinski
            num = (f / 8.0) * (Re - 1000.0) * Pr
            den = 1.0 + 12.7 * (f / 8.0)**0.5 * (Pr**(2.0/3.0) - 1.0)
            
            Nu = num / den
            # Segurança física: Nu turbulento não deve ser menor que o laminar na transição imediata
            return max(Nu, 3.66)

def calc_nusselt_externo_aletas(Re, Pr):
    """
    Calcula Nu para o ar passando pelas aletas (Louvered fins ou similar).
    Correlação aproximada para radiadores compactos.
    """
    if Re <= 0: return 0.0
    # Correlação típica para fluxo cruzado em bancos de tubos/aletas
    # Hilpert ou similar para geometria compacta
    return 0.683 * (Re**0.466) * (Pr**(1/3))

def calc_coeficiente_convectivo_h(Nu, k_fluido, diametro_hidraulico):
    """Calcula 'h' [W/m^2K] a partir do Nusselt."""
    if diametro_hidraulico <= 0: return 0.0
    return (Nu * k_fluido) / diametro_hidraulico

def calc_efetividade_ntu_cruzado(C_min, C_max, UA):
    """
    Calcula a efetividade (epsilon) para Trocador de Calor Fluxo Cruzado (Cross-flow).
    Ambos fluidos 'unmixed' (padrão automotivo).
    """
    if C_min <= 0 or UA <= 0: return 0.0
    
    C_r = C_min / C_max
    NTU = UA / C_min
    
    # Equação FSAE/Incropera para Cross-flow unmixed/unmixed
    parte_exp = np.exp(-C_r * (NTU**0.78)) - 1
    epsilon = 1 - np.exp((1/C_r) * (NTU**0.22) * parte_exp)
    
    # Trava física (não pode trocar mais de 100%)
    return min(max(epsilon, 0.0), 0.999)

def atualizar_temp_massa(massa, cp, calor_input, calor_output, dt, temp_atual):
    """
    Integração de Euler para atualizar temperatura de uma massa térmica.
    T_nova = T_velha + (Saldo_Energia / Capacidade_Termica) * dt
    """
    saldo_q = calor_input - calor_output
    dT = saldo_q / (massa * cp)
    return temp_atual + dT * dt

# ==============================================================================
# 2. CLASSES DE COMPONENTES
# ==============================================================================

class PropriedadesFluido:
    def __init__(self):
        # Mix 50/50 Água-Glicol padrão
        pass 
        
    def get_props(self, T_celsius):
        T = T_celsius + 273.15
        rho = 1061 - 0.6 * (T - 293)
        cp = 3300 + 2.5 * (T - 293)
        k = 0.38 + 0.0005 * (T - 293)
        mu = 0.0035 * np.exp(-0.025 * (T - 293))
        mu = max(mu, 0.0002) # Limite físico
        Pr = (cp * mu) / k
        return rho, cp, k, mu, Pr

class MotorEmrax268:
    def __init__(self):
        # Dados do Emrax 268
        self.massa = 22.0          # kg
        self.cp_material = 900.0   # J/kgK (Mix Alumínio/Cobre)
        self.eff = 0.92            # Eficiência média
        self.temp_nucleo = 25.0    # Estado inicial
        
        # Geometria Interna (Water Jacket)
        self.Dh = 0.010            # Diâmetro hidráulico estimado [m]
        self.area_contato = 0.15   # Área superficial interna [m^2]
        
    def calcular_troca_termica(self, potencia_eletr_kw, T_agua_entrada, m_dot_agua, props_fluido, dt):
        """Executa um passo de simulação do motor."""
        rho, cp_f, k_f, mu, Pr = props_fluido
        
        # 1. Calcular Calor Gerado (Perda Joule + Atrito)
        Q_gerado = (potencia_eletr_kw * 1000) * (1 - self.eff)
        
        # 2. Determinar Coeficiente de Convecção (h) interno
        # Velocidade estimada nos canais internos
        vel_fluido = m_dot_agua / (rho * (np.pi * (self.Dh/2)**2)) 
        
        Re = calc_reynolds(rho, vel_fluido, self.Dh, mu)
        Nu = calc_nusselt_interno(Re, Pr)
        h_interno = calc_coeficiente_convectivo_h(Nu, k_f, self.Dh)
        
        # 3. Calor transferido para a água (Lei de Resfriamento de Newton)
        # Q = h * A * (T_parede - T_fluido)
        # Usamos T_agua_entrada para simplificar o passo (método explícito)
        Q_absorvido_agua = h_interno * self.area_contato * (self.temp_nucleo - T_agua_entrada)
        
        # 4. Atualizar temperatura da massa do motor
        self.temp_nucleo = atualizar_temp_massa(
            self.massa, self.cp_material, Q_gerado, Q_absorvido_agua, dt, self.temp_nucleo
        )
        
        # 5. Calcular temperatura da água saindo do motor
        # Balanço na água: Q = m * cp * (T_out - T_in)
        if m_dot_agua > 0:
            T_agua_saida = T_agua_entrada + Q_absorvido_agua / (m_dot_agua * cp_f)
        else:
            T_agua_saida = T_agua_entrada + 0.1 # Segurança p/ vazão zero
            
        return T_agua_saida, Q_absorvido_agua, h_interno

class RadiadorFSAE:
    def __init__(self):
        # Geometria
        self.width = 0.3
        self.height = 0.3
        self.depth = 0.03
        self.area_frontal = self.width * self.height
        
        # Detalhes construtivos
        self.Dh_tubos = 0.003     # Diâmetro hidr. dos tubos de água
        self.Dh_ar = 0.004        # Diâmetro hidr. das aletas de ar
        
        self.area_troca_agua = 0.8 # m^2 (total interno)
        self.area_troca_ar = 3.5   # m^2 (total aletas)
        self.eta_aleta = 0.85     # Eficiência da aleta
        
    def calcular_performance(self, T_agua_entrada, m_dot_agua, vel_ar, props_fluido, T_amb=25.0):
        """Calcula Q dissipado e temperaturas de saída."""
        rho, cp_f, k_f, mu_f, Pr_f = props_fluido
        
        # --- Lado da Água ---
        # Velocidade aprox nos tubos (Assumindo 40% de área aberta p/ fluxo)
        area_secao_tubos = self.depth * self.width * 0.4 
        if area_secao_tubos > 0:
            vel_agua = (m_dot_agua / rho) / area_secao_tubos
        else: 
            vel_agua = 0
            
        Re_agua = calc_reynolds(rho, vel_agua, self.Dh_tubos, mu_f)
        Nu_agua = calc_nusselt_interno(Re_agua, Pr_f)
        h_agua = calc_coeficiente_convectivo_h(Nu_agua, k_f, self.Dh_tubos)
        
        # --- Lado do Ar ---
        # Propriedades do Ar (Simplificado constante ou função de T)
        rho_ar, cp_ar, k_ar, mu_ar = 1.16, 1007, 0.026, 1.8e-5
        Pr_ar = 0.7
        
        Re_ar = calc_reynolds(rho_ar, vel_ar, self.Dh_ar, mu_ar)
        Nu_ar = calc_nusselt_externo_aletas(Re_ar, Pr_ar)
        h_ar = calc_coeficiente_convectivo_h(Nu_ar, k_ar, self.Dh_ar)
        
        # --- Coeficiente Global (UA) ---
        # Resistência total = R_conv_ar + R_conv_agua (desprezando condução parede)
        if h_ar <= 0 or h_agua <= 0:
            UA = 0
        else:
            R_ar = 1 / (h_ar * self.area_troca_ar * self.eta_aleta)
            R_agua = 1 / (h_agua * self.area_troca_agua)
            UA = 1 / (R_ar + R_agua)
            
        # --- Método e-NTU ---
        m_dot_ar = rho_ar * vel_ar * self.area_frontal
        C_ar = m_dot_ar * cp_ar
        C_agua = m_dot_agua * cp_f
        
        if C_ar <= 0 or C_agua <= 0:
            return T_agua_entrada, T_amb, 0.0
            
        C_min = min(C_ar, C_agua)
        C_max = max(C_ar, C_agua)
        
        epsilon = calc_efetividade_ntu_cruzado(C_min, C_max, UA)
        
        Q_max = C_min * (T_agua_entrada - T_amb)
        Q_real = epsilon * Q_max
        
        # --- Temperaturas de Saída ---
        T_agua_saida = T_agua_entrada - Q_real / C_agua
        T_ar_saida = T_amb + Q_real / C_ar
        
        return T_agua_saida, T_ar_saida, Q_real

class Reservatorio:
    def __init__(self, volume_litros):
        self.vol_m3 = volume_litros / 1000.0
        self.temp_atual = 25.0 # Inicial
        
    def misturar(self, T_entrada, m_dot, props_fluido, dt):
        """
        Mistura o fluido que chega do radiador com o volume armazenado.
        Retorna a nova temperatura do tanque (que irá para o motor).
        """
        rho, cp, _, _, _ = props_fluido
        massa_res = self.vol_m3 * rho
        
        # Balanço de energia no tanque (Considerando perfeitamente misturado)
        # m_res * cp * dT/dt = m_dot * cp * (T_in - T_tank)
        # dT = (m_dot / m_res) * (T_in - T_tank) * dt
        
        if massa_res > 0:
            fator_renovacao = m_dot / massa_res
            dT = fator_renovacao * (T_entrada - self.temp_atual) * dt
            self.temp_atual += dT
            
        return self.temp_atual

# ==============================================================================
# 3. LOOP PRINCIPAL (SIMULAÇÃO)
# ==============================================================================

def rodar_simulacao_modular():
    # A. Instanciar Objetos
    fluido = PropriedadesFluido()
    motor = MotorEmrax268()
    radiador = RadiadorFSAE()
    reservatorio = Reservatorio(volume_litros=2.5)
    
    # B. Definir Tempo e Inputs
    t_total = 1000 # segundos
    dt = 0.1      # alta resolução temporal
    steps = int(t_total/dt)
    time_array = np.linspace(0, t_total, steps)
    
    # Inputs (Hardcoded para teste, substitua por CSV se necessário)
    potencia_kw = np.ones(steps) * 80.0 # 80kW contínuos
    vel_ar_ms = np.ones(steps) * 20.0   # 15 m/s (54 km/h)
    vazao_lpm = 15.0 # constante
    
    # C. Arrays de Resultados
    log_T_motor = []
    log_T_agua_entrada = []
    log_T_ar_saida = [] # <-- O dado que você pediu
    log_Q_dissipado = []

    print("Iniciando Loop Modular...")
    
    # D. Loop Temporal
    for i in range(steps):
        # 1. Obter estado atual do fluido (baseado na temp do reservatório)
        T_in_motor = reservatorio.temp_atual
        props = fluido.get_props(T_in_motor)
        rho = props[0]
        m_dot = (vazao_lpm / 60000.0) * rho
        
        # 2. Computar MOTOR
        T_out_motor, Q_absorvido, h_int = motor.calcular_troca_termica(
            potencia_eletr_kw=potencia_kw[i],
            T_agua_entrada=T_in_motor,
            m_dot_agua=m_dot,
            props_fluido=props,
            dt=dt
        )
        
        # 3. Computar RADIADOR
        # A água sai do motor e entra no radiador
        T_in_rad = T_out_motor 
        
        T_out_rad, T_ar_out, Q_dissipado = radiador.calcular_performance(
            T_agua_entrada=T_in_rad,
            m_dot_agua=m_dot,
            vel_ar=vel_ar_ms[i],
            props_fluido=props,
            T_amb=25.0
        )
        
        # 4. Computar RESERVATÓRIO (Fecha o ciclo)
        reservatorio.misturar(T_out_rad, m_dot, props, dt)
        
        # 5. Log
        log_T_motor.append(motor.temp_nucleo)
        log_T_agua_entrada.append(T_in_motor)
        log_T_ar_saida.append(T_ar_out)
        log_Q_dissipado.append(Q_dissipado)

    # E. Plotagem
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    ax1.set_xlabel('Tempo (s)')
    ax1.set_ylabel('Temperatura (°C)')
    ax1.plot(time_array, log_T_motor, 'r-', label='Núcleo Motor', linewidth=2)
    ax1.plot(time_array, log_T_agua_entrada, 'b-', label='Água (Entrada Motor)')
    
    # Destacando a temperatura do ar de saída (sua solicitação)
    ax1.plot(time_array, log_T_ar_saida, 'g--', label='Ar (Saída Radiador)', alpha=0.7)
    
    ax1.legend(loc='center right')
    ax1.grid(True)
    ax1.set_title('Simulação Modular: Emrax 268 + Radiador FSAE')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    rodar_simulacao_modular()