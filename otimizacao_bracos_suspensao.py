# CÉLULA 1: IMPORTAÇÕES
import numpy as np
from scipy.optimize import differential_evolution



#SISTEMA DE COORDENADAS: X - LONGITUDINAL, Y - LATERAL, Z - VERTICAL



print("Bibliotecas importadas com sucesso.")

# ---
# CÉLULA 2: FUNÇÕES AUXILIARES (Cálculo de Ângulo)
# (Coloque a função 'calculate_angle_between_vectors' da resposta anterior aqui)
def calculate_angle_between_vectors(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    dot_product = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    angle_rad = np.arccos(dot_product)
    return np.degrees(angle_rad)

print("Função auxiliar 'calculate_angle_between_vectors' definida.")

# ---
# CÉLULA 3.A: O SOLVER ESTÁTICO 6-DOF
# (Coloque a função 'solve_statics' que acabei de escrever acima aqui)
def solve_statics(points, F_external, P_external):
    links = {
        'UCA_F': (points['UCA_F'], points['UBJ']),
        'UCA_R': (points['UCA_R'], points['UBJ']),
        'LCA_F': (points['LCA_F'], points['LBJ']),
        'LCA_R': (points['LCA_R'], points['LBJ']),
        'PR':    (points['PR_In'], points['PR_Out']),
        'TR':    (points['TR_In'], points['TR_Out']) # O 6º link de direção
    }
    link_names = list(links.keys())
    A = np.zeros((6, 6))
    b = np.zeros(6)
    
    r_ext = P_external
    M_external = np.cross(r_ext, F_external)
    b[0:3] = -F_external
    b[3:6] = -M_external
    
    for i, name in enumerate(link_names):
        p_chassi, p_manga = links[name]
        v_link = p_manga - p_chassi
        u_link = v_link / np.linalg.norm(v_link)
        r_link = p_manga 
        m_link = np.cross(r_link, u_link)
        
        A[0:3, i] = u_link
        A[3:6, i] = m_link

    try:
        forces_vector = np.linalg.solve(A, b)
        return {
            'UCA_F_force': forces_vector[0], 'UCA_R_force': forces_vector[1],
            'LCA_F_force': forces_vector[2], 'LCA_R_force': forces_vector[3],
            'PR_force':    forces_vector[4], 'TR_force':    forces_vector[5]
        }
    except np.linalg.LinAlgError:
        # Geometria singular
        return {
            'UCA_F_force': 1e9, 'UCA_R_force': 1e9,
            'LCA_F_force': 1e9, 'LCA_R_force': 1e9,
            'PR_force':    1e9, 'TR_force':    1e9
        }

print("Função 'solve_statics' definida.")

# ---
# CÉLULA 3.B: FUNÇÃO DE SIMULAÇÃO (ATUALIZADA)

def run_suspension_simulation(points):
    """
    Função de simulação atualizada.
    Agora chama o solver estático real.
    A cinemática (RC, Cambagem) AINDA É PLACEHOLDER.
    """
    
    # =================================================================
    # 1. CÁLCULOS CINEMÁTICOS (PLACEHOLDER)
    # !!! ATENÇÃO: SUBSTITUA ESTA SEÇÃO !!!
    # (Este é um problema de simulação de movimento, separado da estática)
    rc_height = 0.035
    camber_gain_per_bump = -1.5
    # =================================================================
    
    # =================================================================
    # 2. CÁLCULOS ESTÁTICOS (REAIS)
    # (Substitui os placeholders da image_664922.png)
    
    # Definir a força de carga externa (ex: 2G de frenagem e 1.5G de curva)
    # (Estes são seus parâmetros de entrada FIXOS)
    F_EXT = np.array([0, -2300, 1700])
    
    # Ponto de aplicação (furo da manga)
    P_EXT = (points['UBJ'] + points['LBJ']) / 2.0 
    
    # Chama o solver estático
    force_results = solve_statics(points, F_EXT, P_EXT)
    # =================================================================

    # 3. CÁLCULOS DE ÂNGULO (REAIS)
    try:
        v_LCA_1 = points['LBJ'] - points['LCA_F']; v_LCA_2 = points['LBJ'] - points['LCA_R']
        lca_angle = calculate_angle_between_vectors(v_LCA_1, v_LCA_2)
        
        v_UCA_1 = points['UBJ'] - points['UCA_F']; v_UCA_2 = points['UBJ'] - points['UCA_R']
        uca_angle = calculate_angle_between_vectors(v_UCA_1, v_UCA_2)
        
        v_pushrod = points['PR_In'] - points['PR_Out']
        pushrod_angle_rad = np.arcsin(np.abs(v_pushrod[2]) / np.linalg.norm(v_pushrod))
        pushrod_angle = np.degrees(pushrod_angle_rad)
            
    except Exception as e:
        raise ValueError(f"Geometria inválida para ângulos: {e}")

    # 4. Retornar TODOS os resultados
    # (Combina os dicionários)
    all_results = {
        'rc_height': rc_height,
        'camber_gain': camber_gain_per_bump,
        'lca_angle': lca_angle,
        'uca_angle': uca_angle,
        'pushrod_angle': pushrod_angle,
        # Nós nos importamos com as forças TOTAIS nos braços, não em cada link
        'pushrod_force': force_results['PR_force'],
        'lca_force_total': abs(force_results['LCA_F_force']) + abs(force_results['LCA_R_force']),
        'uca_force_total': abs(force_results['UCA_F_force']) + abs(force_results['UCA_R_force']),
        'tie_rod_force': force_results['TR_force']
    }
    
    return all_results

print("Função 'run_suspension_simulation' atualizada com solver estático.")

# ---
# CÉLULA 4: FUNÇÃO OBJETIVO (ATUALIZADA)

def objective_function(x):
    """
    Função objetivo agora usa 27 variáveis.
    """
    
    # 1. "Desempacotar" o vetor 'x' (AGORA COM 27 VARIÁVEIS)
    points = {
        'UBJ':   x[0:3],   'LBJ':   x[3:6],
        'UCA_F': x[6:9],   'UCA_R': x[9:12],
        'LCA_F': x[12:15], 'LCA_R': x[15:18],
        'PR_In': x[18:21], 'PR_Out':x[21:24],
        'TR_In': x[24:27], 'TR_Out':x[27:30] # <-- ATUALIZADO
    }
    # CORREÇÃO: Os pontos do Tie Rod na manga (TR_Out) e chassi (TR_In)
    # Você tem 10 pontos = 30 variáveis. 
    # Vou redefinir o desempacotamento para 10 pontos:
    points = {
        'UBJ':    x[0:3],   'LBJ':    x[3:6],
        'UCA_F':  x[6:9],   'UCA_R':  x[9:12],
        'LCA_F':  x[12:15], 'LCA_R':  x[15:18],
        'PR_In':  x[18:21], 'PR_Out': x[21:24],
        'TR_In':  x[24:27], # Ponto do chassi (caixa de direção)
        'TR_Out': x[27:30]  # Ponto na manga de eixo
    }

    # 2. Rodar a simulação
    try:
        results = run_suspension_simulation(points)
    except Exception as e:
        return 1e9 # Custo muito alto por falha

    # 3. Calcular o Custo
    
    # --- Alvos (Ajuste conforme necessário) ---
    TARGET_RC_HEIGHT = 0.030
    TARGET_CAMBER_GAIN = -1.2
    TARGET_LCA_ANGLE_MIN = 52; TARGET_LCA_ANGLE_MAX = 52
    TARGET_UCA_ANGLE_MIN = 40; TARGET_UCA_ANGLE_MAX = 40
    TARGET_PUSHROD_ANGLE_MIN = 50; TARGET_PUSHROD_ANGLE_MAX = 60
    
    # --- Pesos (Ajuste conforme necessário) ---
    W_RC = 100.0; W_CAMBER = 50.0
    W_LCA_ANGLE = 10.0; W_UCA_ANGLE = 10.0; W_PUSHROD_ANGLE = 20.0
    
    # Pesos das Forças (para minimizar a magnitude)
    W_PUSHROD_FORCE = 0.000001
    W_LCA_FORCE = 0.000001
    W_UCA_FORCE = 0.000001
    W_TR_FORCE = 0.000001 # Peso para a barra de direção

    # --- Calcular Erros (Intervalo) ---
    err_rc = results['rc_height'] - TARGET_RC_HEIGHT
    err_camber = results['camber_gain'] - TARGET_CAMBER_GAIN
    
    lca_angle = results['lca_angle']
    err_lca_angle = max(0, lca_angle - TARGET_LCA_ANGLE_MAX, TARGET_LCA_ANGLE_MIN - lca_angle)
    
    uca_angle = results['uca_angle']
    err_uca_angle = max(0, uca_angle - TARGET_UCA_ANGLE_MAX, TARGET_UCA_ANGLE_MIN - uca_angle)
    
    pushrod_angle = results['pushrod_angle']
    err_pushrod_angle = max(0, pushrod_angle - TARGET_PUSHROD_ANGLE_MAX, TARGET_PUSHROD_ANGLE_MIN - pushrod_angle)
    
    # --- Calcular Custo das Forças ---
    cost_force_pushrod = results['pushrod_force']
    cost_force_lca = results['lca_force_total']
    cost_force_uca = results['uca_force_total']
    cost_force_tr = results['tie_rod_force']

    # --- Custo Final ---
    cost = (W_RC * (err_rc**2)) + \
           (W_CAMBER * (err_camber**2)) + \
           (W_LCA_ANGLE * (err_lca_angle**2)) + \
           (W_UCA_ANGLE * (err_uca_angle**2)) + \
           (W_PUSHROD_ANGLE * (err_pushrod_angle**2)) + \
           (W_PUSHROD_FORCE * (cost_force_pushrod**2)) + \
           (W_LCA_FORCE * (cost_force_lca**2)) + \
           (W_UCA_FORCE * (cost_force_uca**2)) + \
           (W_TR_FORCE * (cost_force_tr**2))
           
    return cost

print("Função 'objective_function' atualizada para 30 variáveis.")

# ---
# CÉLULA 5: CONFIGURAÇÃO DE LIMITES (ATUALIZADA)

# Lista de 30 tuplas (min, max) para os 10 pontos.
# !!! VOCÊ DEVE PREENCHER OS LIMITES PARA TR_In e TR_Out !!!
bounds = [
    # UBJ (x,y,z)
    (0, 0), (0, 0), (320, 320),
    # LBJ (x,y,z)
    (0, 0), (0, 0), (100, 100),
    # UCA_F (x,y,z)
    (81, 114), (244, 257), (320, 320),
    # UCA_R (x,y,z)
    (-114, -81), (244, 257), (320, 320),
    # LCA_F (x,y,z)
    (109, 160), (277, 300), (100,100),
    # LCA_R (x,y,z)
    (-160, -109), (277, 300), (100, 100),
    # PR_In (x,y,z)
    (0, 0), (0, 0), (0, 0),
    # PR_Out (x,y,z)
    (0, 0), (171, 321), (383, 470),
    
    # --- NOVOS LIMITES OBRIGATÓRIOS ---
    # TR_In (Ponto do chassi / caixa de direção)
    (70, 70), (700, 700), (210, 210),
    # TR_Out (Ponto na manga)
    (70, 70), (0, 0), (210, 210),
]

print(f"Limites (bounds) definidos. Número de variáveis: {len(bounds)}")

# ---
# CÉLULA 6: EXECUTAR A OTIMIZAÇÃO

print("Iniciando a otimização...")

# Use workers=1 para testar a lógica e ver o 'disp=True'
# Use workers=-1 DEPOIS que a simulação real (pesada) estiver implementada
result = differential_evolution(
    func=objective_function,
    bounds=bounds,
    maxiter=1000,      # Aumente para 1000+ para resultados reais
    popsize=60,       # Aumente com o nº de variáveis (30*2 = 60)
    disp=True,        # Mostra o progresso
    workers=1         # Use 1 para depurar e ver o progresso
)

print("\n... Otimização concluída.")

# ---
# CÉLULA 7: IMPRESSÃO DOS RESULTADOS FINAIS (ATUALIZADA)

if result.success:
    print("\n--- SUCESSO: Otimização convergiu! ---")
    print(f"Melhor Custo (Fitness) Final: {result.fun}")
    
    best_x = result.x
    
    # "Desempacota" os 30 pontos
    best_points = {
        'UBJ':    best_x[0:3],   'LBJ':    best_x[3:6],
        'UCA_F':  best_x[6:9],   'UCA_R':  best_x[9:12],
        'LCA_F':  best_x[12:15], 'LCA_R':  best_x[15:18],
        'PR_In':  best_x[18:21], 'PR_Out': best_x[21:24],
        'TR_In':  best_x[24:27], 'TR_Out': best_x[27:30]
    }
    
    print("\n--- Geometria Ótima Encontrada (Coordenadas) ---")
    for point, coords in best_points.items():
        print(f"  {point}: {np.round(coords, 2)}")
        
    print("\n--- Métricas Finais da Geometria Otimizada ---")
    final_metrics = run_suspension_simulation(best_points)
    for metric, value in final_metrics.items():
        print(f"  {metric}: {value:.3f}")
    
else:
    print("\n--- FALHA: A otimização não convergiu ---")
    print(f"Mensagem: {result.message}")