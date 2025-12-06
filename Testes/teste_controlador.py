import sympy as sp
import numpy as np
import control as ctrl
import matplotlib.pyplot as plt

# ==============================================================================
# 1. DEFINIÇÃO SIMBÓLICA DO SISTEMA (Símbolos e Equações)
# ==============================================================================
print("=== 1. DEFINIÇÃO DAS EQUAÇÕES DIFERENCIAIS DO MOTOR ===")

# Definir variáveis no domínio do tempo
t = sp.symbols('t')
isd, isq = sp.symbols('isd isq')    # Estados (Correntes)
vd, vq = sp.symbols('vd vq')        # Entradas (Tensões)

# Definir parâmetros físicos
Rs, Ld, Lq = sp.symbols('Rs Ld Lq')      # Resistência e Indutâncias
we = sp.symbols('we')                    # Velocidade Elétrica (rad/s)
lambda_m = sp.symbols('lambda_m')        # Fluxo Magnético

# --- AS EDOS DO MOTOR (Baseado no seu Motor.py) ---
# d_isd/dt = (vd - Rs*isd + we*Lq*isq) / Ld
f_isd = (vd - Rs * isd + we * Lq * isq) / Ld

# d_isq/dt = (vq - Rs*isq - we*(Ld*isd + lambda_m)) / Lq
f_isq = (vq - Rs * isq - we * (Ld * isd + lambda_m)) / Lq

print("Equação f_isd (dId/dt):")
sp.pprint(f_isd)
print("\nEquação f_isq (dIq/dt):")
sp.pprint(f_isq)

# ==============================================================================
# 2. CÁLCULO DAS MATRIZES A e B (JACOBIANAS)
# ==============================================================================
print("\n=== 2. CÁLCULO DAS MATRIZES A e B ===")

# Vetores de Estado e Entrada
X = sp.Matrix([isd, isq])
U = sp.Matrix([vd, vq])
F = sp.Matrix([f_isd, f_isq])  # Vetor das derivadas

# Calcular Jacobianas
A_sym = F.jacobian(X)
B_sym = F.jacobian(U)

print("Matriz A Simbólica (Dinâmica):")
sp.pprint(A_sym)

print("\nMatriz B Simbólica (Entrada):")
sp.pprint(B_sym)

# ==============================================================================
# 3. CÁLCULO DA FUNÇÃO DE TRANSFERÊNCIA (SIMBÓLICA)
# ==============================================================================
print("\n=== 3. OBTENÇÃO DA FUNÇÃO DE TRANSFERÊNCIA G(s) ===")
print("Fórmula: G(s) = C * (sI - A)^-1 * B")

s = sp.symbols('s')
C = sp.eye(2)  # Matriz Identidade (Saída = Estados)
I = sp.eye(2)

# Inversa de (sI - A)
sI_A = s * I - A_sym
# Simplificação: Para sintonia de PID, assumimos we=0 (desacoplamento perfeito)
# Se we != 0, a FT fica acoplada e complexa.
sI_A_decoupled = sI_A.subs(we, 0) 
inv_sIA = sI_A_decoupled.inv()

# Cálculo da FT (Matriz 2x2)
# G(s) = [ G_id_vd   G_id_vq ]
#        [ G_iq_vd   G_iq_vq ]
G_sym = C * inv_sIA * B_sym

print("\nMatriz de Funções de Transferência (Desacoplada/Sintonia):")
sp.pprint(G_sym)

# ==============================================================================
# 4. SUBSTITUIÇÃO NUMÉRICA (SEUS DADOS)
# ==============================================================================
print("\n=== 4. FUNÇÕES DE TRANSFERÊNCIA NUMÉRICAS ===")

# Parâmetros do seu projeto (Motor.py / main.py)
# Nota: Lq = 0.00172 (conforme Motor.py instanciado no build_defaults)
valores = {
    Rs: 0.0732,
    Ld: 0.000078,
    Lq: 0.000078,
    we: 0.0  # Assumindo motor parado ou desacoplado para sintonia
}

# Substituir valores na Matriz A e B
A_num = np.array(A_sym.subs(valores)).astype(np.float64)
B_num = np.array(B_sym.subs(valores)).astype(np.float64)
C_num = np.array(C).astype(np.float64)
D_num = np.array([[0, 0], [0, 0]]).astype(np.float64)

print("Matriz A Numérica:")
print(A_num)
print("\nMatriz B Numérica:")
print(B_num)

# Criar sistema no Python Control
sys_ss = ctrl.StateSpace(A_num, B_num, C_num, D_num)

# --- CORREÇÃO DO ERRO SLYCOT ---
# Em vez de converter tudo de uma vez (MIMO), extraímos os pares SISO primeiro.
# O comando sys_ss[saida, entrada] cria um subsistema.

# Eixo D: Saída 0 (Id) / Entrada 0 (Vd)
sys_d_ss = sys_ss[0, 0]
tf_d = ctrl.ss2tf(sys_d_ss)

# Eixo Q: Saída 1 (Iq) / Entrada 1 (Vq)
sys_q_ss = sys_ss[1, 1]
tf_q = ctrl.ss2tf(sys_q_ss)

print("\nFunção de Transferência Eixo D (Id / Vd):")
print(tf_d)

print("\nFunção de Transferência Eixo Q (Iq / Vq):")
print(tf_q)

# ==============================================================================
# 5. VERIFICAÇÃO VISUAL (BODE)
# ==============================================================================
plt.figure(figsize=(10, 5))
# Usamos tf_q diretamente, que já é uma TransferFunction válida
ctrl.bode_plot(tf_q, dB=True, Hz=True, title="Diagrama de Bode - Eixo Q (Torque)")
plt.tight_layout()
plt.show()

# ==============================================================================
# 6. PLOTS DE RESPOSTA NO TEMPO
# ==============================================================================
print("\n=== 6. PLOTS DE RESPOSTA NO TEMPO ===")

# Definir tempo de simulação
t_sim = np.linspace(0, 0.1, 1000)  # 100ms

# Resposta ao degrau para ambos os eixos
t_d, y_d = ctrl.step_response(tf_d, t_sim)
t_q, y_q = ctrl.step_response(tf_q, t_sim)

# Plot das respostas ao degrau
plt.figure(figsize=(12, 8))

# Resposta Eixo D
plt.subplot(2, 2, 1)
plt.plot(t_d, y_d, 'b-', linewidth=2)
plt.grid(True, alpha=0.3)
plt.xlabel('Tempo (s)')
plt.ylabel('Corrente Id (A)')
plt.title('Resposta ao Degrau - Eixo D (Id/Vd)')
plt.xlim([0, 0.1])

# Resposta Eixo Q
plt.subplot(2, 2, 2)
plt.plot(t_q, y_q, 'r-', linewidth=2)
plt.grid(True, alpha=0.3)
plt.xlabel('Tempo (s)')
plt.ylabel('Corrente Iq (A)')
plt.title('Resposta ao Degrau - Eixo Q (Iq/Vq)')
plt.xlim([0, 0.1])

# Resposta ao impulso
t_imp_d, y_imp_d = ctrl.impulse_response(tf_d, t_sim)
t_imp_q, y_imp_q = ctrl.impulse_response(tf_q, t_sim)

# Resposta ao Impulso Eixo D
plt.subplot(2, 2, 3)
plt.plot(t_imp_d, y_imp_d, 'g-', linewidth=2)
plt.grid(True, alpha=0.3)
plt.xlabel('Tempo (s)')
plt.ylabel('Corrente Id (A)')
plt.title('Resposta ao Impulso - Eixo D (Id/Vd)')
plt.xlim([0, 0.1])

# Resposta ao Impulso Eixo Q
plt.subplot(2, 2, 4)
plt.plot(t_imp_q, y_imp_q, 'm-', linewidth=2)
plt.grid(True, alpha=0.3)
plt.xlabel('Tempo (s)')
plt.ylabel('Corrente Iq (A)')
plt.title('Resposta ao Impulso - Eixo Q (Iq/Vq)')
plt.xlim([0, 0.1])

plt.tight_layout()
plt.show()


# ==============================================================================
# 7. ANÁLISE DE POLOS E ZEROS
# ==============================================================================
print("\n=== 8. ANÁLISE DE POLOS E ZEROS ===")

# Polos e zeros do sistema
polos_d = tf_d.poles()
zeros_d = tf_d.zeros()

polos_q = tf_q.poles()
zeros_q = tf_q.zeros()

print(f"Polos Eixo D: {polos_d}")
print(f"Zeros Eixo D: {zeros_d}")
print(f"Polos Eixo Q: {polos_q}")
print(f"Zeros Eixo Q: {zeros_q}")

# Plot polos e zeros
plt.figure(figsize=(10, 4))

plt.subplot(1, 2, 1)
ctrl.pzmap(tf_d, plot=True, grid=True, title='Mapa de Polos e Zeros - Eixo D')
plt.subplot(1, 2, 2)
ctrl.pzmap(tf_q, plot=True, grid=True, title='Mapa de Polos e Zeros - Eixo Q')

plt.tight_layout()
plt.show()

# ==============================================================================
# 9. CÁLCULO DE CARACTERÍSTICAS TEMPORAIS
# ==============================================================================
print("\n=== 9. CARACTERÍSTICAS TEMPORAIS ===")

# Informações da resposta ao degrau
info_d = ctrl.step_info(tf_d)
info_q = ctrl.step_info(tf_q)

print("\nCaracterísticas Eixo D:")
print(f"  Tempo de Subida (Rise Time): {info_d['RiseTime']:.6f} s")
print(f"  Tempo de Pico: {info_d['PeakTime']:.6f} s")
print(f"  Sobressinal: {info_d['Overshoot']:.2f} %")
print(f"  Tempo de Acomodação: {info_d['SettlingTime']:.6f} s")

print("\nCaracterísticas Eixo Q:")
print(f"  Tempo de Subida (Rise Time): {info_q['RiseTime']:.6f} s")
print(f"  Tempo de Pico: {info_q['PeakTime']:.6f} s")
print(f"  Sobressinal: {info_q['Overshoot']:.2f} %")
print(f"  Tempo de Acomodação: {info_q['SettlingTime']:.6f} s")

# Constante de tempo teórica (para sistema de primeira ordem)
tau_d = Ld.subs(valores) / Rs.subs(valores)
tau_q = Lq.subs(valores) / Rs.subs(valores)
print(f"\nConstante de tempo teórica Eixo D: {tau_d:.6f} s")
print(f"Constante de tempo teórica Eixo Q: {tau_q:.6f} s")