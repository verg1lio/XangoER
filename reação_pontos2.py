import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# --- Funções de Cálculo de Forças ---

def calcular_vetor_unitario(p1, p2):
    """Calcula o vetor unitário de p2 para p1."""
    v = np.array(p1) - np.array(p2)
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

def calcular_forcas_carga_vertical(geometria, W_vert_N):
    """
    Calcula as forças nos braços de suspensão para o caso de máxima carga vertical.
    """
    BJ_inf = geometria['BJ_inf']
    Top_pushrod = geometria['Top_pushrod']
    Braço_inf_front = geometria['Braço_inf_front']
    Braço_inf_tras = geometria['Braço_inf_tras']
    


    F_vertical = np.array([0, 0, W_vert_N])

    vet_inf_1 = calcular_vetor_unitario(Braço_inf_front, BJ_inf)
    vet_inf_2 = calcular_vetor_unitario(Braço_inf_tras, BJ_inf)
    vet_pushrod = calcular_vetor_unitario(Top_pushrod, BJ_inf)
    

    A = np.column_stack([vet_inf_1, vet_inf_2, vet_pushrod])
    b = -F_vertical
    
    try:
        forcas = np.linalg.solve(A, b)
        T_inf_1, T_inf_2, T_pushrod = forcas
    except np.linalg.LinAlgError:
        print("Aviso: Sistema singular para carga vertical. Verifique a geometria.")
        return None

    return {
        "Braço Inferior Frontal (N)": T_inf_1,
        "Braço Inferior Traseiro (N)": T_inf_2,
        "Pushrod (N)": T_pushrod,
        "Braço Superior Frontal (N)": 0.0,
        "Braço Superior Traseiro (N)": 0.0
    }

def calcular_forcas_frenagem(geometria, W_long_N, W_vert_N):
    """
    Calcula as forças para o caso de máxima frenagem.
    """
    BJ_inf = np.array(geometria['BJ_inf'])
    BJ_sup = np.array(geometria['BJ_sup'])
    Braço_inf_front = geometria['Braço_inf_front']
    Braço_inf_tras = geometria['Braço_inf_tras']
    Braço_sup_front = geometria['Braço_sup_front']
    Braço_sup_tras = geometria['Braço_sup_tras']
    
    forcas_vert = calcular_forcas_carga_vertical(geometria, W_vert_N)
    T_inf_1_vert = forcas_vert["Braço Inferior Frontal (N)"]
    T_inf_2_vert = forcas_vert["Braço Inferior Traseiro (N)"]
    T_pushrod = forcas_vert["Pushrod (N)"]

    h1 = BJ_inf[2]
    h2 = BJ_sup[2] - BJ_inf[2]

    F_top_long = -W_long_N * (h1 / h2)
    F_bottom_long = -W_long_N * ((h1 + h2) / h2)
    print(F_top_long)
    print(F_bottom_long)

    vet_sup_1 = calcular_vetor_unitario(Braço_sup_front, BJ_sup)
    vet_sup_2 = calcular_vetor_unitario(Braço_sup_tras, BJ_sup)
    A_sup = np.column_stack([vet_sup_1, vet_sup_2])
    b_sup = np.array([0, F_top_long, 0])
    T_sup = np.linalg.lstsq(A_sup, b_sup, rcond=None)[0]

    vet_inf_1 = calcular_vetor_unitario(Braço_inf_front, BJ_inf)
    vet_inf_2 = calcular_vetor_unitario(Braço_inf_tras, BJ_inf)
    A_inf = np.column_stack([vet_inf_1, vet_inf_2])
    b_inf = np.array([0, F_bottom_long, 0])
    T_inf_long = np.linalg.lstsq(A_inf, b_inf, rcond=None)[0]
    
    T_inf_1_total = T_inf_1_vert + T_inf_long[0]
    T_inf_2_total = T_inf_2_vert + T_inf_long[1]
    
    return {
        "Braço Inferior Frontal (N)": T_inf_1_total,
        "Braço Inferior Traseiro (N)": T_inf_2_total,
        "Pushrod (N)": T_pushrod,
        "Braço Superior Frontal (N)": T_sup[0],
        "Braço Superior Traseiro (N)": T_sup[1]
    }

def calcular_forcas_curva(geometria, W_lat_N, W_vert_N):
    """
    Calcula as forças para o caso de máxima força em curva.
    """
    BJ_inf = np.array(geometria['BJ_inf'])
    BJ_sup = np.array(geometria['BJ_sup'])
    Braço_inf_front = geometria['Braço_inf_front']
    Braço_inf_tras = geometria['Braço_inf_tras']
    Braço_sup_front = geometria['Braço_sup_front']
    Braço_sup_tras = geometria['Braço_sup_tras']

    forcas_vert = calcular_forcas_carga_vertical(geometria, W_vert_N)
    T_inf_1_vert = forcas_vert["Braço Inferior Frontal (N)"]
    T_inf_2_vert = forcas_vert["Braço Inferior Traseiro (N)"]
    T_pushrod = forcas_vert["Pushrod (N)"]
    
    h1 = BJ_inf[2]
    h2 = BJ_sup[2] - BJ_inf[2]

    F_top_lat = -W_lat_N * (h1 / h2)
    F_bottom_lat = -W_lat_N * ((h1 + h2) / h2)

    vet_sup_1 = calcular_vetor_unitario(Braço_sup_front, BJ_sup)
    vet_sup_2 = calcular_vetor_unitario(Braço_sup_tras, BJ_sup)
    A_sup = np.column_stack([vet_sup_1, vet_sup_2])
    b_sup = np.array([F_top_lat, 0, 0])
    T_sup = np.linalg.lstsq(A_sup, b_sup, rcond=None)[0]

    vet_inf_1 = calcular_vetor_unitario(Braço_inf_front, BJ_inf)
    vet_inf_2 = calcular_vetor_unitario(Braço_inf_tras, BJ_inf)
    A_inf = np.column_stack([vet_inf_1, vet_inf_2])
    b_inf = np.array([F_bottom_lat, 0, 0])
    T_inf_lat = np.linalg.lstsq(A_inf, b_inf, rcond=None)[0]

    T_inf_1_total = T_inf_1_vert + T_inf_lat[0]
    T_inf_2_total = T_inf_2_vert + T_inf_lat[1]
    
    return {
        "Braço Inferior Frontal (N)": T_inf_1_total,
        "Braço Inferior Traseiro (N)": T_inf_2_total,
        "Pushrod (N)": T_pushrod,
        "Braço Superior Frontal (N)": T_sup[0],
        "Braço Superior Traseiro (N)": T_sup[1]
    }

def calcular_forcas_aceleracao(geometria, W_long_N, W_vert_N):
    """
    Calcula as forças para o caso de máxima aceleração.
    """
    BJ_inf = np.array(geometria['BJ_inf'])
    BJ_sup = np.array(geometria['BJ_sup'])
    Braço_inf_front = geometria['Braço_inf_front']
    Braço_inf_tras = geometria['Braço_inf_tras']
    Braço_sup_front = geometria['Braço_sup_front']
    Braço_sup_tras = geometria['Braço_sup_tras']

    forcas_vert = calcular_forcas_carga_vertical(geometria, W_vert_N)
    T_inf_1_vert = forcas_vert["Braço Inferior Frontal (N)"]
    T_inf_2_vert = forcas_vert["Braço Inferior Traseiro (N)"]
    T_pushrod = forcas_vert["Pushrod (N)"]

    F_long_cada_braco = -W_long_N / 2.0

    vet_sup_1 = calcular_vetor_unitario(Braço_sup_front, BJ_sup)
    u_sup_2 = calcular_vetor_unitario(Braço_sup_tras, BJ_sup)
    A_sup = np.column_stack([vet_sup_1, u_sup_2])
    b_sup = np.array([0, F_long_cada_braco, 0])
    T_sup = np.linalg.lstsq(A_sup, b_sup, rcond=None)[0]

    vet_inf_1 = calcular_vetor_unitario(Braço_inf_front, BJ_inf)
    u_inf_2 = calcular_vetor_unitario(Braço_inf_tras, BJ_inf)
    A_inf = np.column_stack([vet_inf_1, u_inf_2])
    b_inf = np.array([0, F_long_cada_braco, 0])
    T_inf_long = np.linalg.lstsq(A_inf, b_inf, rcond=None)[0]

    T_inf_1_total = T_inf_1_vert + T_inf_long[0]
    T_inf_2_total = T_inf_2_vert + T_inf_long[1]
    
    return {
        "Braço Inferior Frontal (N)": T_inf_1_total,
        "Braço Inferior Traseiro (N)": T_inf_2_total,
        "Pushrod (N)": T_pushrod,
        "Braço Superior Frontal (N)": T_sup[0],
        "Braço Superior Traseiro (N)": T_sup[1]
    }



# Define os pontos fixos da geometria
geometria_suspensao = {
    'Braço_sup_front': [-0.302,  1.670,  0.240],
    'Braço_sup_tras': [-0.293,  1.950,  0.250],
    'BJ_sup': [0.0, 1.75, 0.25],
    'Braço_inf_front': [-0.271,  1.665,  0.030],
    'Braço_inf_tras': [-0.271,  1.835,  0.030],
    'BJ_inf': [0.0, 1.75, 0.05]
}

p_inf_front = geometria_suspensao['Braço_inf_front']
p_inf_tras = geometria_suspensao['Braço_inf_tras']
BJ_inf = geometria_suspensao['BJ_inf']

px = (p_inf_front[0] + p_inf_tras[0]) / 2.0
py = (p_inf_front[1] + p_inf_tras[1]) / 2.0
comp_pushrod = 0.400 #tamanho do pushrod em m
pz = np.sqrt(comp_pushrod**2-(px-BJ_inf[0])**2-(py-BJ_inf[1])**2)+BJ_inf[2]


print(f"PONTO Z {px};{py};{pz}")

Top_pushrod = [px, py, pz]

# Adiciona a coordenada calculada ao dicionário de geometria
geometria_suspensao['Top_pushrod'] = Top_pushrod

# Cargas para cada cenário 
max_carga_vertical = 4500 # N
max_frenagem_long = 6300  # N
max_frenagem_vert = 4500  # N
max_curva_lat = 17712       # N
max_curva_vert = 14760      # N
max_aceleracao_long = 5124  # N
max_aceleracao_vert = 3416  # N

# --- Cálculo e Impressão dos Resultados ---

print("--- ANÁLISE DE FORÇAS NA SUSPENSÃO ---\n")

print("CASO 1: Máxima Carga Vertical")
forcas_vert = calcular_forcas_carga_vertical(geometria_suspensao, max_carga_vertical)
for nome, valor in forcas_vert.items():
    print(f"  {nome:<30}: {valor:8.2f} N")
print("-" * 45)

print("\nCASO 2: Máxima Frenagem")
forcas_fren = calcular_forcas_frenagem(geometria_suspensao, max_frenagem_long, max_frenagem_vert)
for nome, valor in forcas_fren.items():
    print(f"  {nome:<30}: {valor:8.2f} N")
print("-" * 45)

print("\nCASO 3: Máxima Força em Curva")
forcas_curva = calcular_forcas_curva(geometria_suspensao, max_curva_lat, max_curva_vert)
for nome, valor in forcas_curva.items():
    print(f"  {nome:<30}: {valor:8.2f} N")
print("-" * 45)

print("\nCASO 4: Máxima Aceleração")
forcas_acel = calcular_forcas_aceleracao(geometria_suspensao, max_aceleracao_long, max_aceleracao_vert)
for nome, valor in forcas_acel.items():
    print(f"  {nome:<30}: {valor:8.2f} N")
print("-" * 45)

# --- Plotagem 3D da Geometria (Visualização) ---

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Geometria da Suspensão - Duplo Wishbone com Pushrod", fontsize=14)

def plot_barra(p1, p2, cor='k', label=''):
    xs, ys, zs = zip(p1, p2)
    ax.plot(xs, ys, zs, color=cor)
    if label:
        ax.text(*p2, label, fontsize=9)

# Plotar a geometria
plot_barra(geometria_suspensao['BJ_inf'], geometria_suspensao['Braço_inf_front'], 'b')
plot_barra(geometria_suspensao['BJ_inf'], geometria_suspensao['Braço_inf_tras'], 'b')
plot_barra(geometria_suspensao['BJ_inf'], geometria_suspensao['Top_pushrod'], 'r')
plot_barra(geometria_suspensao['BJ_sup'], geometria_suspensao['Braço_sup_front'], 'g')
plot_barra(geometria_suspensao['BJ_sup'], geometria_suspensao['Braço_sup_tras'], 'g')
plot_barra(geometria_suspensao['BJ_sup'], geometria_suspensao['BJ_inf'], 'grey', 'Upright') # Linha do upright

# Plotar os pontos
for nome, pt in geometria_suspensao.items():
    ax.scatter(*pt, s=40, label=nome)
ax.legend()

ax.set_xlabel("Eixo X (m)")
ax.set_ylabel("Eixo Y (m)")
ax.set_zlabel("Eixo Z (m)")
ax.set_box_aspect([2, 2, 2]) # Aspecto igual para melhor visualização
ax.view_init(elev=20, azim=130)

plt.tight_layout()
plt.show()