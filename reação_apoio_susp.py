# %matplotlib widget

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# -------------------------------
# Função de cálculo 3D
def calcular_forcas_duplo_wishbone_pushrod(
    Braço_sup_front, Braço_sup_tras, BJ_sup,
    Braço_inf_front, Braço_inf_tras, BJ_inf,
    Top_pushrod, F_total
):
    def vetor(p1, p2):
        v = np.array(p1) - np.array(p2)
        return v / np.linalg.norm(v)

    vet_inf_1 = vetor(Braço_inf_front, BJ_inf)
    vet_inf_2 = vetor(Braço_inf_tras, BJ_inf)
    vet_pushrod = vetor(Top_pushrod, BJ_inf)

    u_sup_1 = vetor(Braço_sup_front, BJ_sup)
    u_sup_2 = vetor(Braço_sup_tras, BJ_sup)

    A_inf = np.column_stack([vet_inf_1, vet_inf_2, vet_pushrod])
    b_inf = np.array(F_total)
    T_inf = np.linalg.solve(A_inf, b_inf)
    T_inf_1, T_inf_2, T_pushrod = T_inf

    h_contato_ate_bj_inferior = BJ_inf[2]     # Altura do centro de contato até o ball joint inferior [m]
    h_bj_inferior_ate_bj_superior = BJ_sup[2]-BJ_inf[2]

    h1 =float(h_contato_ate_bj_inferior)
    h2 =float(h_bj_inferior_ate_bj_superior)
    F_sup_parcial = F_total[1] * (h1 / h2)
    F_lat_primaria_inferior = F_total[1] * ((h1 + h2) / h2)

    F_sup = [F_sup_parcial, 0, 0]
    A_sup = np.column_stack([u_sup_1, u_sup_2])
    b_sup = -np.array(F_sup)
    T_sup, _, _, _ = np.linalg.lstsq(A_sup, b_sup, rcond=None)
    T_sup_1, T_sup_2 = T_sup

    return {
        "inferior": {
            "T_frontal": T_inf_1,
            "T_traseiro": T_inf_2,
            "T_pushrod": T_pushrod,
            "direcoes": [vet_inf_1, vet_inf_2, vet_pushrod]
        },
        "superior": {
            "T_frontal": T_sup_1,
            "T_traseiro": T_sup_2,
            "direcoes": [u_sup_1, u_sup_2]
        }
    }


# -------------------------------
# Geometria do sistema (exemplo)
Braço_sup_front = [-0.181,  0.000,  0.360]
Braço_sup_tras = [-0.285,  0.555,  0.270]
BJ_sup = [0.0, 0.2, 0.3]
# Braço_sup_front = [0.3, 0.5, 0.6]
# Braço_sup_tras = [-0.3, 0.5, 0.6]
# BJ_sup = [0.0, 0.2, 0.05]

Braço_inf_front = [-0.181,  0.000,  0.050]
Braço_inf_tras = [-0.285,  0.495,  0.045]
BJ_inf = [0.0, 0.2, 0.01]
# Braço_inf_front = [0.3, 0.5, 0.1]
# Braço_inf_tras = [-0.3, 0.5, 0.1]
# BJ_inf = [0.0, 0.2, 0.0]

px = (Braço_inf_front[0]+Braço_inf_tras[0])/2
py = (Braço_inf_front[1]+Braço_inf_tras[1])/2
alfa = 43   #angulo entre o pushrod e o plano do solo
alfa_rad=np.radians(alfa)
pz = np.sqrt(px**2+py**2)*np.tan(alfa_rad)
print(px, py, pz)
Top_pushrod = [px, py, pz ]
F_total = [4440,5000,0]


# -------------------------------
# Cálculo das forças
res = calcular_forcas_duplo_wishbone_pushrod(
    Braço_sup_front, Braço_sup_tras, BJ_sup,
    Braço_inf_front, Braço_inf_tras, BJ_inf,
    Top_pushrod, F_total
)

# -------------------------------
# Plotagem 3D com anotações de força
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Modelo 3D - Duplo Wishbone com Pushrod", fontsize=14)

def plot_barra_f(p1, p2, color='k', label='', force_val=None, lw=2, linestyle='-'):
    xs, ys, zs = zip(p1, p2)
    ax.plot(xs, ys, zs, color=color, linewidth=lw, linestyle=linestyle)
    if label:
        ax.text(*p2, label, fontsize=8)
    if force_val is not None:
        mid = [(a + b)/2 for a, b in zip(p1, p2)]
        ax.text(*mid, f"{force_val:.1f} N", fontsize=8, color=color)

# Forças calculadas
T_inf = res["inferior"]
T_sup = res["superior"]

print("\n🔹 Braço Inferior:")
for nome, valor in T_inf.items():
    if nome.startswith("T_"):
        print(f"{nome}: {valor:.2f} N")

print("\n🔹 Braço Superior:")
for nome, valor in T_sup.items():
    if nome.startswith("T_"):
        print(f"{nome}: {valor:.2f} N")
# Inferior
plot_barra_f(BJ_inf, Braço_inf_front, 'b', 'Inf-Frontal', T_inf["T_frontal"])
plot_barra_f(BJ_inf, Braço_inf_tras, 'b', 'Inf-Traseiro', T_inf["T_traseiro"])
plot_barra_f(BJ_inf, Top_pushrod, 'r', 'Pushrod', T_inf["T_pushrod"])

# Superior
plot_barra_f(BJ_sup, Braço_sup_front, 'g', 'Sup-Frontal', T_sup["T_frontal"])
plot_barra_f(BJ_sup, Braço_sup_tras, 'g', 'Sup-Traseiro', T_sup["T_traseiro"])

# Vetor da força da roda
escala = 0.0003
Fvec = np.array(BJ_inf) + escala * np.array(F_total)
plot_barra_f(BJ_inf, Fvec.tolist(), 'm', 'Força da roda', None, lw=3)

# Pontos de interesse
for pt, nome in zip(
    [Braço_sup_front, Braço_sup_tras, BJ_sup, Braço_inf_front, Braço_inf_tras, BJ_inf, Top_pushrod],
    ['Braço_sup_front', 'Braço_sup_tras', 'BJ_sup', 'Braço_inf_front', 'Braço_inf_tras', 'BJ_inf', 'Pushrod_top']
):
    ax.scatter(*pt, s=30)
    ax.text(*pt, nome, fontsize=0.5)

# Ajustes visuais
ax.set_xlabel("X [m]")
ax.set_ylabel("Y [m]")
ax.set_zlabel("Z [m]")
ax.set_xlim([-0.5, 0.5])
ax.set_ylim([0, 0.7])
ax.set_zlim([0, 0.7])
ax.view_init(elev=20, azim=130)

plt.tight_layout()
plt.show()
