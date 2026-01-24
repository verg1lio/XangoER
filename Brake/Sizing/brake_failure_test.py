import numpy as np
import matplotlib.pyplot as plt

# =====================
# FUNÇÕES DE FALHA
# =====================

def falha_tresca(s1, s2, s3, Sy):
    return max(abs(s1 - s2), abs(s2 - s3), abs(s3 - s1)) >= Sy


def falha_von_mises(s1, s2, s3, Sy):
    sigma_vm = np.sqrt(0.5*((s1-s2)**2 + (s2-s3)**2 + (s3-s1)**2))
    return sigma_vm >= Sy


# =====================
# CONTORNOS
# =====================

def contorno_tresca(Sy):
    v = [
        ( Sy,  0),
        ( Sy, Sy),
        ( 0,  Sy),
        (-Sy/2,  Sy/2),
        (-Sy,  0),
        (-Sy, -Sy),
        ( 0, -Sy),
        ( Sy/2, -Sy/2),
        ( Sy,  0),
    ]
    s1 = np.array([p[0] for p in v])
    s2 = np.array([p[1] for p in v])
    return s1, s2


def contorno_von_mises(Sy, n=800):
    s1 = np.linspace(-1.2*Sy, 1.2*Sy, n)

    a = 1
    b = s1
    c = s1**2 - Sy**2

    disc = b**2 - 4*a*c
    mask = disc >= 0

    s1v = s1[mask]
    d = np.sqrt(disc[mask])

    s2p = (b[mask] + d) / 2
    s2m = (b[mask] - d) / 2

    return (
        np.concatenate([s1v, s1v[::-1]]),
        np.concatenate([s2p, s2m[::-1]])
    )


# =====================
# MATERIAIS
# =====================

materiais = [
    {"nome": "Aço SAE 1020", "Sy": 250e6},
    {"nome": "Aço Inox 304", "Sy": 215e6},
    {"nome": "Aço 4140", "Sy": 415e6},
]

# =====================
# ESTADOS DE TENSÃO
# =====================

estados_tensao = [
    {"label": "Ponto A (freio leve)", "s": (200e6, 50e6, 0)},
    {"label": "Ponto B (freio severo)", "s": (300e6, 0e6, 0)},
    {"label": "Ponto C (compressão radial)", "s": (150e6, -100e6, 0)},
]

# =====================
# LOOP POR MATERIAL
# =====================

for mat in materiais:

    Sy = mat["Sy"]
    nome = mat["nome"]

    plt.figure(figsize=(6,6))

    # Contornos
    s1_t, s2_t = contorno_tresca(Sy)
    s1_vm, s2_vm = contorno_von_mises(Sy)

    plt.plot(s1_t/1e6, s2_t/1e6, label="Tresca", linewidth=2)
    plt.plot(s1_vm/1e6, s2_vm/1e6, '--', label="Von Mises", linewidth=2)

    # Pontos
    for p in estados_tensao:
        s1, s2, s3 = p["s"]

        ft = falha_tresca(s1, s2, s3, Sy)
        fv = falha_von_mises(s1, s2, s3, Sy)

        cor = 'green'
        if ft and fv:
            cor = 'red'
        elif ft or fv:
            cor = 'orange'

        plt.plot(s1/1e6, s2/1e6, 'o', color=cor)
        plt.text(
            s1/1e6*1.02,
            s2/1e6*1.02,
            f"{p['label']}\nT:{'F' if ft else 'OK'} V:{'F' if fv else 'OK'}",
            fontsize=8
        )

    plt.axhline(0, color='black', linewidth=0.8)
    plt.axvline(0, color='black', linewidth=0.8)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.gca().set_aspect('equal', 'box')

    plt.xlabel("σ₁ [MPa]")
    plt.ylabel("σ₂ [MPa]")
    plt.title(f"Critérios de Escoamento — σ₃ = 0\n{nome}")
    plt.legend()
    plt.show()
