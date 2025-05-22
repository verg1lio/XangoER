import numpy as np

def calcular_forcas_bracos(Fx, Fy, theta_s_deg, theta_i_deg):

    theta_s = np.radians(theta_s_deg)
    theta_i = np.radians(theta_i_deg)
#-------
    A = np.array([
        [np.cos(theta_s), np.cos(theta_i)],
        [np.sin(theta_s), np.sin(theta_i)]
    ])
    b = np.array([Fx, Fy])

#------------
    T = np.linalg.solve(A, b)
    Ts, Ti = T
    return Ts, Ti
#------------
def reacoes_em_apoios(T, L, a, theta_deg):
    theta = np.radians(theta_deg)
    s = np.array([np.cos(theta), np.sin(theta)])

    
    RA = T * (L - a) / L
    RB = T * a / L

    RA_vec = RA * s
    RB_vec = RB * s
    return RA_vec, RB_vec

#------------
Fx = 600     # força lateral [N]
Fy = -300    # força vertical [N] (peso para baixo)
L = 0.3      # comprimento do braço [m]
a = 0.15     # ponto médio [m]

theta_s = 30  # braço superior [°]
theta_i = 10  # braço inferior [°]

# --------------------------
Ts, Ti = calcular_forcas_bracos(Fx, Fy, theta_s, theta_i)
print(f"Força axial no braço superior (Ts): {Ts:.2f} N")
print(f"Força axial no braço inferior (Ti): {Ti:.2f} N")

# --------------------------
RA_s, RB_s = reacoes_em_apoios(Ts, L, a, theta_s)  # braço superior
RA_i, RB_i = reacoes_em_apoios(Ti, L, a, theta_i)  # braço inferior

# --------------------------
print("--- Reações nos pontos de fixação ---")
print(f"Ponto A (braço sup. interno): {RA_s.round(2)} N")
print(f"Ponto B (braço sup. externo): {RB_s.round(2)} N")
print(f"Ponto C (braço inf. interno): {RA_i.round(2)} N")
print(f"Ponto D (braço inf. externo): {RB_i.round(2)} N")

"""
O código calcula as forças nos braços de suspensão (superior e inferior) 
e as reações nos pontos de fixação com o chassi.

Na função calcular_forcas_bracos, usamos os ângulos dos braços para montar
um sistema de equações baseado no equilíbrio das forças lateral (Fx) e vertical (Fy).
A solução do sistema fornece as forças axiais: Tu (braço superior) e Tl (braço inferior),
que indicam se os braços estão em tração ou compressão.

Depois, na função reacoes_em_apoios, essas forças axiais são distribuídas entre os dois 
pontos de fixação de cada braço (A e B no superior, C e D no inferior). A divisão é 
proporcional à posição da carga ao longo do braço, garantindo o equilíbrio de forças e momentos.
As reações são retornadas como vetores com direção e intensidade.

Na parte principal do código, definimos as forças e geometrias do sistema, calculamos Tu e Tl,
e depois obtemos as reações nos apoios. O resultado final mostra as forças que cada ponto de 
ligação transmite ao chassi, o que é essencial para o dimensionamento mecânico da suspensão.
"""