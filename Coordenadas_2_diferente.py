import numpy as np

def energia_deformacao(epsilon, sigma, E, volume):
    return np.trapz(np.trapz((E * epsilon ** 2) / 2, axis=-1), axis=-1) * volume

def energia_potencial_forcas_externas(forcas, posicoes):
    return np.sum(forcas * posicoes)

def energia_potencial_total(U, W):
    return -U - W

# Definindo os parâmetros
epsilon = np.array([[0.1, 0.2], [0.3, 0.4]])  # Deformação de engenharia
sigma = np.array([[1.0, 2.0], [3.0, 4.0]])  # Tensão normal
E = 10.0  # Módulo de elasticidade
volume = 2.0  # Volume do elemento
forcas = np.array([1, 2])  # Forças aplicadas
posicoes = np.array([0.1, 0.2])  # Posições independentes

# Calculando as energias
U = energia_deformacao(epsilon, sigma, E, volume)
W = energia_potencial_forcas_externas(forcas, posicoes)
Pi = energia_potencial_total(U, W)

# Imprimindo os resultados
print("Energia de deformação:", U)
print("Energia potencial das forças externas:", W)
print("Energia potencial total:", Pi)
