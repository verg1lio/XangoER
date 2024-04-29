import numpy as np

# Definição das coordenadas dos nós
# (substitua os valores conforme necessário)
coordenadas = {
    1: [0, 0, 0],  # Nó 1: coordenadas x, y e z
    2: [1, 0, 0],  # Nó 2: coordenadas x, y e z
    3: [1, 1, 0],  # Nó 3: coordenadas x, y e z
    4: [0, 1, 0]   # Nó 4: coordenadas x, y e z
}

# Função para calcular a matriz de rigidez do elemento no sistema local
def calcular_matriz_rigidez_local(AE, L):
    k = AE / L * np.array([[1, -1],
                           [-1, 1]])
    return k

# Função para calcular a matriz de transformação
def calcular_matriz_transformacao(theta):
    T = np.array([[np.cos(theta), np.sin(theta)],
                  [-np.sin(theta), np.cos(theta)]])
    return T

# Função para calcular a matriz de rigidez do elemento no sistema global
def calcular_matriz_rigidez_global(K_local, T):
    K_global = np.dot(np.dot(T.T, K_local), T)
    return K_global

# AE e L para o elemento (substitua os valores conforme necessário)
AE = 1
L = 1

# Calcular a matriz de rigidez local
K_local = calcular_matriz_rigidez_local(AE, L)

# Calcular a matriz de transformação (supondo theta = 0 para simplificar)
theta = 0
T = calcular_matriz_transformacao(theta)

# Calcular a matriz de rigidez global
K_global = calcular_matriz_rigidez_global(K_local, T)

print("Matriz de rigidez global:")
print(K_global)
