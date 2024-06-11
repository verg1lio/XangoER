import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class BodyAndFrameFEM:
    def __init__(self, A, E, L, I, G, a, B, K, y, rho, k, g, J, num_elements):
        self.A = A    # Área da seção transversal
        self.E = E    # Módulo de elasticidade longitudinal
        self.L = L    # Comprimento da barra
        self.I = I    # Momento de inércia
        self.G = G    # Módulo de elasticidade transversal
        self.a = a    # Largura da seção transversal
        self.B = B    # Altura da seção transversal
        self.K = K    # Rigidez da mola
        self.y = y    # Deslocamento
        self.rho = rho # Massa específica do material
        self.k = k    # Coeficiente de amortecimento
        self.g = g    # Aceleração devido à gravidade
        self.J = J    # Momento polar de inércia
        self.num_elements = num_elements # Número de elementos finitos

    def analyze(self, F):
        deformacao_axial = self.deformacao_axial(F)
        flexao = self.flexao(F)
        torsao = self.torsao(F)

        return deformacao_axial, flexao, torsao

    def deformacao_axial(self, F):
        deformacao = (F / (self.A * self.E)) * (self.L / self.num_elements)
        return np.ones(self.num_elements + 1) * deformacao

    def flexao(self, F):
        num_nodes = self.num_elements + 1
        L_e = self.L / self.num_elements
        
        # Montagem da matriz de rigidez global
        K_global = lil_matrix((2 * num_nodes, 2 * num_nodes))
        for i in range(self.num_elements):
            K_e = (self.E * self.I / L_e**3) * np.array([
                [12, 6*L_e, -12, 6*L_e],
                [6*L_e, 4*L_e**2, -6*L_e, 2*L_e**2],
                [-12, -6*L_e, 12, -6*L_e],
                [6*L_e, 2*L_e**2, -6*L_e, 4*L_e**2]
            ])
            dofs = [2*i, 2*i+1, 2*i+2, 2*i+3]
            for ii in range(4):
                for jj in range(4):
                    K_global[dofs[ii], dofs[jj]] += K_e[ii, jj]
        
        # Aplicação das condições de contorno (engaste nas extremidades)
        K_global[0, :] = 0
        K_global[:, 0] = 0
        K_global[1, :] = 0
        K_global[:, 1] = 0
        K_global[0, 0] = 1
        K_global[1, 1] = 1

        K_global[-2, :] = 0
        K_global[:, -2] = 0
        K_global[-1, :] = 0
        K_global[:, -1] = 0
        K_global[-2, -2] = 1
        K_global[-1, -1] = 1

        # Montagem do vetor de forças global
        F_global = np.zeros(2 * num_nodes)
        for i in range(len(F)):
            F_global[2 * i] += F[i]

        # Ajustando as condições de contorno para permitir rotação nas extremidades
        K_global[0, 0] = 1
        K_global[1, 1] = 1e10  # Simula uma rigidez muito alta, mas permite rotação
        K_global[-2, -2] = 1e10  # Simula uma rigidez muito alta, mas permite rotação
        K_global[-1, -1] = 1

        # Resolução do sistema linear
        v = spsolve(K_global.tocsr(), F_global)

        return v[::2]  # Apenas deslocamentos verticais

    def torsao(self, F):
        torsao = (F * self.L) / (self.G * self.J)
        return np.ones(self.num_elements + 1) * torsao

    def plot_results(self, F):
        deformacao_axial, flexao, torsao = self.analyze(F)

        x = np.linspace(0, self.L, self.num_elements + 1)

        plt.figure(figsize=(18, 6))

        # Plot deformação axial
        plt.subplot(1, 3, 1)
        plt.plot(x, deformacao_axial, label='Deformação Axial')
        plt.xlabel('Posição ao longo da barra (m)')
        plt.ylabel('Deformação Axial')
        plt.title('Deformação Axial ao longo da barra')
        plt.legend()

        # Plot flexão
        plt.subplot(1, 3, 2)
        plt.plot(x, flexao, label='Flexão', color='orange')
        plt.xlabel('Posição ao longo da barra (m)')
        plt.ylabel('Deflexão (m)')
        plt.title('Deflexão devido à Flexão ao longo da barra')
        plt.legend()

        # Plot torção
        plt.subplot(1, 3, 3)
        plt.plot(x, torsao, label='Torção', color='green')
        plt.xlabel('Posição ao longo da barra (m)')
        plt.ylabel('Torção')
        plt.title('Torção ao longo da barra')
        plt.legend()

        plt.tight_layout()
        plt.show()

# Exemplo de uso:
A = 0.01    # m^2
E = 210e9   # Pa
L = 2.0     # m
I = 1.6667e-5 # m^4
G = 81.2e9  # Pa
a = 0.1     # m
B = 0.1     # m
K = 1000    # N/m
y = 0.01    # m
rho = 7850  # kg/m^3
k = 0.01    # Ns/m
g = 9.81    # m/s^2
J = 1e-6    # m^4 (momento polar de inércia)
num_elements = 10  # Número de elementos finitos

beam = BodyAndFrameFEM(A, E, L, I, G, a, B, K, y, rho, k, g, J, num_elements)
F = np.linspace(0, 4000, num_elements + 1)  # Exemplo de forças aplicadas

beam.plot_results(F)
