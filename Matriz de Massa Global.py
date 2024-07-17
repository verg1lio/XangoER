import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def modal_analysis(K, M, num_modes=20):
    # Adicionar um valor pequeno na diagonal da matriz de massa para garantir que seja positiva definida
    M += np.eye(M.shape[0]) * 1e-8
    
    unsorted_eigenvalues, unsorted_eigenvectors = eigh(K, M)
    unsorted_frequencies = np.sqrt(unsorted_eigenvalues) / (2 * np.pi)

    sorted_indices = np.argsort(unsorted_frequencies)
    top_indices = sorted_indices[:num_modes]

    eigenvalues = np.array(unsorted_eigenvalues)[top_indices]
    eigenvectors = np.array(unsorted_eigenvectors)[top_indices]
    frequencies = np.array(unsorted_frequencies)[top_indices]

    return eigenvalues, eigenvectors, frequencies

class Estrutura:
    def __init__(self, E, I, rho, A, nodes, elements):
        self.E = E
        self.I = I
        self.rho = rho
        self.A = A
        self.nodes = nodes
        self.elements = elements
        self.num_dofs_per_node = 6
        self.num_nodes = len(nodes)
        self.num_dofs = self.num_nodes * self.num_dofs_per_node
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))
        self.M_global = np.zeros((self.num_dofs, self.num_dofs))
    
    def connect_matrix(self):
        num_connect = len(self.elements)
        CM = np.zeros((num_connect, 3), dtype=int)
        for i, (no1, no2) in enumerate(self.elements, start=0):
            CM[i][0] = i + 1
            CM[i][1] = no1
            CM[i][2] = no2
        print("\n Conexão   1º Nó   2º Nó")
        print(CM)
        
    def node_loc_matrix(self, node_tags, node_coord):
        num_nodes = len(node_tags)
        node_loc_matrix = np.zeros((num_nodes, 4), dtype=float)
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i]
            node_loc_matrix[i][1] = x
            node_loc_matrix[i][2] = y
            node_loc_matrix[i][3] = z
        print("\n   Nó   x   y   z")
        print(node_loc_matrix)
    
    def matriz_rigidez_global(self):
        for node1, node2 in self.elements:
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            L_e = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])
            coef = 2 * self.E * self.I / L_e**3
            K_e = coef * np.array([
                [6, -3 * L_e, -6, -3 * L_e],
                [3 * L_e, 2 * L_e**2, -3 * L_e, L_e**2],
                [-6, -3 * L_e, 6, 3 * L_e],
                [-3 * L_e, L_e**2, 3 * L_e, 2 * L_e**2]
            ])
            dofs_node1 = [node1 * self.num_dofs_per_node + i for i in range(4)]
            dofs_node2 = [node2 * self.num_dofs_per_node + i for i in range(4)]
            dofs = dofs_node1 + dofs_node2
            for i in range(len(dofs)):
                for j in range(len(dofs)):
                    self.K_global[dofs[i], dofs[j]] += K_e[i % 4, j % 4] if i < 4 and j < 4 else 0
        return self.K_global
    
    def bar_mass_matrix(self):
        self.M_global.fill(0)  # Resetting M_global to avoid overlap with beam matrix
        for node1, node2 in self.elements:
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            L_e = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])
            m_e = (self.rho * self.A * L_e / 6) * np.array([
                [2, 1, 0, 0, 0, 0, 1, 0],
                [1, 2, 0, 0, 0, 0, 1, 0],
                [0, 0, 2, 1, 0, 0, 0, 1],
                [0, 0, 1, 2, 0, 0, 0, 1],
                [0, 0, 0, 0, 2, 1, 0, 0],
                [0, 0, 0, 0, 1, 2, 0, 0],
                [1, 1, 0, 0, 0, 0, 2, 0],
                [0, 0, 1, 1, 0, 0, 0, 2],
            ])
            dofs_node1 = [node1 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
            dofs_node2 = [node2 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
            dofs = dofs_node1 + dofs_node2
            for i in range(len(dofs)):
                for j in range(len(dofs)):
                    self.M_global[dofs[i], dofs[j]] += m_e[i % 8, j % 8] if i < 8 and j < 8 else 0
        return self.M_global
    
    def beam_mass_matrix(self):
        self.M_global.fill(0)  # Resetting M_global to avoid overlap with bar matrix
        for node1, node2 in self.elements:
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            L_e = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])
            m_e = (self.rho * self.A * L_e / 420) * np.array([
                [156, 22 * L_e, 54, -13 * L_e],
                [22 * L_e, 4 * L_e**2, 13 * L_e, -3 * L_e**2],
                [54, 13 * L_e, 156, -22 * L_e],
                [-13 * L_e, -3 * L_e**2, -22 * L_e, 4 * L_e**2]
            ])
            dofs_node1 = [node1 * self.num_dofs_per_node + i for i in range(4)]
            dofs_node2 = [node2 * self.num_dofs_per_node + i for i in range(4)]
            dofs = dofs_node1 + dofs_node2
            for i in range(4):
                for j in range(4):
                    self.M_global[dofs[i], dofs[j]] += m_e[i, j]
        return self.M_global
    
    def plot_structure(self):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Adicionando os pontos
        for i, (x, y, z) in enumerate(self.nodes):
            ax.scatter(x, y, z, color='b', s=50)
            ax.text(x, y, z, f'{i}', color='red')

        # Adicionando as linhas
        for node1, node2 in self.elements:
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            ax.plot([x1, x2], [y1, y2], [z1, z2], color='k')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.title('Estrutura 3D')
        plt.show()

# Definindo parâmetros da estrutura
E = 210e9  # Módulo de elasticidade (Pa)
I = 8.1e-6  # Momento de inércia (m^4)
rho = 7850  # Densidade (kg/m^3)
A = 0.01  # Área da seção transversal (m^2)

# Definindo nós e elementos (ajustado)
nodes = [
    (0, 0, 0),
    (0, 0.375, 0),
    (0, 0.700, 0),
    (1.500, 0.375, 0),
    (1.500, 0, 0),
    (1.500, 0.700, 0)
]
elements = [
    (0, 1),
    (1, 2),
    (4, 3),
    (3, 5),
    (1, 3)
]

# Criando objeto da estrutura
estrutura = Estrutura(E, I, rho, A, nodes, elements)
estrutura.connect_matrix()
estrutura.node_loc_matrix([0, 1, 2, 3, 4, 5], nodes)
K_global = estrutura.matriz_rigidez_global()
M_global = estrutura.bar_mass_matrix()  # Ou estrutura.beam_mass_matrix(), dependendo do tipo de análise desejada

# Exibindo a matriz de massa global
print("\n Matriz de Massa Global:")
print(M_global)

# Realizando a análise modal
autovalores, autovetores, frequências = modal_analysis(K_global, M_global)

# Exibindo os resultados
print("\n Autovalores:")
print(autovalores)
print("\n Autovetores:")
print(autovetores)
print("\n Frequências:")
print(frequências)

# Plotando a estrutura
estrutura.plot_structure()
