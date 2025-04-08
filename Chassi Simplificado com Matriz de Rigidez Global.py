import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Estrutura:
    def __init__(self, E, I, nodes, elements):
        self.E = E
        self.I = I
        self.nodes = nodes
        self.elements = elements
        self.num_dofs_per_node = 6
        self.num_nodes = len(nodes)
        self.num_dofs = self.num_nodes * self.num_dofs_per_node
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))
    
    def connect_matrix(self):
        # Número de conexões
        num_connect = len(self.elements)

        # Gerando a matriz de conectividade: (Para cada linha, teremos: [índice a coneção, 1º nó, 2º nó])
        CM = np.zeros((num_connect, 3), dtype=int)

        # Preenchendo a matriz:
        for i, (no1, no2) in enumerate(self.elements, start=0):
            CM[i][0] = i + 1
            CM[i][1] = no1
            CM[i][2] = no2

        print("\n Conexão   1º Nó   2º Nó")
        print(CM)
        
    def node_loc_matrix(self, node_tags, node_coord):
        # Número de Nós
        num_nodes = len(node_tags)

        # Gerando uma matrix número de nos x 4: (Para cada linha, teremos: [índice do nó, x, y, z])
        node_loc_matrix = np.zeros((num_nodes, 4), dtype=float)  # Na primeira coluna, teremos os índices, nas seguintes, x y e z.

        # Preenchendo a matriz de zeros
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i]  # Número do nó na primeira coluna     
            node_loc_matrix[i][1] = x  # Coordenada x na segunda coluna
            node_loc_matrix[i][2] = y  # Coordenada y na terceira coluna
            node_loc_matrix[i][3] = z  # Coordenada z na quarta coluna

        print("\n   Nó   x   y   z")
        print(node_loc_matrix)

    def matriz_rigidez_global(self):
        for element in self.elements:
            node1, node2 = element
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            L_e = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            K_e = np.zeros((12, 12))
            coef = (2 * self.E * self.I / L_e**3)
            K_e[:4, :4] = coef * np.array([
                [6, -3 * L_e, -6, -3 * L_e],
                [3 * L_e, 2 * L_e**2, -3 * L_e, L_e**2],
                [-6, -3 * L_e, 6, 3 * L_e],
                [-3 * L_e, L_e**2, 3 * L_e, 2 * L_e**2]])
            K_e[6:10, 6:10] = K_e[:4, :4]
            dofs_node1 = [node1 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
            dofs_node2 = [node2 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
            dofs = dofs_node1 + dofs_node2
            for i in range(len(dofs)):
                for j in range(len(dofs)):
                    self.K_global[dofs[i], dofs[j]] += K_e[i, j]
        return self.K_global

# Exemplo de uso
E = 210e9  # Módulo de elasticidade em Pa
I = 8.33e-6  # Momento de inércia em m^4

# Coordenadas dos nós (x, y, z)
nodes = [
    (0, 0, 0),
    (0, 375, 0),
    (0, 700, 0),
    (1500, 375, 0),
    (1500, 0, 0),
    (1500, 700, 0)
]

# Conectividade dos elementos (índices dos nós)
elements = [
    (0, 1),  # Elemento entre os nós 0 e 1
    (1, 2),  # Elemento entre os nós 1 e 2
    (4, 3),  # Elemento entre os nós 4 e 3
    (3, 5),  # Elemento entre os nós 3 e 5
    (1, 3)   # Elemento entre os nós 1 e 3
]

# Criar a estrutura e montar a matriz de rigidez global
estrutura = Estrutura(E, I, nodes, elements)
K_global = estrutura.matriz_rigidez_global()
# Gerar as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
estrutura.node_loc_matrix(node_tags, nodes)
estrutura.connect_matrix()

print("\n Matriz de Rigidez Global")
print(K_global)

# Convertendo a matriz de conectividade para um array numpy
connectivity_array = np.array(elements)

# Plotando o gráfico 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Adicionando os pontos
for i, (x, y, z) in enumerate(nodes):
    ax.scatter(x, y, z, color='b', s=100)
    ax.text(x, y, z, f'  {i}', color='black', fontsize=12)

# Adicionando as linhas de ligação entre os nós
for node1, node2 in elements:
    x = [nodes[node1][0], nodes[node2][0]]
    y = [nodes[node1][1], nodes[node2][1]]
    z = [nodes[node1][2], nodes[node2][2]]
    ax.plot(x, y, z, marker='o')

# Configurações adicionais
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Estrutura 3D')

plt.show()