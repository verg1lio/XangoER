import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


# Coordenadas dos nós
nodes = [
    (0, 0, 0),
    (0, 0.375, 0),
    (0, 0.700, 0),
    (1.500, 0.375, 0),
    (1.500, 0, 0),
    (1.500, 0.700, 0)
]

# Conexões entre os nós (elementos)
elements = [
    (0, 1),  # Elemento entre os nós 0 e 1
    (1, 2),  # Elemento entre os nós 1 e 2
    (4, 3),  # Elemento entre os nós 4 e 3
    (3, 5),  # Elemento entre os nós 3 e 5
    (1, 3)   # Elemento entre os nós 1 e 3
]

# Inicializar uma lista para armazenar as conexões
connections = []

for i, element in enumerate(elements):
    node_start, node_end = element
    connections.append([ i+1, node_start, node_end])

# Converter a lista em um array numpy
connections_matrix = np.array(connections)

#print("Matriz de conexões:")
#print(connections_matrix)

nodes =np.array(nodes)

#Plotando a estrutura    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plotar os nós
ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], c='b')   
        
# Conectando os nós
for element in elements:
    node_start, node_end = element
    ax.plot([nodes[node_start, 0], nodes[node_end, 0]], # X
            [nodes[node_start, 1], nodes[node_end, 1]], # Y
            [nodes[node_start, 2], nodes[node_end, 2]]) # Z
        
# Numerar os nós
for i, node in enumerate(nodes):
     ax.text(node[0], node[1], node[2], str(i), color='black')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title("Estrutura 3D")
plt.show()
