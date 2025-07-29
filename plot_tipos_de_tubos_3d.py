import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D

# Funções principais

def calcular_comprimento(node1, node2):
    return np.linalg.norm(node2 - node1)

def calcular_comprimentos_e_quantidades_por_tipo(nodes, elements):
    resultados = {}
    for idx1, idx2, tipo in elements:
        node1 = nodes[idx1]
        node2 = nodes[idx2]
        comprimento = calcular_comprimento(node1, node2)

        if tipo not in resultados:
            resultados[tipo] = {'comprimento_total': 0.0, 'quantidade': 0}

        resultados[tipo]['comprimento_total'] += comprimento
        resultados[tipo]['quantidade'] += 1

    return resultados

def set_axes_equal(ax):
    # Deixa a escala dos eixos x, y e z iguais
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max(x_range, y_range, z_range)

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)

    ax.set_xlim3d([x_middle - max_range / 2, x_middle + max_range / 2])
    ax.set_ylim3d([y_middle - max_range / 2, y_middle + max_range / 2])
    ax.set_zlim3d([z_middle - max_range / 2, z_middle + max_range / 2])

def plotar_chassi_3d(nodes, elements):
    # Definindo cores para tipos diferentes de tubo
    cores_tubos = {
        'Tubo A': 'red',
        'Tubo B': 'blue',
        'Tubo C': 'green',
        'Tubo D': 'orange'
    }

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plotar os nós
    for i, (x, y, z) in enumerate(nodes):
        ax.scatter(x, y, z, color='black', s=30)
        ax.text(x, y, z, f'{i}', color='black', fontsize=9)

    # Plotar os tubos
    for idx1, idx2, tipo in elements:
        cor = cores_tubos.get(tipo, 'gray')
        x = [nodes[idx1][0], nodes[idx2][0]]
        y = [nodes[idx1][1], nodes[idx2][1]]
        z = [nodes[idx1][2], nodes[idx2][2]]
        ax.plot(x, y, z, color=cor, linewidth=2)

    # Ajustar escala dos eixos para manter proporção
    set_axes_equal(ax)

    # Adicionar legenda
    legend_elements = [
        Line2D([0], [0], color=cor, lw=3, label=tipo)
        for tipo, cor in cores_tubos.items()
    ]
    ax.legend(handles=legend_elements, title="Tipos de Tubo")

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('Visualização 3D da Estrutura do Chassi')
    plt.tight_layout()
    plt.show()

# Dados dos nodes e elements

"""nodes = np.array([
    [1.92, 0.0, 0.0], [1.92, 0.64, 0.0], [1.92, 0.0, 0.48], [1.92, 0.64, 0.48],
    [1.77, 0.0, 0.21], [1.77, 0.64, 0.21], [1.92, 0.0, 0.09], [1.92, 0.64, 0.09],
    [1.5, 0.0, 0.03], [1.5, 0.64, 0.03], [1.14, 0.08, 0.03], [1.14, 0.56, 0.03],
    [1.14, 0.0, 0.09], [1.14, 0.64, 0.09], [1.23, 0.0, 0.36], [1.23, 0.64, 0.36],
    [1.14, 0.03, 0.72], [1.14, 0.6, 0.72], [0.63, 0.0, 0.54], [0.63, 0.64, 0.54],
    [0.69, 0.0, 0.24], [0.69, 0.64, 0.24], [0.69, 0.0, 0.0], [0.69, 0.64, 0.0],
    [0.45, 0.0, 0.21], [0.45, 0.64, 0.21], [0.24, 0.0, 0.09], [0.24, 0.64, 0.09],
    [0.0, 0.16, 0.21], [0.0, 0.48, 0.21], [0.0, 0.16, 0.09], [0.0, 0.48, 0.09],
    [0.0, 0.16, 0.42], [0.0, 0.48, 0.42], [0.33, 0.03, 0.66], [0.33, 0.6, 0.66],
    [0.57, 0.03, 1.2], [0.57, 0.6, 1.2], [0.54, 0.32, 1.35], [1.14, 0.32, 0.78]
])

elements = [
    (0, 1, 'Tubo B'), (0, 6, 'Tubo B'), (6, 2, 'Tubo B'), (1, 7, 'Tubo B'), (7, 3, 'Tubo B'),
    (2, 3, 'Tubo C'), (4, 0, 'Tubo C'), (4, 2, 'Tubo C'), (5, 1, 'Tubo C'), (5, 3, 'Tubo C'),
    (4, 5, 'Tubo A'), (6, 7, 'Tubo A'), (0, 8, 'Tubo A'), (1, 9, 'Tubo A'), (4, 8, 'Tubo A'),
    (5, 9, 'Tubo D'), (8, 9, 'Tubo D'), (10, 8, 'Tubo D'), (10, 4, 'Tubo D'), (11, 9, 'Tubo D'),
    (11, 5, 'Tubo A'), (10, 11, 'Tubo A'), (12, 10, 'Tubo A'), (12, 4, 'Tubo A'), (13, 11, 'Tubo A'),
    (13, 5, 'Tubo A'), (14, 12, 'Tubo A'), (14, 4, 'Tubo A'), (15, 13, 'Tubo A'), (15, 5, 'Tubo A'),
    (16, 14, 'Tubo A'), (16, 4, 'Tubo A'), (17, 15, 'Tubo A'), (17, 5, 'Tubo A'), (2, 16, 'Tubo A'),
    (3, 17, 'Tubo A'), (16, 18, 'Tubo A'), (17, 19, 'Tubo A'), (20, 18, 'Tubo A'), (20, 16, 'Tubo A'),
    (20, 14, 'Tubo A'), (20, 10, 'Tubo A'), (21, 19, 'Tubo A'), (21, 17, 'Tubo A'), (21, 15, 'Tubo A'),
    (21, 11, 'Tubo A'), (22, 10, 'Tubo A'), (22, 20, 'Tubo A'), (23, 11, 'Tubo A'), (23, 21, 'Tubo A'),
    (22, 23, 'Tubo A'), (24, 18, 'Tubo A'), (24, 20, 'Tubo A'), (24, 22, 'Tubo A'), (25, 19, 'Tubo A'),
    (25, 21, 'Tubo A'), (25, 23, 'Tubo A'), (26, 22, 'Tubo A'), (26, 24, 'Tubo A'), (27, 23, 'Tubo A'),
    (27, 25, 'Tubo A'), (26, 27, 'Tubo A'), (28, 30, 'Tubo A'), (28, 32, 'Tubo A'), (29, 31, 'Tubo A'),
    (29, 33, 'Tubo A'), (30, 26, 'Tubo A'), (31, 27, 'Tubo A'), (30, 31, 'Tubo A'), (28, 24, 'Tubo A'),
    (29, 25, 'Tubo A'), (32, 24, 'Tubo A'), (32, 18, 'Tubo A'), (33, 25, 'Tubo A'), (33, 19, 'Tubo A'),
    (32, 33, 'Tubo A'), (34, 18, 'Tubo A'), (34, 32, 'Tubo A'), (35, 19, 'Tubo A'), (35, 33, 'Tubo A'),
    (34, 35, 'Tubo A'), (36, 34, 'Tubo A'), (36, 18, 'Tubo A'), (37, 35, 'Tubo A'), (37, 19, 'Tubo A'),
    (36, 38, 'Tubo A'), (37, 38, 'Tubo A'), (16, 39, 'Tubo A'), (17, 39, 'Tubo A')
]
""" 
df_nodes = pd.read_csv("solucao_otimizada_nodes.csv")
nodes = df_nodes[['x', 'y', 'z']].to_numpy()

 #Carregar elementos
df_elements = pd.read_csv("solucao_otimizada_elements.csv")
elements = [(row.node_i, row.node_j, row.perfil) for row in df_elements.itertuples()]

# Executando


resultados = calcular_comprimentos_e_quantidades_por_tipo(nodes, elements)

for tipo, dados in resultados.items():
    print(f"{tipo}:")
    print(f"  - Comprimento total: {dados['comprimento_total']:.3f} metros")
    print(f"  - Quantidade de tubos: {dados['quantidade']}")

plotar_chassi_3d(nodes, elements)






