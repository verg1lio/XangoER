import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# NÓS E CONEXÕES 

nodes = np.array([
    # LADO ESQUERDO 
    [-0.1600, -0.1500, 0.3500],   # 0
    [-0.1600, -0.1500, 0.0000],   # 1
    [-0.2450,  0.2750, 0.2200],   # 2
    [-0.2450,  0.6050, 0.0000],   # 3
    [-0.2450,  0.2550, 0.0000],   # 4
    [-0.2450,  0.5850, 0.2200],   # 5
    [-0.2100,  0.5550, 0.5000],   # 6
    [-0.2900,  1.3700, 0.2500],   # 7
    [-0.2680,  1.3500, 0.0000],   # 8
    [-0.3020,  1.6700, 0.2400],   # 9
    [-0.2710,  1.6650, 0.0300],   #10
    [-0.2710,  1.9500, 0.0300],   #11
    [-0.2710,  2.2300, 0.2500],   #12
    [-0.2710,  2.2300, 0.0300],   #13
    [-0.1700,  1.4000, 0.9650],   #14
    [-0.2710,  1.9500, 0.2500],   #15
    [ 0.0000,  1.4100, 1.1050],   #16
    [ 0.0000,  0.5550, 0.5500],   #17
    [-0.1850,  0.2025, 0.4250],   #18
    # LADO DIREITO 
    [ 0.1600, -0.1500, 0.3500],   #19
    [ 0.1600, -0.1500, 0.0000],   #20
    [ 0.2450,  0.2750, 0.2200],   #21
    [ 0.2450,  0.6050, 0.0000],   #22
    [ 0.2450,  0.2550, 0.0000],   #23
    [ 0.2450,  0.5850, 0.2200],   #24
    [ 0.2100,  0.5550, 0.5000],   #25
    [ 0.2900,  1.3700, 0.2500],   #26
    [ 0.2680,  1.3500, 0.0000],   #27
    [ 0.3020,  1.6700, 0.2400],   #28
    [ 0.2710,  1.6650, 0.0300],   #29
    [ 0.2710,  1.9500, 0.0300],   #30
    [ 0.2710,  2.2300, 0.2500],   #31
    [ 0.2710,  2.2300, 0.0300],   #32
    [ 0.1700,  1.4000, 0.9650],   #33
    [ 0.2710,  1.9500, 0.2500],   #34
    [ 0.0000,  1.4100, 1.1050],   #35
    [ 0.0000,  0.5550, 0.5500],   #36
    [ 0.1850,  0.2025, 0.4250],   #37
    #NÓS CENTRAIS
    [-0.2710, 2.0900, 0.2500],    #38
    [ 0.2710, 2.0900, 0.2500],    #39
])

connections = [
    # LADO ESQUERDO
    (0,1),(0,2),(0,6),(1,2),(1,4),(2,3),(2,4),(2,5),(2,6),
    (3,4),(3,5),(3,8),(5,6),(5,7),(5,8),(6,7),(6,17),(7,8),
    (7,9),(7,14),(8,9),(8,10),(9,10),(9,11),(9,15),(10,11),
    (11,12),(11,13),(11,15),(12,13),(12,15),(12,14),(14,16),
    # LADO DIREITO
    (19,20),(19,21),(19,25),(20,21),(20,23),(21,22),(21,23),(21,24),(21,25),
    (22,23),(22,24),(22,27),(24,25),(24,26),(24,27),(25,26),(25,36),(26,27),
    (26,28),(26,33),(27,28),(27,29),(28,29),(28,30),(28,34),(29,30),
    (30,31),(30,32),(30,34),(31,32),(31,34),(31,33),(33,35),
    # BARRAS TRANSVERSAIS
    (18,37),(12,31),(13,32),(0,19),(1,20),(11,30),(8,27),
    # BARRA CENTRAL
    (38,39),
]

# PROJEÇÃO 3D DE UM NÓ SOBRE UMA RETA

def project_node_to_line(node_idx, P0_idx, P1_idx):
    """Projeta node_idx na reta P0 → P1 (3D)."""
    P0, P1, P = nodes[P0_idx], nodes[P1_idx], nodes[node_idx]
    v = P1 - P0
    w = P - P0
    t = np.dot(w, v) / np.dot(v, v)
    nodes[node_idx] = P0 + t * v

# nós 18 e 37 projetados nas diagonais 0-6 e 19-25
project_node_to_line(18, 0, 6)
project_node_to_line(37, 19, 25)

# EXPORTAÇÃO PARA STEP

import cadquery as cq
from cadquery import exporters

nodes_mm = nodes * 1000  # m → mm

edges = [
    cq.Edge.makeLine(cq.Vector(*nodes_mm[a]), cq.Vector(*nodes_mm[b]))
    for a, b in connections
]

chassi = cq.Compound.makeCompound(edges)
exporters.export(chassi, "chassi_final_linhas.step")
print("Arquivo STEP gerado: chassi_final_linhas.step")

# PLOTAGEM

fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], c='red', s=45)
for i, (x, y, z) in enumerate(nodes):
    ax.text(x, y, z, f'{i}', fontsize=9)

for a, b in connections:
    ax.plot([nodes[a,0], nodes[b,0]],
            [nodes[a,1], nodes[b,1]],
            [nodes[a,2], nodes[b,2]], color='black', linewidth=2)

ax.set_box_aspect([1, 3, 2])
plt.show()