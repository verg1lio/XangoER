import cadquery as cq
from cadquery import exporters
import numpy as np
import matplotlib.pyplot as plt

# Escalas para controlar coordenadas
i, j, k = 0.3, 0.4, 0.3

#nós
nodes = np.array ([
[64*i, 0*j, 0*k] , #nó 0
[64*i, 16*j, 0*k] , #nó 1
[64*i, 0*j, 16*k] , #nó 2
[64*i, 16*j, 16*k] ,#nó 3
[59*i, 0*j, 7*k] ,#nó 4
[59*i, 16*j, 7*k] ,#nó 5
[64*i, 0*j, 3*k] ,#nó 6
[64*i, 16*j, 3*k] , #nó 7
[50*i, 0*j, 1*k] ,#nó 8
[50*i, 16*j, 1*k] , #nó 9
[38*i, 2*j, 1*k] ,#nó 10
[38*i, 14*j, 1*k] ,#nó 11
[38*i, 0*j, 3*k] , #nó 12
[38*i, 16*j, 3*k] , #nó 13
[38*i, 0*j, 12*k] , #nó 14
[41*i, 16*j, 12*k] , #nó 56
[38*i, 1*j, 24*k] , #nó 16
[38*i, 15*j, 24*k] , #nó 17
[21*i, 0*j, 18*k] ,#nó 18
[21*i, 16*j, 18*k] ,#nó 19
[23*i, 0*j, 8*k] , #20
[23*i, 16*j, 8*k] ,#21
[23*i, 0*j, 0*k] ,#22
[23*i, 16*j, 0*k] , #23
[15*i, 0*j, 7*k] ,#24
[15*i, 16*j, 7*k] ,#25
[8*i, 0*j, 3*k] ,#26
[8*i, 16*j, 3*k] ,#27
[0*i, 4*j, 7*k] ,#28
[0*i, 12*j, 7*k] ,#29
[0*i, 4*j, 3*k] ,#30
[0*i, 12*j, 3*k] ,#31
[0*i, 4*j, 14*k],#32
[0*i, 12*j, 14*k] ,#33
[11*i, 1*j, 22*k] ,#34
[11*i, 15*j, 22*k] , #35
[19*i, 1*j, 40*k] ,#36
[19*i, 15*j, 40*k] , #37
[18*i, 8*j, 45*k] ,#38
[38*i, 8*j, 26*k]#39
])
# Conexões entre os nós
elements = [
(0,1), 
(0,2),
(1,3),
(2,3),
(4,0),
(4,2),
(5,1),
(5,3),
(4,5),
(6,7),
(0,8),
(1,9),
(4,8),
(5,9),
(8,9),
(10,8),
(10,4),
(11,9),
(11,5),
(10,11),
(12,10),
(12,4),
(13,11),
(13,5),
(12,14),
(14,4),
(15,13),
(15,5),
(16,14),
(16,4),
(17,15),
(17,5),
(2,16),
(3,17),
(16,18),
(17,19),
(17,21),
(20,18),
(20,16),
(20,14),
(20,10),
(21,19),
(22,10),
(22,26),
(36,38),
(37,38),
(36,18),
(37,19),
(23,11),
(27,23),
(16,39),
(17,39),
(36,34),
(37,35),
(30,28),
(28,32),
(31,29),
(29,31),
(30,26),
(31,27),
(36,34),
(36,34),
(30,31),
(34,18),
(35,19),
(35,33),
(34,35),
(24,20),
(26,27),
(28,24),
(32,33),
(34,32),
(28,24),
(32,18),
(33,19),
(32,24),
(33,25),
(23,21),
(25,21),
(25,23),
(29,25),
(22,20),
(24,18),
(24,22),
(25,19),
(29,33),
(11,21),
(21,15),
(26,24),
(27,25)
]

# Parâmetros das vigas tubulares
beam_outer_radius = 0.1
beam_wall_thickness = 0.05

def tubo_oco_entre(p1, p2, outer_radius, wall_thickness):
    vec = p2 - p1
    length = np.linalg.norm(vec)
    if length == 0:
        return None

    outer_cyl = cq.Workplane("XY").circle(outer_radius).extrude(length)
    inner_cyl = cq.Workplane("XY").circle(outer_radius - wall_thickness).extrude(length)
    tubo = outer_cyl.cut(inner_cyl)

    z_axis = np.array([0, 0, 1])
    vec_normalized = vec / length

    axis = np.cross(z_axis, vec_normalized)
    axis_length = np.linalg.norm(axis)

    if axis_length < 1e-6:
        # Vetor alinhado com eixo Z, sem rotação necessária
        tubo = tubo.translate(tuple(p1))
    else:
        axis = axis / axis_length
        angle = np.degrees(np.arccos(np.dot(z_axis, vec_normalized)))
        tubo = tubo.rotate((0, 0, 0), tuple(axis), angle).translate(tuple(p1))

    return tubo

# Construir estrutura 
estrutura = cq.Workplane("XY")
for start_idx, end_idx in elements:
    p1 = nodes[start_idx]
    p2 = nodes[end_idx]
    tubo = tubo_oco_entre(p1, p2, beam_outer_radius, beam_wall_thickness)
    if tubo:
        estrutura = estrutura.union(tubo)

# Exportar STEP
exporters.export(estrutura, "retangulo_4_vigas_custom.step")
print("✅ STEP file 'retangulo_4_vigas_custom.step' generated!")

# Visualização 3D 
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Mostrar nós 
for idx, (x, y, z) in enumerate(nodes):
    ax.scatter(x, y, z, c='red', s=60)
    ax.text(x, y, z, f'{idx}', color='black', fontsize=10,
            verticalalignment='bottom', horizontalalignment='right')

# Mostrar vigas 
for start_idx, end_idx in elements:
    p1 = nodes[start_idx]
    p2 = nodes[end_idx]
    xs, ys, zs = zip(p1, p2)
    ax.plot(xs, ys, zs, color='blue', linewidth=5, alpha=0.7)

ax.set_xlabel("X [mm]")
ax.set_ylabel("Y [mm]")
ax.set_zlabel("Z [mm]")
ax.set_title("Rectangular Frame (4 Tubular Beams)")
ax.set_box_aspect([1, 1, 0.5])

plt.tight_layout()
plt.show()