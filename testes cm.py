import numpy as np
import csv
import math


nodes = np.array([
    [0, 0.0, 0.0],    # 64*
    [0, 20, 5],   # 64*i
    
])

elements = [(0, 1, 'Tubo A')]

def read_tube_properties(csv_file):
    props = {}

    with open(csv_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        for row in reader:
            name = row['Tube']
            density = float(row['Densidade(kg/m^3)'])
            thickness = float(row['Espessura(m)'])
            shape = row['Shape']

            # Área da seção transversal
            if shape == 'Circular':
                D = float(row['Diametro(m)'])
                r_ext = D / 2
                r_int = r_ext - thickness
                area = math.pi * (r_ext**2 - r_int**2)

            elif shape == 'Square':
                side = float(row['Side(m)'])
                area = side**2 - (side - 2*thickness)**2

            else:
                raise ValueError(f"Forma não suportada: {shape}")

            props[name] = {
                'rho': density,
                'area': area
            }

    return props

def center_of_mass_truss(nodes, elements, tube_props):
    total_mass = 0.0
    weighted_sum = np.zeros(3)

    for i, j, tube in elements:
        p1 = nodes[i]
        p2 = nodes[j]

        length = np.linalg.norm(p2 - p1)
        centroid = 0.5 * (p1 + p2)

        rho = tube_props[tube]['rho']
        area = tube_props[tube]['area']
        mass = rho * area * length

        weighted_sum += mass * centroid
        total_mass += mass

    return weighted_sum / total_mass

# Ler propriedades do CSV
tube_props = read_tube_properties("tubos_atualizado.csv")

# Calcular centro de massa
cm = center_of_mass_truss(nodes, elements, tube_props)

print("Centro de massa da estrutura:")
print(f"x = {cm[0]:.6f}")
print(f"y = {cm[1]:.6f}")
print(f"z = {cm[2]:.6f}")


