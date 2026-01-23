import math
import matplotlib.pyplot as plt

def calcular_momento_inercia_polar(diametro):
    raio = diametro / 2
    return (math.pi * (raio**4)) / 2

def calcular_momento_inercia_area(diametro):
    return (math.pi * (diametro**4)) / 64

def calcular_tensao_torcao(torque, diametro, kt_torcao=1.0):
    raio = diametro / 2
    j = calcular_momento_inercia_polar(diametro)
    if j == 0:
        return float('inf')
    tensao_torcao = (torque * raio) / j
    return tensao_torcao * kt_torcao

def calcular_tensao_flexao(momento_fletor, diametro, kt_flexao=1.0):
    c = diametro / 2
    i = calcular_momento_inercia_area(diametro)
    if i == 0:
        return float('inf')
    tensao_flexao = (momento_fletor * c) / i
    return tensao_flexao * kt_flexao

def calcular_tensao_von_mises(tensao_normal, tensao_cisalhamento):
    return math.sqrt(tensao_normal**2 + 3 * tensao_cisalhamento**2)

def plotar_tensao_von_mises(secoes_semieixo, torque_max, momento_fletor_max, material_propriedades):
    tensoes_von_mises = []
    diametros = []

    for secao in secoes_semieixo:
        diametro = secao['diametro'] / 1000  # mm para m
        kt_torcao = secao.get('kt_torcao', 1.0)
        kt_flexao = secao.get('kt_flexao', 1.0)

        tensao_torcao = calcular_tensao_torcao(torque_max, diametro, kt_torcao)
        tensao_flexao = calcular_tensao_flexao(momento_fletor_max, diametro, kt_flexao)
        tensao_von_mises = calcular_tensao_von_mises(tensao_flexao, tensao_torcao)

        tensoes_von_mises.append(tensao_von_mises / 1e6)  # convertendo para MPa
        diametros.append(secao['diametro'])

    # Gráfico
    plt.figure(figsize=(8, 5))
    plt.plot(diametros, tensoes_von_mises, marker='o', color='b', label='Tensão de von Mises (MPa)')
    plt.axhline(y=material_propriedades['limite_escoamento'] / 1e6, color='r', linestyle='--', label='Limite de Escoamento (MPa)')
    plt.axhline(y=material_propriedades['limite_resistencia'] / 1e6, color='g', linestyle='--', label='Limite de Resistência (MPa)')

    plt.title('Tensão Equivalente de von Mises ao longo do Semieixo')
    plt.xlabel('Diâmetro da Seção (mm)')
    plt.ylabel('Tensão von Mises (MPa)')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    material_propriedades = {  # Aço 4340
        'limite_escoamento': 1000e6,
        'limite_resistencia': 1200e6
    }

    torque_max_semieixo = 500  # Nm
    momento_fletor_max_semieixo = 300  # Nm

    secoes_semieixo = [
        {'diametro': 30, 'kt_flexao': 1.2, 'kt_torcao': 1.1},
        {'diametro': 20, 'kt_flexao': 2.0, 'kt_torcao': 1.8},
    ]

    # Gráfico de tensão equivalente de von Mises
    plotar_tensao_von_mises(secoes_semieixo, torque_max_semieixo, momento_fletor_max_semieixo, material_propriedades)
