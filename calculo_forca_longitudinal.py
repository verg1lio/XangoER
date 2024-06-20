import math
import matplotlib.pyplot as plt

n = 4          # numero_de_pistoes
pi = 3.14                 #valor de pi
area_do_pistao = 285.07  # Área do pistão Acm em mm^2
diametro_dos_pistoes = 30  
area_efetiva_dos_pistoes = pi * ( (diametro_dos_pistoes / 2) ** 2)
nova_area_efetiva_dos_pistoes = area_efetiva_dos_pistoes * n   # Área efetiva levando em conta o numero de pistões
forca_pedal = 880

forca_longitudinal = forca_pedal * nova_area_efetiva_dos_pistoes / area_do_pistao


print(f'Força Longitudinal {forca_longitudinal:.2f}, quando aplicada uma força no pedal de {forca_pedal}N.')


''' variando 200N
    =>Força aplicada no pedal : [245, 445, 665, 845, 1045]
    =>Força longitudinal : [2428.77, 4411.44, 6394.11, 8376.78, 10359.45] '''

''' variando 120N
    =>Força aplicada no pedal : [120, 240, 360, 480, 500, 620, 740, 860, 980]
    =>Força longitudinal : [1189, 2379.21, 3568.81, 3568.81, 4956.68, 6146.28, 7335.88,  8525.48, 9715.09] '''

forca_pedal = [120, 240, 360, 480, 500, 620, 740, 860, 980]  # Exemplo de força no pedal 
forca_longitudinal = [1189, 2379.21, 3568.81, 3568.81, 4956.68, 6146.28, 7335.88,  8525.48, 9715.09]  # Exemplo de força longitudinal no pneu 
# Criando o gráfico de dispersão
plt.figure(figsize=(8, 6))
plt.plot(forca_pedal, forca_longitudinal, color='blue', marker='o', label='Força no pedal vs Força longitudinal no pneu')

# Adicionando rótulos e título
plt.title('Relação entre Força no Pedal e Força Longitudinal no Pneu')
plt.xlabel('Força no Pedal')
plt.ylabel('Força Longitudinal no Pneu')
plt.grid(True)
plt.legend()

# Exibindo o gráfico
plt.tight_layout()
plt.show()

