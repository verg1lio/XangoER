import numpy as np

# Penalidade de retidão vertical para um elemento definido entre dois nós
def penalidade_retidão_vertical(nodes, no1, no2):
    pen = 0 
    componente_x_vetor = nodes[(no1)][0] - nodes[(no2)][0]                          
    componente_y_vetor = nodes[(no1)][1] - nodes[(no2)][1]   
    componente_z_vetor = nodes[(no1)][2] - nodes[(no2)][2]
    vetor_no1_no2 = (componente_x_vetor, componente_y_vetor, componente_z_vetor)                          # vetor entre os nós
    modulo_vetor = np.linalg.norm(vetor_no1_no2)                    # comprimento do vetor
    if modulo_vetor == 0:
        return 0  # evita divisão por zero

    vetor_vertical = np.array([0, 0, 1])                             # vetor unitário vertical
    produto_escalar_no1_no2_vertical = np.dot(vetor_no1_no2, vetor_vertical)
    cos_theta_no1_no2_vertical = produto_escalar_no1_no2_vertical / modulo_vetor
    theta_no1_no2 = np.degrees(np.acos(cos_theta_no1_no2_vertical))
    if theta_no1_no2 > 1e-6 :
        pen += ((theta_no1_no2) ** 2) * 1e6

    return pen        


# Exemplo de dados
nodes = np.array([[0, 0, 1], [1, 1, 1]])  # dois nós na mesma altura (z = 1) vetor horizontal
pen = penalidade_retidão_vertical(nodes, 0, 1)
print("Penalidade de retidão vertical:", pen)

"""
para fins de teste aconselho alterar as cordenadas de no1 um de forma a formar vetores horizontais, oblíquo e vertical com no2

exemplos de cordenadas que formam respectivamente esses vetores [0. 0, 1] , [0, 0, 2] , [1, 1, 9]
"""

