import numpy as np

# Dados fornecidos
matriz_dados = [
    {"rpm": 375, "ptc": 9.645583738270778, "trq": 245.61},
    {"rpm": 750, "ptc": 19.252680965147455, "trq": 245.12},
    {"rpm": 1125, "ptc": 28.866061704088473, "trq": 245.01},
    {"rpm": 1500, "ptc": 38.46766085790885, "trq": 244.88},
    {"rpm": 1875, "ptc": 47.988359793900806, "trq": 244.39},
    {"rpm": 2250, "ptc": 57.54597436327078, "trq": 244.22},
    {"rpm": 2625, "ptc": 67.11497779825737, "trq": 244.14},
    {"rpm": 3000, "ptc": 76.65570542895443, "trq": 243.99},
    {"rpm": 3375, "ptc": 86.08922063505362, "trq": 243.57},
    {"rpm": 3750, "ptc": 95.60756325402146, "trq": 243.45},
    {"rpm": 4125, "ptc": 105.02144248491958, "trq": 243.11},
    {"rpm": 4500, "ptc": 114.40861678954424, "trq": 242.77},
    {"rpm": 4875, "ptc": 123.84566647117964, "trq": 242.58},
    {"rpm": 5250, "ptc": 133.2403024463807, "trq": 242.34},
    {"rpm": 5625, "ptc": 142.61019709282843, "trq": 242.09},
    {"rpm": 6000, "ptc": 152.06727546916892, "trq": 242.01},
    {"rpm": 6375, "ptc": 161.53809902815016, "trq": 241.96},
    {"rpm": 6750, "ptc": 170.60206518096516, "trq": 241.34},
    {"rpm": 7125, "ptc": 178.6025469168901, "trq": 239.36},
    {"rpm": 7500, "ptc": 186.73026977211796, "trq": 237.74},
    {"rpm": 7875, "ptc": 194.7307515080429, "trq": 236.12},
    {"rpm": 8250, "ptc": 203.42477588806972, "trq": 235.45},
    {"rpm": 8625, "ptc": 210.73839121146113, "trq": 233.31},
    {"rpm": 9000, "ptc": 217.41265918230565, "trq": 230.67},
    {"rpm": 9375, "ptc": 221.6508880697051, "trq": 225.76},
    {"rpm": 9750, "ptc": 224.86019185656838, "trq": 220.22},
    {"rpm": 10125, "ptc": 227.2738459282842, "trq": 214.34},
    {"rpm": 10500, "ptc": 228.22501256702415, "trq": 207.55},
    {"rpm": 10875, "ptc": 229.2806425938338, "trq": 201.32},
    {"rpm": 11250, "ptc": 226.73660564678286, "trq": 192.45},
    {"rpm": 11625, "ptc": 219.9531616538204, "trq": 180.67},
    {"rpm": 12000, "ptc": 218.2263739946381, "trq": 173.65},
    {"rpm": 12375, "ptc": 209.6627324899464, "trq": 161.78},
    {"rpm": 12750, "ptc": 198.2173152647453, "trq": 148.45},
    {"rpm": 13125, "ptc": 189.38112642426276, "trq": 137.78},
    {"rpm": 13500, "ptc": 179.38170241286863, "trq": 126.88},
    {"rpm": 13875, "ptc": 170.95276369805632, "trq": 117.65},
    {"rpm": 14250, "ptc": 163.32104557640753, "trq": 109.44},
    {"rpm": 14625, "ptc": 152.94618171916892, "trq": 99.86},
    {"rpm": 15000, "ptc": 137.40470006702415, "trq": 87.47}
]

# Função para calcular o torque no semieixo com uma pequena perda constante
def calcular_torque_motor_semieixo(rpm, trq_motor):
    # Suponha que o torque no semieixo seja uma fração do torque do motor
    trq_semieixo = trq_motor * redp * red1 * 1.65 * 1   # Variação entre 94% e 96%
    return trq_semieixo

# Função para calcular a eficiência baseada na potência
def calcular_eficiencia_potencia(rpm, trq_motor, trq_semieixo):
    # Converte RPM para rad/s
    omega_motor = rpm * (2 * np.pi / 60)
    omega_semieixo = omega_motor  # Assumindo que a relação de transmissão é 1:1 para simplificação
    
    # Potência de entrada e saída
    ptc_motor = omega_motor * trq_motor
    ptc_semieixo = omega_semieixo * trq_semieixo
    
    # Eficiência como a razão entre a potência de saída e a potência de entrada
    eficiencia = (trq_semieixo / trq_motor) * 10
    return eficiencia

# Calcula o torque do motor, torque no semieixo e eficiência para cada ponto da matriz de dados
for dado in matriz_dados:
    rpm = dado["rpm"]
    trq_motor = dado["trq"]
    redp = 2.12  # redução primária
    red1 = 2.76  # redução da marcha única
    
    trq_semieixo = calcular_torque_motor_semieixo(rpm, trq_motor)
    eficiencia = calcular_eficiencia_potencia(rpm, trq_motor, trq_semieixo)
    
    # Output dos resultados
    print(f"RPM: {rpm} rpm")
    print(f"Torque do Motor: {trq_motor:.2f} Nm")
    print(f"Torque no Semieixo: {trq_semieixo:.2f} Nm")
    print(f"Eficiência: {eficiencia:.2f}%")
    print()
