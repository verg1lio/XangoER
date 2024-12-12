import numpy as np
import matplotlib.pyplot as plt

class Pista:
    def __init__(self, gps_coords, largura_pista):
        self.gps_coords = gps_coords
        self.largura_pista = largura_pista

    def calcular_segmentos(self):
        segmentos = []
        for i in range(len(self.gps_coords) - 1):
            ponto1 = self.gps_coords[i]
            ponto2 = self.gps_coords[i + 1]
            distancia = np.sqrt((ponto2[0] - ponto1[0])**2 + (ponto2[1] - ponto1[1])**2)
            segmentos.append(distancia)
        return segmentos

    def plotar_pista(self):
        latitudes = [coord[0] for coord in self.gps_coords]
        longitudes = [coord[1] for coord in self.gps_coords]
        plt.plot(longitudes, latitudes, marker="o")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Mapa da Pista")
        plt.grid(True)
        plt.show()

# Função para simular a pista
def simulate_track(segments, max_acceleration, max_deceleration, time_step=1):
    """
    Simula a aceleração, desaceleração e velocidade constante ao longo de uma pista.

    Parâmetros:
    - segments: Lista de dicionários, cada um contendo 'length' (comprimento do segmento)
                e 'speed' (velocidade alvo para o segmento).
    - max_acceleration: Aceleração máxima permitida (positiva).
    - max_deceleration: Desaceleração máxima permitida (positiva).
    - time_step: Intervalo de tempo para a simulação (em segundos).

    Retorna:
    - time_vector: Vetor de tempo total.
    - accel_vector: Vetor de aceleração em cada instante.
    """
    time_vector = []
    accel_vector = []

    current_speed = 0
    current_time = 0

    for segment in segments:
        segment_length = segment["length"]
        target_speed = segment["speed"]

        # Acelerar ou desacelerar até a velocidade alvo
        delta_speed = target_speed - current_speed
        acceleration = max_acceleration if delta_speed > 0 else -max_deceleration

        accel_time = abs(delta_speed) / abs(acceleration)
        accel_distance = 0.5 * abs(acceleration) * accel_time**2

        if accel_distance > segment_length:
            accel_time = np.sqrt(2 * segment_length / abs(acceleration))
            accel_distance = segment_length

        accel_steps = int(accel_time / time_step)
        for step in range(accel_steps):
            accel_vector.append(acceleration)
            time_vector.append(current_time)
            current_time += time_step

        # Atualiza a velocidade atual após aceleração
        current_speed = target_speed

        # Manter velocidade constante no restante do segmento
        remaining_distance = segment_length - accel_distance
        if remaining_distance > 0:
            cruise_time = remaining_distance / current_speed
            cruise_steps = int(cruise_time / time_step)
            for step in range(cruise_steps):
                accel_vector.append(0)  # Velocidade constante implica aceleração zero
                time_vector.append(current_time)
                current_time += time_step

    return time_vector, accel_vector

gps_coords = [
    (-22.73895, -47.53324),
    (-22.7387, -47.53312),
    (-22.7377, -47.53264),
    (-22.7376, -47.53259),
    (-22.73752, -47.53251),
    (-22.73746, -47.53245),
    (-22.73741, -47.53236),
    (-22.73739, -47.53226),
    (-22.73739, -47.53213),
    (-22.73742, -47.53202),
    (-22.73746, -47.53192),
    (-22.73762, -47.53171),
    (-22.73762, -47.5317),
    (-22.73776, -47.53157),
    (-22.73803, -47.53131),
    (-22.73816, -47.53118),
    (-22.73829, -47.53106),
    (-22.73842, -47.53098),
    (-22.73855, -47.53093),
    (-22.73863, -47.53091),
    (-22.73865, -47.5309),
    (-22.73875, -47.5309),
    (-22.73887, -47.53093),
    (-22.73902, -47.53097),
    (-22.73912, -47.53104),
    (-22.73914, -47.53106),
    (-22.7392, -47.53113),
    (-22.73927, -47.53122),
    (-22.73935, -47.53139),
    (-22.7394, -47.53195),
    (-22.73942, -47.53221),
    (-22.73945, -47.53247),
    (-22.73944, -47.53259),
    (-22.74001, -47.53296),
    (-22.74006, -47.53299),
    (-22.74068, -47.53339),
    (-22.7407, -47.5334),
    (-22.74073, -47.53344),
    (-22.74075, -47.53354),
    (-22.74076, -47.53362),
    (-22.74076, -47.53372),
    (-22.74075, -47.53374),
    (-22.74074, -47.53376),
    (-22.7407, -47.53383),
    (-22.74066, -47.53388),
    (-22.74061, -47.53393),
    (-22.74057, -47.53394),
    (-22.7405, -47.53395),
    (-22.74043, -47.53395),
    (-22.74034, -47.53393),
    (-22.7403, -47.53391),
    (-22.73969, -47.5336),
    (-22.73895,-47.53324)
]

pista = Pista(gps_coords, largura_pista=10)
pista.plotar_pista()

segments = [
    {"type": "straight", "length": 300, "speed": 30},
    {"type": "curve", "length": 100, "speed": 15},
    {"type": "straight", "length": 400, "speed": 35},
    {"type": "curve", "length": 150, "speed": 20},
    {"type": "straight", "length": 250, "speed": 25},
]

# Parâmetros do carro
max_acceleration = 5.5  # m/s^2
max_deceleration = -4.5  # m/s^2

# Simular a pista
time_vector, accel_vector = simulate_track(segments, max_acceleration, max_deceleration, time_step=1)

# Substituir valores no modelo de temperatura
vector = np.array(accel_vector)
i = j = k = tempcounter = otimecounter = laps = 0
timevector = []
Temp = []

Ti = 32  # Temperatura inicial (°C)
Tinf = 30  # Temperatura ambiente (°C)
m = 22 * 10**-3  # Massa do cobre (kg)
c = 385  # Calor específico do cobre (J/kg-K)
h_bar = 10  # Coeficiente de transferência de calor (W/m^2-K)
As = 0.791 * 0.838 * (0.0254)**2 * 2
qbat = 0.3  # Potência dissipada pela bateria (W)
Tss = qbat / (h_bar * As) + Tinf  # Temperatura de equilíbrio

tt = 2 * m * c / (h_bar * As)  # Constante de tempo

while laps <= 30:  # Executar o código para 30 voltas
    while i < len(vector):
        if i > 0:
            Ti_new = Temp[-1]

        if vector[i] > 0:  # Acelerando
            t_on_accel = abs(vector[i])
            timecounter = 0
            while j < t_on_accel:
                if i == 0 and laps == 0:
                    Temp.append((Ti - Tss) * np.exp(-timecounter / tt) + Tss)
                else:
                    Temp.append((Ti_new - Tss) * np.exp(-timecounter / tt) + Tss)

                timevector.append(otimecounter)  # Atualiza aqui
                tempcounter += 1
                otimecounter += 1  # Passo ajustado para simular minutos reais
                timecounter += 1
                j += 1
            j = 0
       
        else:  # Desacelerando
            t_on_deccel = abs(vector[i])
            timecounter = 0
            while k < t_on_deccel:
                if i == 0 and laps == 0:
                    Temp.append((Ti - Tinf) * np.exp(-h_bar * As / (m * c) * timecounter) + Tinf)
                else:
                    Temp.append((Ti_new - Tinf) * np.exp(-h_bar * As / (m * c) * timecounter) + Tinf)

                timevector.append(otimecounter)  # Atualiza aqui
                tempcounter += 1
                otimecounter += 1  # Passo ajustado para simular minutos reais
                timecounter += 1
                k += 1
            k = 0
            
        i += 1
    i = 0
    laps += 1

# Plotar o resultado da temperatura
timevector_min = np.array(timevector) / 60
plt.figure(figsize=(10, 6))
plt.plot(timevector_min, Temp, label="Temperatura da bateria")
plt.xlabel("Tempo (min)")
plt.ylabel("Temperatura (°C)")
plt.title("Evolução da Temperatura da Bateria")
plt.grid(True)
plt.legend()
plt.show()

# Gerar o gráfico para diferentes valores de qbat
qbat_values = [0.2, 0.3, 0.4, 0.5, 0.6]  # Diferentes valores de qbat
time = np.arange(0, 5400, 1)  # Tempo de 0 a 90 minutos (5400 segundos)

plt.figure(figsize=(10, 6))
for qbat in qbat_values:
    Tss = qbat / (h_bar * As) + Tinf  # Temperatura em regime permanente
    tt = 2 * m * c / (h_bar * As)  # Constante de tempo
    Temp_qbat = (Ti - Tss) * np.exp(-time / tt) + Tss  # Equação para a temperatura
    plt.plot(time / 60, Temp_qbat, label=f"qbat = {qbat} W")  # Converter tempo para minutos

# Configurar o gráfico de qbat
plt.xlabel('Tempo (min)')
plt.ylabel('Temperatura (°C)')
plt.title('Curvas Analíticas para Diferentes Valores de qbat')
plt.legend()
plt.grid(True)
plt.show()