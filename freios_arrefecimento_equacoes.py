import numpy as np
import matplotlib.pyplot as plt

# Dicionário contendo as propriedades dos materiais
materiais = {
    'aluminio': {'densidade': 2600, 'calor_especifico': 950},
    'cobre': {'densidade': 8960, 'calor_especifico': 385},
    'aço': {'densidade': 7850, 'calor_especifico': 500},
    'titânio': {'densidade': 4500, 'calor_especifico': 522},
    'latão': {'densidade': 8530, 'calor_especifico': 380},
    'polietileno': {'densidade': 950, 'calor_especifico': 1900}
}

# Parâmetros do veículo e do disco de freio
massa_veiculo = 285  # kg
velocidade_inicial = 32 * 1000 / 3600  # Convertendo de km/h para m/s (32 km/h)
momento_inercia_disco = 0.0000502687  # m^4
diametro_disco = 0.15  # m
espessura_disco = 0.004  # m
condutividade_termica_aco = 169200  # W/mK
tempo_resfriamento = 90 / 3600  # Convertendo de segundos para horas

# Cálculo da energia gerada pela frenagem
energia_frenagem = (massa_veiculo * velocidade_inicial**2) / 2  # Nm

# Cálculo da potência média gerada na frenagem
aceleracao = 9.81  # m/s² (considerando desaceleração constante)
tempo_frenagem = velocidade_inicial / aceleracao
tempo = np.linspace(0, tempo_frenagem, 100)
potencia_frenagem = (massa_veiculo * aceleracao * velocidade_inicial) * (1 - (tempo / tempo_frenagem))

# Cálculo do fluxo de calor inicial no disco
area_varrida = np.pi * (diametro_disco / 2)**2
fluxo_calor_inicial = potencia_frenagem[0] / area_varrida  # Nm/h/m²

# Cálculo do ganho de temperatura para cada material em sucessivas frenagens
numero_frenagens = np.arange(1, 31)  # Número de frenagens para o gráfico

# Dicionário para armazenar o ganho de temperatura de cada material
ganho_temperatura = {}

plt.figure(figsize=(12, 8))
for material, propriedades in materiais.items():
    densidade = propriedades['densidade']
    calor_especifico = propriedades['calor_especifico']
    delta_temperatura = (fluxo_calor_inicial * tempo_resfriamento) / (densidade * calor_especifico * (area_varrida * espessura_disco))
    
    temperatura_sucessivas = [
        (1 - np.exp(-n * condutividade_termica_aco * area_varrida * tempo_resfriamento / (densidade * calor_especifico * area_varrida * espessura_disco))) * delta_temperatura
        for n in numero_frenagens
    ]
    
    ganho_temperatura[material] = temperatura_sucessivas  # Armazenando os ganhos de temperatura por material
    
    plt.plot(numero_frenagens, temperatura_sucessivas, label=material, marker='o', markersize=3)  # Adicionando bolinhas

# Configuração do gráfico
plt.xlabel("Número de Frenagens")
plt.ylabel("Ganho de Temperatura (K)")
plt.title("Ganho de Temperatura x Número de Frenagens")
plt.legend()
plt.grid(True)
plt.show()

# Cálculo do resfriamento com a Lei de Resfriamento de Newton
k = 0.1  # Constante de resfriamento (ajustável)
temperaturas_resfriamento = []
temperatura_ambiente = 25  # Temperatura ambiente em graus Celsius

# Simulação do resfriamento para cada material
tempos_resfriamento = np.linspace(0, 3600, 1000)  # 1 hora em segundos
for material in materiais.keys():
    temperatura_maxima = ganho_temperatura[material][-1] + temperatura_ambiente  # T é a temperatura máxima
    temperatura_resfriamento = [
        temperatura_ambiente + (temperatura_maxima - temperatura_ambiente) * np.exp(-k * t)  # Lei de resfriamento de Newton
        for t in tempos_resfriamento
    ]
    temperaturas_resfriamento.append(temperatura_resfriamento)

# Gráfico de resfriamento
plt.figure(figsize=(12, 8))
for idx, (material, temp_resfriamento) in enumerate(zip(materiais.keys(), temperaturas_resfriamento)):
    plt.plot(tempos_resfriamento, temp_resfriamento, label=f'Resfriamento {material}', marker='x', markersize=3)

# Configuração do gráfico de resfriamento
plt.xlabel("Tempo (s)")
plt.ylabel("Temperatura (°C)")
plt.title("Resfriamento do Disco de Freio ao Longo do Tempo")
plt.legend()
plt.grid(True)
plt.show()
