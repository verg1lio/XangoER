import math

class DiscoFreio:
    def __init__(self):
        # Parâmetros para calcular a energia térmica dissipada pelo disco
        self.k_energia = 1.5  # Fator de correção para massas em rotação em adm
        self.peso_veiculo_energia = 250  # Peso do veículo em N
        self.gravidade = 9.81  # Aceleração da gravidade em m/s²
        self.velocidade_inicial_energia = 0  # Velocidade inicial em m/s
        self.velocidade_final_energia = 10.8  # Velocidade final em m/s

        # Parâmetros para calcular a potência média de frenagem
        self.k_potencia = 1.3  # Fator de correção em adm
        self.massa_veiculo = 250  # Massa do veículo em kg
        self.desaceleracao = 0.8  # Desaceleração em m/s²
        self.velocidade_inicial_potencia = 10.8  # Velocidade inicial em m/s

        # Parâmetros para calcular o fluxo de calor que entra nos dois lados do disco
        self.k_fluxo = 1.3  # Fator de correção para massas em rotação
        self.deslizamento = 0.10  # Deslizamento do pneu (10%)
        self.velocidade_inicial_fluxo = 10.8  # Velocidade inicial do veículo em m/s
        self.desaceleracao_fluxo = 7.85  # Desaceleração do veículo em m/s^2
        self.peso_veiculo_fluxo = 2452.5  # Peso do veículo em N

        # Parâmetros para calcular a temperatura máxima atingida pelo disco
        self.pi = math.pi
        self.Dext = 0.160  # Diâmetro externo do disco [m]
        self.Dint = 0.130  # Diâmetro interno do disco [m]
        self.q0_linha = 329905.414  # Fluxo de calor que entra nos dois lados do disco [W/m²]
        self.V0 = 10.8  # Velocidade inicial do veículo [m/s]
        self.a = 0.8 * 9.81  # Desaceleração do veículo [m/s²]
        self.rho_Disc = 8027  # Densidade do material do disco [kg/m³]
        self.c_Disc = 502  # Calor específico do material do disco [J/kg.K]
        self.k_Disc = 174.465  # Condutividade térmica do material do disco [J/h.K.m]
        self.T0 = 300  # Temperatura ambiente [K]

        # Parâmetros para calcular o número de Reynolds
        self.raio_disco = 0.085  # Raio do disco [m]
        self.viscosidade_ar = 15.89 * 10**-6  # Viscosidade cinemática do ar [m²/s]

        # Parâmetros para calcular o coeficiente de convecção
        self.kar = 94.68  # Condutividade térmica do ar [J/h.K.m]
        self.diametro_disco = 0.17  # Diâmetro do disco [m]

        # Parâmetros para calcular a área total do disco e o volume do disco
        self.D_disc = 0.17  # Diâmetro do disco
        self.D_cubo = 0.06  # Diâmetro do cubo
        self.L_disc = 0.005  # Comprimento do disco

        # Parâmetros para calcular a temperatura máxima para n freadas
        self.Delta_T = 45.76 - 27  # Variação da temperatura máxima [K]
        self.A_Disc = 0.0198706  # Área do disco [m²]
        self.h_Disc = 162126.012  # Coeficiente de convecção do disco [J/h.K.m²]
        self.V_Disc = 9.9353e-5  # Volume do disco [m³]
        self.t_c = 88 / 3600  # Tempo de resfriamento do disco [h]
        self.na = 1440  # Número de freadas
        self.rho_Disc_temp = 8027  # Densidade do material do disco [kg/m³]
        self.c_Disc_temp = 502  # Calor específico do material do disco [J/kg.K]
        self.T0_temp = 300  # Temperatura ambiente [K]

    def calcular_energia_termica(self):
        return self.k_energia * self.peso_veiculo_energia * self.gravidade * (self.velocidade_inicial_energia ** 2 - self.velocidade_final_energia ** 2) / 2

    def calcular_potencia_media_frenangem(self):
        return self.k_potencia * self.massa_veiculo * (self.desaceleracao * self.gravidade) * self.velocidade_inicial_potencia / 2

    def calcular_fluxo_calor(self):
        return self.k_fluxo * (1 - self.deslizamento) * self.velocidade_inicial_fluxo * self.desaceleracao_fluxo * self.peso_veiculo_fluxo * (3600) / (778 * 3.412)

    def calcular_temperatura_maxima_disco(self):
        Ts = self.V0 / self.a
        Tmax = math.sqrt(5 / 18) * self.q0_linha * math.sqrt(Ts) / math.sqrt(self.rho_Disc * self.c_Disc * self.k_Disc) + self.T0
        return Tmax

    def calcular_numero_reynolds(self):
        return (self.velocidade_inicial_fluxo / self.raio_disco) * (self.raio_disco ** 2) / self.viscosidade_ar

    def calcular_coeficiente_conveccao(self):
        return 0.7 * self.kar / self.diametro_disco * self.calcular_numero_reynolds() ** 0.55

    def calcular_area_total_disco(self):
        area_total = (self.pi / 4) * (self.D_disc ** 2 - self.D_cubo ** 2)
        return area_total

    def calcular_volume_disco(self):
        volume = self.calcular_area_total_disco() * self.L_disc
        return volume

    def calcular_temperatura_maxima_freadas(self):
        numerator = self.Delta_T * (1 - math.exp(-self.na * self.h_Disc * self.A_Disc * self.t_c))
        denominator = (self.rho_Disc_temp * self.c_Disc_temp * self.V_Disc) / self.Delta_T * (1 - math.exp(-self.h_Disc * self.A_Disc * self.t_c))
        T_max_na = self.T0_temp + numerator / denominator
        return T_max_na

# Instanciar o objeto DiscoFreio
disco = DiscoFreio()

# Calcular e imprimir os resultados
print("Energia Térmica Dissipada pelo Disco:", disco.calcular_energia_termica(), "J")
print("Potência Média de Frenagem:", disco.calcular_potencia_media_frenangem(), "W")
print("Fluxo de Calor que Entra nos Dois Lados do Disco:", disco.calcular_fluxo_calor(), "W/h")
print("Temperatura Máxima Atingida pelo Disco:", disco.calcular_temperatura_maxima_disco(), "K")
print("Número de Reynolds:", disco.calcular_numero_reynolds())
print("Coeficiente de Convecção do Disco:", disco.calcular_coeficiente_conveccao(), "J/h.K.m²")
print("Área total do disco:", disco.calcular_area_total_disco(), "m²")
print("Volume do disco:", disco.calcular_volume_disco(), "m³")
print("Temperatura máxima para", disco.na, "freadas:", disco.calcular_temperatura_maxima_freadas(), "K")

import matplotlib.pyplot as plt

class DiscoFreioPerfurado:
    def __init__(self, potencia_media_frenagem, coef_convec, area_total_disco,
                 fator_freio, temp_maxima_freadas, temp_ambiente, temp_limite):
        self.potencia_media_frenagem = potencia_media_frenagem
        self.coef_convec = coef_convec
        self.area_total_disco = area_total_disco
        self.fator_freio = fator_freio
        self.temp_maxima_freadas = temp_maxima_freadas
        self.temp_ambiente = temp_ambiente
        self.temp_limite = temp_limite
        self.temperatura = temp_ambiente

    def simular_frenagens(self, num_frenagens, tempo_frenagem):
        temperaturas = []
        tempos = []
        for _ in range(num_frenagens):
            tempo = 0
            while tempo < tempo_frenagem:
                temperatura_media = self.temp_maxima_freadas - ((self.temp_maxima_freadas - self.temp_ambiente) * self.fator_freio)
                taxa_resfriamento = -self.coef_convec * self.area_total_disco * (self.temperatura - self.temp_ambiente)
                if self.temperatura > self.temp_limite:
                    taxa_resfriamento *= 2
                self.temperatura += taxa_resfriamento * tempo_frenagem
                if self.temperatura < self.temp_ambiente:
                    self.temperatura = self.temp_ambiente
                tempo += 1
                tempos.append(tempo)
                temperaturas.append(self.temperatura)
            self.temperatura = self.temp_maxima_freadas
        return tempos, temperaturas

# Parâmetros do disco de freio
potencia_media_frenagem = 13773.24  # W
coef_convec = 162126.012  # J/(h*K*m^2)
area_total_disco = 0.019870573  # m^2
fator_freio = 0.9
temp_maxima_freadas = 307.69  # Kelvin
temp_ambiente = 300  # Kelvin
temp_limite = 350  # Kelvin

# Instanciar o objeto DiscoFreio com os parâmetros fornecidos
disco = DiscoFreioPerfurado(potencia_media_frenagem, coef_convec, area_total_disco,
                   fator_freio, temp_maxima_freadas, temp_ambiente, temp_limite)

# Simular frenagens e obter dados de temperatura
num_frenagens = 1440
tempo_frenagem = 10  # segundos
tempos, temperaturas = disco.simular_frenagens(num_frenagens, tempo_frenagem)

# Plotar o gráfico da temperatura do disco de freio ao longo do tempo
plt.plot(tempos, temperaturas)
plt.title('Temperatura do Disco de Freio ao Longo do Tempo')
plt.xlabel('Tempo (s)')
plt.ylabel('Temperatura (K)')
plt.grid(True)
plt.show()

import matplotlib.pyplot as plt
import numpy as np

# Parâmetros do disco
diametro_disco = 0.17  # metros
espessura_disco = 0.08  # metros
numero_perfuracoes = 12
diametro_perfuracao = 0.02  # metros

# Função para criar as coordenadas de um círculo
def criar_circulo(raio, num_pontos=100):
    t = np.linspace(0, 2 * np.pi, num_pontos)
    x = raio * np.cos(t)
    y = raio * np.sin(t)
    return x, y

# Função para criar as coordenadas de um disco perfurado
def criar_disco_perfurado(diametro_disco, espessura_disco, numero_perfuracoes, diametro_perfuracao):
    fig, ax = plt.subplots()
    
    # Disco externo
    disco_exterior = plt.Circle((0, 0), diametro_disco/2, color='purple', alpha=0.5)
    ax.add_artist(disco_exterior)
    
    # Perfurações
    angulo = 2 * np.pi / numero_perfuracoes
    raio_perfuracao = diametro_perfuracao / 2
    for i in range(numero_perfuracoes):
        x_offset = (diametro_disco/2 - raio_perfuracao) * np.cos(i * angulo)
        y_offset = (diametro_disco/2 - raio_perfuracao) * np.sin(i * angulo)
        perfuracao = plt.Circle((x_offset, y_offset), raio_perfuracao, color='white')
        ax.add_artist(perfuracao)
    
    # Configurações do gráfico
    ax.set_xlim(-diametro_disco/2, diametro_disco/2)
    ax.set_ylim(-diametro_disco/2, diametro_disco/2)
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('X (metros)')
    ax.set_ylabel('Y (metros)')
    ax.set_title('Disco Perfurado')
    ax.grid(True)
    ax.axhline(0, color='black',linewidth=0.5)
    ax.axvline(0, color='black',linewidth=0.5)
    
    # Exibir
    plt.show()

# Criar o disco perfurado
criar_disco_perfurado(diametro_disco, espessura_disco, numero_perfuracoes, diametro_perfuracao)
    
