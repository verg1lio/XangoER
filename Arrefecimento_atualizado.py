import numpy as np
import matplotlib.pyplot as plt
import math

# Dicionário contendo as propriedades dos materiais
materiais = {
    'aluminio': {'densidade': 2600, 'calor_especifico': 950},  # Densidade em kg/m^3, Calor específico em J/(kg·K)
    'cobre': {'densidade': 8960, 'calor_especifico': 385},
    'aço': {'densidade': 7850, 'calor_especifico': 500},
    'titânio': {'densidade': 4500, 'calor_especifico': 522},
    'latão': {'densidade': 8530, 'calor_especifico': 380},     
    'polietileno': {'densidade': 950, 'calor_especifico': 1900}  
}

class Arrefecimento:
    def __init__(self, velocidade, area, comprimento, calor_gerado,tipo_geometria='placa'):
        self.rho = 1.184                                        # Massa específica do ar (kg/m^3)
        self.c = 1007                                           # Calor específico do ar (J/(kg*K))
        self.v = velocidade                                     # Velocidade do ar (m/s)
        self.mi = 1.849e-5                                      # Viscosidade dinâmica do ar (kg/m*s)
        self.ni = self.mi / self.rho                            # Viscosidade cinemática do ar (m²/s)
        self.k_ar = 0.02551                                     # Condutividade térmica do ar (W/m*K)
        self.tf = 25                                            # Temperatura do ar (°C)
        self.Pr = self.mi * self.c / self.k_ar                  # Número de Prandtl
        self.a_res = area                                       # Área de resfriamento (m²)
        self.L_c = comprimento                                  # Comprimento característico (m)
        self.calor_gerado = calor_gerado                        # Calor gerado pelo objeto (W)
        self.tipo_geometria = tipo_geometria                    # Geometria do objeto ('placa' ou 'cilindro')

    def reynolds(self):
        return self.rho * self.v * self.L_c / self.mi
     
    def calc_Nu(self):
        Rey = self.reynolds()
        
        if self.tipo_geometria == 'placa':
            if Rey < 200000:  # Fluxo laminar
                Nu = 0.332 * Rey**0.5 * self.Pr**(1/3)
            else:  # Fluxo Turbulento
                X = 200000 * self.ni / self.v  # Ponto de transição para turbulento
                Rey_X = X * self.v / self.ni
                A = 0.037 * Rey_X**0.8 - 0.664 * Rey_X**0.5
                Nu = 0.037 * (Rey**0.8 - A) * self.Pr**(1/3)
        elif self.tipo_geometria == 'cilindro':
            Nu = 0.3 + (0.62 * Rey**0.5 * self.Pr**(1/3)) / (1 + (0.4/self.Pr)**(2/3))**0.25
        else:
            raise ValueError("Tipo de geometria não suportada. Escolha 'placa' ou 'cilindro'.")
        
        return Nu
    
    def calor_conveccao(self, temp_sup):
        nu = self.calc_Nu()
        h = nu * self.k_ar / self.L_c
        return h * self.a_res * (temp_sup - self.tf)

    def troca_termica(self, tempo_simulacao, temp_objeto, c_objeto, m_objeto):
        temp_atual = temp_objeto
        temperaturas = []
        tempos = []
        for i in range(tempo_simulacao):
            q_conveccao = self.calor_conveccao(temp_atual)
            q_total = self.calor_gerado - q_conveccao
            temp_atual = (q_total / (c_objeto * m_objeto)) + temp_atual
            tempos.append(i)
            temperaturas.append(temp_atual)
        return tempos, temperaturas

    def plotar_grafico(self, temp_objeto, c_objeto, m_objeto, tempo_simulacao,objeto):
        tempos, temperaturas = self.troca_termica(tempo_simulacao, temp_objeto, c_objeto, m_objeto)
        plt.plot(tempos, temperaturas)
        plt.title(f'Variação da Temperatura ao Longo do Tempo ({objeto})')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Temperatura (°C)')
        plt.grid(True)
        plt.show()

# Função para calcular o volume da camada externa de um cilindro
def volume_camada_cilindro(raio_externo, espessura, comprimento):
    raio_interno = raio_externo - espessura
    return math.pi * comprimento * (raio_externo**2 - raio_interno**2)

# Parâmetros de simulação
velocidade = 5
calor_gerado = 10000 / 3600  # Potência convertida para J/s
temp_objeto = 25  # Temperatura inicial do objeto (°C)
tempo_simulacao = 6000  # Tempo total da simulação (s)
material = 'aluminio'  # Escolha o material ('aluminio', 'cobre', 'aço', 'titânio', 'latão', 'polietileno')
rho = materiais[material]['densidade']
calor_especifico = materiais[material]['calor_especifico']

# Dimensões da placa
L = 0.60  # Comprimento (m)
W = 0.22  # Largura (m)
espessura_caixa = 0.003  # Espessura (m)
area_caixa = 5 * (L * W) + 0.3  # Área total de troca térmica (m²)
volume_caixa = 6 * L * W * espessura_caixa  # Volume da placa (m³)
massa_caixa = rho * volume_caixa  # Massa da placa (kg)

# Simulação para a placa
placa = Arrefecimento(velocidade, area_caixa, L, calor_gerado, tipo_geometria='placa')
placa.plotar_grafico(temp_objeto, calor_especifico, massa_caixa, tempo_simulacao,"bateria")

# Dimensões do motor
comprimento_motor = 0.355  # Comprimento do motor (m)
diametro_motor = 0.20  # Diâmetro do motor (m)
raio_externo_motor = diametro_motor / 2  # Raio externo do motor (m)
espessura_motor = 0.003  # Considerando apenas os 5 mm mais externos

# Área e volume da camada externa do motor
area_motor = math.pi * diametro_motor * comprimento_motor  # Área externa do cilindro (m²)
volume_motor = volume_camada_cilindro(raio_externo_motor, espessura_motor, comprimento_motor)
massa_motor = volume_motor * rho  # Massa do motor considerando 3 mm externos

# Simulação para o motor
motor = Arrefecimento(velocidade, area_motor, diametro_motor, calor_gerado, tipo_geometria='cilindro')
motor.plotar_grafico(temp_objeto, calor_especifico, massa_motor, tempo_simulacao,"motor")