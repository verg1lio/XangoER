import numpy as np
import matplotlib.pyplot as plt
import math

# Dicionário contendo as propriedades dos materiais
materiais = {
    'aluminio': {'densidade': 2600, 'calor_especifico': 950},
    'cobre': {'densidade': 8960, 'calor_especifico': 385},
    'aço': {'densidade': 7850, 'calor_especifico': 500},
    'titânio': {'densidade': 4500, 'calor_especifico': 522},
    'latão': {'densidade': 8530, 'calor_especifico': 380},
    'polietileno': {'densidade': 950, 'calor_especifico': 1900}
}

class Arrefecimento:
    def __init__(self, velocidade, area, comprimento, calor_gerado, tipo_geometria='placa'):
        self.rho = 1.184
        self.c = 1007
        self.v = velocidade
        self.mi = 1.849e-5
        self.ni = self.mi / self.rho
        self.k_ar = 0.02551
        self.tf = 25
        self.Pr = self.mi * self.c / self.k_ar
        self.a_res = area
        self.L_c = comprimento                                  #Comprimento caracteristico
        self.calor_gerado = calor_gerado            
        self.tipo_geometria = tipo_geometria
        self.tipos_aletas = ['retangular', 'triangular']
        self.alturas_aletas = [0.005, 0.01]                     # Altura de aletas variadas (menores que 0.01m)
        self.espacos_entre_aletas = [0.01, 0.03]                # Distâncias entre as aletas
        self.espessura_aleta = 0.003

    def geometria_aleta(self, tipo, altura_aleta, comprimento_aleta, espessura_aleta, rho):
        if tipo == 'retangular':
            area = 2 * altura_aleta * comprimento_aleta  # Área da seção transversal
            volume = altura_aleta * comprimento_aleta * espessura_aleta
        elif tipo == 'triangular':
            area = 2 * math.sqrt(altura_aleta**2 + (espessura_aleta/2)**2) * comprimento_aleta
            volume = (espessura_aleta * altura_aleta / 2) * comprimento_aleta
        else:
            raise ValueError("Tipo de geometria não suportado. Use 'retangular' ou 'triangular'.")
        return area, volume * rho

    def numero_aletas(self, espaco):
        if self.tipo_geometria == 'placa':
            n1 = math.floor(self.L_c / (espaco + self.espessura_aleta))                         #Calculo do numero de aletas em uma face da caixa de baterias
            num_aletas = 4 * n1                                                                 #Calculo para o numero total das aletas, considerando 4 faces 
        elif self.tipo_geometria == 'cilindro':
            num_aletas = math.floor((self.L_c * math.pi) / (espaco + self.espessura_aleta))     #Calculo para o numero de aletas no motor
        else:
            raise ValueError("Tipo de geometria não suportada.")
        return num_aletas

    def reynolds(self):
        return self.rho * self.v * self.L_c / self.mi

    def calc_Nu(self):
        Rey = self.reynolds()
        if self.tipo_geometria == 'placa':
            if Rey < 200000:
                Nu = 0.332 * Rey**0.5 * self.Pr**(1/3)
            else:
                X = 200000 * self.ni / self.v
                Rey_X = X * self.v / self.ni
                A = 0.037 * Rey_X**0.8 - 0.664 * Rey_X**0.5
                Nu = 0.037 * (Rey**0.8 - A) * self.Pr**(1/3)
        elif self.tipo_geometria == 'cilindro':
            Nu = 0.3 + (0.62 * Rey**0.5 * self.Pr**(1/3)) / (1 + (0.4/self.Pr)**(2/3))**0.25
        else:
            raise ValueError("Tipo de geometria não suportada.")
        return Nu

    def arrefecimento(self, tempo, a_res, temp_inicial, c_objeto, m_objeto):
        h = (self.calc_Nu() * self.k_ar) / self.L_c
        T_final = self.tf + (self.calor_gerado / (h * a_res))
        tau = (m_objeto * c_objeto) / (h * a_res)
        return T_final + (temp_inicial - T_final) * np.exp(-tempo / tau)

    def troca_termica(self, tempo_simulacao, temp_inicial, a_res, c_objeto, m_objeto):
        tempos = np.linspace(0, tempo_simulacao, tempo_simulacao)
        temperaturas = self.arrefecimento(tempos, a_res, temp_inicial, c_objeto, m_objeto)
        return tempos, temperaturas

    def plotar_grafico(self, temp_objeto, a_res, c_objeto, m_objeto, tempo_simulacao, legenda):
        tempos, temperaturas = self.troca_termica(tempo_simulacao, temp_objeto, a_res, c_objeto, m_objeto)
        plt.plot(tempos, temperaturas, label=legenda)

    def plotar_graficos_aletas(self, temp_inicial, c_objeto, m_objeto, tempo_simulacao, objeto):
        plt.figure(figsize=(12, 8))
        # Plot do objeto sem aletas
        self.plotar_grafico(temp_inicial, self.a_res, c_objeto, m_objeto, tempo_simulacao, 'Sem aletas')
        # Plots com diferentes configurações de aletas
        for tipo in self.tipos_aletas:
            for altura in self.alturas_aletas:
                for espaco in self.espacos_entre_aletas:
                    area_aleta, massa_aleta = self.geometria_aleta(tipo, altura, self.L_c, self.espessura_aleta, rho)
                    num_aletas = self.numero_aletas(espaco)
                    area_com_aletada = self.a_res + num_aletas * area_aleta
                    massa_com_aletas = m_objeto + num_aletas * massa_aleta
                    legenda = f'Aleta {tipo}, Altura: {altura*1000:.1f}mm, Espaço: {espaco*1000:.1f}mm'
                    self.plotar_grafico(temp_inicial, area_com_aletada, c_objeto, massa_com_aletas, tempo_simulacao, legenda)

        plt.title(f'Variação da Temperatura ao Longo do Tempo ({objeto})')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Temperatura (°C)')
        plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        plt.grid(True)
        plt.tight_layout()
        plt.show()

# Função para calcular o volume da camada externa de um cilindro
def volume_camada_cilindro(raio_externo, espessura, comprimento):
    raio_interno = raio_externo - espessura
    return math.pi * comprimento * (raio_externo**2 - raio_interno**2)

# Parâmetros de simulação
velocidade = 10
calor_gerado = 100000 / 3600         # Potência convertida para J/s
temp_objeto = 60                    # Temperatura inicial do objeto (°C)
tempo_simulacao = 2500              # Tempo total da simulação (s)
material = 'aluminio'
rho = materiais[material]['densidade']
calor_especifico = materiais[material]['calor_especifico']

# Dimensões da placa e do motor
L = 0.60
W = 0.22
espessura_caixa = 0.003
area_caixa = 5 * (L * W)
volume_caixa = 6 * L * W * espessura_caixa
massa_caixa = rho * volume_caixa

# Simulação para a placa
bateria = Arrefecimento(velocidade, area_caixa, L, calor_gerado, 'placa')
bateria.plotar_graficos_aletas(temp_objeto, calor_especifico, massa_caixa, tempo_simulacao, "bateria")

comprimento_motor = 0.355
diametro_motor = 0.20
raio_externo_motor = diametro_motor / 2
espessura_motor = 0.003
area_motor = math.pi * diametro_motor * comprimento_motor
volume_motor = volume_camada_cilindro(raio_externo_motor, espessura_motor, comprimento_motor)
massa_motor = volume_motor * rho

# Simulação para o motor
motor = Arrefecimento(velocidade, area_motor, comprimento_motor, calor_gerado, 'cilindro')
motor.plotar_graficos_aletas(temp_objeto, calor_especifico, massa_motor, tempo_simulacao, "motor")

