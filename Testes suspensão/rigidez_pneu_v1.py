import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

class RigidezDeformacao:
    def __init__(self, x, F):
        self.x = x
        self.F = F
        self.k = None
        self.b = None
    
    def calcular_rigidez(self):
        # Calculando a rigidez (k) e a constante (b) usando regressão linear
        slope, intercept, r_value, p_value, std_err = linregress(self.x, self.F)
        self.k = slope  # A inclinação da reta é a rigidez
        self.b = intercept  # O intercepto é a constante
        return self.k, self.b
    

    def plotar_grafico(self):
        
        # Gerando os valores de F usando a fórmula F = kx + b
        F_calculado = self.k * self.x + self.b

        # Criando o gráfico
        plt.figure(figsize=(8, 6))
        plt.plot(self.x, self.F, 'o', label='Dados experimentais')
        plt.plot(self.x, F_calculado, color='b')

        # Adicionando títulos e rótulos
        plt.title('Gráfico de Carga Aplicada vs Deformação')
        plt.xlabel('Deformação (m)')
        plt.ylabel('Carga Aplicada (N)')
        plt.legend()
        plt.grid(True)
        plt.show()

# Dados de entrada: carga aplicada (F) e deformação (x)
x = np.array([0.0, 2.2, 3.6, 5.7, 7.2, 8.8, 9.7, 10.3, 11.3, 12.5, 13.6, 14.4, 15.3, 16.7, 17.7, 18.9, 20.1, 21.2, 22.7, 23.9, 25.2])  # Exemplo de deformação em metros
F = np.array([-29.41995, 323.28045, 558.96405, 921.0981, 1214.5456, 1579.17965, 1875.12015, 1893.33045, 2149.96135, 2435.2792, 2689.9101, 2817.92255, 3058.2738, 3312.9047, 3567.5356, 3778.49725, 4033.12815, 4287.75905, 4542.38995, 4785.8328, 5029.27565])  # Exemplo de carga em Newtons

