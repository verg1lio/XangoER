import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio
import math


pio.renderers.default = "browser"

class BatteryPack:
    def __init__(self, tipo_celula, n_serie, n_paralelo, soc_inicial=1.0):
        self.parametros = self.definir_celula(tipo_celula)
        if self.parametros is None:
            raise ValueError(f"Tipo de célula desconhecido: {tipo_celula}")

        self.E0 = self.parametros['E0']
        self.K = self.parametros['K']
        self.Q = self.parametros['Q']
        self.A = self.parametros['A']
        self.B = self.parametros['B']
        self.R = self.parametros['R']
        self.peso = self.parametros['peso']
        self.volume = self.parametros['volume']

        self.diametro_celula = self.parametros['diametro_cm']
        self.altura_celula = self.parametros.get('altura_cm',
                                                  (self.volume * 1000) / (math.pi * (self.diametro_celula / 2)**2) if self.diametro_celula > 0 else 6.5)


        self.n_serie = n_serie
        self.n_paralelo = n_paralelo

        self.soc = soc_inicial
        self.carga_total = self.Q * n_paralelo

        self.Iast = 0.0
        self.tempo_acumulado = 0.0

        self.tensao_hist = []
        self.corrente_hist = []
        self.soc_hist = []
        self.tempo = []

    def definir_celula(self, tipo):
        catalogo = {
            'Li-ion': {
                'E0': 3.7,
                'K': 0.005,
                'Q': 3.0,
                'A': 0.1,
                'B': 0.01,
                'R': 0.02,
                'peso': 0.045,
                'volume': 0.0165,
                'diametro_cm': 1.8,
                'altura_cm': 6.5
            },
            'LiFePO4': {
                'E0': 3.2,
                'K': 0.003,
                'Q': 2.8,
                'A': 0.05,
                'B': 0.015,
                'R': 0.008,
                'peso': 0.07,
                'volume': 0.028,
                'diametro_cm': 2.6,
                'altura_cm': 6.5
            }
        }
        return catalogo.get(tipo)

    def calcular_tensao(self, corrente, dt):
        carga_utilizada = (1 - self.soc) * self.carga_total
        Iast = self.Iast
        self.tempo_acumulado += dt

        if corrente >= 0:
            try:
                Vt = (self.E0
                      - self.R * corrente
                      - self.K * (self.carga_total / (self.carga_total - carga_utilizada)) * (carga_utilizada + Iast / 3600)
                      + self.A * np.exp(-self.B * carga_utilizada))
            except ZeroDivisionError:
                Vt = 0
        else:
            corrente_abs = abs(corrente)
            try:
                Vt = (self.E0
                      - self.R * corrente_abs
                      - self.K * (self.carga_total / (corrente_abs * self.tempo_acumulado - 0.1 * self.carga_total)) * (Iast / 3600)
                      - self.K * ((self.carga_total / (self.carga_total - carga_utilizada)) * carga_utilizada)
                      + self.A * np.exp(-self.B * carga_utilizada))
            except ZeroDivisionError:
                Vt = 0

        tensao_total = Vt * self.n_serie
        return tensao_total

    def atualizar_soc(self, corrente, dt):
        delta_ah = corrente * dt / 3600
        self.soc -= delta_ah / self.carga_total
        self.soc = max(0, min(self.soc, 1))
        self.Iast += corrente * dt
        return self.soc

    def simular(self, tempo_total, dt, demanda_potencia):
        print("\n--- Iniciando Simulação do Comportamento Elétrico ---")
        for t in range(0, tempo_total, dt):
            potencia = demanda_potencia(t)
            tensao_estimativa = self.calcular_tensao(0.1, dt)
            corrente = potencia / max(tensao_estimativa, 0.1)
            corrente = max(-self.Q * self.n_paralelo * 2, min(corrente, self.Q * self.n_paralelo * 2))
            tensao = self.calcular_tensao(corrente, dt)
            self.atualizar_soc(corrente, dt)
            self.tempo.append(t)
            self.tensao_hist.append(tensao)
            self.corrente_hist.append(corrente)
            self.soc_hist.append(self.soc)
        self.plot_resultados()
        print("--- Simulação Elétrica Concluída e Gráficos Gerados ---")

    def plot_resultados(self):
        tempo = self.tempo
        plt.figure(figsize=(10, 8))
        plt.subplot(4, 1, 1)
        plt.plot(tempo, self.tensao_hist, color='blue')
        plt.ylabel('Tensão (V)')
        plt.grid()
        plt.title('Simulação - Modelo Modified Shepherd')
        plt.subplot(4, 1, 2)
        corrente = np.array(self.corrente_hist)
        plt.plot(tempo, corrente, color='red')
        plt.ylabel('Corrente (A)')
        plt.grid()
        plt.axhline(0, color='black', linewidth=0.8)
        plt.subplot(4, 1, 3)
        plt.plot(tempo, self.soc_hist, color='green')
        plt.ylabel('SoC (%)')
        plt.grid()
        plt.subplot(4, 1, 4)
        potencia = np.array(self.corrente_hist) * np.array(self.tensao_hist)
        plt.plot(tempo, potencia, color='purple')
        plt.ylabel('Potência (W)')
        plt.xlabel('Tempo (s)')
        plt.grid()
        plt.axhline(0, color='black', linewidth=0.8)
        plt.tight_layout()
        plt.show()

    def calcular_peso_total(self):
        return self.n_serie * self.n_paralelo * self.peso

    def calcular_volume_total(self):
        return self.n_serie * self.n_paralelo * self.volume

    def calcular_tensao_nominal(self):
        return self.E0 * self.n_serie

    def plot_estrutura_3d(self):
        print("\n--- Gerando Visualização 3D da Estrutura da Bateria ---")
        radius = self.diametro_celula / 2
        height = self.altura_celula

        data = []
        offset_x = self.diametro_celula * 1.5
        offset_y = self.altura_celula * 1.2

        for s in range(self.n_serie):
            for p in range(self.n_paralelo):
                center_x = p * offset_x
                center_y = s * offset_y
                center_z = height / 2

                for z_coord in [center_z - height / 2, center_z + height / 2]:
                    theta = np.linspace(0, 2 * np.pi, 50)
                    x_circle = center_x + radius * np.cos(theta)
                    y_circle = center_y + radius * np.sin(theta)
                    data.append(go.Scatter3d(x=x_circle, y=y_circle, z=z_coord * np.ones_like(theta),
                                            mode='lines', line=dict(color='gray'), showlegend=False))

                phi = np.linspace(0, 2 * np.pi, 20)
                z = np.linspace(center_z - height / 2, center_z + height / 2, 2)
                x_cyl = radius * np.outer(np.cos(phi), np.ones_like(z)) + center_x
                y_cyl = radius * np.outer(np.sin(phi), np.ones_like(z)) + center_y
                z_cyl = np.outer(np.ones_like(phi), z)

                data.append(go.Surface(x=x_cyl, y=y_cyl, z=z_cyl,
                                       colorscale=[[0, 'red'], [1, 'red']],
                                       showscale=False, opacity=0.7, showlegend=False))

        for p in range(self.n_paralelo):
            for s in range(self.n_serie - 1):
                center1_x = p * offset_x
                center1_y = s * offset_y
                center1_z = height / 2

                center2_x = p * offset_x
                center2_y = (s + 1) * offset_y
                center2_z = height / 2

                z_conn1 = center1_z + height / 2
                z_conn2 = center2_z - height / 2

                data.append(go.Scatter3d(
                    x=[center1_x, center2_x],
                    y=[center1_y + radius, center2_y - radius],
                    z=[z_conn1, z_conn2],
                    mode='lines',
                    line=dict(color='orange', width=5),
                    name='Conexão Série',
                    showlegend=True if s == 0 and p == 0 else False
                ))

        for s in range(self.n_serie):
            for p in range(self.n_paralelo - 1):
                center1_x = p * offset_x
                center1_y = s * offset_y
                center1_z_mid = height / 2

                center2_x = (p + 1) * offset_x
                center2_y = s * offset_y
                center2_z_mid = height / 2

                data.append(go.Scatter3d(
                    x=[center1_x + radius, center2_x - radius],
                    y=[center1_y, center2_y],
                    z=[center1_z_mid, center2_z_mid],
                    mode='lines',
                    line=dict(color='blue', width=5),
                    name='Conexão Paralelo',
                    showlegend=True if s == 0 and p == 0 else False
                ))

        fig = go.Figure(data=data)

        fig.update_layout(
            scene=dict(
                xaxis_title='Largura (cm)',
                yaxis_title='Profundidade (cm)',
                zaxis_title='Altura (cm)',
                aspectmode='data',
            ),
            title=f'Estrutura 3D do Pack de Baterias ({self.n_serie}S{self.n_paralelo}P)',
            showlegend=True
        )

        # No Colab, fig.show() sozinho já funciona com o renderer 'colab' configurado
        fig.show()
        print("--- Visualização 3D incorporada na saída da célula do Colab. ---")

# 🔋 Exemplo de uso
if __name__ == "__main__":
    def demanda(t):
        if t < 300:
            return 5000
        elif t < 500:
            return -3000
        else:
            return 4000

    bateria_pequena = BatteryPack(tipo_celula='Li-ion', n_serie=5, n_paralelo=3, soc_inicial=0.8)
    print(f"\n--- Detalhes da Bateria Pequena ({bateria_pequena.n_serie}S{bateria_pequena.n_paralelo}P) ---")
    print(f"Tensão nominal: {bateria_pequena.calcular_tensao_nominal():.2f} V")
    print(f"Peso total: {bateria_pequena.calcular_peso_total():.2f} kg")
    print(f"Volume total: {bateria_pequena.calcular_volume_total():.2f} litros")

    bateria_pequena.plot_estrutura_3d()

    # Opcional: Descomente a linha abaixo para simular o comportamento elétrico
    # bateria_pequena.simular(tempo_total=800, dt=1, demanda_potencia=demanda)


    bateria_grande = BatteryPack(tipo_celula='LiFePO4', n_serie=10, n_paralelo=5, soc_inicial=1.0)
    print(f"\n--- Detalhes da Bateria Grande ({bateria_grande.n_serie}S{bateria_grande.n_paralelo}P) ---")
    print(f"Tensão nominal: {bateria_grande.calcular_tensao_nominal():.2f} V")
    print(f"Peso total: {bateria_grande.calcular_peso_total():.2f} kg")
    print(f"Volume total: {bateria_grande.calcular_volume_total():.2f} litros")

    bateria_grande.plot_estrutura_3d()

    # Opcional: Descomente a linha abaixo para simular o comportamento elétrico
    # bateria_grande.simular(tempo_total=800, dt=1, demanda_potencia=demanda)