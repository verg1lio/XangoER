import numpy as np
import matplotlib.pyplot as plt


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

        self.n_serie = n_serie
        self.n_paralelo = n_paralelo

        self.soc = soc_inicial
        self.carga_total = self.Q * n_paralelo  # Ah

        # Variáveis acumuladas
        self.Iast = 0.0  # Corrente acumulada (Coulombs)
        self.tempo_acumulado = 0.0  # Tempo total (s)

        # Histórico
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
                'volume': 0.025
            },
            'LiFePO4': {
                'E0': 3.2,
                'K': 0.003,
                'Q': 2.8,
                'A': 0.05,
                'B': 0.015,
                'R': 0.008,
                'peso': 0.07,
                'volume': 0.03
            }
        }
        return catalogo.get(tipo)

    def calcular_tensao(self, corrente, dt):
        carga_utilizada = (1 - self.soc) * self.carga_total

        Iast = self.Iast  # Corrente acumulada em Ampere-segundo

        # Acumula tempo total
        self.tempo_acumulado += dt

        if corrente >= 0:  # Descarga
            try:
                Vt = (self.E0
                      - self.R * corrente
                      - self.K * (self.carga_total / (self.carga_total - carga_utilizada)) * (carga_utilizada + Iast / 3600)
                      + self.A * np.exp(-self.B * carga_utilizada))
            except ZeroDivisionError:
                Vt = 0
        else:  # Carga
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

        # Atualiza corrente acumulada (em Coulombs)
        self.Iast += corrente * dt
        return self.soc

    def simular(self, tempo_total, dt, demanda_potencia):
        for t in range(0, tempo_total, dt):
            potencia = demanda_potencia(t)

            tensao_estimativa = self.calcular_tensao(0.1, dt)
            corrente = potencia / max(tensao_estimativa, 0.1)

            tensao = self.calcular_tensao(corrente, dt)
            self.atualizar_soc(corrente, dt)

            self.tempo.append(t)
            self.tensao_hist.append(tensao)
            self.corrente_hist.append(corrente)
            self.soc_hist.append(self.soc)

            print(f"Tempo: {t}s | Potência: {potencia:.2f}W | Corrente: {corrente:.2f}A | "
                  f"Tensão: {tensao:.2f}V | SoC: {self.soc:.2%}")

        self.plot_resultados()

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
        plt.bar(tempo, corrente, color=['red' if i > 0 else 'green' for i in corrente], width=1)
        plt.ylabel('Corrente (A)')
        plt.grid()
        plt.axhline(0, color='black', linewidth=0.8)

        plt.subplot(4, 1, 3)
        plt.plot(tempo, self.soc_hist, color='green')
        plt.ylabel('SoC (%)')
        plt.grid()

        plt.subplot(4, 1, 4)
        potencia = np.array(self.corrente_hist) * np.array(self.tensao_hist)
        plt.bar(tempo, potencia, color=['red' if p > 0 else 'green' for p in potencia], width=1)
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


# 🔋 Exemplo de uso
if __name__ == "__main__":
    def demanda(t):
        if t < 300:
            return 5000
        elif t < 500:
            return 5000
        else:
            return 5000

    bateria = BatteryPack(tipo_celula='Li-ion', n_serie=50, n_paralelo=100, soc_inicial=1)
    bateria.simular(tempo_total=13500, dt=1, demanda_potencia=demanda)

    print(f"Tensão nominal: {bateria.calcular_tensao_nominal():.2f} V")
    print(f"Peso total: {bateria.calcular_peso_total():.2f} kg")
    print(f"Volume total: {bateria.calcular_volume_total():.2f} litros")
