import numpy as np
import matplotlib.pyplot as plt

class Bateria:
    """
    Classe que modela o comportamento de um pack de baterias.
    Inclui modelagem de tensão, corrente, temperatura, SoC e análise temporal.
    """
    
    def __init__(self, tipo_celula, num_celulas_serie, num_celulas_paralelo, soh_inicial=1, soc_inicial=1, resistencia_inicial=0.02, temperatura_ambiente=25):
        self.parametros = self.definir_tipo_celula(tipo_celula)
        if not self.parametros:
            raise ValueError(f"Tipo de célula desconhecido: {tipo_celula}")
        
        self.tensao_celula = self.parametros['tensao']
        self.capacidade_celula = self.parametros['capacidade']
        self.resistencia_inicial = self.parametros['resistencia']
        
        self.num_celulas_serie = num_celulas_serie
        self.num_celulas_paralelo = num_celulas_paralelo
        
        self.tensao_nominal = self.calcular_tensao_total()
        self.capacidade_total = self.calcular_capacidade_total()
        self.soh = soh_inicial
        self.soc = soc_inicial
        
        self.resistencia_interna = self.resistencia_inicial
        self.temperatura_ambiente = temperatura_ambiente
        self.temperatura_pack = temperatura_ambiente
        
        # Histórico
        self.tempo = []
        self.tensao_hist = []
        self.corrente_hist = []
        self.temperatura_hist = []
        self.soc_hist = []

    def definir_tipo_celula(self, tipo_celula):
        tipos = {
            'Li-ion': {'tensao': 3.7, 'capacidade': 3.0, 'resistencia': 0.02},
            'LiFePO4': {'tensao': 3.2, 'capacidade': 2.8, 'resistencia': 0.007},
        }
        return tipos.get(tipo_celula)

    def calcular_tensao_total(self):
        return self.num_celulas_serie * self.tensao_celula

    def calcular_capacidade_total(self):
        return self.num_celulas_paralelo * self.capacidade_celula

    def calcular_tensao(self, corrente):
        queda_tensao = corrente * self.resistencia_interna
        return max(0, self.tensao_nominal - queda_tensao)

    def calcular_corrente(self, potencia):
        if self.tensao_nominal == 0:
            return 0
        return potencia / self.tensao_nominal

    def calcular_temperatura(self, corrente):
        aquecimento = (corrente ** 2) * self.resistencia_interna
        self.temperatura_pack = self.temperatura_ambiente + aquecimento
        return self.temperatura_pack

    def calcular_soc(self, corrente, delta_tempo):
        energia_consumida = corrente * delta_tempo / 3600
        self.soc -= energia_consumida / self.capacidade_total
        self.soc = max(0, min(self.soc, 1))
        return self.soc

    def simular(self, tempo_total, delta_tempo, demanda_func):
        for t in range(0, tempo_total, delta_tempo):
            potencia_aplicada = demanda_func(t)
            corrente = self.calcular_corrente(potencia_aplicada)
            tensao = self.calcular_tensao(corrente)
            temperatura = self.calcular_temperatura(corrente)
            soc = self.calcular_soc(corrente, delta_tempo)
            
            print(f"Tempo: {t}s | Potência: {potencia_aplicada}W | Corrente: {corrente:.2f}A | "
                f"Tensão: {tensao:.2f}V | SoC: {soc:.2%} | Temp: {temperatura:.2f}°C")

            self.tempo.append(t)
            self.tensao_hist.append(tensao)
            self.corrente_hist.append(corrente)
            self.temperatura_hist.append(temperatura)
            self.soc_hist.append(soc)
        
        self.plot_resultados()
        self.analisar_fft()

    def plot_resultados(self):
        plt.figure(figsize=(12, 8))

        plt.subplot(2, 2, 1)
        plt.plot(self.tempo, self.tensao_hist)
        plt.title("Tensão ao longo do tempo")
        plt.xlabel("Tempo (s)")
        plt.ylabel("Tensão (V)")
        plt.grid(True)

        plt.subplot(2, 2, 2)
        plt.plot(self.tempo, self.corrente_hist, color='orange')
        plt.title("Corrente ao longo do tempo")
        plt.xlabel("Tempo (s)")
        plt.ylabel("Corrente (A)")
        plt.grid(True)

        plt.subplot(2, 2, 3)
        plt.plot(self.tempo, self.temperatura_hist, color='red')
        plt.title("Temperatura ao longo do tempo")
        plt.xlabel("Tempo (s)")
        plt.ylabel("Temperatura (°C)")
        plt.grid(True)

        plt.subplot(2, 2, 4)
        plt.plot(self.tempo, self.soc_hist, color='purple')
        plt.title("Estado de Carga (SoC) ao longo do tempo")
        plt.xlabel("Tempo (s)")
        plt.ylabel("SoC (%)")
        plt.grid(True)

        plt.tight_layout()
        plt.show()

    def analisar_fft(self):
        """Realiza a análise de FFT para tensão, corrente, temperatura e SoC e plota os espectros de frequência."""
        
        def fft_plot(signal_data, label, color):
            """Função auxiliar para gerar o gráfico de FFT."""
            signal = np.array(signal_data)
            freq_amostragem = 1 / (self.tempo[1] - self.tempo[0])  # Hz
            N = len(signal)
            fft_signal = np.fft.fft(signal)
            freq = np.fft.fftfreq(N, d=1/freq_amostragem)

            plt.figure(figsize=(8, 4))
            plt.plot(freq[:N//2], np.abs(fft_signal[:N//2]), label=label, color=color)
            plt.title(f"Espectro de Frequência - {label}")
            plt.xlabel("Frequência (Hz)")
            plt.ylabel("Magnitude")
            plt.grid(True)
            plt.show()

        fft_plot(self.tensao_hist, 'Tensão', 'blue')
        fft_plot(self.corrente_hist, 'Corrente', 'orange')

    # Exemplo de uso
    def exemplo():
        def demanda_motor(t):
            return 500

        bateria = Bateria(tipo_celula='Li-ion', num_celulas_serie=108, num_celulas_paralelo=6)
        bateria.simular(tempo_total=600, delta_tempo=1, demanda_func=demanda_motor)

Bateria.exemplo()
