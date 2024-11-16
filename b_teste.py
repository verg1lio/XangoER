import matplotlib.pyplot as plt 

class Bateria:
    def __init__(self, tipo_celula, num_celulas_serie, num_celulas_paralelo, soh_inicial=1, soc_inicial=1, resistencia_inicial=0.02, temperatura_ambiente=25):
        # Definir os parâmetros da célula com base no tipo
        parametros = self.definir_tipo_celula(tipo_celula)
        self.tensao_celula = parametros['tensao']
        self.capacidade_celula = parametros['capacidade']
        self.resistencia_inicial = resistencia_inicial
        
        # Definir o pack de baterias
        self.num_celulas_serie = num_celulas_serie
        self.num_celulas_paralelo = num_celulas_paralelo
        
        self.tensao_nominal = self.calcular_tensao_total()
        self.capacidade_total = self.calcular_capacidade_total()
        self.soh = soh_inicial  # Estado de Saúde inicial
        self.soc = soc_inicial  # Estado de Carga inicial
        
        # Parâmetros iniciais de operação
        self.resistencia_interna = self.resistencia_inicial
        self.temperatura_ambiente = temperatura_ambiente
        self.temperatura_pack = temperatura_ambiente
        
        # Para armazenamento de dados de plot
        self.tempo = []
        self.tensao_hist = []
        self.corrente_hist = []
        self.potencia_hist = []
        self.temperatura_hist = []
        self.soc_hist = []
        self.soh_hist = []
        self.resistencia_hist = []

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
        tensao = self.tensao_nominal - queda_tensao
        self.tensao_hist.append(tensao)
        return tensao

    def calcular_corrente(self, potencia):
        tensao = self.calcular_tensao(potencia / self.capacidade_total)
        corrente = potencia / tensao
        self.corrente_hist.append(corrente)
        return corrente

    def calcular_potencia(self, corrente):
        tensao = self.calcular_tensao(corrente)
        potencia = tensao * corrente
        self.potencia_hist.append(potencia)
        return potencia

    def calcular_temperatura(self, corrente):
        aquecimento = (corrente ** 2) * self.resistencia_interna
        self.temperatura_pack = self.temperatura_ambiente + aquecimento
        self.temperatura_hist.append(self.temperatura_pack)
        return self.temperatura_pack

    def calcular_soc(self, corrente, delta_tempo):
        energia_consumida = corrente * delta_tempo / 3600  # Energia em Ah
        self.soc -= (energia_consumida / self.capacidade_total)
        self.soc = max(0, min(self.soc, 1))  # Limita o SoC entre 0 e 1
        self.soc_hist.append(self.soc)
        return self.soc

    def calcular_soh(self, ciclos, fator_degradacao):
        self.soh *= (1 - ciclos * fator_degradacao)
        self.soh_hist.append(self.soh)
        return max(0, self.soh)

    def calcular_resistencia_interna(self):
        fator_soh = 1 + (1 - self.soh) * 0.2
        fator_temperatura = 1 + (self.temperatura_pack - 25) * 0.01
        self.resistencia_interna = self.resistencia_inicial * fator_soh * fator_temperatura
        self.resistencia_hist.append(self.resistencia_interna)
        return self.resistencia_interna

    def plot_resultados(self):
        plt.figure(figsize=(12, 8))

        # Plot de Tensão
        plt.subplot(2, 3, 1)
        plt.plot(self.tempo, self.tensao_hist, label="Tensão (V)")
        plt.xlabel("Tempo (s)")
        plt.ylabel("Tensão (V)")
        plt.title("Tensão vs. Tempo")
        plt.grid(True)

        # Plot de Corrente
        plt.subplot(2, 3, 2)
        plt.plot(self.tempo, self.corrente_hist, label="Corrente (A)", color='orange')
        plt.xlabel("Tempo (s)")
        plt.ylabel("Corrente (A)")
        plt.title("Corrente vs. Tempo")
        plt.grid(True)

        # Plot de Potência
        plt.subplot(2, 3, 3)
        plt.plot(self.tempo, self.potencia_hist, label="Potência (W)", color='green')
        plt.xlabel("Tempo (s)")
        plt.ylabel("Potência (W)")
        plt.title("Potência vs. Tempo")
        plt.grid(True)

        # Plot de Temperatura
        plt.subplot(2, 3, 4)
        plt.plot(self.tempo, self.temperatura_hist, label="Temperatura (°C)", color='red')
        plt.xlabel("Tempo (s)")
        plt.ylabel("Temperatura (°C)")
        plt.title("Temperatura vs. Tempo")
        plt.grid(True)

        # Plot de SoC
        plt.subplot(2, 3, 5)
        plt.plot(self.tempo, self.soc_hist, label="SoC (%)", color='purple')
        plt.xlabel("Tempo (s)")
        plt.ylabel("SoC (%)")
        plt.title("Estado de Carga vs. Tempo")
        plt.grid(True)

        # Plot de SoH
        plt.subplot(2, 3, 6)
        plt.plot(self.tempo, self.soh_hist, label="SoH (%)", color='brown')
        plt.xlabel("Tempo (s)")
        plt.ylabel("SoH (%)")
        plt.title("Estado de Saúde vs. Ciclos")
        plt.grid(True)

        plt.tight_layout()
        plt.show()
    
    def print_resultados(self):
        print("Resultados da Simulação:")
        print(f"Tempo (s): {self.tempo}")
        print(f"Tensão (V): {self.tensao_hist}")
        print(f"Corrente (A): {self.corrente_hist}")
        print(f"Potência (W): {self.potencia_hist}")
        print(f"Temperatura (°C): {self.temperatura_hist}")
        print(f"SoC (%): {self.soc_hist}")
        print(f"SoH (%): {self.soh_hist}")
        print(f"Resistência (Ω): {self.resistencia_hist}")

    # Função para simular por determinado tempo
    def simular(self, tempo_total, delta_tempo, potencia_aplicada):
        for t in range(0, tempo_total, delta_tempo):
            corrente = self.calcular_corrente(potencia_aplicada)
            potencia = self.calcular_potencia(corrente)  
            temperatura = self.calcular_temperatura(corrente)
            soc = self.calcular_soc(corrente, delta_tempo)
            resistencia = self.calcular_resistencia_interna()
            
            # Adiciona o tempo atual e todos os históricos no mesmo ponto do tempo
            self.tempo.append(t)
            self.tensao_hist.append(self.calcular_tensao(corrente))  # Adiciona a tensão ao histórico aqui
            self.corrente_hist.append(corrente)
            self.potencia_hist.append(potencia)
            self.temperatura_hist.append(temperatura)
            self.soc_hist.append(soc)
            self.soh_hist.append(self.soh)  # Se a SoH for calculada a cada iteração, adicione aqui.
            self.resistencia_hist.append(resistencia)
        
        self.plot_resultados()  # Gera gráficos ao final da simulação

    @staticmethod
    def exemplo():
        bateria_18650 = Bateria(tipo_celula='Li-ion', num_celulas_serie=108, num_celulas_paralelo=6)
        potencia = 5000  # Potência em watts (W)
        bateria_18650.simular(tempo_total=3600, delta_tempo=10, potencia_aplicada=potencia)
        bateria_18650.print_resultados()

# Executa o exemplo
Bateria.exemplo()
