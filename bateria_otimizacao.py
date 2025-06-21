import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio
import math

pio.renderers.default = "browser" # Use "browser" if running locally

class BatteryPack:
    def __init__(self, tipo_celula, n_serie, n_paralelo, soc_inicial=1.0):
        self.parametros = self.definir_celula(tipo_celula)
        if self.parametros is None:
            raise ValueError(f"Tipo de célula desconhecida: {tipo_celula}")

        self.E0 = self.parametros['E0']
        self.K = self.parametros['K']
        self.Q = self.parametros['Q']
        self.A = self.parametros['A']
        self.B = self.parametros['B']
        self.R = self.parametros['R']
        self.peso = self.parametros['peso']
        self.volume_celula_litros = self.parametros['volume'] # Volume individual da célula em litros

        # Dimensões para a representação do cubo
        self.lado_cubo_x = self.parametros['dim_x_cm']
        self.lado_cubo_y = self.parametros['dim_y_cm']
        self.lado_cubo_z = self.parametros['dim_z_cm']

        self.n_serie = n_serie
        self.n_paralelo = n_paralelo

        self.soc = soc_inicial
        self.carga_total = self.Q * n_paralelo  # Ah

        self.Iast = 0.0
        self.tempo_acumulado = 0.0

        self.tensao_hist = []
        self.corrente_hist = []
        self.soc_hist = []
        self.tempo = []

    def definir_celula(self, tipo):
        # Assumindo que as células são retangulares ou cilíndricas aproximadas por um cuboide
        # e fornecendo dimensões de um cuboide para melhor cálculo de volume de empacotamento.
        # Dimensões em cm. Volume da célula em litros para o cálculo original.
        catalogo = {
            'Li-ion': {
                'E0': 3.7, 'K': 0.005, 'Q': 3.0, 'A': 0.1, 'B': 0.01, 'R': 0.02,
                'peso': 0.045, # kg
                'volume': 0.0165, # litros
                'dim_x_cm': 1.8, # Largura do cuboide (aprox. diâmetro 18650)
                'dim_y_cm': 1.8, # Profundidade do cuboide (aprox. diâmetro 18650)
                'dim_z_cm': 6.5  # Altura do cuboide (altura 18650)
            },
            'LiFePO4': {
                'E0': 3.2, 'K': 0.003, 'Q': 2.8, 'A': 0.05, 'B': 0.015, 'R': 0.008,
                'peso': 0.07, # kg
                'volume': 0.028, # litros
                'dim_x_cm': 2.6, # Largura do cuboide (aprox. diâmetro 26650)
                'dim_y_cm': 2.6, # Profundidade do cuboide (aprox. diâmetro 26650)
                'dim_z_cm': 6.5  # Altura do cuboide (altura 26650)
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

    def calcular_volume_total_litros(self):
        return self.n_serie * self.n_paralelo * self.volume_celula_litros

    def calcular_tensao_nominal(self):
        return self.E0 * self.n_serie

    def calcular_volume_empacotamento_cm3(self):
        num_parallel_in_x = math.ceil(self.n_paralelo**0.5)
        num_parallel_in_y = math.ceil(self.n_paralelo / num_parallel_in_x)

        total_width_cm = num_parallel_in_x * (self.lado_cubo_x * 1.2)
        total_depth_cm = num_parallel_in_y * (self.lado_cubo_y * 1.2)
        total_height_cm = self.n_serie * (self.lado_cubo_z * 1.2)

        return total_width_cm * total_depth_cm * total_height_cm

    def plot_estrutura_3d(self):
        print("\n--- Gerando Visualização 3D da Estrutura da Bateria (Células como Cubos) ---")
        lx, ly, lz = self.lado_cubo_x, self.lado_cubo_y, self.lado_cubo_z

        data = []
        offset_x_spacing = lx * 0.2
        offset_y_spacing = ly * 0.2
        offset_z_spacing = lz * 0.2

        num_parallel_in_x = math.ceil(self.n_paralelo**0.5)
        num_parallel_in_y = math.ceil(self.n_paralelo / num_parallel_in_x)

        _cell_positions = {}

        for s in range(self.n_serie):
            for p_row in range(num_parallel_in_y):
                for p_col in range(num_parallel_in_x):
                    current_p_index = p_row * num_parallel_in_x + p_col
                    if current_p_index >= self.n_paralelo:
                        continue

                    x_start = p_col * (lx + offset_x_spacing)
                    y_start = p_row * (ly + offset_y_spacing)
                    z_start = s * (lz + offset_z_spacing)

                    data.append(go.Mesh3d(
                        x=[x_start, x_start + lx, x_start + lx, x_start, x_start, x_start + lx, x_start + lx, x_start],
                        y=[y_start, y_start, y_start + ly, y_start + ly, y_start, y_start, y_start + ly, y_start + ly],
                        z=[z_start, z_start, z_start, z_start, z_start + lz, z_start + lz, z_start + lz, z_start + lz],
                        i = [7, 0, 0, 0, 4, 4, 2, 6, 6, 4, 0, 1],
                        j = [3, 4, 1, 2, 5, 6, 5, 5, 2, 2, 6, 7],
                        k = [0, 7, 2, 3, 6, 7, 1, 1, 3, 3, 5, 5],
                        color='red', opacity=0.7, showscale=False, showlegend=False
                    ))

                    _cell_positions[(s, current_p_index)] = {
                        'x_start': x_start, 'y_start': y_start, 'z_start': z_start,
                        'lx': lx, 'ly': ly, 'lz': lz
                    }

        for p_idx in range(self.n_paralelo):
            for s in range(self.n_serie - 1):
                pos1 = _cell_positions.get((s, p_idx))
                pos2 = _cell_positions.get((s + 1, p_idx))

                if pos1 and pos2:
                    x1_conn = pos1['x_start'] + pos1['lx'] / 2
                    y1_conn = pos1['y_start'] + pos1['ly'] / 2
                    z1_conn = pos1['z_start'] + pos1['lz']

                    x2_conn = pos2['x_start'] + pos2['lx'] / 2
                    y2_conn = pos2['y_start'] + pos2['ly'] / 2
                    z2_conn = pos2['z_start']

                    data.append(go.Scatter3d(
                        x=[x1_conn, x2_conn],
                        y=[y1_conn, y2_conn],
                        z=[z1_conn, z2_conn],
                        mode='lines',
                        line=dict(color='orange', width=5),
                        name='Conexão Série',
                        showlegend=True if s == 0 and p_idx == 0 else False
                    ))

        p_rows = num_parallel_in_y
        p_cols = num_parallel_in_x

        for s in range(self.n_serie):
            for p_row in range(p_rows):
                for p_col in range(p_cols):
                    current_p_index = p_row * p_cols + p_col
                    if current_p_index >= self.n_paralelo:
                        continue

                    pos1 = _cell_positions.get((s, current_p_index))

                    if p_col < p_cols - 1:
                        next_p_index = p_row * p_cols + (p_col + 1)
                        if next_p_index < self.n_paralelo:
                            pos2 = _cell_positions.get((s, next_p_index))
                            if pos1 and pos2:
                                x1_conn = pos1['x_start'] + pos1['lx']
                                y1_conn = pos1['y_start'] + pos1['ly'] / 2
                                z1_conn = pos1['z_start'] + pos1['lz'] / 2

                                x2_conn = pos2['x_start']
                                y2_conn = pos2['y_start'] + pos2['ly'] / 2
                                z2_conn = pos2['z_start'] + pos2['lz'] / 2

                                data.append(go.Scatter3d(
                                    x=[x1_conn, x2_conn],
                                    y=[y1_conn, y2_conn],
                                    z=[z1_conn, z2_conn],
                                    mode='lines',
                                    line=dict(color='blue', width=5),
                                    name='Conexão Paralelo',
                                    showlegend=True if s == 0 and current_p_index == 0 else False
                                ))

                    if p_row < p_rows - 1:
                        next_p_index = (p_row + 1) * p_cols + p_col
                        if next_p_index < self.n_paralelo:
                            pos2 = _cell_positions.get((s, next_p_index))
                            if pos1 and pos2:
                                x1_conn = pos1['x_start'] + pos1['lx'] / 2
                                y1_conn = pos1['y_start'] + pos1['ly']
                                z1_conn = pos1['z_start'] + pos1['lz'] / 2

                                x2_conn = pos2['x_start'] + pos2['lx'] / 2
                                y2_conn = pos2['y_start']
                                z2_conn = pos2['z_start'] + pos2['lz'] / 2

                                data.append(go.Scatter3d(
                                    x=[x1_conn, x2_conn],
                                    y=[y1_conn, y2_conn],
                                    z=[z1_conn, z2_conn],
                                    mode='lines',
                                    line=dict(color='blue', width=5),
                                    name='Conexão Paralelo',
                                    showlegend=False
                                ))


        fig = go.Figure(data=data)

        fig.update_layout(
            scene=dict(
                xaxis_title='Largura Total (cm)',
                yaxis_title='Profundidade Total (cm)',
                zaxis_title='Altura Total (cm)',
                aspectmode='data',
            ),
            title=f'Estrutura 3D do Pack de Baterias ({self.n_serie}S{self.n_paralelo}P)',
            showlegend=True
        )

        try:
            fig.show()
            print("--- Visualização 3D incorporada na saída da célula do Colab. ---")
        except Exception as e:
            print(f"Não foi possível exibir o gráfico no Colab: {e}")
            print("Tentando salvar o gráfico como um arquivo HTML...")
            file_name = f'estrutura_bateria_{self.n_serie}S{self.n_paralelo}P.html'
            pio.write_html(fig, file=file_name, auto_open=True)
            print(f"Gráfico salvo como {file_name}. Por favor, baixe e abra este arquivo no seu navegador.")


def find_optimal_battery_pack_by_volume(max_volume_cm3, min_voltage_V=None, max_voltage_V=None):
    """
    Encontra a bateria com maior ENERGIA (kWh) que se encaixa
    em um volume de empacotamento máximo especificado,
    e opcionalmente atende a uma tensão nominal mínima e máxima.

    Args:
        max_volume_cm3 (float): O volume máximo disponível em centímetros cúbicos.
        min_voltage_V (float, optional): Tensão nominal mínima requerida em Volts. Defaults to None.
        max_voltage_V (float, optional): Tensão nominal máxima permitida em Volts. Defaults to None.

    Returns:
        tuple: (optimal_pack_details_dict, optimal_pack_object)
                Retorna um dicionário com os detalhes da melhor bateria e o objeto BatteryPack.
                Retorna (None, None) se nenhuma combinação se encaixar.
    """
    print(f"\n--- Buscando a Bateria Ideal para Volume Máximo de {max_volume_cm3:.2f} cm³ ---")
    if min_voltage_V:
        print(f"    e Tensão Nominal Mínima de {min_voltage_V:.2f} V")
    if max_voltage_V:
        print(f"    e Tensão Nominal Máxima de {max_voltage_V:.2f} V ---")


    optimal_energy_kwh = -1.0 # Objetivo: Maximizar energia, não só Ah
    optimal_pack_details = None
    optimal_pack_object = None

    cell_types = ['Li-ion', 'LiFePO4']

    # Definir os limites de busca para n_serie e n_paralelo
    # Para FSAE, n_serie pode ir até ~160 para 600V (Li-ion) ou ~188 (LiFePO4)
    # n_paralelo pode variar de 1 a 20 ou mais, dependendo da capacidade da célula e do Ah total desejado
    # Ajuste esses ranges para explorar um espaço de design relevante para o seu caso.
    max_n_serie = 200 # Um pouco acima do necessário para 600V para ter margem
    max_n_paralelo = 30 # Depende da corrente e capacidade desejada

    for cell_type in cell_types:
        for ns in range(1, max_n_serie + 1):
            for np_val in range(1, max_n_paralelo + 1):
                try:
                    current_pack = BatteryPack(tipo_celula=cell_type, n_serie=ns, n_paralelo=np_val)
                    pack_volume_cm3 = current_pack.calcular_volume_empacotamento_cm3()
                    current_nominal_voltage_V = current_pack.calcular_tensao_nominal()
                    current_charge_capacity_ah = current_pack.carga_total
                    
                    # CÁLCULO DA ENERGIA (kWh) - NOVO OBJETIVO
                    current_energy_kwh = (current_nominal_voltage_V * current_charge_capacity_ah) / 1000  

                    # 1. Verificar restrição de volume
                    if pack_volume_cm3 <= max_volume_cm3:
                        # 2. Verificar restrição de tensão mínima (se especificada)
                        if min_voltage_V is None or current_nominal_voltage_V >= min_voltage_V:
                            # 3. Verificar restrição de tensão máxima (NOVA CONDIÇÃO)
                            if max_voltage_V is None or current_nominal_voltage_V <= max_voltage_V:
                                # 4. Se todas as restrições forem atendidas, verificar se é o melhor em termos de ENERGIA
                                if current_energy_kwh > optimal_energy_kwh:
                                    optimal_energy_kwh = current_energy_kwh
                                    optimal_pack_details = {
                                        'tipo_celula': cell_type,
                                        'n_serie': ns,
                                        'n_paralelo': np_val,
                                        'tensao_nominal_V': current_nominal_voltage_V,
                                        'carga_total_Ah': current_charge_capacity_ah,
                                        'energia_total_kWh': current_energy_kwh, # Adicionar energia aos detalhes
                                        'peso_total_kg': current_pack.calcular_peso_total(),
                                        'volume_empacotamento_cm3': pack_volume_cm3
                                    }
                                    optimal_pack_object = current_pack
                except ValueError:
                    # Ignorar tipos de célula desconhecidos (já tratado na classe)
                    pass
                except Exception as e:
                    # Capturar outras exceções que possam ocorrer durante a criação do pack
                    # print(f"Erro ao criar pack {cell_type}, {ns}S{np_val}P: {e}") # Descomente para depurar erros específicos
                    pass

    if optimal_pack_details:
        print("\n--- Bateria Otimizada Encontrada ---")
        print(f"Tipo de Célula: {optimal_pack_details['tipo_celula']}")
        print(f"Configuração: {optimal_pack_details['n_serie']}S{optimal_pack_details['n_paralelo']}P")
        print(f"Tensão Nominal: {optimal_pack_details['tensao_nominal_V']:.2f} V")
        print(f"Carga Total: {optimal_pack_details['carga_total_Ah']:.2f} Ah")
        print(f"Energia Total: {optimal_pack_details['energia_total_kWh']:.2f} kWh")
        print(f"Peso Total: {optimal_pack_details['peso_total_kg']:.2f} kg")
        print(f"Volume de Empacotamento: {optimal_pack_details['volume_empacotamento_cm3']:.2f} cm³")
        return optimal_pack_details, optimal_pack_object
    else:
        print(f"Nenhuma combinação de bateria encontrada que se encaixe no volume de {max_volume_cm3:.2f} cm³.")
        if min_voltage_V:
            print(f"Considerando a tensão nominal mínima de {min_voltage_V:.2f} V.")
        if max_voltage_V:
            print(f"Considerando a tensão nominal máxima de {max_voltage_V:.2f} V.")
        return None, None

# --- Exemplo de uso para encontrar a bateria ideal ---
if __name__ == "__main__":
    # Defina o volume máximo disponível para a bateria (em litros)
    max_volume_litros = 10.0 # Exemplo: 10 litros de espaço disponível
    max_volume_cm3 = max_volume_litros * 1000 # Convertendo para cm³

    # Adicione a tensão nominal mínima e MÁXIMA desejada para FSAE
    min_fsae_voltage = 300.0 # Tensão mínima requerida em Volts
    max_fsae_voltage = 600.0 # Tensão máxima permitida em Volts

    optimal_pack_details, optimal_pack_obj = find_optimal_battery_pack_by_volume(
        max_volume_cm3,
        min_voltage_V=min_fsae_voltage, # Passar a restrição de tensão mínima
        max_voltage_V=max_fsae_voltage  # Passar a restrição de tensão máxima aqui
    )

    if optimal_pack_obj:
        # Se uma bateria otimizada for encontrada, visualize sua estrutura 3D
        optimal_pack_obj.plot_estrutura_3d()

