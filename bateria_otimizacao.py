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

        # Para armazenar o layout otimizado para plotagem
        self._plot_layout_factors = None

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

    def calcular_dimensoes_empacotamento_cm(self):
        """
        Calcula as dimensões totais de empacotamento do pack em cm,
        tentando distribuir as células em série para equilibrar as dimensões.
        Armazena os fatores de layout para a plotagem 3D.
        Retorna (largura_total_cm, profundidade_total_cm, altura_total_cm).
        """
        lx, ly, lz = self.lado_cubo_x, self.lado_cubo_y, self.lado_cubo_z
        spacing_factor = 1.1 # Margem de 10% para espaçamento e invólucro

        # Base para um único bloco paralelo (n_paralelo células em um plano)
        num_parallel_in_x_base = math.ceil(self.n_paralelo**0.5)
        num_parallel_in_y_base = math.ceil(self.n_paralelo / num_parallel_in_x_base)

        width_base_dim = num_parallel_in_x_base * (lx * spacing_factor)
        depth_base_dim = num_parallel_in_y_base * (ly * spacing_factor)
        height_base_dim = lz * spacing_factor # Altura de um "bloco" paralelo

        best_total_width = float('inf')
        best_total_depth = float('inf')
        best_total_height = float('inf')
        min_aspect_ratio_diff = float('inf') 
        
        best_s_w_factor = 1
        best_s_d_factor = 1
        best_s_h_factor = self.n_serie

        max_distribute_factor = min(self.n_serie, 10) 

        for s_w_candidate in range(1, max_distribute_factor + 1):
            remaining_s_for_d_h = self.n_serie / s_w_candidate
            if remaining_s_for_d_h != int(remaining_s_for_d_h):
                continue
            remaining_s_for_d_h = int(remaining_s_for_d_h)

            for s_d_candidate in range(1, min(remaining_s_for_d_h, max_distribute_factor) + 1):
                s_h_candidate_exact = remaining_s_for_d_h / s_d_candidate
                if s_h_candidate_exact != int(s_h_candidate_exact):
                    continue
                s_h_candidate = int(s_h_candidate_exact)

                current_width = s_w_candidate * width_base_dim
                current_depth = s_d_candidate * depth_base_dim
                current_height = s_h_candidate * height_base_dim

                dims = sorted([current_width, current_depth, current_height])
                
                if dims[0] <= 0:
                    ar_diff = float('inf')
                else:
                    ar_diff = max(dims) / min(dims) 

                if ar_diff < min_aspect_ratio_diff:
                    min_aspect_ratio_diff = ar_diff
                    best_total_width, best_total_depth, best_total_height = current_width, current_depth, current_height
                    best_s_w_factor = s_w_candidate
                    best_s_d_factor = s_d_candidate
                    best_s_h_factor = s_h_candidate

        self._plot_layout_factors = {
            'num_parallel_in_x_base': num_parallel_in_x_base,
            'num_parallel_in_y_base': num_parallel_in_y_base,
            's_w_factor': best_s_w_factor,
            's_d_factor': best_s_d_factor,
            's_h_factor': best_s_h_factor
        }
        
        return best_total_width, best_total_depth, best_total_height


    def plot_estrutura_3d(self):
        print("\n--- Gerando Visualização 3D da Estrutura da Bateria (Células como Cubos) ---")
        lx, ly, lz = self.lado_cubo_x, self.lado_cubo_y, self.lado_cubo_z
        
        spacing_factor = 1.1
        offset_x_spacing = lx * 0.1
        offset_y_spacing = ly * 0.1
        offset_z_spacing = lz * 0.1

        data = []
        _cell_positions = {}
        
        if not self._plot_layout_factors:
            _,_,_ = self.calcular_dimensoes_empacotamento_cm()
        
        num_parallel_in_x_base = self._plot_layout_factors['num_parallel_in_x_base']
        num_parallel_in_y_base = self._plot_layout_factors['num_parallel_in_y_base']
        s_w_factor = self._plot_layout_factors['s_w_factor']
        s_d_factor = self._plot_layout_factors['s_d_factor']
        s_h_factor = self._plot_layout_factors['s_h_factor']

        current_s_idx = 0
        for s_z_block in range(s_h_factor):
            for s_y_block in range(s_d_factor):
                for s_x_block in range(s_w_factor):
                    if current_s_idx >= self.n_serie:
                        break
                    
                    block_x_offset = s_x_block * (num_parallel_in_x_base * (lx + offset_x_spacing))
                    block_y_offset = s_y_block * (num_parallel_in_y_base * (ly + offset_y_spacing))
                    block_z_offset = s_z_block * (lz + offset_z_spacing)

                    for p_row in range(num_parallel_in_y_base):
                        for p_col in range(num_parallel_in_x_base):
                            current_p_index = p_row * num_parallel_in_x_base + p_col
                            if current_p_index >= self.n_paralelo:
                                continue

                            x_start = block_x_offset + p_col * (lx + offset_x_spacing)
                            y_start = block_y_offset + p_row * (ly + offset_y_spacing)
                            z_start = block_z_offset 

                            data.append(go.Mesh3d(
                                x=[x_start, x_start + lx, x_start + lx, x_start, x_start, x_start + lx, x_start + lx, x_start],
                                y=[y_start, y_start, y_start + ly, y_start + ly, y_start, y_start, y_start + ly, y_start + ly],
                                z=[z_start, z_start, z_start, z_start, z_start + lz, z_start + lz, z_start + lz, z_start + lz],
                                i = [7, 0, 0, 0, 4, 4, 2, 6, 6, 4, 0, 1],
                                j = [3, 4, 1, 2, 5, 6, 5, 5, 2, 2, 6, 7],
                                k = [0, 7, 2, 3, 6, 7, 1, 1, 3, 3, 5, 5],
                                color='red', opacity=0.7, showscale=False, showlegend=False
                            ))
                            _cell_positions[(current_s_idx, current_p_index)] = {
                                'x_start': x_start, 'y_start': y_start, 'z_start': z_start,
                                'lx': lx, 'ly': ly, 'lz': lz
                            }
                    current_s_idx += 1 

        series_conn_points = []
        for s_idx in range(self.n_serie - 1):
            for p_idx in range(self.n_paralelo):
                pos1 = _cell_positions.get((s_idx, p_idx))
                pos2 = _cell_positions.get((s_idx + 1, p_idx))
                if pos1 and pos2:
                    series_conn_points.append(
                        (pos1['x_start'] + pos1['lx']/2, pos1['y_start'] + pos1['ly']/2, pos1['z_start'] + pos1['lz'])
                    )
                    series_conn_points.append(
                        (pos2['x_start'] + pos2['lx']/2, pos2['y_start'] + pos2['ly']/2, pos2['z_start'])
                    )
        
        if series_conn_points:
            data.append(go.Scatter3d(
                x=[p[0] for p in series_conn_points],
                y=[p[1] for p in series_conn_points],
                z=[p[2] for p in series_conn_points],
                mode='lines',
                line=dict(color='orange', width=5),
                name='Conexão Série',
                showlegend=True
            ))

        parallel_conn_points = []
        for s_idx in range(self.n_serie):
            for p_row in range(num_parallel_in_y_base):
                for p_col in range(num_parallel_in_x_base):
                    current_p_index_in_base_block = p_row * num_parallel_in_x_base + p_col
                    if current_p_index_in_base_block >= self.n_paralelo:
                        continue

                    pos1 = _cell_positions.get((s_idx, current_p_index_in_base_block))

                    if pos1:
                        if p_col < num_parallel_in_x_base - 1:
                            next_p_index_in_base_block_h = p_row * num_parallel_in_x_base + (p_col + 1)
                            if next_p_index_in_base_block_h < self.n_paralelo:
                                pos2_h = _cell_positions.get((s_idx, next_p_index_in_base_block_h))
                                if pos2_h:
                                    parallel_conn_points.append(
                                        (pos1['x_start'] + pos1['lx'], pos1['y_start'] + pos1['ly']/2, pos1['z_start'] + pos1['lz']/2)
                                    )
                                    parallel_conn_points.append(
                                        (pos2_h['x_start'], pos2_h['y_start'] + pos2_h['ly']/2, pos2_h['z_start'] + pos2_h['lz']/2)
                                    )
                        if p_row < num_parallel_in_y_base - 1:
                            next_p_index_in_base_block_v = (p_row + 1) * num_parallel_in_x_base + p_col
                            if next_p_index_in_base_block_v < self.n_paralelo:
                                pos2_v = _cell_positions.get((s_idx, next_p_index_in_base_block_v))
                                if pos2_v:
                                    parallel_conn_points.append(
                                        (pos1['x_start'] + pos1['lx']/2, pos1['y_start'] + pos1['ly'], pos1['z_start'] + pos1['lz']/2)
                                    )
                                    parallel_conn_points.append(
                                        (pos2_v['x_start'] + pos2_v['lx']/2, pos2_v['y_start'], pos2_v['z_start'] + pos2_v['lz']/2)
                                    )
        
        if parallel_conn_points:
            data.append(go.Scatter3d(
                x=[p[0] for p in parallel_conn_points],
                y=[p[1] for p in parallel_conn_points],
                z=[p[2] for p in parallel_conn_points],
                mode='lines',
                line=dict(color='blue', width=5),
                name='Conexão Paralelo',
                showlegend=True
            ))

        fig = go.Figure(data=data)

        total_width, total_depth, total_height = self.calcular_dimensoes_empacotamento_cm()
        x_title = f'Largura Total ({total_width:.2f} cm)'
        y_title = f'Profundidade Total ({total_depth:.2f} cm)'
        z_title = f'Altura Total ({total_height:.2f} cm)'

        fig.update_layout(
            scene=dict(
                xaxis_title=x_title,
                yaxis_title=y_title,
                zaxis_title=z_title,
                aspectmode='data',
            ),
            title=f'Estrutura 3D do Pack de Baterias ({self.n_serie}S{self.n_paralelo}P) - Layout Otimizado',
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


def find_optimal_battery_pack_by_volume(
    max_width_cm, max_depth_cm, max_height_cm,
    min_voltage_V=None, max_voltage_V=None, max_energy_kwh=None # NOVO PARÂMETRO
):
    """
    Encontra a bateria com maior ENERGIA (kWh) que se encaixa
    em um volume de empacotamento máximo especificado,
    e opcionalmente atende a uma tensão nominal mínima e máxima, e energia máxima.

    Args:
        max_width_cm (float): A largura máxima disponível em centímetros.
        max_depth_cm (float): A profundidade máxima disponível em centímetros.
        max_height_cm (float): A altura máxima disponível em centímetros.
        min_voltage_V (float, optional): Tensão nominal mínima requerida em Volts. Defaults to None.
        max_voltage_V (float, optional): Tensão nominal máxima permitida em Volts. Defaults to None.
        max_energy_kwh (float, optional): Energia máxima permitida em kWh. Defaults to None. # NOVO

    Returns:
        tuple: (optimal_pack_details_dict, optimal_pack_object)
                Retorna um dicionário com os detalhes da melhor bateria e o objeto BatteryPack.
                Retorna (None, None) se nenhuma combinação se encaixe.
    """
    print(f"\n--- Buscando a Bateria Ideal para Dimensões Máximas: {max_width_cm:.2f}x{max_depth_cm:.2f}x{max_height_cm:.2f} cm ---")
    if min_voltage_V:
        print(f"    e Tensão Nominal Mínima de {min_voltage_V:.2f} V")
    if max_voltage_V:
        print(f"    e Tensão Nominal Máxima de {max_voltage_V:.2f} V")
    if max_energy_kwh: # NOVO PRINT
        print(f"    e Energia Máxima de {max_energy_kwh:.2f} kWh ---")
    else:
        print("---")


    optimal_energy_kwh = -1.0 # Objetivo: Maximizar energia, não só Ah
    optimal_pack_details = None
    optimal_pack_object = None
    optimal_volume_utilization = -1.0 # Novo: Acompanhar o quão bem o volume é utilizado

    cell_types = ['Li-ion', 'LiFePO4']

    max_n_serie = 200 
    max_n_paralelo = 30 

    # Calcula o volume máximo disponível
    max_total_volume_cm3 = max_width_cm * max_depth_cm * max_height_cm

    for cell_type in cell_types:
        for ns in range(1, max_n_serie + 1):
            for np_val in range(1, max_n_paralelo + 1):
                try:
                    current_pack = BatteryPack(tipo_celula=cell_type, n_serie=ns, n_paralelo=np_val)
                    
                    pack_width_cm, pack_depth_cm, pack_height_cm = current_pack.calcular_dimensoes_empacotamento_cm()
                    current_pack_volume_cm3 = pack_width_cm * pack_depth_cm * pack_height_cm 

                    current_nominal_voltage_V = current_pack.calcular_tensao_nominal()
                    current_charge_capacity_ah = current_pack.carga_total
                    
                    current_energy_kwh = (current_nominal_voltage_V * current_charge_capacity_ah) / 1000  

                    # 1. Verificar restrições de dimensão (largura, profundidade, altura)
                    if (pack_width_cm <= max_width_cm and
                        pack_depth_cm <= max_depth_cm and
                        pack_height_cm <= max_height_cm):
                        
                        # 2. Verificar restrição de tensão mínima (se especificada)
                        if min_voltage_V is None or current_nominal_voltage_V >= min_voltage_V:
                            # 3. Verificar restrição de tensão máxima (NOVA CONDIÇÃO)
                            if max_voltage_V is None or current_nominal_voltage_V <= max_voltage_V:
                                # 4. Verificar restrição de energia máxima (NOVA CONDIÇÃO)
                                if max_energy_kwh is None or current_energy_kwh <= max_energy_kwh:
                                    
                                    # Calcular a utilização do volume para a bateria atual
                                    current_volume_utilization = current_pack_volume_cm3 / max_total_volume_cm3
                                    
                                    # 5. Priorizar a maior energia, e em caso de empate, o maior aproveitamento de volume
                                    if (current_energy_kwh > optimal_energy_kwh or
                                        (abs(current_energy_kwh - optimal_energy_kwh) < 1e-9 and # Considera energias muito próximas como empate
                                         current_volume_utilization > optimal_volume_utilization)):

                                        optimal_energy_kwh = current_energy_kwh
                                        optimal_volume_utilization = current_volume_utilization
                                        optimal_pack_details = {
                                            'tipo_celula': cell_type,
                                            'n_serie': ns,
                                            'n_paralelo': np_val,
                                            'tensao_nominal_V': current_nominal_voltage_V,
                                            'carga_total_Ah': current_charge_capacity_ah,
                                            'energia_total_kWh': current_energy_kwh,
                                            'peso_total_kg': current_pack.calcular_peso_total(),
                                            'volume_empacotamento_cm3': current_pack_volume_cm3,
                                            'largura_cm': pack_width_cm,
                                            'profundidade_cm': pack_depth_cm,
                                            'altura_cm': pack_height_cm,
                                            'utilizacao_volume_%': current_volume_utilization * 100 # Adicionado
                                        }
                                        optimal_pack_object = current_pack
                except ValueError:
                    pass
                except Exception as e:
                    # print(f"Erro ao processar ns={ns}, np={np_val}: {e}") # Descomente para depurar erros
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
        print(f"Dimensões (L x P x A): {optimal_pack_details['largura_cm']:.2f} cm x {optimal_pack_details['profundidade_cm']:.2f} cm x {optimal_pack_details['altura_cm']:.2f} cm")
        print(f"Utilização do Volume Disponível: {optimal_pack_details['utilizacao_volume_%']:.2f} %") # Novo
        return optimal_pack_details, optimal_pack_object
    else:
        print(f"Nenhuma combinação de bateria encontrada que se encaixe nas dimensões máximas de {max_width_cm:.2f}x{max_depth_cm:.2f}x{max_height_cm:.2f} cm³.")
        if min_voltage_V:
            print(f"Considerando a tensão nominal mínima de {min_voltage_V:.2f} V.")
        if max_voltage_V:
            print(f"Considerando a tensão nominal máxima de {max_voltage_V:.2f} V.")
        if max_energy_kwh:
            print(f"Considerando a energia máxima de {max_energy_kwh:.2f} kWh.")
        return None, None

# --- Exemplo de uso para encontrar a bateria ideal ---
if __name__ == "__main__":
    # Defina as dimensões máximas disponíveis para a bateria (em metros)
    max_width_m = 0.3 # 50 cm
    max_depth_m = 0.3 # 30 cm
    max_height_m = 0.5 # 30 cm

    # Converte para centímetros para uso no código
    max_width_cm = max_width_m * 100
    max_depth_cm = max_depth_m * 100
    max_height_cm = max_height_m * 100

    # Adicione a tensão nominal mínima e MÁXIMA desejada para FSAE
    min_fsae_voltage = 220.0 # Tensão mínima requerida em Volts
    max_fsae_voltage = 600.0 # Tensão máxima permitida em Volts

    optimal_pack_details, optimal_pack_obj = find_optimal_battery_pack_by_volume(
        max_width_cm=max_width_cm,
        max_depth_cm=max_depth_cm,
        max_height_cm=max_height_cm,
        min_voltage_V=min_fsae_voltage,
        max_voltage_V=max_fsae_voltage,
    )

    if optimal_pack_obj:
        optimal_pack_obj.plot_estrutura_3d()
