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
        # Este spacing_factor afeta os CÁLCULOS DE DIMENSÃO do pack.
        # Deve ser mantido baixo para refletir um empacotamento real.
        spacing_factor = 1.1 

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


    def plot_estrutura_3d(self, title_suffix=""):
        print(f"\n--- Gerando Visualização 3D da Estrutura da Bateria (Células como Cubos) {title_suffix} ---")
        lx, ly, lz = self.lado_cubo_x, self.lado_cubo_y, self.lado_cubo_z
        
    
        visual_spacing_factor = 1.5 # Ajuste este valor para mais ou menos espaço
        offset_x_spacing_visual = lx * (visual_spacing_factor - 1)
        offset_y_spacing_visual = ly * (visual_spacing_factor - 1)
        offset_z_spacing_visual = lz * (visual_spacing_factor - 1)

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
                    
                    # Usando offset_x/y/z_spacing_visual para a plotagem
                    block_x_offset = s_x_block * (num_parallel_in_x_base * (lx + offset_x_spacing_visual))
                    block_y_offset = s_y_block * (num_parallel_in_y_base * (ly + offset_y_spacing_visual))
                    block_z_offset = s_z_block * (lz + offset_z_spacing_visual)

                    for p_row in range(num_parallel_in_y_base):
                        for p_col in range(num_parallel_in_x_base):
                            current_p_index = p_row * num_parallel_in_x_base + p_col
                            if current_p_index >= self.n_paralelo:
                                continue

                            x_start = block_x_offset + p_col * (lx + offset_x_spacing_visual)
                            y_start = block_y_offset + p_row * (ly + offset_y_spacing_visual)
                            z_start = block_z_offset 

                            data.append(go.Mesh3d(
                                x=[x_start, x_start + lx, x_start + lx, x_start, x_start, x_start + lx, x_start + lx, x_start],
                                y=[y_start, y_start, y_start + ly, y_start + ly, y_start, y_start, y_start + ly, y_start + ly],
                                z=[z_start, z_start, z_start, z_start, z_start + lz, z_start + lz, z_start + lz, z_start + lz],
                                i = [7, 0, 0, 0, 4, 4, 2, 6, 6, 4, 0, 1],
                                j = [3, 4, 1, 2, 5, 6, 5, 5, 2, 2, 6, 7],
                                k = [0, 7, 2, 3, 6, 7, 1, 1, 3, 3, 5, 5],
                                color='red', opacity=0.4, showscale=False, showlegend=False,
                                lighting=dict(ambient=0.4, diffuse=0.8, specular=0.2, roughness=0.5, fresnel=0.0), # Melhorar iluminação
                                lightposition=dict(x=0, y=0, z=10000) # Posição da luz
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
                    # Conecta o topo de uma célula ao fundo da próxima na série
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
                line=dict(color='darkorange', width=10), # Maior espessura para destaque
                name='Conexão Série',
                showlegend=True
            ))

        parallel_conn_points = []
        # Conexões paralelas para células dentro da mesma camada "paralela"
        for s_idx in range(self.n_serie): # Para cada "string" em série
            # Conexões ao longo do eixo X (largura) para células na mesma linha
            for p_row in range(num_parallel_in_y_base):
                for p_col in range(num_parallel_in_x_base - 1): # Vai até a penúltima para conectar com a próxima
                    current_p_index = p_row * num_parallel_in_x_base + p_col
                    next_p_index_h = p_row * num_parallel_in_x_base + (p_col + 1)
                    if current_p_index < self.n_paralelo and next_p_index_h < self.n_paralelo:
                        pos1 = _cell_positions.get((s_idx, current_p_index))
                        pos2_h = _cell_positions.get((s_idx, next_p_index_h))
                        if pos1 and pos2_h:
                            # Conecta o meio da face direita de uma célula com o meio da face esquerda da próxima
                            parallel_conn_points.append(
                                (pos1['x_start'] + pos1['lx'], pos1['y_start'] + pos1['ly']/2, pos1['z_start'] + pos1['lz']/2)
                            )
                            parallel_conn_points.append(
                                (pos2_h['x_start'], pos2_h['y_start'] + pos2_h['ly']/2, pos2_h['z_start'] + pos2_h['lz']/2)
                            )
            # Conexões ao longo do eixo Y (profundidade) para células na mesma coluna
            for p_col in range(num_parallel_in_x_base):
                for p_row in range(num_parallel_in_y_base - 1): # Vai até a penúltima para conectar com a próxima
                    current_p_index = p_row * num_parallel_in_x_base + p_col
                    next_p_index_v = (p_row + 1) * num_parallel_in_x_base + p_col
                    if current_p_index < self.n_paralelo and next_p_index_v < self.n_paralelo:
                        pos1 = _cell_positions.get((s_idx, current_p_index))
                        pos2_v = _cell_positions.get((s_idx, next_p_index_v))
                        if pos1 and pos2_v:
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
                line=dict(color='blue', width=10), # Maior espessura
                name='Conexão Paralelo (Interna à Chapa)',
                showlegend=True
            ))

        fig = go.Figure(data=data)

        # As dimensões totais no título ainda refletem as dimensões REAIS de empacotamento
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
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1.5) 
                ),
                xaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
                yaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
                zaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
            ),
            title=f'Estrutura 3D do Pack de Baterias ({self.n_serie}S{self.n_paralelo}P) - Layout Otimizado {title_suffix}',
            showlegend=True,
            margin=dict(l=0, r=0, b=0, t=50)
        )

        try:
            fig.show()
            print("--- Visualização 3D incorporada na saída da célula do Colab. ---")
        except Exception as e:
            print(f"Não foi possível exibir o gráfico no Colab: {e}")
            print("Tentando salvar o gráfico como um arquivo HTML...")
            file_name = f'estrutura_bateria_{self.n_serie}S{self.n_paralelo}P{title_suffix.replace(" ", "_").replace("(", "").replace(")", "")}.html'
            pio.write_html(fig, file=file_name, auto_open=True)
            print(f"Gráfico salvo como {file_name}. Por favor, baixe e abra este arquivo no seu navegador.")


def plot_pack_3_chapas_paralelas(optimal_pack_obj_chapa, num_chapas=3):
    """
    Gera uma visualização 3D de N chapas de bateria idênticas
    dispostas em paralelo ao longo do eixo de profundidade (Y).
    Inclui as células individuais e as conexões paralelas entre as chapas.

    Args:
        optimal_pack_obj_chapa (BatteryPack): O objeto BatteryPack otimizado para uma única chapa.
        num_chapas (int): O número de chapas a serem visualizadas. Padrão para 3.
    """
    print(f"\n--- Gerando Visualização 3D do Pack Completo ( {num_chapas} Chapas em Paralelo) ---")

    chapa_width, chapa_depth, chapa_height = optimal_pack_obj_chapa.calcular_dimensoes_empacotamento_cm()

    full_pack_data = []
    
    all_chapa_positive_terminals = [] 
    all_chapa_negative_terminals = [] 

    lx, ly, lz = optimal_pack_obj_chapa.lado_cubo_x, optimal_pack_obj_chapa.lado_cubo_y, optimal_pack_obj_chapa.lado_cubo_z
    
    # Aumentando o espaçamento para FINS DE VISUALIZAÇÃO APENAS
    visual_spacing_factor = 1.5 
    offset_x_spacing_visual = lx * (visual_spacing_factor - 1)
    offset_y_spacing_visual = ly * (visual_spacing_factor - 1)
    offset_z_spacing_visual = lz * (visual_spacing_factor - 1)
    
    # Garantir que _plot_layout_factors esteja preenchido
    optimal_pack_obj_chapa.calcular_dimensoes_empacotamento_cm()

    num_parallel_in_x_base = optimal_pack_obj_chapa._plot_layout_factors['num_parallel_in_x_base']
    num_parallel_in_y_base = optimal_pack_obj_chapa._plot_layout_factors['num_parallel_in_y_base']
    s_w_factor = optimal_pack_obj_chapa._plot_layout_factors['s_w_factor']
    s_d_factor = optimal_pack_obj_chapa._plot_layout_factors['s_d_factor']
    s_h_factor = optimal_pack_obj_chapa._plot_layout_factors['s_h_factor']

    spacing_between_chapas_cm = 5 # Espaçamento entre as chapas para visualização (ajustável)

    for i in range(num_chapas):
        chapa_offset_y = i * (chapa_depth + spacing_between_chapas_cm) 

        current_s_idx = 0
        current_chapa_cell_positions = {} 

        for s_z_block in range(s_h_factor):
            for s_y_block in range(s_d_factor):
                for s_x_block in range(s_w_factor):
                    if current_s_idx >= optimal_pack_obj_chapa.n_serie:
                        break
                    
                    block_x_offset = s_x_block * (num_parallel_in_x_base * (lx + offset_x_spacing_visual))
                    block_y_offset = s_y_block * (num_parallel_in_y_base * (ly + offset_y_spacing_visual))
                    block_z_offset = s_z_block * (lz + offset_z_spacing_visual)

                    for p_row in range(num_parallel_in_y_base):
                        for p_col in range(num_parallel_in_x_base):
                            current_p_index = p_row * num_parallel_in_x_base + p_col
                            if current_p_index >= optimal_pack_obj_chapa.n_paralelo:
                                continue

                            x_start = block_x_offset + p_col * (lx + offset_x_spacing_visual)
                            y_start = chapa_offset_y + block_y_offset + p_row * (ly + offset_y_spacing_visual)
                            z_start = block_z_offset 

                            full_pack_data.append(go.Mesh3d(
                                x=[x_start, x_start + lx, x_start + lx, x_start, x_start, x_start + lx, x_start + lx, x_start],
                                y=[y_start, y_start, y_start + ly, y_start + ly, y_start, y_start, y_start + ly, y_start + ly],
                                z=[z_start, z_start, z_start, z_start, z_start + lz, z_start + lz, z_start + lz, z_start + lz],
                                i = [7, 0, 0, 0, 4, 4, 2, 6, 6, 4, 0, 1],
                                j = [3, 4, 1, 2, 5, 6, 5, 5, 2, 2, 6, 7],
                                k = [0, 7, 2, 3, 6, 7, 1, 1, 3, 3, 5, 5],
                                color='red', opacity=0.4, showscale=False, showlegend=False,
                                lighting=dict(ambient=0.4, diffuse=0.8, specular=0.2, roughness=0.5, fresnel=0.0), # Melhorar iluminação
                                lightposition=dict(x=0, y=0, z=10000) # Posição da luz
                            ))
                            current_chapa_cell_positions[(current_s_idx, current_p_index)] = {
                                'x_start': x_start, 'y_start': y_start, 'z_start': z_start,
                                'lx': lx, 'ly': ly, 'lz': lz
                            }
                    current_s_idx += 1 
        
        # Conexões em série DENTRO de CADA CHAPA replicada
        series_conn_points_chapa = []
        for s_idx in range(optimal_pack_obj_chapa.n_serie - 1):
            for p_idx in range(optimal_pack_obj_chapa.n_paralelo):
                pos1 = current_chapa_cell_positions.get((s_idx, p_idx))
                pos2 = current_chapa_cell_positions.get((s_idx + 1, p_idx))
                if pos1 and pos2:
                    series_conn_points_chapa.append(
                        (pos1['x_start'] + pos1['lx']/2, pos1['y_start'] + pos1['ly']/2, pos1['z_start'] + pos1['lz'])
                    )
                    series_conn_points_chapa.append(
                        (pos2['x_start'] + pos2['lx']/2, pos2['y_start'] + pos2['ly']/2, pos2['z_start'])
                    )
        if series_conn_points_chapa:
            full_pack_data.append(go.Scatter3d(
                x=[p[0] for p in series_conn_points_chapa],
                y=[p[1] for p in series_conn_points_chapa],
                z=[p[2] for p in series_conn_points_chapa],
                mode='lines',
                line=dict(color='darkorange', width=10), 
                name=f'Conexão Série (Chapa {i+1})', 
                showlegend=False 
            ))

        # Conexões paralelas INTERNAS DENTRO de CADA CHAPA replicada
        parallel_conn_points_chapa_internal = []
        for s_idx in range(optimal_pack_obj_chapa.n_serie):
            # Conexões ao longo do eixo X (largura)
            for p_row in range(num_parallel_in_y_base):
                for p_col in range(num_parallel_in_x_base - 1): 
                    current_p_index = p_row * num_parallel_in_x_base + p_col
                    next_p_index_h = p_row * num_parallel_in_x_base + (p_col + 1)
                    if current_p_index < optimal_pack_obj_chapa.n_paralelo and next_p_index_h < optimal_pack_obj_chapa.n_paralelo:
                        pos1 = current_chapa_cell_positions.get((s_idx, current_p_index))
                        pos2_h = current_chapa_cell_positions.get((s_idx, next_p_index_h))
                        if pos1 and pos2_h:
                            parallel_conn_points_chapa_internal.append(
                                (pos1['x_start'] + pos1['lx'], pos1['y_start'] + pos1['ly']/2, pos1['z_start'] + pos1['lz']/2)
                            )
                            parallel_conn_points_chapa_internal.append(
                                (pos2_h['x_start'], pos2_h['y_start'] + pos2_h['ly']/2, pos2_h['z_start'] + pos2_h['lz']/2)
                            )
            # Conexões ao longo do eixo Y (profundidade dentro da chapa)
            for p_col in range(num_parallel_in_x_base):
                for p_row in range(num_parallel_in_y_base - 1):
                    current_p_index = p_row * num_parallel_in_x_base + p_col
                    next_p_index_v = (p_row + 1) * num_parallel_in_x_base + p_col
                    if current_p_index < optimal_pack_obj_chapa.n_paralelo and next_p_index_v < optimal_pack_obj_chapa.n_paralelo:
                        pos1 = current_chapa_cell_positions.get((s_idx, current_p_index))
                        pos2_v = current_chapa_cell_positions.get((s_idx, next_p_index_v))
                        if pos1 and pos2_v:
                            parallel_conn_points_chapa_internal.append(
                                (pos1['x_start'] + pos1['lx']/2, pos1['y_start'] + pos1['ly'], pos1['z_start'] + pos1['lz']/2)
                            )
                            parallel_conn_points_chapa_internal.append(
                                (pos2_v['x_start'] + pos2_v['lx']/2, pos2_v['y_start'], pos2_v['z_start'] + pos2_v['lz']/2)
                            )
        if parallel_conn_points_chapa_internal:
            full_pack_data.append(go.Scatter3d(
                x=[p[0] for p in parallel_conn_points_chapa_internal],
                y=[p[1] for p in parallel_conn_points_chapa_internal],
                z=[p[2] for p in parallel_conn_points_chapa_internal],
                mode='lines',
                line=dict(color='blue', width=10), 
                name=f'Conexão Paralelo Interna (Chapa {i+1})',
                showlegend=False 
            ))


        # Encontrar os terminais das chapas para as conexões paralelas EXTERNAS
        term_x_coord = chapa_width / 2
        term_z_coord = chapa_height / 2

        # Ponto frontal central para os barramentos (negativo)
        all_chapa_negative_terminals.append((term_x_coord, chapa_offset_y, term_z_coord))
        # Ponto traseiro central para os barramentos (positivo)
        all_chapa_positive_terminals.append((term_x_coord, chapa_offset_y + chapa_depth, term_z_coord))

    # Adicionar as conexões paralelas (barramentos) entre as chapas
    if len(all_chapa_positive_terminals) > 1:
        full_pack_data.append(go.Scatter3d(
            x=[p[0] for p in all_chapa_positive_terminals],
            y=[p[1] for p in all_chapa_positive_terminals],
            z=[p[2] for p in all_chapa_positive_terminals],
            mode='lines',
            line=dict(color='lime', width=12), 
            name='Barramento Positivo (Pack Final)',
            showlegend=True
        ))

    if len(all_chapa_negative_terminals) > 1:
        full_pack_data.append(go.Scatter3d(
            x=[p[0] for p in all_chapa_negative_terminals],
            y=[p[1] for p in all_chapa_negative_terminals],
            z=[p[2] for p in all_chapa_negative_terminals],
            mode='lines',
            line=dict(color='cyan', width=12), 
            name='Barramento Negativo (Pack Final)',
            showlegend=True
        ))
    
    fig = go.Figure(data=full_pack_data)

    total_width_pack = chapa_width
    total_depth_pack = chapa_depth * num_chapas + (num_chapas - 1) * spacing_between_chapas_cm 
    total_height_pack = chapa_height

    x_title = f'Largura Total ({total_width_pack:.2f} cm)'
    y_title = f'Profundidade Total ({total_depth_pack:.2f} cm)'
    z_title = f'Altura Total ({total_height_pack:.2f} cm)'

    fig.update_layout(
        scene=dict(
            xaxis_title=x_title,
            yaxis_title=y_title,
            zaxis_title=z_title,
            aspectmode='data',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5) 
            ),
            xaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
            yaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
            zaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
        ),
        title=f'Visualização 3D do Pack de {num_chapas} Chapas de Bateria em Paralelo',
        showlegend=True,
        margin=dict(l=0, r=0, b=0, t=50),
        scene_camera=dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=1.5, y=1.5, z=1.5) 
        )
    )

    try:
        fig.show()
        print("--- Visualização 3D incorporada na saída da célula do Colab. ---")
    except Exception as e:
        print(f"Não foi possível exibir o gráfico no Colab: {e}")
        print("Tentando salvar o gráfico como um arquivo HTML...")
        file_name = f'estrutura_bateria_pack_final_{num_chapas}_chapas.html'
        pio.write_html(fig, file=file_name, auto_open=True)
        print(f"Gráfico salvo como {file_name}. Por favor, baixe e abra este arquivo no seu navegador.")


def find_optimal_battery_pack_by_volume(
    max_width_cm, max_depth_cm, max_height_cm,
    min_voltage_V=None, max_voltage_V=None, max_energy_kwh=None, max_weight_kg=None
):
    """
    Encontra a bateria com maior ENERGIA (kWh) que se encaixa
    em um volume de empacotamento máximo especificado,
    e opcionalmente atende a uma tensão nominal mínima e máxima, energia máxima e peso máximo.

    Args:
        max_width_cm (float): A largura máxima disponível em centímetros.
        max_depth_cm (float): A profundidade máxima disponível em centímetros.
        max_height_cm (float): A altura máxima disponível em centímetros.
        min_voltage_V (float, optional): Tensão nominal mínima requerida em Volts. Defaults to None.
        max_voltage_V (float, optional): Tensão nominal máxima permitida em Volts. Defaults to None.
        max_energy_kwh (float, optional): Energia máxima permitida em kWh. Defaults to None.
        max_weight_kg (float, optional): Peso máximo permitido em kg. Defaults to None.

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
    if max_energy_kwh:
        print(f"    e Energia Máxima de {max_energy_kwh:.2f} kWh")
    if max_weight_kg:
        print(f"    e Peso Máximo de {max_weight_kg:.2f} kg ---")
    else:
        print("---")


    optimal_energy_kwh = -1.0 # Objetivo: Maximizar energia, não só Ah
    optimal_pack_details = None
    optimal_pack_object = None
    optimal_volume_utilization = -1.0 # Novo: Acompanhar o quão bem o volume é utilizado

    cell_types = ['Li-ion', 'LiFePO4']

    max_n_serie = 200 
    max_n_paralelo = 30 

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
                    
                    peso_total_kg = current_pack.calcular_peso_total()


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
                                    # 5. Verificar restrição de peso máximo (NOVA CONDIÇÃO)
                                    if max_weight_kg is None or peso_total_kg <= max_weight_kg:
                                    
                                        # Calcular a utilização do volume para a bateria atual
                                        current_volume_utilization = current_pack_volume_cm3 / max_total_volume_cm3
                                        
                                        # 6. Priorizar a maior energia, e em caso de empate, o maior aproveitamento de volume
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

    if optimal_pack_details:
        print("\n--- Bateria Otimizada para UMA Chapa Encontrada ---")
        print(f"Tipo de Célula: {optimal_pack_details['tipo_celula']}")
        print(f"Configuração por Chapa: {optimal_pack_details['n_serie']}S{optimal_pack_details['n_paralelo']}P")
        print(f"Tensão Nominal por Chapa: {optimal_pack_details['tensao_nominal_V']:.2f} V")
        print(f"Carga Total por Chapa: {optimal_pack_details['carga_total_Ah']:.2f} Ah")
        print(f"Energia Total por Chapa: {optimal_pack_details['energia_total_kWh']:.2f} kWh")
        print(f"Peso Total por Chapa: {optimal_pack_details['peso_total_kg']:.2f} kg")
        print(f"Volume de Empacotamento por Chapa: {optimal_pack_details['volume_empacotamento_cm3']:.2f} cm³")
        print(f"Dimensões por Chapa (L x P x A): {optimal_pack_details['largura_cm']:.2f} cm x {optimal_pack_details['profundidade_cm']:.2f} cm x {optimal_pack_details['altura_cm']:.2f} cm")
        print(f"Utilização do Volume Disponível por Chapa: {optimal_pack_details['utilizacao_volume_%']:.2f} %") # Novo
        return optimal_pack_details, optimal_pack_object
    else:
        print(f"Nenhuma combinação de bateria encontrada que se encaixe nas dimensões máximas de {max_width_cm:.2f}x{max_depth_cm:.2f}x{max_height_cm:.2f} cm³ para UMA Chapa.")
        if min_voltage_V:
            print(f"Considerando a tensão nominal mínima de {min_voltage_V:.2f} V.")
        if max_voltage_V:
            print(f"Considerando a tensão nominal máxima de {max_voltage_V:.2f} V.")
        if max_energy_kwh:
            print(f"Considerando a energia máxima de {max_energy_kwh:.2f} kWh.")
        if max_weight_kg:
            print(f"Considerando o peso máximo de {max_weight_kg:.2f} kg.")
        return None, None

# --- Exemplo de uso para encontrar a bateria ideal ---
if __name__ == "__main__":

    max_width_m_chapa = 0.3 
    max_depth_m_chapa = 0.1   
    max_height_m_chapa = 0.5

    # Converte para centímetros para uso no código
    max_width_cm_chapa = max_width_m_chapa * 100
    max_depth_cm_chapa = max_depth_m_chapa * 100
    max_height_cm_chapa = max_height_m_chapa * 100

    # Restrição de peso para o PACK TOTAL (50 kg)
    # Dividimos pela quantidade de chapas para que a otimização de UMA chapa respeite a proporção
    max_total_weight_kg = 50
    max_weight_kg_chapa_limit = max_total_weight_kg / 3 

    # Restrições de tensão nominal para o PACK COMPLETO (Série é mantida por chapa)
    min_fsae_voltage = 280.0 
    max_fsae_voltage = 600.0 

    # --- Otimiza para UMA ÚNICA CHAPA de bateria ---
    optimal_pack_details_chapa, optimal_pack_obj_chapa = find_optimal_battery_pack_by_volume(
        max_width_cm=max_width_cm_chapa,
        max_depth_cm=max_depth_cm_chapa,
        max_height_cm=max_height_cm_chapa,
        min_voltage_V=min_fsae_voltage, 
        max_voltage_V=max_fsae_voltage, 
        max_weight_kg=max_weight_kg_chapa_limit # Passando o limite de peso por chapa
    )

    if optimal_pack_details_chapa:
        # --- Cálculos para o PACK FINAL com 3 CHAPAS ---
        num_chapas = 3

        final_n_serie = optimal_pack_details_chapa['n_serie']
        # O Np do pack final é a soma dos Np de cada chapa em paralelo
        final_n_paralelo_total_pack = optimal_pack_details_chapa['n_paralelo'] * num_chapas 

        final_pack_voltage = optimal_pack_details_chapa['tensao_nominal_V'] # Tensão se mantém (chapas em paralelo)
        final_charge_capacity_ah = optimal_pack_details_chapa['carga_total_Ah'] * num_chapas
        final_energy_kwh = optimal_pack_details_chapa['energia_total_kWh'] * num_chapas
        final_weight_kg = optimal_pack_details_chapa['peso_total_kg'] * num_chapas

        # Dimensões do pack final: profundidade é 3x a da chapa individual + espaçamentos
        # Consideramos um pequeno espaçamento entre as chapas para visualização
        spacing_between_chapas_cm_visual = 5 
        final_width_cm = optimal_pack_details_chapa['largura_cm']
        final_depth_cm_visual = optimal_pack_details_chapa['profundidade_cm'] * num_chapas + (num_chapas - 1) * spacing_between_chapas_cm_visual
        final_height_cm = optimal_pack_details_chapa['altura_cm']
        
        # Para o cálculo de volume, usamos a profundidade real (sem espaçamento visual extra entre as chapas)
        final_depth_cm_real = optimal_pack_details_chapa['profundidade_cm'] * num_chapas
        final_volume_cm3_real = final_width_cm * final_depth_cm_real * final_height_cm

        # O volume máximo disponível para as 3 chapas deve ser o volume total do "contêiner" (30x30x50cm)
        max_total_volume_for_3_chapas_cm3 = max_width_cm_chapa * (max_depth_cm_chapa * num_chapas) * max_height_cm_chapa

        print("\n" + "="*50)
        print("--- DETALHES DO PACK FINAL (3 CHAPAS IDÊNTICAS) ---")
        print("="*50)
        print(f"Tipo de Célula: {optimal_pack_details_chapa['tipo_celula']}")
        print(f"Configuração do Pack Completo: {final_n_serie}S{final_n_paralelo_total_pack}P") 
        print(f"Tensão Nominal do Pack Completo: {final_pack_voltage:.2f} V")
        print(f"Carga Total do Pack Completo: {final_charge_capacity_ah:.2f} Ah")
        print(f"Energia Total do Pack Completo: {final_energy_kwh:.2f} kWh")
        print(f"Peso Total do Pack Completo: {final_weight_kg:.2f} kg")
        print(f"Volume de Empacotamento Real do Pack Completo: {final_volume_cm3_real:.2f} cm³") 
        print(f"Dimensões do Pack Completo (L x P x A): {final_width_cm:.2f} cm x {final_depth_cm_real:.2f} cm x {final_height_cm:.2f} cm")
        print("="*50 + "\n")

        # Visualiza a estrutura 3D de UMA ÚNICA CHAPA
        optimal_pack_obj_chapa.plot_estrutura_3d(title_suffix="(Uma Chapa)")

        # Visualiza as TRÊS CHAPAS JUNTAS
        plot_pack_3_chapas_paralelas(optimal_pack_obj_chapa, num_chapas)

    else:
        print("\nNão foi possível encontrar uma configuração de bateria para uma chapa que atenda aos critérios.")