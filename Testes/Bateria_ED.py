import numpy as np
import math
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor
import time
import os
import sys # Importado para verificar se está no ambiente do Google Colab

# --- IMPORTAÇÕES NECESSÁRIAS PARA PLOTAGEM 3D ---
try:
    import plotly.graph_objects as go
    import plotly.io as pio
    PLOTLY_DISPONIVEL = True
    
    # --- NOVO: CONFIGURAÇÃO LIMPA DO PLOTLY ---
    # Define o modo de renderização para abrir no navegador, a menos que esteja no Colab.
    if 'google.colab' in sys.modules:
        pio.renderers.default = 'colab' # Padrão para visualização no output do Colab
    else:
        # Tenta forçar a abertura em um novo browser.
        # Se você estiver em um terminal sem ambiente gráfico (headless),
        # pode ser necessário mudar para 'browser' ou 'json'.
        pio.renderers.default = 'browser'
    # ------------------------------------------

except ImportError:
    print("Atenção: A biblioteca Plotly não está instalada. Execute 'pip install plotly' para a visualização 3D funcionar.")
    class DummyPlotly:
        def Figure(self, data): return None
        def Mesh3d(self, **kwargs): return None
        def Scatter3d(self, **kwargs): return None
        def update_layout(self, **kwargs): return
        def show(self): return
        def write_html(self, fig, file, auto_open=False): pass 
    go = DummyPlotly()
    pio = DummyPlotly()
    PLOTLY_DISPONIVEL = False
# ----------------------------------------


# --- CÓDIGO BatteryPack INTEGRADO COM PLOTAGEM ---
class BatteryPack:
    def __init__(self, tipo_celula, n_serie, n_paralelo, soc_inicial=1.0):
        self.parametros = self.definir_celula(tipo_celula)
        if self.parametros is None:
            raise ValueError(f"Tipo de célula desconhecida: {tipo_celula}")
        self.E0 = self.parametros['E0']
        self.Q = self.parametros['Q']
        self.peso = self.parametros['peso']
        self.lado_cubo_x = self.parametros['dim_x_cm']
        self.lado_cubo_y = self.parametros['dim_y_cm']
        self.lado_cubo_z = self.parametros['dim_z_cm']
        self.n_serie = n_serie
        self.n_paralelo = n_paralelo
        self.carga_total = self.Q * n_paralelo 
        self._plot_layout_factors = None

    def definir_celula(self, tipo):
        catalogo = {
            'Li-ion': {
                'E0': 3.7, 'K': 0.005, 'Q': 3.0, 'A': 0.1, 'B': 0.01, 'R': 0.02,
                'peso': 0.045, 'volume': 0.0165, 'dim_x_cm': 1.8, 'dim_y_cm': 1.8, 'dim_z_cm': 6.5 
            },
            'LiFePO4': {
                'E0': 3.2, 'K': 0.003, 'Q': 2.8, 'A': 0.05, 'B': 0.015, 'R': 0.008,
                'peso': 0.07, 'volume': 0.028, 'dim_x_cm': 2.6, 'dim_y_cm': 2.6, 'dim_z_cm': 6.5 
            }
        }
        return catalogo.get(tipo)
        
    def calcular_peso_total(self):
        return self.n_serie * self.n_paralelo * self.peso

    def calcular_tensao_nominal(self):
        return self.E0 * self.n_serie

    def calcular_dimensoes_empacotamento_cm(self):
        lx, ly, lz = self.lado_cubo_x, self.lado_cubo_y, self.lado_cubo_z
        spacing_factor = 1.1 
          
        num_parallel_in_x_base = math.ceil(self.n_paralelo**0.5)
        num_parallel_in_y_base = math.ceil(self.n_paralelo / num_parallel_in_x_base)

        width_base_dim = num_parallel_in_x_base * (lx * spacing_factor)
        depth_base_dim = num_parallel_in_y_base * (ly * spacing_factor)
        height_base_dim = lz * spacing_factor 

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
        """
        Gera uma visualização 3D da estrutura do pack de baterias usando Plotly.
        Usa o renderer configurado para exibição limpa (browser ou colab).
        """
        if not PLOTLY_DISPONIVEL:
             return
             
        # Removido print de início para evitar mensagens no terminal

        lx, ly, lz = self.lado_cubo_x, self.lado_cubo_y, self.lado_cubo_z
        
        visual_spacing_factor = 1.5
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
                                lighting=dict(ambient=0.4, diffuse=0.8, specular=0.2, roughness=0.5, fresnel=0.0),
                                lightposition=dict(x=0, y=0, z=10000)
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
                line=dict(color='darkorange', width=10), 
                name='Conexão Série',
                showlegend=True
            ))

        parallel_conn_points = []
        for s_idx in range(self.n_serie): 
            for p_row in range(num_parallel_in_y_base):
                for p_col in range(num_parallel_in_x_base - 1):
                    current_p_index = p_row * num_parallel_in_x_base + p_col
                    next_p_index_h = p_row * num_parallel_in_x_base + (p_col + 1)
                    if current_p_index < self.n_paralelo and next_p_index_h < self.n_paralelo:
                        pos1 = _cell_positions.get((s_idx, current_p_index))
                        pos2_h = _cell_positions.get((s_idx, next_p_index_h))
                        if pos1 and pos2_h:
                            parallel_conn_points.append(
                                (pos1['x_start'] + pos1['lx'], pos1['y_start'] + pos1['ly']/2, pos1['z_start'] + pos1['lz']/2)
                            )
                            parallel_conn_points.append(
                                (pos2_h['x_start'], pos2_h['y_start'] + pos2_h['ly']/2, pos2_h['z_start'] + pos2_h['lz']/2)
                            )
            for p_col in range(num_parallel_in_x_base):
                for p_row in range(num_parallel_in_y_base - 1):
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
                line=dict(color='blue', width=10),
                name='Conexão Paralelo (Interna à Chapa)',
                showlegend=True
            ))

        fig = go.Figure(data=data)

        total_width, total_depth, total_height = self.calcular_dimensoes_empacotamento_cm()
        x_title = f'Largura Total ({total_width:.2f} cm)'
        y_title = f'Profundidade Total ({total_depth:.2f} cm)'
        z_title = f'Altura Total ({total_height:.2f} cm)'

        fig.update_layout(
            scene=dict(
                xaxis_title=x_title, yaxis_title=y_title, zaxis_title=z_title,
                aspectmode='data', camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)),
                xaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
                yaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
                zaxis=dict(showbackground=True, backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(200,200,200,0.5)"),
            ),
            title=f'Estrutura 3D do Pack de Baterias ({self.n_serie}S{self.n_paralelo}P) - Layout Otimizado {title_suffix}',
            showlegend=True,
            margin=dict(l=0, r=0, b=0, t=50)
        )

        try:
            # fig.show() respeita o pio.renderers.default definido acima
            fig.show()
        except Exception as e:
            # Caso fig.show() falhe (ambiente headless, por exemplo)
            # Mantemos o fallback para salvar o HTML, mas sem abrir.
            print(f"A exibição interativa falhou ({e}). O gráfico foi salvo como HTML.")
            file_name = f'estrutura_bateria_{self.n_serie}S{self.n_paralelo}P{title_suffix.replace(" ", "_").replace("(", "").replace(")", "")}.html'
            pio.write_html(fig, file=file_name, auto_open=False) 


    def plot_pack_3_chapas_paralelas(self, num_chapas=3): pass


# --- ESTRUTURA GLOBAL E FUNÇÕES AUXILIARES ---

global OPTIMIZER
OPTIMIZER = None

def init_worker(optimizer_instance):
    global OPTIMIZER
    OPTIMIZER = optimizer_instance

# ... (Funções penalidade_bateria, decode_and_evaluate e evaluate - Inalteradas)
def penalidade_bateria(V_nom, carga_Ah, energia_kWh, massa_kg,
                         max_w, max_d, max_h, V_min_req, V_max_perm, max_mass_perm,
                         real_w, real_d, real_h):
    
    penalty_V = 0.0
    V_MIN_RESTRICAO = V_min_req 
    V_MAX_RESTRICAO = V_max_perm

    if V_nom < V_MIN_RESTRICAO:
        penalty_V += 1e6 * (V_MIN_RESTRICAO - V_nom) / V_MIN_RESTRICAO
    if V_nom > V_MAX_RESTRICAO:
        penalty_V += 1e6 * (V_nom - V_MAX_RESTRICAO) / V_MAX_RESTRICAO

    penalty_dim = 0.0
    if real_w > max_w:
        penalty_dim += 1e6 * ((real_w - max_w) / max_w)**2
    if real_d > max_d:
        penalty_dim += 1e6 * ((real_d - max_d) / max_d)**2
    if real_h > max_h:
        penalty_dim += 1e6 * ((real_h - max_h) / max_h)**2

    penalty_mass = 0.0
    if max_mass_perm is not None and massa_kg > max_mass_perm:
        penalty_mass += 1e6 * ((massa_kg - max_mass_perm) / max_mass_perm)**2
        
    if penalty_V > 1e-4 or penalty_dim > 1e-4 or penalty_mass > 1e-4:
        return (penalty_V + penalty_dim + penalty_mass), massa_kg, V_nom, energia_kWh, real_w * real_d * real_h

    fitness_energia = 1.0 / (energia_kWh + 1e-6)

    volume_max = max_w * max_d * max_h
    volume_real = real_w * real_d * real_h
    
    volume_utilization_penalty = 1.0 * (volume_max - volume_real) / (volume_max + 1e-6)
    
    W_ENERGIA = 1000.0
    W_UTILIZACAO = 10.0
    
    fitness_total = W_ENERGIA * fitness_energia + W_UTILIZACAO * volume_utilization_penalty
    
    return fitness_total, massa_kg, V_nom, energia_kWh, volume_real

def decode_and_evaluate(x: np.ndarray):
    global OPTIMIZER
    
    cell_type_idx = int(np.clip(np.round(x[0]), 0, 1))
    n_serie = int(np.clip(np.round(x[1]), OPTIMIZER.min_ns, OPTIMIZER.max_ns))
    n_paralelo = int(np.clip(np.round(x[2]), OPTIMIZER.min_np, OPTIMIZER.max_np))
    
    cell_types = ['Li-ion', 'LiFePO4']
    tipo_celula = cell_types[cell_type_idx]
    
    try:
        pack = BatteryPack(tipo_celula, n_serie, n_paralelo)
        
        pack_w_cm, pack_d_cm, pack_h_cm = pack.calcular_dimensoes_empacotamento_cm()
        
        massa_kg = pack.calcular_peso_total()
        V_nom = pack.calcular_tensao_nominal()
        carga_Ah = pack.carga_total
        energia_kWh = (V_nom * carga_Ah) / 1000
        
        fitness, mass, V, E, Vol = penalidade_bateria(
            V_nom, carga_Ah, energia_kWh, massa_kg,
            OPTIMIZER.max_w_cm, OPTIMIZER.max_d_cm, OPTIMIZER.max_h_cm,
            OPTIMIZER.min_v, OPTIMIZER.max_v, OPTIMIZER.max_mass,
            pack_w_cm, pack_d_cm, pack_h_cm
        )
        
        details = {
            'tipo_celula': tipo_celula,
            'n_serie': n_serie,
            'n_paralelo': n_paralelo,
            'largura_cm': pack_w_cm,
            'profundidade_cm': pack_d_cm,
            'altura_cm': pack_h_cm,
            'tensao_nominal_V': V_nom,
            'carga_total_Ah': carga_Ah,
            'energia_total_kWh': energia_kWh,
            'massa_total': mass
        }
        
        return fitness, mass, V, E, Vol, details

    except Exception:
        return 1e12, 1e12, 0.0, 0.0, 0.0, None

def evaluate(x):
    global OPTIMIZER
    
    if OPTIMIZER is None:
        return 1e12, 1e12, 0.0, 0.0, 0.0, None
        
    cell_type_idx = int(np.clip(np.round(x[0]), 0, 1))
    n_serie = int(np.clip(np.round(x[1]), OPTIMIZER.min_ns, OPTIMIZER.max_ns))
    n_paralelo = int(np.clip(np.round(x[2]), OPTIMIZER.min_np, OPTIMIZER.max_np))
    key = (cell_type_idx, n_serie, n_paralelo)
    
    if key in OPTIMIZER.global_decoded_cache:
        return OPTIMIZER.global_decoded_cache[key]
    
    result = decode_and_evaluate(x)
    
    if result[0] < 1e11:
          OPTIMIZER.global_decoded_cache[key] = result
        
    return result

# --- CLASSE DE OTIMIZAÇÃO (BatteryDEOptimizer) ---

class BatteryDEOptimizer: 
    
    def __init__(
        self,
        max_w_cm: float, max_d_cm: float, max_h_cm: float,
        min_v: float, max_v: float, max_mass: float = None,
        min_ns: int = 1, max_ns: int = 200, min_np: int = 1, max_np: int = 200,
        pop_size: int = 50, F: float = 0.6, CR: float = 0.8,
        max_generations: int = 500, use_parallel: bool = True, n_workers: int = None
    ):
        self.max_w_cm = max_w_cm
        self.max_d_cm = max_d_cm
        self.max_h_cm = max_h_cm
        self.min_v = min_v
        self.max_v = max_v
        self.max_mass = max_mass
        
        self.min_ns = min_ns
        self.max_ns = max_ns
        self.min_np = min_np
        self.max_np = max_np

        self.pop_size = pop_size
        self.F = F
        self.CR = CR
        self.max_gens = max_generations

        self.dim = 3 
        self.bounds = np.array([
            [0, 1.999], 
            [self.min_ns, self.max_ns], 
            [self.min_np, self.max_np] 
        ])

        self.use_parallel = use_parallel
        self.n_workers = max(1, os.cpu_count() if n_workers is None else n_workers)
        self.global_decoded_cache = {}
        
    def initialize_individual(self) -> np.ndarray:
        x = np.random.uniform(self.bounds[:, 0], self.bounds[:, 1])
        return x

    def initialize_population(self):
        pop = []
        while len(pop) < self.pop_size:
            x = self.initialize_individual()
            pop.append(x)
        pop = np.array(pop)
        
        results = self.parallel_evaluate(pop)
        fitness = np.array([res[0] for res in results])
        
        return pop, fitness

    def mutate_and_crossover(self, pop):
        new_pop = pop.copy()
        for i in range(self.pop_size):
            idxs = [idx for idx in range(self.pop_size) if idx != i]
            a, b, c = np.random.choice(idxs, 3, replace=False)
            
            v = pop[a] + self.F * (pop[b] - pop[c])
            
            j_rand = np.random.randint(self.dim)
            u = np.where((np.random.rand(self.dim) < self.CR) | (np.arange(self.dim) == j_rand), v, pop[i])
            
            u = np.clip(u, self.bounds[:, 0], self.bounds[:, 1])
            
            new_pop[i] = u
        return new_pop
    
    def parallel_evaluate(self, pop):
        if self.use_parallel:
            with ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=init_worker,
                initargs=(self,)
            ) as executor:
                results = list(executor.map(evaluate, pop))
        else:
            results = [evaluate(x) for x in pop]
            
        return results 
        
    def optimize(self):
        print("Otimização de Bateria Iniciada (Differential Evolution)")
        
        global OPTIMIZER
        OPTIMIZER = self
        
        pop, fit = self.initialize_population()
        
        best_idx = np.argmin(fit)
        best_result = evaluate(pop[best_idx])
        
        history = {
            'best_fit': [best_result[0]], 'avg_fit': [np.mean(fit)], 'std_dev': [np.std(fit)],
            'best_E': [best_result[3]], 'best_V': [best_result[2]], 'best_Mass': [best_result[1]]
        }

        start_time = time.time()
        
        executor = None
        if self.use_parallel:
            executor = ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=init_worker,
                initargs=(self,)
            )

        for gen in range(1, self.max_gens + 1):
            
            cand_pop = self.mutate_and_crossover(pop)
            
            if executor:
                results = list(executor.map(evaluate, cand_pop))
            else:
                results = [evaluate(x) for x in cand_pop]
                
            new_pop = pop.copy()
            new_fit = fit.copy()

            for i, res in enumerate(results):
                f = res[0]
                if f <= fit[i]:
                    new_pop[i] = cand_pop[i]
                    new_fit[i] = f

            pop, fit = new_pop, new_fit
            
            current_best_idx = np.argmin(fit)
            current_best_fit = fit[current_best_idx]
            current_best_result = evaluate(pop[current_best_idx])
            
            history['best_fit'].append(current_best_fit)
            history['avg_fit'].append(np.mean(fit))
            history['std_dev'].append(np.std(fit))
            history['best_E'].append(current_best_result[3])
            history['best_V'].append(current_best_result[2])
            history['best_Mass'].append(current_best_result[1])

            # MENSAGEM LIMPA DE PROGRESSO
            print(
                f"Gen {gen}/{self.max_gens} | Best F: {current_best_fit:.4e} | E: {current_best_result[3]:.2f} kWh | V: {current_best_result[2]:.1f} V | Mass: {current_best_result[1]:.2f} kg",
                end='\r'
            )
        
        if executor:
            executor.shutdown()
            
        duration = time.time() - start_time
        
        final_best_idx = np.argmin(fit)
        final_x = pop[final_best_idx]
        final_fit, final_mass, final_V, final_E, final_Vol, final_details = evaluate(final_x)
        
        final_pack_obj = BatteryPack(
            final_details['tipo_celula'],
            final_details['n_serie'],
            final_details['n_paralelo']
        )
        
        return final_pack_obj, final_details, history, duration


if __name__ == "__main__":
    
    # Dimensões máximas para o PACK COMPLETO (cm)
    MAX_W_PACK_CM = 30.0
    MAX_D_PACK_CM = 30.0 
    MAX_H_PACK_CM = 50.0

    # Restrições de tensão nominal (para o PACK COMPLETO)
    MIN_FSAE_VOLTAGE = 220.0
    MAX_FSAE_VOLTAGE = 600.0
    
    # Restrição de peso para o PACK COMPLETO (50 kg)
    MAX_MASS_PACK_KG = 50.0 
    
    # Configurações do DE
    MAX_GEN = 100
    POP_SIZE = 50
    
    # Limites para Ns e Np 
    MAX_NS = 200
    MAX_NP = 190
    
    # --- EXECUÇÃO DA OTIMIZAÇÃO ---
    otimizador = BatteryDEOptimizer(
        max_w_cm=MAX_W_PACK_CM, max_d_cm=MAX_D_PACK_CM, max_h_cm=MAX_H_PACK_CM,
        min_v=MIN_FSAE_VOLTAGE, max_v=MAX_FSAE_VOLTAGE, max_mass=MAX_MASS_PACK_KG,
        max_ns=MAX_NS, max_np=MAX_NP,
        pop_size=POP_SIZE, max_generations=MAX_GEN, 
        use_parallel=True, n_workers=os.cpu_count()
    )
    
    optimal_pack_obj, optimal_pack_details, history, duration = otimizador.optimize()

    # --- SAÍDA DOS RESULTADOS OTIMIZADOS NO TERMINAL (Pack Único) ---
    
    # Limpa a linha de progresso antes de imprimir a saída final
    print(" "*100, end='\r') 
    print("\n" + "="*70)
    print("--- RESULTADO OTIMIZADO PARA O PACK DE BATERIA COMPLETO ---")
    print(f"Duração da Otimização: {duration:.2f} segundos ({MAX_GEN * POP_SIZE} avaliações teóricas)")
    print("="*70)
    print(f"Melhor Fitness (Menor Custo): {history['best_fit'][-1]:.4e}")
    print(f"Tipo de Célula Otimizada: {optimal_pack_details['tipo_celula']}")
    print(f"Configuração Otimizada (S x P): {optimal_pack_details['n_serie']}S{optimal_pack_details['n_paralelo']}P")
    print(f"Tensão Nominal: {optimal_pack_details['tensao_nominal_V']:.2f} V")
    print(f"Energia Total: {optimal_pack_details['energia_total_kWh']:.2f} kWh")
    print(f"Massa Total: {optimal_pack_details['massa_total']:.2f} kg (Limite: {MAX_MASS_PACK_KG:.2f} kg)")
    print(f"Dimensões (L x P x A): {optimal_pack_details['largura_cm']:.2f} x {optimal_pack_details['profundidade_cm']:.2f} x {optimal_pack_details['altura_cm']:.2f} cm")
    print("="*70)

    # --- VISUALIZAÇÃO 3D DO RESULTADO OTIMIZADO ---
    # Isto irá abrir a aba do navegador (ou visualizar no Colab/Jupyter)
    optimal_pack_obj.plot_estrutura_3d(title_suffix=f"({optimal_pack_details['tipo_celula']})")