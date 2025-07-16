import numpy as np
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback
from scipy.spatial.distance import pdist
from rascunho_propriedadestubos import Estrutura #Substituir pelo nome correto da main
import matplotlib.pyplot as plt
import os
from datetime import datetime
import time

from rascunho_propriedadestubos import K_global

def init_worker(optimizer_instance):
    """
    Inicializador de cada worker: atribui a instância global OPTIMIZER
    """
    global OPTIMIZER
    OPTIMIZER = optimizer_instance

class ChassisDEOptimizer:   
    """
    Otimizador de geometria de chassi tubular usando Differential Evolution (DE).

    Atributos:
    - base_nodes: np.ndarray, coordenadas iniciais dos nós de um lado (x>=0).
    - base_connections: list of tuple, pares de índices representando arestas base.
    - mandatory_indices: índices de nós que têm deslocamento limitado a radius_mand.
    - pop_size, F, CR, max_gens: parâmetros do DE (tamanho pop., taxa de mutação, crossover, gerações).
    - radius_mand, radius_opt: raios máximos de deslocamento para nós mandatórios e opcionais.
    - tipos_tubos: lista de strings com perfis de tubo permitidos.
    """

    def __init__(
        self,
        base_nodes: np.ndarray,
        base_connections: list,
        mandatory_indices: list,
        pop_size: int = 2,
        F: float = 0.5,
        CR: float = 0.9,
        max_generations: int = 5,
        radius_mand: float = 0.025,
        radius_opt: float = 0.05,
        use_parallel: bool = True,
        n_workers: int = None,
    ):
        """
        Inicializa o otimizador.

        Entradas:
        - base_nodes: array (n,3) de floats.
        - base_connections: lista de tuplas (i,j).
        - mandatory_indices: lista de inteiros.
        - pop_size, F, CR, max_generations: parâmetros DE.
        - radius_mand, radius_opt: floats definindo limites de deslocamento.

        Retorno:
        - Nenhum (configura atributos internos).
        """
        self.base_nodes = base_nodes.copy()
        self.n = base_nodes.shape[0]
        self.n_tubes = len(base_connections)
        self.base_connections = base_connections
        self.mandatory = set(mandatory_indices)
        self.radius_mand = radius_mand
        self.radius_opt = radius_opt

        self.tipos_tubos = ['Tubo A', 'Tubo B', 'Tubo C', 'Tubo D']

        self.pop_size = pop_size
        self.F = F
        self.CR = CR
        self.max_gens = max_generations

        self.dim_coords = 3 * self.n
        self.dim_tubes = self.n_tubes
        self.dim = self.dim_coords + self.dim_tubes

        self.use_parallel = use_parallel
        self.n_workers = max(1, os.cpu_count() if n_workers is None else n_workers)

        self.is_central = np.isclose(self.base_nodes[:, 0], 0.0)
        self.central_alignment_map = self._create_alignment_map()

        self.global_decoded_cache = {}
        
    def _create_alignment_map(self,exceptions=[15,20]):
        """
        Cria um dicionário mapeando cada nó central ao nó lateral que ele deve espelhar.
        exceptions não são mapeados.
        """
        alignment_map = {}
        central_nodes = [i for i in range(self.n) if self.is_central[i]]
        for i, j in self.base_connections:
            if i in central_nodes and not self.is_central[j] and i not in exceptions:
                alignment_map.setdefault(i, []).append(j)
            elif j in central_nodes and not self.is_central[i] and j not in exceptions:
                alignment_map.setdefault(j, []).append(i)
        return alignment_map

    def apply_central_alignment(self, coords: np.ndarray) -> np.ndarray:
        """
        Força nós centrais a terem o mesmo Y e Z dos nós laterais conectados.
        Para múltiplas conexões, usa a média das coordenadas.
        """
        aligned_coords = coords.copy()
        
        for central_idx, lateral_indices in self.central_alignment_map.items():
            aligned_coords[central_idx, 1] = np.mean(coords[lateral_indices, 1])
            aligned_coords[central_idx, 2] = np.mean(coords[lateral_indices, 2])
            aligned_coords[central_idx, 0] = 0.0
            
        return aligned_coords
    
    def enforce_bounds(self, coords: np.ndarray) -> np.ndarray:
        """
        Aplica limites de deslocamento e arredonda as coordenadas.

        Entrada:
        - coords: array (n,3) de floats (proposta de deslocamento).
        Saída:
        - adjusted: array (n,3) de floats, cada deslocamento limitado a radius_mand ou radius_opt
          e arredondado para 3 casas decimais.
        """
        adjusted = coords.copy()
        for i in range(self.n):
            orig = self.base_nodes[i]

            if self.is_central[i]:
                adjusted[i, 0] = 0.0
                delta = coords[i, 1:] - orig[1:]
                dist = np.linalg.norm(delta)
                r = self.radius_mand if i in self.mandatory else self.radius_opt
                if dist > r:
                    adjusted[i, 1:] = orig[1:] + (delta / dist) * r
            else:
                delta = coords[i] - orig
                dist = np.linalg.norm(delta)
                r = self.radius_mand if i in self.mandatory else self.radius_opt
                if dist > r:
                    adjusted[i] = orig + (delta / dist) * r
        adjusted[self.is_central, 0] = 0.0
        return np.round(adjusted, 3)
        
    def find_new_index(self, old_index):
        """
        Calcula o indice do nó após passar pela otimização
        Retorna o novo indice do nó e do seu espelhado correspondente
        Caso seja um nó central retorna apenas o novo indice do nó
        
        Entradas:
        - old_index: Índice do nó

        Retorno:
        - new_index: Índice do nó após passar pela otimização
        - mirrored_index: índice do nó espelhado correspondente
        - new_central_index: novo index se for um nó central
        """
        new_index = old_index*2
        mirrored_index = old_index*2+1

        for i in range(len(self.base_nodes)):
            if self.is_central[i]:
                first_central=i
                break
        
        new_central_index=(first_central*2)+(old_index-first_central)
        
        return new_central_index if self.is_central[old_index] else new_index,mirrored_index
    
    def validate_min_distance(self, coords: np.ndarray, min_dist: float = 0.05) -> bool:
        """
        Verifica se todas as distâncias entre pares de nós são maiores ou iguais a uma distância mínima.
        Vetorizado usando pdist.

        Parâmetros:
        - coords: np.ndarray (N, 3), coordenadas dos nós.
        - min_dist: float, distância mínima permitida.

        Retorno:
        - bool: True se válido, False se algum par violar a distância mínima.
        """
        return np.all(pdist(coords) >= min_dist)
 
    def decode_individual(self, x: np.ndarray):
        """
        Converte um vetor genotípico em nós completos e lista de elementos.

        Entrada:
        - x: array (dim_coords+dim_tubes,) de floats.
        Saída:
        - nodes_full: array (N,3) de floats com nós de ambos os lados.
        - elements: lista de tuplas (i, j, perfil), com índices nos nodes_full.

        Processos:
        1. Extrai coords e tube_vars do vetor x.
        2. Aplica enforce_bounds às coords.
        3. Para nós com base x≈0, inclui apenas um nó (central).
           Para outros, inclui coord e seu espelho.
        4. Usa mapeamento para gerar conexões simétricas com tipo de tubo.
        """
        # Separa coords e tube_vars
        coords = x[:self.dim_coords].reshape((self.n, 3))
        coords = self.enforce_bounds(coords)
        coords = self.apply_central_alignment(coords)
        tube_vars = x[self.dim_coords:]

        full_nodes = []
        mapping = {}
        for i, coord in enumerate(coords):
            if self.is_central[i]:
                idx = len(full_nodes)
                full_nodes.append(coord)
                mapping[i] = [idx]
            else:
                idx1 = len(full_nodes)
                full_nodes.append(coord)
                mirrored = coord.copy()
                mirrored[0] *= -1
                idx2 = len(full_nodes)
                full_nodes.append(mirrored)
                mapping[i] = [idx1, idx2]

        nodes_full = np.array(full_nodes)
        elements = []
        for idx_conn, (i, j) in enumerate(self.base_connections):
            perfil = self.tipos_tubos[int(np.clip(np.floor(tube_vars[idx_conn]), 0, len(self.tipos_tubos) - 1))]
            ids_i, ids_j = mapping[i], mapping[j]
            if len(ids_i) == 1 and len(ids_j) == 1:
                elements.append((ids_i[0], ids_j[0], perfil))
            elif len(ids_i) == 2 and len(ids_j) == 2:
                elements.append((ids_i[0], ids_j[0], perfil))
                elements.append((ids_i[1], ids_j[1], perfil))
            else:
                cent, lats = (ids_i[0], ids_j) if len(ids_i) == 1 else (ids_j[0], ids_i)
                for lat in lats:
                    elements.append((cent, lat, perfil))

        return nodes_full, elements

    def decode_population(self, pop):
        decoded = []
        for x in pop:
            key = hash(x.tobytes())
            if key not in self.global_decoded_cache:
                nodes, elems = self.decode_individual(x)
                if not self.validate_min_distance(nodes):
                    continue
                self.global_decoded_cache[key] = (nodes, elems)
            decoded.append(self.global_decoded_cache[key])
        return decoded 
   
    def initialize_individual(self) -> np.ndarray:
        """
        Gera um indivíduo válido.

        Saída:
        - x: array (dim,) de floats, contendo coords arredondadas e tube_vars iniciais.
        """
        while True:
            deltas = np.random.normal(size=(self.n, 3))
            deltas[self.is_central, 0] = 0.0
            norms = np.linalg.norm(deltas, axis=1, keepdims=True)
            radii = (np.random.rand(self.n, 1) ** (1/3)) * self.radius_opt
            coords = self.base_nodes + (deltas / norms) * radii
            coords = self.enforce_bounds(coords)
            coords = self.apply_central_alignment(coords)
            tube_vars = np.random.uniform(0, len(self.tipos_tubos), size=(self.n_tubes,))
            x = np.concatenate([coords.reshape(-1), tube_vars])
            nodes, _ = self.decode_individual(x)
            if self.validate_min_distance(nodes):
                return x

    def initialize_population(self):
        pop = []
        seen = set()
        while len(pop) < self.pop_size:
            x = self.initialize_individual()
            key = hash(x.tobytes())
            if key in seen:
                continue
            seen.add(key)
            pop.append(x)
        pop = np.array(pop)
        decoded = self.decode_population(pop)
        nodes_list, elems_list = zip(*decoded)
        results = self.parallel_evaluate(nodes_list, elems_list)
        fitness = np.array([res[0] for res in results])
        return pop, fitness

    def mutate_and_crossover(self, pop):
        new_pop = pop.copy()
        for i in range(self.pop_size):
            a, b, c = np.random.choice(np.delete(np.arange(self.pop_size), i), 3, replace=False)
            v = pop[a] + self.F * (pop[b] - pop[c])
            v[0:self.n * 3].reshape(self.n, 3)[self.is_central, 0] = 0.0
            j_rand = np.random.randint(self.dim)
            u = np.where((np.random.rand(self.dim) < self.CR) | (np.arange(self.dim) == j_rand), v, pop[i])
            if np.allclose(u, pop[i], atol=1e-4): continue
            key = hash(u.tobytes())
            if key not in self.global_decoded_cache:
                try:
                    nodes, _ = self.decode_individual(u)
                    if not self.validate_min_distance(nodes): continue
                    self.global_decoded_cache[key] = (nodes, _)
                except: continue
            new_pop[i] = u
        return new_pop
    
    def parallel_evaluate(self, nodes_list, elements_list):
        if self.use_parallel:
            with ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=init_worker,
                initargs=(self,)
            ) as executor:
                results = list(executor.map(evaluate, nodes_list, elements_list))
        else:
            results = [evaluate(nodes_list[i], elements_list[i]) for i in range(len(nodes_list))]
        return results  
      
    def optimize(self):
        """
        Executa o loop principal de otimização, com avaliação paralela.

        Retorna:
        - best_solution: tupla (nodes, elements)
        - best_cost: float
        - best_mass, best_KT, best_KF: métricas físicas
        - history: dicionário com histórico da evolução
        - duration: float (tempo total em segundos)
        """
        print("Otimização Iniciada")
        pop, fit = self.initialize_population()

        history = {
            'best_fit': [], 'avg_fit': [], 'std_dev': [],
            'best_mass': [], 'avg_mass': [],
            'best_KT': [], 'avg_KT': [],
            'best_KF': [], 'avg_KF': []
        }

        convergence_count = 0
        convergence_threshold = 3
        start_time = time.time()
        # Abre pool de processos apenas uma vez
        executor = None
        if self.use_parallel:
            executor = ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=init_worker,
                initargs=(self,)
            )

        # Loop de gerações
        for gen in range(1, self.max_gens + 1):
            t_gen = time.time()
            
            t1 = time.time()
            cand_pop = self.mutate_and_crossover(pop)
            t2 = time.time()

            decoded = self.decode_population(cand_pop)
            t3 = time.time()

            nodes_list, elems_list = zip(*decoded)
            t4 = time.time()
            # Avaliação (paralela ou sequencial)
            if executor:
                results = list(executor.map(evaluate, nodes_list, elems_list))
            else:
                results = [evaluate(n, e) for n, e in zip(nodes_list, elems_list)]
            t5 = time.time()

            t6 = time.time()
            new_pop = pop.copy()
            new_fit = fit.copy()

            new_mass = np.full(self.pop_size, np.nan)
            new_KT   = np.full(self.pop_size, np.nan)
            new_KF   = np.full(self.pop_size, np.nan)

            for i, (f, m, kt, kf) in enumerate(results):
                # Armazena métricas para análise estatística
                new_mass[i] = m
                new_KT[i]   = kt
                new_KF[i]   = kf

                # Substituição por seleção de DE
                if f <= fit[i]:
                    new_pop[i] = cand_pop[i]
                    new_fit[i] = f

            pop, fit = new_pop, new_fit

            best_idx = np.argmin(fit)
            best_fit = fit[best_idx]
            avg_fit  = np.mean(fit)
            std_dev  = float(np.std(pop))

            history['best_fit'].append(best_fit)
            history['avg_fit'].append(avg_fit)
            history['std_dev'].append(std_dev)
            history['best_mass'].append(new_mass[best_idx])
            history['avg_mass'].append(np.nanmean(new_mass))
            history['best_KT'].append(new_KT[best_idx])
            history['avg_KT'].append(np.nanmean(new_KT))
            history['best_KF'].append(new_KF[best_idx])
            history['avg_KF'].append(np.nanmean(new_KF))

            if std_dev < 0.02:
                convergence_count += 1
            else:
                convergence_count = 0

            t7 = time.time()

            print(
                f"[Gen {gen}] mutate: {t2 - t1:.2f}s | decode: {t3 - t2:.2f}s | eval: {t5 - t4:.2f}s | stats: {t7 - t6:.2f}s | total: {t7 - t_gen:.2f}s"
            )

            status = (f"Gen {gen}/{self.max_gens} — Best={best_fit:.4e} Std={std_dev:.4f}"
                    if convergence_count < convergence_threshold
                    else f"Convergência após {gen} gerações")
            print(status, end='\r')

            if convergence_count >= convergence_threshold:
                print()
                break

        # Fecha pool de processos
        if executor:
            executor.shutdown()

        duration = time.time() - start_time

        best_idx       = np.argmin(fit)
        best_solution  = self.decode_individual(pop[best_idx])
        best_cost      = fit[best_idx]
        best_mass      = history['best_mass'][-1]
        best_KT        = history['best_KT'][-1]
        best_KF        = history['best_KF'][-1]

        return best_solution, best_cost, best_mass, best_KT, best_KF, history, duration

    def plotar(self, individuo, save_path=None):
        """
        Plota a geometria do chassi evoluído em 3D.

        Entrada:
        - individuo: tupla (nodes, elements)
        - save_path: caminho opcional para salvar a imagem

        Saída:
        - exibe o plot e/ou salva a imagem se o caminho for fornecido.
        """
        nodes, elements = individuo
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')

        num_nodes = len(nodes)
        elements_valid = [(i, j, t) for i, j, t in elements if 0 <= i < num_nodes and 0 <= j < num_nodes]

        xs, ys, zs = zip(*nodes)
        ax.scatter(ys, xs, zs, s=25, c='black')
        for i, j, t in elements_valid:
            ni, nj = nodes[i], nodes[j]
            ax.plot([ni[1], nj[1]], [ni[0], nj[0]], [ni[2], nj[2]], 'b-')

        ax.set_xlabel('Y')
        ax.set_ylabel('X')
        ax.set_zlabel('Z')
        ax.set_box_aspect([3, 1, 2])
        plt.title("Chassi Evoluído")

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Visualização 3D salva em: {save_path}")

        plt.show()

    def plot_metrics(self, history, save_path=None, show=True):
        """
        Gera gráfico de evolução das métricas: massa, rigidez torcional e flexional.

        Entradas:
        - history: dict com as listas 'best_mass', 'best_KT', 'best_KF'.
        - save_path: caminho opcional para salvar o PNG.
        - show: bool, se True exibe o plot.
        """
        plt.figure(figsize=(10, 6))
        gens = range(len(history['best_mass']))
        
        # Configura eixo primário (massa)
        ax1 = plt.gca()
        ax1.plot(gens, history['best_mass'], 'b-', linewidth=2, label='Massa (kg)')
        ax1.set_xlabel('Geração')
        ax1.set_ylabel('Massa (kg)', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        
        # Configura eixo secundário (rigidezes)
        ax2 = ax1.twinx()
        ax2.plot(gens, history['best_KT'], 'r-', linewidth=2, label='Rigidez Torcional (KT)')
        ax2.plot(gens, history['best_KF'], 'g-', linewidth=2, label='Rigidez Flexional (KF)')
        ax2.set_ylabel('Rigidez (N·m/rad ou N/m)', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        
        plt.title('Evolução das Métricas do Melhor Indivíduo')
        
        # Unificar legendas
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
        
        plt.grid(True)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            print(f"Gráfico de métricas salvo em: {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()

    def plot_convergence(self, history, save_path=None, show=True):
        """
        Gera gráfico de convergência: fitness vs geração e std_dev.

        Entradas:
        - history: dict com 'best_fit', 'avg_fit', 'std_dev'.
        - save_path: caminho opcional para salvar o PNG.
        - show: bool, se True exibe o plot.
        """
        plt.figure(figsize=(10, 6))
        
        # Gráfico duplo eixo Y
        ax1 = plt.gca()
        ax2 = ax1.twinx()
        
        # Curva de fitness
        ax1.plot(history['best_fit'], 'b-', linewidth=2, label='Melhor Fitness')
        ax1.plot(history['avg_fit'], 'b--', alpha=0.7, label='Fitness Médio')
        ax1.set_ylabel('Fitness', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        ax1.set_yscale('log')
        
        # Curva de desvio padrão
        ax2.plot(history['std_dev'], 'r-', linewidth=2, label='Desvio Padrão')
        ax2.axhline(y=0.2, color='g', linestyle='--', label='Limite Convergência')
        ax2.set_ylabel('Desvio Padrão', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        
        plt.title('Progresso da Otimização')
        plt.xlabel('Geração')
        plt.grid(True)
        
        # Unificar legendas
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            print(f"Gráfico de convergência salvo em: {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()

    def plot_tradeoff(self, history, axises=('mass','KT'), save_path=None, show=True):
        x_key = {'mass':'best_mass','KT':'best_KT','KF':'best_KF'}[axises[0]]
        y_key = {'mass':'best_mass','KT':'best_KT','KF':'best_KF'}[axises[1]]
        x, y = history[x_key], history[y_key]
        plt.figure(figsize=(8,6))
        sc = plt.scatter(x, y, c=range(len(x)), cmap='viridis', s=50)
        plt.colorbar(sc, label='Geração')
        plt.xlabel(axises[0].upper())
        plt.ylabel(axises[1].upper())
        plt.title(f"Trade-off: {axises[0]} × {axises[1]}")
        plt.grid(True)
        if save_path: plt.savefig(save_path, dpi=300)
        if show: plt.show()
        else: plt.close()

    def plot_fitness_vs_metric(self, history, metric='mass', save_path=None, show=True):
        mkey = f"best_{metric}"
        x, y = history[mkey], history['best_fit']
        plt.figure(figsize=(8,6))
        sc = plt.scatter(x, y, c=range(len(x)), cmap='plasma', s=40)
        plt.colorbar(sc, label='Geração')
        plt.xlabel(metric.upper())
        plt.ylabel('Fitness')
        plt.title(f"Fitness × {metric.upper()}")
        plt.grid(True)
        if save_path: plt.savefig(save_path, dpi=300)
        if show: plt.show()
        else: plt.close()

    def plot_all(self, history, results_dir, show=False):
        os.makedirs(results_dir, exist_ok=True)
        self.plot_convergence(history, save_path=os.path.join(results_dir,'convergencia.png'), show=show)
        self.plot_metrics(history, save_path=os.path.join(results_dir,'evolucao_metricas.png'), show=show)
        #self.plot_tradeoff(history, axises=('mass','KT'), save_path=os.path.join(results_dir,'tradeoff_mass_KT.png'), show=show)
        #self.plot_tradeoff(history, axises=('mass','KF'), save_path=os.path.join(results_dir,'tradeoff_mass_KF.png'), show=show)
        #self.plot_fitness_vs_metric(history,'mass', save_path=os.path.join(results_dir,'fitness_vs_mass.png'), show=show)
        #self.plot_fitness_vs_metric(history,'KT',   save_path=os.path.join(results_dir,'fitness_vs_KT.png'), show=show)
        #self.plot_fitness_vs_metric(history,'KF',   save_path=os.path.join(results_dir,'fitness_vs_KF.png'), show=show)
        print(f"Todos os gráficos salvos em: {results_dir}")

    def save_solution(self, nodes, elements, file_path, fitness=None, mass=None, KT=None, KF=None, duration=None):
        """
        Salva relatório completo em TXT, incluindo duração e valores do melhor indivíduo.
        """
        with open(file_path, 'w') as f:
            f.write("SOLUÇÃO FINAL DO CHASSI\n")
            f.write("="*60 + "\n\n")
            f.write(f"{self.max_gens} Gerações, {self.pop_size} individuos por população \n")
            if duration is not None:
                hrs, rem = divmod(duration, 3600)
                mins, secs = divmod(rem, 60)
                f.write(f"Duração da Otimização: {int(hrs)}h {int(mins)}m {secs:.1f}s\n")
            if fitness is not None:
                f.write(f"Melhor Fitness: {fitness:.6e}\n")
            if mass is not None:
                f.write(f"Massa: {mass:.3f} kg\n")
            if KT is not None and KF is not None:
                f.write(f"Rigidez Torcional (KT): {KT:.3e} N·m/rad\n")
                f.write(f"Rigidez Flexional (KF): {KF:.3e} N/m\n")
            f.write("\nNÓS (X, Y, Z):\n")
            for i, node in enumerate(nodes):
                f.write(f"{i:3d}: {node[0]:6.3f}, {node[1]:6.3f}, {node[2]:6.3f}\n")
            f.write("\nELEMENTOS (Conexões):\n")
            for i, (ni, nj, tp) in enumerate(elements):
                f.write(f"{i:3d}: Nó {ni} - Nó {nj} | Perfil: {tp}\n")
        print(f"Relatório completo salvo em: {file_path}")

    def export_solution(self, nodes: np.ndarray, elements: list, directory: str, filename: str = "melhor_solucao"):
        """
        Exporta a solução para arquivos CSV dentro do diretório especificado.
        
        Args:
            nodes: array (N,3) com coordenadas dos nós
            elements: lista de tuplas (i, j, perfil)
            directory: diretório onde os arquivos serão salvos
            filename: nome base para os arquivos (sem extensão)
        """
        # Garante que o diretório existe
        os.makedirs(directory, exist_ok=True)
        
        # Cria caminhos completos
        nodes_path = os.path.join(directory, f"{filename}_nodes.csv")
        elements_path = os.path.join(directory, f"{filename}_elements.csv")
        
        # Exportar nós
        with open(nodes_path, 'w') as f_nodes:
            f_nodes.write("index,x,y,z\n")
            for i, node in enumerate(nodes):
                f_nodes.write(f"{i},{node[0]:.6f},{node[1]:.6f},{node[2]:.6f}\n")
        
        # Exportar elementos
        with open(elements_path, 'w') as f_elements:
            f_elements.write("node_i,node_j,perfil\n")
            for elem in elements:
                f_elements.write(f"{elem[0]},{elem[1]},{elem[2]}\n")
        
        print(f"Solucao exportada para:")
        print(f" - Nodes: {nodes_path}")
        print(f" - Elements: {elements_path}")

def evaluate(nodes,elements) -> float:
    """
    Avalia o custo de um indivíduo.

    Entrada:
    - nodes
    - elements
    Saída:
    - float, valor de penalidade (fitness).

    Processo:
    1. Decodifica para nodes e elements.
    2. Checa distância mínima.
    3. Monta e analisa pela classe Estrutura (FEA).
    4. Calcula penalidade via penalidade_chassi.
    """
    try:
        t0 = time.perf_counter()

        penalty = 0
        # Penalidade em relação ao FH e FHB
        fronthoop_node = nodes[ChassisDEOptimizer.find_new_index(20)]       #declara o nó do fronthoop com indice novo
        fhb_node = nodes[ChassisDEOptimizer.find_new_index(5)[0]]           #declara o nó de um front hoop bracing com indice novo
        dist_fh_fhb = fronthoop_node[2] - fhb_node[2]                       #declara a distância no eixo z entre esses dois nós 
        if dist_fh_fhb > 0.05:                                              #condição retirada do regulamento
            penalty += ((dist_fh_fhb - 0.05) ** 2 ) * 1e6                       #aplicação de penalidade

        # Penalidade em relação ao MH e o MHB
        mainhoop_node = nodes[ChassisDEOptimizer.find_new_index(15)]                            #declara o nó do main hoop com indice novo
        mhb_node = nodes[ChassisDEOptimizer.find_new_index(14)[0]]                              #declara o nó do main hoop bracing com indice novo não espelhado
        deltax_mh_mhb = mainhoop_node[0] - mhb_node[0]                                          #diferença das coordenadas "x" em ambos os nós
        deltay_mh_mhb = mainhoop_node[1] - mhb_node[1]                                          #diferença das coordenadas "y" em ambos os nós
        deltaz_mh_mhb = mainhoop_node[2] - mhb_node[2]                                          #diferença das coordenadas "x" em ambos os nós
        dist_mh_mhb = np.sqrt(deltax_mh_mhb ** 2 + deltay_mh_mhb ** 2 + deltaz_mh_mhb ** 2)     #declara a distância entre os dois nós      
        if dist_mh_mhb > 0.16:                                                                  #condição retirada do regulamento
            penalty += ((dist_mh_mhb - 0.16) ** 2) * 1e6                                        #aplicação de penalidade
        
        # Penalidade ângulo entre o Main Hoop e o Main Hoop Bracing
        x_porcao_mh = nodes[ChassisDEOptimizer.find_new_index(14)][0] - nodes[ChassisDEOptimizer.find_new_index(6)][0]                                          #coordenada x do vetor formado pelos nós do elemento da porção do mainhoop analisada
        y_porcao_mh = nodes[ChassisDEOptimizer.find_new_index(14)][1] - nodes[ChassisDEOptimizer.find_new_index(6)][1]                                          #coordenada y do vetor formado pelos nós do elemento da porção do mainhoop analisada
        z_porcao_mh = nodes[ChassisDEOptimizer.find_new_index(14)][2] - nodes[ChassisDEOptimizer.find_new_index(6)][2]                                          #coordenada z do vetor formado pelos nós do elemento da porção do mainhoop analisada
        x_mhb = nodes[ChassisDEOptimizer.find_new_index(14)][0] - nodes[ChassisDEOptimizer.find_new_index(16)][0]                                               #coordenada x do vetor formado pelos nós do elemento do Main Hoop Bracing
        y_mhb = nodes[ChassisDEOptimizer.find_new_index(14)][1] - nodes[ChassisDEOptimizer.find_new_index(16)][1]                                               #coordenada y do vetor formado pelos nós do elemento do Main Hoop Bracing
        z_mhb = nodes[ChassisDEOptimizer.find_new_index(14)][2] - nodes[ChassisDEOptimizer.find_new_index(16)][2]                                               #coordenada z do vetor formado pelos nós do elemento do Main Hoop Bracing
        vetor_porcao_mh = (x_porcao_mh, y_porcao_mh, z_porcao_mh)                                                                                               #vetor formado pelos nós do elemento da porção do mainhoop analisada
        vetor_mhb = (x_mhb, y_mhb, z_mhb)                                                                                                                       #vetor formado pelos nós do elemento do Main Hoop Bracing
        modulo_vetor_porcao_mh = np.sqrt(vetor_porcao_mh[0] ** 2 + vetor_porcao_mh[1] ** 2 + vetor_porcao_mh[2] ** 2 )                                          #módulo do vetor formado pelos nós do elemento da porção do mainhoop analisada
        modulo_vetor_mhb = np.sqrt(vetor_mhb[0] ** 2 + vetor_mhb[1] ** 2 + vetor_mhb[2] ** 2 )                                                                  #módulo do vetor formado pelos nós do elemento do Main Hoop Bracing
        produto_escalar_mh_porcao_e_mhb = (vetor_porcao_mh[0] * vetor_mhb[0]) + (vetor_porcao_mh[1] * vetor_mhb[1]) + (vetor_porcao_mh[2] * vetor_mhb[2])       #produto escalar entre os dois vetores criados
        cos_theta_mh_mhb = produto_escalar_mh_porcao_e_mhb / (modulo_vetor_porcao_mh * modulo_vetor_mhb)                                                        #valor do cosseno do ângulo formado pelos dois vetores
        theta_mh_mhb = np.degrees(np.acos(cos_theta_mh_mhb))                                                                                                    #valor do ângulo formado pelos dois vetores
        if theta_mh_mhb < 30:                                                                                                                                   #condição retirada do regulamento
            penalty += ((theta_mh_mhb - 30) ** 2) * 1e6                                                                                                         #aplicação da penalidade
        

        # Instanciamento da estrutura
        estrutura = Estrutura(elements, nodes)
        t1 = time.perf_counter()

        estrutura.matrizes_global()
        t2 = time.perf_counter()

        fixed = list(range(6))
        Fg = np.zeros(estrutura.num_dofs)
        for i in [7, 8, 21, 22]:
            Fg[i * 6 + 2] = (65 * 9.81) / 4

        _, _, frequencies = estrutura.modal_analysis()
        t3 = time.perf_counter()

        disp = estrutura.static_analysis(Fg, fixed)
        t4 = time.perf_counter()

        strain = estrutura.compute_strain(disp)
        stresses = estrutura.compute_stress(strain)
        von = estrutura.compute_von_mises(stresses)
        t5 = time.perf_counter()

        massa = estrutura.mass()
        LFnode,RFnode = otimizador.find_new_index(3) #Left front suspension node index (3)
        LRnode,RRnode = otimizador.find_new_index(11) #Left rear suspension node index (11)
        LCnode,RCnode = otimizador.find_new_index(7) #Left central node index (7)
        KT, KF = estrutura.compute_Kf_Kt(K_global,LFnode,LRnode,RFnode,RRnode,LCnode,RCnode)
        t6 = time.perf_counter()

        penalty += penalidade_chassi(KT, KF, massa, von, frequencies)
        t7 = time.perf_counter()

        print(f"""[EVALUATE]
        - Instância:     {t1 - t0:.4f}s
        - Montagem K/M:  {t2 - t1:.4f}s
        - Modal:         {t3 - t2:.4f}s
        - Estática:      {t4 - t3:.4f}s
        - Tensão/Von:    {t5 - t4:.4f}s
        - Massa/Rigidez: {t6 - t5:.4f}s
        - Penalidade:    {t7 - t6:.4f}s
        - Total:         {t7 - t0:.4f}s
        """)

        return penalty, massa, KT, KF

    except Exception:
        traceback.print_exc()
        return 1e9, 1e9, 0.0, 0.0

def penalidade_chassi(KT, KF, massa, tensoes, frequencias):
    """
    Calcula penalidade total do chassi.

    Entradas:
    - KT, KF: floats de rigidezes.
    - massa: float.
    - tensoes: array de floats.
    - frequencias: array de frequencias naturais.
    Saída:
    - penalidade_total: float.
    """

    # Limites e parâmetros
    KT_min = 1e7                # Rigidez torcional mínima (N·m/rad)
    KF_min = 1e6                # Rigidez flexão mínima (N/m)
    massa_ideal = 23            # Massa alvo (kg)
    K_mola = 5e5                # Constante da mola do amortecedor (N/m)
    tensao_adm = 250e6          # Tensão admissível do material (Pa)
    alpha = 0.5                 # Fator de sensibilidade exponencial
    beta = 10                   # Fator de escala logarítmica
    freq_motor = 4800           # Rotação do motor dada em (RPM) e convertida para (Hz)
    penalidade_total = 0

    # 1. Rigidez Torcional (Função Exponencial)
    if KT < KT_min:
        deficit = (KT_min - KT) / KT_min
        # Penalidade cresce exponencialmente com o déficit
        penalidade_total += np.exp(alpha * deficit) - 1

    # 2. Rigidez em Flexão (Função Logarítmica)
    if KF < KF_min:
        deficit = (KF_min - KF) / KF_min
        # Penalidade logarítmica: suave para pequenas violações, forte para grandes
        penalidade_total += beta * np.log(1 + deficit)

    # 3. Massa (Função Híbrida)
    if massa > massa_ideal:
        excesso = (massa - massa_ideal) / massa_ideal
        # Combina resposta linear inicial com crescimento exponencial
        penalidade_total += excesso + np.exp(alpha * excesso) - 1

    # 4. Compatibilidade com Mola (Lógica Aprimorada)
    ratio_KT = K_mola / KT if KT > 0 else float('inf')
    ratio_KF = K_mola / KF if KF > 0 else float('inf')
    
    if ratio_KT > 25 or ratio_KF > 25:
        # Penalidade proporcional ao nível de incompatibilidade
        violacao = max(ratio_KT/25, ratio_KF/25) - 1
        penalidade_total += 100 * violacao**2

    # 5. Tensões (Abordagem Baseada em Risco)
    tensao_max = max(tensoes)
    if tensao_max > tensao_adm:
        # Penalidade exponencial para tensões acima do admissível
        excesso = (tensao_max - tensao_adm) / tensao_adm
        penalidade_total += np.exp(5 * excesso) - 1
    
    # Penalidade por distribuição desigual de tensões (logarítmica)
    razao_tensoes = np.ptp(tensoes) / np.mean(tensoes) if np.mean(tensoes) > 0 else 0
    penalidade_total += np.log(1 + razao_tensoes)

    # 7. Frequências Naturais em Zona de Ressonância com o Motor (Função exponencial)
    freq_crit_min = (0.95*freq_motor)/60        #Frequência crítica mínima convertida para Hz
    freq_crit_max = (1.05*freq_motor)/60        #Frequência crítica máxima convertida para Hz

    for f in frequencias:
        # Penalidade exponencial mais severa para frequências próximas a do motor
        if freq_crit_min <= f <= freq_crit_max:
            severidade = np.exp(alpha * (1 - abs((f - freq_motor) / freq_motor)))
            penalidade_total += 100 * severidade 
    return penalidade_total * 100  # Fator de escala global

if __name__ == "__main__":
    nodes = np.array([[-0.181,  0.000,  0.360],             #00
    [-0.181,  0.000,  0.050],                               #01
    [-0.280,  0.275,  0.240],                               #02
    [-0.285,  0.495,  0.045],                               #03
    [-0.285,  0.555,  0.270],                               #04
    [-0.212,  0.555,  0.550],                               #05
    [-0.293,  1.370,  0.250],                               #06
    [-0.268,  1.350,  0.000],                               #07
    [-0.268,  1.495,  0.015],                               #08
    [-0.302,  1.670,  0.240],                               #09
    [-0.271,  1.665,  0.030],                               #10
    [-0.271,  1.835,  0.070],                               #11
    [-0.183,  2.015,  0.285],                               #12
    [-0.183,  2.015,  0.060],                               #13
    [-0.170,  1.400,  0.965],                               #14
    [ 0.000,  1.410,  1.105],                               #15
    [-0.293,  1.950,  0.250],                               #16
    [0.000,  0.000,  0.360],                                #conection between the 2 sides
    [0.000,  0.000,  0.050],                                #conection between the 2 sides
    [0.000,  0.495,  0.045],                                #conection between the 2 sides
    [0.000,  0.555,  0.550],                                #conection between the 2 sides
    [0.000,  1.350,  0.000],                                #conection between the 2 sides
    [0.000,  1.495,  0.015],                                #conection between the 2 sides
    [0.000,  1.665,  0.030],                                #conection between the 2 sides
    [0.000,  1.835,  0.070],                                #conection between the 2 sides
    [0.000,  2.015,  0.285],                                #conection between the 2 sides
    [0.000,  2.015,  0.060],                                #conection between the 2 sides
    [0.000,  1.370,  0.250]])                               #conection between the 2 sides 
                                
    connections = [(0,1)  ,(0,2)  ,(1,2)  ,(0,5)  ,(2,5)  ,(2,4)  ,(2,3)  ,(1,3)  ,
                (3,7)  ,(3,4)  ,(4,7)  ,(4,6)  ,(5,6)  ,(7,6)  ,(7,8)  ,(6,8)  ,
                (6,9)  ,(8,9)  ,(8,10) ,(9,10) ,(9,11) ,(10,11),(11,16),(16,13),
                (11,12),(11,13),(12,13),(12,16),(16,14),(14,15),(14,6) ,(9,16) ,
                (4,5)  ,(0,17) ,(1,18) ,(3,19) ,(5,20) ,(7,21) ,(8,22) ,(10,23),
                (11,24),(12,25),(6,27),(13,26)]

    indices = [0,1,3,4,5,6,7,8,11,16,14,15,20]

    # Criar diretório para resultados
    timestamp = datetime.now().strftime("%Y-%m-%d %H%M")
    max_gen = 5
    pop_size = 5
    otimizador = ChassisDEOptimizer(
        base_nodes=nodes,
        base_connections=connections,
        mandatory_indices=indices,
        pop_size=pop_size,
        max_generations=max_gen,
        use_parallel=True,     # ou False, se quiser forçar modo serial
        n_workers=4            # defina quantos núcleos usar
    )
    
    results_dir = f"Resultados_Otimizacao_LOCAL__{timestamp}_{max_gen}GEN_{pop_size}POP"
    os.makedirs(results_dir, exist_ok=True)

    best_indiv, best_cost, best_mass, best_KT, best_KF, history,duration = otimizador.optimize()
    print(f"\nRESULTADOS FINAIS:")
    print(f"Melhor custo: {best_cost:.6e}")
    print(f"Massa: {best_mass:.2f} kg")
    print(f"Rigidez Torcional (KT): {best_KT:.2e} N·m/rad")
    print(f"Rigidez Flexional (KF): {best_KF:.2e} N/m")
    nodes_final, elements_final = best_indiv

    # Salvar solução em arquivo TXT
    solution_path = os.path.join(results_dir, "solucao_final.txt")
    otimizador.save_solution(
        nodes_final, elements_final, solution_path,
        fitness=best_cost, mass=best_mass, KT=best_KT, KF=best_KF, duration=duration)
    
    # Geração de todos os gráficos de uma única vez
    otimizador.plot_all(history, results_dir, show=False)
    # visualização 3D final
    otimizador.plotar(best_indiv, save_path=os.path.join(results_dir,'visualizacao_3d.png'))

    # Exportar para arquivos
    otimizador.export_solution(
        nodes_final, 
        elements_final, 
        directory=results_dir,
        filename="solucao_otimizada")
    
