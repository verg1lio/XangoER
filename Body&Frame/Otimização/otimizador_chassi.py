import numpy as np
from Estrutura_Tipada import Estrutura
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import os
from datetime import datetime
import time
from itertools import repeat

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
        fixed_nodes: list = None,
        pop_size: int = 50,
        F: float = 0.6,
        CR: float = 0.8,
        max_generations: int = 200,
        radius_mand: float = 0.03,
        radius_opt: float = 0.06,
        use_parallel: bool = True,
        n_workers: int = None,
        tubos_fixos: dict = None,
        connect_with_mirror: list = None,
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

        self.tipos_tubos = ['Tubo A', 'Tubo B', 'Tubo C', 'Tubo D', 'Tubo E']
        self.tubos_fixos = tubos_fixos or {}
        self.connect_with_mirror = connect_with_mirror or []
        self.n_extra_espelhos = len(self.connect_with_mirror)
        self.total_tubes = self.n_tubes + self.n_extra_espelhos

        self.fixed_indices = set(fixed_nodes) if fixed_nodes else set()

        self.pop_size = pop_size
        self.F = F
        self.CR = CR
        self.max_gens = max_generations

        self.dim_coords = 3 * self.n
        self.dim_tubes = self.n_tubes
        self.dim = self.dim_coords + self.total_tubes

        self.use_parallel = use_parallel
        self.n_workers = max(1, os.cpu_count() if n_workers is None else n_workers)

        self.is_central = np.isclose(self.base_nodes[:, 0], 0.0)
        self.central_alignment_map = self._create_alignment_map()

        self.global_decoded_cache = {}
        
    def _create_alignment_map(self,exceptions=[]):
        """
        Gera um mapeamento hierárquico associando nós centrais aos seus nós periféricos conectados.

        Este método varre a lista de conexões base (`base_connections`) para identificar
        relacionamentos diretos onde um nó é classificado como central (`is_central`) e o
        outro não. O resultado agrupa os nós periféricos sob seus respectivos nós centrais.

        Args:
            exceptions (list[int], optional): Lista de índices de nós centrais que devem ser
                excluídos do mapeamento, mesmo que possuam conexões válidas.
                O padrão é [18, 19].

        Returns:
            dict[int, list[int]]: Um dicionário representando o mapa de alinhamento, onde:
                - As chaves (int) são os índices dos nós centrais.
                - Os valores (list[int]) são listas contendo os índices dos nós não-centrais
                conectados a cada chave.
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
        Aplica um alinhamento geométrico aos nós centrais, ajustando suas coordenadas com base na média dos vizinhos laterais.

        Este método reposiciona os nós identificados no mapa de alinhamento (`self.central_alignment_map`).
        Para cada nó central, as coordenadas Y e Z são redefinidas como a média aritmética das coordenadas
        correspondentes dos nós laterais conectados. A coordenada X é explicitamente forçada para 0.0,
        alinhando o nó ao plano de simetria central.

        Args:
            coords (np.ndarray): Matriz de coordenadas de entrada com dimensão (N, 3),
                onde as colunas representam as coordenadas X, Y e Z, respectivamente.

        Returns:
            np.ndarray: Uma nova matriz de coordenadas (cópia da original) contendo as
            posições atualizadas dos nós centrais, mantendo os nós não listados no mapa inalterados.
        """
        aligned_coords = coords.copy()
        
        for central_idx, lateral_indices in self.central_alignment_map.items():
            aligned_coords[central_idx, 1] = np.mean(coords[lateral_indices, 1])
            aligned_coords[central_idx, 2] = np.mean(coords[lateral_indices, 2])
            aligned_coords[central_idx, 0] = 0.0
            
        return aligned_coords
    
    def enforce_bounds(self, coords: np.ndarray) -> np.ndarray:
        """
            Aplica restrições geométricas de deslocamento máximo e fixação de nós, retornando coordenadas válidas.

            Este método valida as coordenadas propostas comparando-as com as posições originais (`self.base_nodes`).
            A aplicação de limites segue uma lógica hierárquica:

            1. Nós Fixos (`fixed_indices`): São ignorados pelo cálculo de raio e rigidamente redefinidos para sua posição original.
            2. Nós Centrais (`is_central`): O deslocamento é calculado e limitado apenas no plano 2D (eixos Y e Z), forçando o eixo X em 0.0.
            3. Demais Nós: O deslocamento é calculado e limitado no espaço 3D (esférico).

            Se o deslocamento (delta) exceder o raio permitido (`radius_mand` para nós obrigatórios ou `radius_opt` para opcionais),
            o vetor é normalizado e reescalado para o limite do raio.

            Args:
                coords (np.ndarray): Matriz (n, 3) contendo as coordenadas propostas para a nova configuração.

            Returns:
                np.ndarray: Uma nova matriz (n, 3) com as coordenadas ajustadas, garantindo que nenhum nó
                viole seus limites de deslocamento ou restrições de plano, arredondada para 3 casas decimais.
            """
        adjusted = coords.copy()
        for i in range(self.n):
            # Se o nó for fixo, reseta para a posição original
            if i in self.fixed_indices:
                adjusted[i] = self.base_nodes[i]
                continue # Pula para o próximo nó, ignorando raio e validação

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

    def validate_min_distance(self, coords: np.ndarray, min_dist: float = 0.12) -> bool:
        """
            Valida a restrição de proximidade física entre todos os pares únicos de nós.

            Este método utiliza `scipy.spatial.distance.pdist` para calcular eficientemente a
            distância euclidiana entre todos os pares possíveis de coordenadas (cálculo vetorizado).
            Sua função é garantir a integridade estrutural, evitando que nós fiquem excessivamente
            próximos ou sobrepostos.

            Args:
                coords (np.ndarray): Matriz de coordenadas (N, 3) contendo as posições atuais dos nós.
                min_dist (float, optional): O limiar mínimo de separação euclidiana permitido entre
                    quaisquer dois nós. O padrão é 0.12.

            Returns:
                bool: Retorna `True` se a distância entre *todos* os pares de nós for maior ou
                igual a `min_dist`. Retorna `False` se houver pelo menos uma violação (colisão).
            """
        return np.all(pdist(coords) >= min_dist)
 
    def decode_individual(self, x: np.ndarray):
        """
            Decodifica o vetor de decisão (genótipo) para gerar a geometria completa e a topologia da estrutura (fenótipo).

            Este método realiza três operações principais:
            1. **Processamento Geométrico:** Extrai as coordenadas, aplica restrições de limites (`enforce_bounds`)
            e alinhamento (`apply_central_alignment`).
            2. **Expansão de Simetria:** Gera a estrutura completa espelhando nós não-centrais no eixo X (inversão de sinal).
            Cria um mapeamento entre os índices originais e os índices expandidos (nós reais e virtuais/espelhados).
            
            3. **Atribuição de Elementos:** Define a conectividade e os perfis dos tubos.
            - Converte variáveis contínuas em índices discretos para seleção de `tipos_tubos`.
            - Gerencia conexões entre Central-Central, Lateral-Lateral (duas conexões) e Central-Lateral.
            - Adiciona conexões de reforço transversais (`connect_with_mirror`) entre um nó e seu par espelhado.

            Args:
                x (np.ndarray): Vetor 1D de floats contendo todas as variáveis de otimização.
                    A primeira parte (`:dim_coords`) corresponde às coordenadas, e a segunda
                    parte (`dim_coords:`) corresponde às variáveis de seleção de tubos.

            Returns:
                tuple: Uma tupla contendo:
                    - **nodes_full (np.ndarray):** Matriz (M, 3) com as coordenadas de todos os nós da estrutura completa (original + espelhada).
                    - **elements (list[tuple]):** Lista de tuplas onde cada tupla representa um elemento finito no formato `(nó_origem, nó_destino, perfil_do_tubo)`.
            """
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
            if idx_conn in self.tubos_fixos:
                perfil = self.tubos_fixos[idx_conn]
            else:
                tval = tube_vars[idx_conn]
                t_int = int(np.clip(np.floor(tval), 0, len(self.tipos_tubos) - 1))
                perfil = self.tipos_tubos[t_int]

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
        # Conexões espelhadas adicionais (entre cada nó e seu espelho)
        for i_extra, node_idx in enumerate(self.connect_with_mirror):
            if node_idx in mapping and len(mapping[node_idx]) == 2:
                id_a, id_b = mapping[node_idx]
                key_fixo = f"espelho_{node_idx}"
                if key_fixo in self.tubos_fixos:
                    perfil = self.tubos_fixos[key_fixo]
                else:
                    tval = tube_vars[self.n_tubes + i_extra]
                    t_int = int(np.clip(np.floor(tval), 0, len(self.tipos_tubos) - 1))
                    perfil = self.tipos_tubos[t_int]
                elements.append((id_a, id_b, perfil))

        return nodes_full, elements

    def decode_population(self, pop):
        """
        Decodifica um lote de indivíduos (população), utilizando memoização para performance e filtragem de validade geométrica.

        Este método itera sobre a população de vetores de decisão. Para evitar cálculos redundantes em otimizações iterativas,
        ele utiliza um hash dos bytes do vetor (`x.tobytes()`) como chave para consultar ou atualizar um cache global
        (`self.global_decoded_cache`).

        O método também atua como um filtro de viabilidade: após decodificar (`decode_individual`), ele valida a geometria
        (`validate_min_distance`). Se o indivíduo violar a distância mínima, ele é descartado e não incluído no retorno.

        Args:
            pop (iterable[np.ndarray]): Uma coleção (geralmente uma lista ou matriz numpy) de vetores de decisão (genótipos) a serem processados.

        Returns:
            list[tuple]: Uma lista contendo apenas os indivíduos válidos e decodificados. Cada item é uma tupla `(nodes, elements)`.
            
            Nota: O comprimento da lista retornada pode ser menor que o da população de entrada (`pop`), pois indivíduos geometricamente inválidos são ignorados.
        """
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
            Gera de maneira aleatoria um candidato inicial (indivíduo) que satisfaz todas as restrições geométricas.

            Este método emprega uma estratégia de 'amostragem por rejeição' para garantir a validade inicial:
            1. Perturbação Geométrica Gera deslocamentos aleatórios a partir dos nós base.
            - Aplica uma distribuição uniforme volumétrica (usando a raiz cúbica `**(1/3)`) para evitar
                agrupamento de pontos no centro da esfera de busca.
            - Zera o deslocamento em X para nós centrais.
            2. Aplicação de Restrições Executa `enforce_bounds` e `apply_central_alignment` para corrigir a geometria.
            3. Inicialização de Variáveis: Sorteia valores uniformes para os tipos de tubos.
            4. Validação: Decodifica o indivíduo e verifica colisões (`validate_min_distance`).
            
            O processo repete indefinidamente (loop `while True`) até que um indivíduo válido seja gerado.

            Returns:
                np.ndarray: Um vetor 1D (genótipo) concatenado, contendo as coordenadas linearizadas
                seguidas pelas variáveis de seleção de tubos.
            """
        while True:
            deltas = np.random.normal(size=(self.n, 3))
            deltas[self.is_central, 0] = 0.0
            norms = np.linalg.norm(deltas, axis=1, keepdims=True)
            radii = (np.random.rand(self.n, 1) ** (1/3)) * self.radius_opt
            coords = self.base_nodes + (deltas / norms) * radii
            coords = self.enforce_bounds(coords)
            coords = self.apply_central_alignment(coords)
            tube_vars = np.random.uniform(0, len(self.tipos_tubos), size=(self.total_tubes,))
            x = np.concatenate([coords.reshape(-1), tube_vars])
            nodes, _ = self.decode_individual(x)
            if self.validate_min_distance(nodes):
                return x

    def initialize_population(self):
        """
            Gera, decodifica e avalia a população inicial do algoritmo de otimização.

            Este método constrói a primeira geração de indivíduos garantindo a diversidade genotípica :
            utiliza um conjunto (`set`) de hashes para assegurar que não existam indivíduos duplicados (clones)
            na população inicial.

            O fluxo de execução é:
            1. Geração: Cria indivíduos iterativamente via `initialize_individual` até atingir `self.pop_size`.
            2. Decodificação em Lote: Converte toda a população de vetores (genótipos) em estruturas físicas (fenótipos).
            3. Avaliação Paralela: Submete as estruturas à simulação física/análise estrutural (`parallel_evaluate`).
            4. Extração de Fitness: Isola os valores da função objetivo dos resultados da avaliação.

            Returns:
                tuple: Uma tupla contendo dois arrays numpy:
                    - pop (np.ndarray):Matriz (pop_size, dimensão_genoma) com os vetores de decisão.
                    - fitness (np.ndarray):Vetor (pop_size,) com os valores de aptidão (score) correspondentes a cada indivíduo.
            """
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
        """
            Gera a população de vetores de teste (Trial Vectors) aplicando operadores de Mutação e Cruzamento da Evolução Diferencial.

            Este método implementa a estratégia **DE/rand/1/bin**:
            1. Mutação: Para cada indivíduo alvo (i), seleciona aleatoriamente três vetores distintos (a, b, c) da população
            e cria um vetor mutante ($V$) seguindo a fórmula:
            V = X_a + F * (X_b - X_c)
            Onde F (`self.F`) é o fator de escala diferencial.
            
            2. Restrição Geométrica: Força explicitamente que as coordenadas X dos nós centrais no vetor mutante sejam 0.0,
            mantendo a consistência com as restrições do problema antes mesmo do cruzamento.

            3.Cruzamento (Crossover Binomial): Combina o vetor mutante (V) com o pai (X_i) para formar o vetor de teste (U),
            baseado na taxa de cruzamento CR (`self.CR`). Garante que pelo menos um gene venha do mutante (`j_rand`).

            4. Validação:
            - Verifica se o vetor de teste é muito similar ao pai (evita estagnação).
            - Decodifica e verifica colisões geométricas (`validate_min_distance`).
            - Se inválido, o pai é mantido na posição correspondente em `new_pop`.

            Args:
                pop (np.ndarray): População atual de vetores de decisão (Genótipos).

            Returns:
                np.ndarray: Uma nova matriz contendo os vetores de teste (candidatos) para os índices onde
                a geração foi bem-sucedida e válida, ou os vetores originais onde a geração falhou.
            """
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
        """
            Organiza a avaliação computacional da população, distribuindo as simulações entre múltiplos processos (CPU-bound).

            Este método gerencia a execução da função de análise estrutural (`evaluate`) para uma lista de candidatos.
            Ele suporta dois modos de operação baseados na flag `self.use_parallel`:
            
            1. Modo Paralelo (Multiprocessing):
            - Utiliza `ProcessPoolExecutor` para criar um pool de trabalhadores.
            - Inicializa cada processo operário com `init_worker` (útil para carregar dados estáticos pesados apenas uma vez por processo, evitando overhead de serialização).
            - Mapeia a função `evaluate` através das listas de nós e elementos de forma assíncrona.
            
            2. Modo Serial:
            - Executa as avaliações sequencialmente em um loop (útil para depuração/debugging ou quando o overhead do paralelismo não compensa).

            Args:
                nodes_list (list[np.ndarray]): Lista contendo as matrizes de coordenadas (nós) para cada indivíduo da população.
                elements_list (list[list[tuple]]): Lista correspondente contendo a topologia (conexões e perfis) de cada indivíduo.
            Returns:
                list: Uma lista contendo os resultados brutos retornados pela função `evaluate` para cada indivíduo.
                A ordem dos resultados corresponde à ordem das listas de entrada.
            """
        if self.use_parallel:
            with ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=init_worker,
                initargs=(self,)
            ) as executor:
                results = list(executor.map(evaluate, nodes_list, elements_list, repeat(self.base_nodes)))

        else:
            results = [evaluate(nodes_list[i], elements_list[i], self.base_nodes, tubos_SAE) for i in range(len(nodes_list))]
        return results  
      
    def optimize(self):
        """
            Executa o ciclo de vida completo do algoritmo de Otimização por Evolução Diferencial (DE).

            Este método atua como o motor principal da classe, coordenando a evolução da população ao longo das gerações.
            O fluxo de execução segue as seguintes etapas:
            
            1. Inicialização Cria a população inicial e instancia o pool de processos (`ProcessPoolExecutor`) 
            para paralelismo persistente (evitando o overhead de recriar processos a cada geração).
            2.Loop Evolutivo Para cada geração:
            - Variação Gera candidatos via `mutate_and_crossover`.
            - Decodificação Transforma genótipos em fenomênos estruturais (`decode_population`).
            - Avaliação Distribui a análise física (FEA) entre os workers (`parallel_evaluate`).
            - Seleção Aplica a seleção determinística/gulosa: o filho substitui o pai apenas se possuir 
                fitness menor ou igual (minimização).
            

        [Image of Differential Evolution flow chart]

            3. Monitoramento Registra métricas estatísticas (Média, Melhores, Desvio Padrão) e métricas 
            físicas (Massa, Rigidez KT/KF) no histórico.
            4. Convergência Verifica se o desvio padrão da aptidão da população caiu abaixo do limiar (0.02) 
            por um número definido de gerações (`convergence_threshold`), interrompendo o loop prematuramente se a população estagnar.

            Returns:
                tuple: Uma tupla contendo os resultados finais e metadados da execução:
                    - best_solution (tuple): A topologia vencedora no formato `(nodes, elements)`.
                    - best_cost (float): O valor da função objetivo da melhor solução.
                    - best_mass (float): Massa da melhor estrutura encontrada.
                    - best_KT (float): Rigidez torcional da melhor estrutura.
                    - best_KF (float): Rigidez à flexão da melhor estrutura.
                    - history (dict): Dicionário contendo a evolução temporal de todas as métricas rastreadas (listas de floats).
                    - duration (float): Tempo total de execução em segundos.
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
                results = list(executor.map(evaluate, nodes_list, elems_list, repeat(self.base_nodes)))
            else:
                results = [
                            evaluate(n, e, self.base_nodes, tubos_SAE) 
                            for n, e in zip(nodes_list, elems_list)
                        ]
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
            std_dev = float(np.std(fit))

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
                f"[Gen {gen}] mutate: {t2 - t1:.2f}s | decode: {t3 - t2:.2f}s | eval: {t5 - t4:.2f}s | stats: {t7 - t6:.2f}s | total: {t7 - t_gen:.2f}s", end="\r")

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
        - save_path: caminho para salvar a imagem

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
        Gera gráficos de evolução das métricas em painéis separados para maior clareza.
        - Painel 1: Massa vs. Rigidez Torcional (KT)
        - Painel 2: Massa vs. Rigidez Flexional (KF)

        Entradas:
        - history: dict com as listas 'best_mass', 'best_KT', 'best_KF'.
        - save_path: caminho para salvar o PNG.
        - show: bool, se True exibe o plot.
        """
        # Cria uma figura com 2 subplots (2 linhas, 1 coluna), compartilhando o eixo X
        fig, (ax1, ax3) = plt.subplots(
            nrows=2, 
            ncols=1, 
            figsize=(12, 10), 
            sharex=True
        )
        
        gens = range(len(history['best_mass']))
        
        # --- GRÁFICO 1: Massa vs. Rigidez Torcional (KT) ---
        
        ax1.set_title('Evolução: Massa vs. Rigidez Torcional (KT)', fontsize=14)
        
        # Eixo esquerdo (ax1) para a Massa
        color_mass = 'b'
        ax1.set_ylabel('Massa (kg)', color=color_mass, fontsize=12)
        line1 = ax1.plot(gens, history['best_mass'], color=color_mass, linestyle='-', linewidth=2, label='Massa (kg)')
        ax1.tick_params(axis='y', labelcolor=color_mass)
        ax1.grid(True, linestyle='--', alpha=0.6)
        
        # Eixo direito (ax2) para a Rigidez Torcional (KT)
        ax2 = ax1.twinx()
        color_kt = 'r'
        ax2.set_ylabel('Rigidez Torcional (N·m/rad)', color=color_kt, fontsize=12)
        line2 = ax2.plot(gens, history['best_KT'], color=color_kt, linestyle='-', linewidth=2, label='Rigidez Torcional (KT)')
        ax2.tick_params(axis='y', labelcolor=color_kt)
        
        # Junta as legendas do primeiro gráfico
        lines = line1 + line2
        labels = [l.get_label() for l in lines]
        ax1.legend(lines, labels, loc='best')

        # --- GRÁFICO 2: Massa vs. Rigidez Flexional (KF) ---
        
        ax3.set_title('Evolução: Massa vs. Rigidez Flexional (KF)', fontsize=14)
        
        # Eixo esquerdo (ax3) para a Massa
        ax3.set_ylabel('Massa (kg)', color=color_mass, fontsize=12)
        line3 = ax3.plot(gens, history['best_mass'], color=color_mass, linestyle='-', linewidth=2, label='Massa (kg)')
        ax3.tick_params(axis='y', labelcolor=color_mass)
        ax3.grid(True, linestyle='--', alpha=0.6)

        # Eixo direito (ax4) para a Rigidez Flexional (KF)
        ax4 = ax3.twinx()
        color_kf = 'g'
        ax4.set_ylabel('Rigidez Flexional (N/m)', color=color_kf, fontsize=12)
        line4 = ax4.plot(gens, history['best_KF'], color=color_kf, linestyle='-', linewidth=2, label='Rigidez Flexional (KF)')
        ax4.tick_params(axis='y', labelcolor=color_kf)
        
        # Junta as legendas do segundo gráfico
        lines_2 = line3 + line4
        labels_2 = [l.get_label() for l in lines_2]
        ax3.legend(lines_2, labels_2, loc='best')

        # Configurações finais da figura
        plt.xlabel('Geração', fontsize=12)
        fig.tight_layout(pad=3.0) # Adiciona espaçamento para evitar sobreposição de títulos
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
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
        - save_path: caminho para salvar o PNG.
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

def find_new_index(old_index, nodes):
    """
        Calcula o mapeamento de índices entre o modelo base (reduzido) e o modelo final (expandido/espelhado).

        Esta função determina a nova posição de um nó na lista global `nodes_full` baseada na lógica de expansão de simetria:
        1. **Nós Laterais:** São duplicados. O índice original $i$ se torna $2i$ (original) e $2i+1$ (espelho).
        2. **Nós Centrais:** São únicos e deslocados. Eles são posicionados após todos os pares de nós laterais terem sido alocados.

        **Nota Importante:** A lógica assume que o array `nodes` de entrada está ordenado de forma que todos os nós laterais
        aparecem antes dos nós centrais.

        Args:
            old_index (int): O índice do nó na lista de nós base (antes da expansão).
            nodes (list | np.ndarray): Matriz (N, 3) contendo as coordenadas dos nós base. Usada para identificar
                quais nós são centrais (X ≈ 0.0).

        Returns:
            int | tuple: 
                - Se o nó for **central**: Retorna um único `int` (`new_index`).
                - Se o nó for **lateral**: Retorna uma tupla `(new_index, mirrored_index)`.

        Raises:
            IndexError: Se `old_index` estiver fora dos limites do array `nodes`.
        """
    
    if not isinstance(nodes, np.ndarray):
        nodes = np.array(nodes)
        
    # is_central será um array de 19 booleanos
    is_central = np.isclose(nodes[:, 0], 0.0)
    
    first_central = -1 # Inicializa
    for i in range(len(is_central)):
        if is_central[i]:
            first_central = i
            break
            
    if first_central == -1:
         first_central = len(is_central)

    # Validação de índice
    if old_index >= len(is_central) or old_index < 0:
        raise IndexError(f"old_index {old_index} está fora dos limites para 'nodes' com tamanho {len(is_central)}")

    if is_central[old_index]:
        # O novo índice começa DEPOIS que todos os não-centrais foram duplicados
        # first_central = 17
        new_index = (first_central * 2) + (old_index - first_central)
        return new_index
    else:
        # Nós não-centrais (índices 0 a 16)
        new_index = old_index * 2
        mirrored_index = new_index + 1
        return new_index, mirrored_index
  
def penalidades_geometricas(nodes, elements, base_nodes):
    '''
    Define e aplica as penalidades geométricas nos indivíduos do processo de otimização

    Entradas:
    - nodes: matriz de localização dos nós
    - elements: matriz de conexão dos nós ou matriz dos elementos

    Saída: 
    - penalidade: valor total da soma das penalidades geométricas acumuladas
    '''

    penalidade = 0

    # Penalidade em relação ao FH e FHB
    def penalidade_fh_fhb(nodes, base_nodes):
        '''
        Segue a regra F.6.3.4 do regulamento da SAE-2025 que limita a distância máxima entre o ponto mais alto do Front Hoop e o ponto que o conecta ao Front Hoop Bacing
        Distância máxima: 0,05 m
        
        Entrada: 
        - nodes: matriz de localização dos nós

        Saída:
        - pen: valor desta penalidade acumulada, especificamente
        '''
        pen = 0
        idx = find_new_index(19, base_nodes)
        fronthoop_node = nodes[idx]                          #declara o nó do fronthoop com indice novo
        fhb_node = nodes[find_new_index(5, base_nodes)[0]]                              #declara o nó de um front hoop bracing com indice novo
        dist_fh_fhb = abs(fronthoop_node[2] - fhb_node[2])                         #declara a distância no eixo z entre esses dois nós 
        if dist_fh_fhb > 0.05 or fronthoop_node[2] < fhb_node[2]:                                                     #condição retirada do regulamento
            pen += ((dist_fh_fhb - 0.05) ** 2 ) * 1e2                              #aplicação de penalidade

        return pen

    # Penalidade em relação ao MH e o MHB
    def penalidade_mh_mhb(nodes, base_nodes):
        '''
        Segue a regra F.5.9.4 do regulamento da SAE-2025 que limita a distância máxima entre o ponto mais alto do Main Hoop e o ponto que o conecta ao(s) Main Hoop Bacing(s).
        Distância máxima: 0,16 m

        Entrada: 
        - nodes: matriz de localização dos nós

        Saída:
        - pen: valor desta penalidade acumulada, especificamente
        '''
        pen = 0
        mainhoop_node = nodes[find_new_index(18, base_nodes)]                                        #declara o nó do main hoop com indice novo
        mhb_node = nodes[find_new_index(14, base_nodes)[0]]                                          #declara o nó do main hoop bracing com indice novo não espelhado
        deltax_mh_mhb = mainhoop_node[0] - mhb_node[0]                                          #diferença das coordenadas "x" em ambos os nós
        deltay_mh_mhb = mainhoop_node[1] - mhb_node[1]                                          #diferença das coordenadas "y" em ambos os nós
        deltaz_mh_mhb = mainhoop_node[2] - mhb_node[2]                                          #diferença das coordenadas "z" em ambos os nós
        dist_mh_mhb = np.sqrt(deltax_mh_mhb ** 2 + deltay_mh_mhb ** 2 + deltaz_mh_mhb ** 2)     #declara a distância entre os dois nós      
        if dist_mh_mhb > 0.16 or dist_mh_mhb < 0:                                                                  #condição retirada do regulamento
            pen += ((dist_mh_mhb - 0.16) ** 2) * 1e2                                            #aplicação de penalidade
        
        return pen
    
    # Penalidade ângulo entre o Main Hoop e o Main Hoop Bracing
    def penalidade_angulo_mh_mhb(nodes, base_nodes):
        '''
        Segue a regra F.5.9.5 do regulamento da SAE-2025 que limita o ângulo mínimo entre o Main Hoop e o(s) Main Hoop Bracing(s)
        Ângulo mínimo: 30°

        Entrada: 
        - nodes: matriz de localização dos nós

        Saída:
        - pen: valor desta penalidade acumulada, especificamente
        '''
        pen = 0
        x_porcao_mh = nodes[find_new_index(14, base_nodes)[0]][0] - nodes[find_new_index(7, base_nodes)[0]][0]                                                                  #coordenada x do vetor formado pelos nós do elemento da porção do mainhoop analisada
        y_porcao_mh = nodes[find_new_index(14, base_nodes)[0]][1] - nodes[find_new_index(7, base_nodes)[0]][1]                                                                  #coordenada y do vetor formado pelos nós do elemento da porção do mainhoop analisada
        z_porcao_mh = nodes[find_new_index(14, base_nodes)[0]][2] - nodes[find_new_index(7, base_nodes)[0]][2]                                                                  #coordenada z do vetor formado pelos nós do elemento da porção do mainhoop analisada
        x_mhb = nodes[find_new_index(14, base_nodes)[0]][0] - nodes[find_new_index(15, base_nodes)[0]][0]                                                                       #coordenada x do vetor formado pelos nós do elemento do Main Hoop Bracing
        y_mhb = nodes[find_new_index(14, base_nodes)[0]][1] - nodes[find_new_index(15, base_nodes)[0]][1]                                                                       #coordenada y do vetor formado pelos nós do elemento do Main Hoop Bracing
        z_mhb = nodes[find_new_index(14, base_nodes)[0]][2] - nodes[find_new_index(15, base_nodes)[0]][2]                                                                       #coordenada z do vetor formado pelos nós do elemento do Main Hoop Bracing
        vetor_porcao_mh = (x_porcao_mh, y_porcao_mh, z_porcao_mh)                                                                                               #vetor formado pelos nós do elemento da porção do mainhoop analisada
        vetor_mhb = (x_mhb, y_mhb, z_mhb)                                                                                                                       #vetor formado pelos nós do elemento do Main Hoop Bracing
        modulo_vetor_porcao_mh = np.sqrt(vetor_porcao_mh[0] ** 2 + vetor_porcao_mh[1] ** 2 + vetor_porcao_mh[2] ** 2 )                                          #módulo do vetor formado pelos nós do elemento da porção do mainhoop analisada
        modulo_vetor_mhb = np.sqrt(vetor_mhb[0] ** 2 + vetor_mhb[1] ** 2 + vetor_mhb[2] ** 2 )                                                                  #módulo do vetor formado pelos nós do elemento do Main Hoop Bracing
        produto_escalar_mh_porcao_e_mhb = (vetor_porcao_mh[0] * vetor_mhb[0]) + (vetor_porcao_mh[1] * vetor_mhb[1]) + (vetor_porcao_mh[2] * vetor_mhb[2])       #produto escalar entre os dois vetores criados
        cos_theta_mh_mhb = produto_escalar_mh_porcao_e_mhb / (modulo_vetor_porcao_mh * modulo_vetor_mhb)                                                        #valor do cosseno do ângulo formado pelos dois vetores
        theta_mh_mhb = np.degrees(np.arccos(cos_theta_mh_mhb))                                                                                                    #valor do ângulo formado pelos dois vetores
        if theta_mh_mhb < 30:                                                                                                                                   #condição retirada do regulamento
            pen += ((theta_mh_mhb - 30) ** 2) * 1e2                                                                                                             #aplicação da penalidade
    
        return pen

    # Penalidade ângulo com a vertical da parte do Front Hoop que fica acima da Upper Side Impact Structure
    def penalidade_angulo_vetical_fh(nodes, base_nodes):
        '''
        Segue a regra F.5.7.6 do regulamento da SAE-2025 que limita o ângulo entre a parte do Front Hoop acima da Upper Side Impact Structure e o eixo da vertical 
        Ângulo máximo: 20°

        Entrada: 
        - nodes: matriz de localização dos nós

        Saída:
        - pen: valor desta penalidade acumulada, especificamente
        '''
        pen = 0
        x_porcao_fh = nodes[find_new_index(5, base_nodes)[0]][0] - nodes[find_new_index(4, base_nodes)[0]][0]                                          #coordenada x do vetor formado pelos nós do elemento da porção do mainhoop analisada
        y_porcao_fh = nodes[find_new_index(5, base_nodes)[0]][1] - nodes[find_new_index(4, base_nodes)[0]][1]                                          #coordenada y do vetor formado pelos nós do elemento da porção do mainhoop analisada
        z_porcao_fh = nodes[find_new_index(5, base_nodes)[0]][2] - nodes[find_new_index(4, base_nodes)[0]][2]                                          #coordenada z do vetor formado pelos nós do elemento da porção do mainhoop analisada
        vetor_porcao_fh = (x_porcao_fh, y_porcao_fh, z_porcao_fh)                                                                      #vetor formado pelos nós que formam a porção do front hoop analisada
        modulo_vetor_porcao_fh = np.sqrt(vetor_porcao_fh[0] ** 2 + vetor_porcao_fh[1] ** 2 + vetor_porcao_fh[2] **2)                   #módulo do vetor formado pelos nós que formam a porção do front hoop analisada
        produto_escalar_porcao_fh_e_vertical = vetor_porcao_fh[2]                                                                      #produto escalar entre o vetor formado pelos nós que formam a porção do front hoop analisada e o versor da vertical
        cos_theta_fh_porcao_vertical = produto_escalar_porcao_fh_e_vertical / modulo_vetor_porcao_fh                                   #cosseno ângulo formado entre a porção do front hoop analisada e a vertical
        theta_fh_porcao_vertical = np.degrees(np.arccos(cos_theta_fh_porcao_vertical))                                                   #cálculo do ângulo através do cosseno
        if theta_fh_porcao_vertical > 20:                                                                                              #condição retirada do regulamento
            pen += ((theta_fh_porcao_vertical - 20) ** 2) * 1e2                                                                       #aplicação da penalidade
    
        return pen

    # Penalidade ângulo com a vertical da parte do Main Hoop que fica acima do ponto que o conecta ao Upper Side Impact Tube 
    def penalidade_angulo_vertical_mh(nodes, base_nodes):
        '''
        Segue a regra F.5.8.3 do regulamento da SAE-2025 que limita o ângulo entre a parte do Main Hoop que fica acima do ponto que o conecta com o Upper Side Impact Tube com o eixo da vertical
        Ângulo máximo: 10° 

        Entrada: 
        - nodes: matriz de localização dos nós

        Saída:
        - pen: valor desta penalidade acumulada, especificamente
        '''
        pen = 0
        x_porcao_mh = nodes[find_new_index(14, base_nodes)[0]][0] - nodes[find_new_index(6, base_nodes)[0]][0]                                          #coordenada x do vetor formado pelos nós do elemento da porção do mainhoop analisada
        y_porcao_mh = nodes[find_new_index(14, base_nodes)[0]][1] - nodes[find_new_index(6, base_nodes)[0]][1]                                          #coordenada y do vetor formado pelos nós do elemento da porção do mainhoop analisada
        z_porcao_mh = nodes[find_new_index(14, base_nodes)[0]][2] - nodes[find_new_index(6, base_nodes)[0]][2]                                          #coordenada z do vetor formado pelos nós do elemento da porção do mainhoop analisada
        vetor_porcao_mh = (x_porcao_mh, y_porcao_mh, z_porcao_mh)                                                                       #vetor formado pelos nós do elemento da porção do mainhoop analisada
        modulo_vetor_porcao_mh = np.sqrt(vetor_porcao_mh[0] ** 2 + vetor_porcao_mh[1] ** 2 + vetor_porcao_mh[2] ** 2 )                  #módulo do vetor formado pelos nós do elemento da porção do mainhoop analisada
        produto_escalar_porcao_mh_e_vertical = vetor_porcao_mh[2]                                                                       #produto escalar do vetor formado pelo elemento da porção do Main Hoop trabalhada com o versor da vertical
        cos_theta_mh_porcao_vertical = produto_escalar_porcao_mh_e_vertical / modulo_vetor_porcao_mh                                    #cosseno do ângulo formado entre este vetor mencionado e a vertical
        theta_mh_porcao_vertical = np.degrees(np.arccos(cos_theta_mh_porcao_vertical))                                                    #ângulo formado entre este vetor mencionado e a vertical
        if theta_mh_porcao_vertical > 10:                                                                                               #condição retirada do regulamento
            pen += ((theta_mh_porcao_vertical -10) ** 2) * 1e2                                                                          #aplicação da penalidade

        return pen

    # Penalidade para posição do Upper Side Impact Member
    def penalidade_posicao_upper_SI_member(nodes, base_nodes):
        '''
        Segue a regra F.6.4.4 do regulamento da SAE-2025 que define uma zona permitida para a posição do Upper Side Impact Member
        Zona: 265 mm a 320 mm da parte frontal do Lower Side Impact Member que é ligado ao Front Hoop

        Entrada: 
        - nodes: matriz de localização dos nós

        Saída:
        - pen: valor desta penalidade acumulada, especificamente
        '''
        pen = 0
        upper_SI_node_front_z = nodes[find_new_index(4, base_nodes)[0]][2]
        upper_SI_node_back_z = nodes[find_new_index(6, base_nodes)[0]][2]
        lowest_point_SI_z =  nodes[find_new_index(7, base_nodes)[0]][2]
        dist_lowestpoint_front_usim = upper_SI_node_front_z - lowest_point_SI_z
        dist_lowestpoint_back_usim = upper_SI_node_back_z - lowest_point_SI_z

        if dist_lowestpoint_front_usim * 1000 not in range(265, 321) and dist_lowestpoint_back_usim * 1000 not in range(265, 321):
            pen += ((dist_lowestpoint_front_usim) ** 2 + (dist_lowestpoint_back_usim) ** 2) * 1e2


        return pen
    
    # Penalidade de retidão vertical para um elemento definido entre dois nós
    def penalidade_retidão_vertical (nodes, base_nodes, no1, no2):
        '''
        Esta penalidade tem como função garantir a retidão de qualquer elemento com a vertical. Ou seja, dados os nós das extremidades do elemento, caso o ângulo formado por ele e o eixo da vertical 
        seja diferente de 0°, uma penalidade alta será gerada

        Entrada: 
        - nodes: matriz de localização dos nós

        Saída:
        - pen: valor desta penalidade acumulada, especificamente
        '''
        pen = 0 
        componente_x_vetor = nodes[find_new_index(no1, base_nodes)[0]][0] - nodes[find_new_index(no2, base_nodes)[0]][0]                          
        componente_y_vetor = nodes[find_new_index(no1, base_nodes)[0]][1] - nodes[find_new_index(no2, base_nodes)[0]][1]   
        componente_z_vetor = nodes[find_new_index(no1, base_nodes)[0]][2] - nodes[find_new_index(no2, base_nodes)[0]][2]
        vetor_no1_no2 = (componente_x_vetor, componente_y_vetor, componente_z_vetor)
        modulo_vetor_no1_no2 = np.linalg.norm(vetor_no1_no2)
        vertical = np.array([0, 0, 1])                                                                                             #vetor vertical unitario
        produto_escalar_no1_no2_vertical = np.dot(vetor_no1_no2, vertical)
        cos_theta_no1_no2_vertical = produto_escalar_no1_no2_vertical / modulo_vetor_no1_no2
        theta_no1_no2 = np.degrees(np.arccos(cos_theta_no1_no2_vertical))
        if theta_no1_no2 > 1e-6 :
            pen += ((theta_no1_no2) ** 2) * 1e2
        
        return pen
    
    #def penalidade_gabarito_F_z(nodes):
    #    '''
    #    Esta penalidade tem como garantir que o gabarito com dimensões informadas pelo regulamento da SAE caiba dentro do chassi
#
    #    Entrada: 
    #    - nodes: matriz de localização dos nós
#
    #    Saída:
    #    - pen: valor desta penalidade acumulada, especificamente
    #    '''
    #    pen = 0
    #    # Lista de pares de nós e seus respectivos limites mínimos
    #    requisitos = [
    #    ((1, 18), 0.20),
    #    ((7, 21), 0.275),
    #    ((8, 22), 0.175),
    #    ]
#
    #    for (n1, n2), largura_min in requisitos:
    #        largura = abs(nodes[find_new_index(n1, nodes)[0]][0] -
    #                  nodes[find_new_index(n2, nodes)[0]][0])
    #    
    #        if largura < largura_min:
    #            pen += ((largura - largura_min) ** 2) * 1e6
#
    #    return pen

    # Controle das penalidades
    penalidade += penalidade_fh_fhb(nodes, base_nodes)
    penalidade += penalidade_mh_mhb(nodes, base_nodes)
    penalidade += penalidade_angulo_vetical_fh(nodes, base_nodes)
    penalidade += penalidade_angulo_vertical_mh(nodes, base_nodes)
    penalidade += penalidade_angulo_mh_mhb(nodes, base_nodes)
    penalidade += penalidade_posicao_upper_SI_member(nodes, base_nodes)
    penalidade += penalidade_retidão_vertical(nodes, base_nodes, 0, 1)
    penalidade += penalidade_retidão_vertical(nodes, base_nodes, 12, 13)
    #penalidade += penalidade_gabarito_F_z(nodes)

    return penalidade

def penalidades_tipo_tubo(nodes, elements, mapeamento_chassi, tubos_SAE):
    """
    Verifica a conformidade normativa dos tubos do chassi em relação ao regulamento (ex: Formula SAE).

    Esta função cruza a geometria otimizada com um mapa de requisitos regionais. Para cada elemento da estrutura,
    ela identifica a qual subsistema ele pertence (ex: Santo Antônio Principal, Impacto Lateral) e verifica se
    o perfil de tubo escolhido pelo otimizador atende aos requisitos mínimos de geometria e inércia exigidos
    para aquela zona específica.

    O processo de validação segue a lógica:
    1. Identificação Regional: Para cada elemento (conexão entre nó A e B), busca sua classificação no `mapeamento_chassi`.    
    2. Recuperação de Regra: Identifica o "Tubo Padrão SAE" exigido para aquela região.
    3. Cálculo de Propriedades: Calcula Área (A), Diâmetro Externo (d), Espessura (e) e Momento de Inércia (I)
       tanto para o tubo proposto (otimizado) quanto para o tubo de referência (SAE).
    4. Comparação: Aplica uma penalidade fixa ("Hard Penalty") se o tubo proposto for inferior ao exigido em
       qualquer uma das quatro métricas físicas.

    Args:
        nodes (np.ndarray): Matriz de coordenadas dos nós.
        elements (list): Lista de elementos contendo a conectividade e o identificador do perfil (ex: 'Tubo B').
        mapeamento_chassi (dict): Dicionário que mapeia regiões do chassi (chaves) para listas de conexões de nós (valores)
            e o tipo de tubo mandatório (último item da lista).
        tubos_SAE (dict): Banco de dados com as propriedades físicas mínimas dos tubos normatizados.

    Returns:
        int: Valor acumulado da penalidade. Cada violação de regra adiciona 10 pontos ao custo total.
    """
    
    if isinstance(mapeamento_chassi, str):
        print("CRITICAL ERROR: 'mapeamento_chassi' is a string!")
        print(f"Value content: {mapeamento_chassi}")
    penalidade = 0

    estrutura = Estrutura(elements, nodes)
    type_tube_sae = ''
    type_tube_otm = ''
    for element in elements:
        elemento = (element[0], element[1])          # nós do elemento
        type_tube_otm = element[2]                   # tipo de tubo otimizado (ex: 'Tubo B')

        for classification in mapeamento_chassi:
            if elemento in mapeamento_chassi[classification] or tuple(reversed(elemento)) in mapeamento_chassi[classification]:
                # tipo mínimo exigido pela SAE
                type_tube_sae = mapeamento_chassi[classification][-1]
                d_sae, e_sae, A_sae, I_sae = tubos_SAE[type_tube_sae]

                # propriedades do tubo otimizado
                props_otm = estrutura.obter_propriedades(type_tube_otm)
                d_otm, e_otm = props_otm[2], props_otm[3]
                A_otm = estrutura.area_seccao_transversal(d_otm, e_otm, props_otm[6])
                I_otm = estrutura.momento_inercia_area_e_polar(d_otm, e_otm, props_otm[6])[0]

                # comparação com limites mínimos da SAE
                if A_otm < A_sae or d_otm < d_sae or e_otm < e_sae or I_otm < I_sae:
                    penalidade += 10
                
    return penalidade

def evaluate(nodes, elements, base_nodes, tubos_SAE) -> float:
    """
        Executa a simulação física completa (FEA) de um indivíduo e calcula sua aptidão (fitness).

        Esta função atua como uma interface entre o algoritmo de otimização e o solver de Elementos Finitos (`Estrutura`).
        Ela submete a geometria proposta a casos de carga estáticos e calcula métricas de desempenho.

        O fluxo de avaliação consiste em:
        1. Instanciação do Modelo Cria o objeto `Estrutura` com a topologia fornecida.
        2. Montagem de Matrizes Calcula as matrizes globais de Rigidez (K) e Massa (M).
        3. Definição de Cargas (Caso Estático)
        - Identifica os nós de suporte de carga (ex: fixação do banco/piloto) usando `find_new_index` 
        - Aplica uma carga vertical distribuída 
        4. Solução do Sistema Resolve $K \cdot u = F$ para obter os deslocamentos.
        5. Pós-Processamento
        - Calcula tensões (Von Mises).
        - Calcula a massa total.
        - Calcula Rigidez à Torção (K_T) e Flexão (K_F) baseada nos deslocamentos dos nós de suspensão .
        
        6. Cálculo de Penalidade: Soma as penalidades de desempenho (massa excessiva, baixa rigidez) e 
        violações de restrições (diâmetro de tubos inválidos, falha de tensão).

        Args:
            nodes (np.ndarray): Matriz de coordenadas (N, 3) do indivíduo decodificado.
            elements (list[tuple]): Lista de conectividade e perfis dos tubos.
            base_nodes (np.ndarray): Geometria de referência para mapeamento de índices.
            tubos_SAE (dict/list): Banco de dados de propriedades materiais e seções transversais.

        Returns:
            tuple: Uma tupla contendo:
                - penalty (float): O valor final da função objetivo (fitness) a ser minimizado.
                - massa (float): Massa total da estrutura (kg).
                - KT (float): Rigidez Torcional (Nm/deg).
                - KF (float): Rigidez à Flexão (N/mm).

        Raises:
            Exception: Captura qualquer erro numérico do solver (ex: matriz singular) e retorna
            valores de penalidade extremos (`1e9`) para descartar o indivíduo falho.
        """
    try:
        t0 = time.perf_counter()

        penalty = 0

        # Instanciamento da estrutura
        estrutura = Estrutura(elements, nodes)
        t1 = time.perf_counter()

        k_global,m_global = estrutura.matrizes_global()
        t2 = time.perf_counter()

        #_, _, frequencies = estrutura.modal_analysis()
        frequencies=[0,0,0]
        t3 = time.perf_counter()

        fixed = list(range(6))
        Fg = np.zeros(estrutura.num_dofs)
        pack_sup1,pack_sup2 = find_new_index(7,nodes)
        pack_sup3,pack_sup4 = find_new_index(10,nodes)

        for i in [pack_sup1,pack_sup2,pack_sup3,pack_sup4]:
            Fg[i * 6 + 2] = (65 * 9.81) / 4
            
        disp = estrutura.static_analysis(Fg, fixed)
        t4 = time.perf_counter()

        strain = estrutura.compute_strain(disp)
        esforcos = estrutura.compute_stress(strain)
        stresses = estrutura.calcular_tensoes_reais(esforcos)
        von = estrutura.compute_von_mises(stresses)
        t5 = time.perf_counter()

        massa = estrutura.mass()
        LFnode,RFnode = find_new_index(3,nodes)                                     #Left front suspension node index (3)
        LRnode,RRnode = find_new_index(11,nodes)                                    #Left rear suspension node index (11)
        LCnode,RCnode = find_new_index(7,nodes)                                     #Left central node index (7)
        KT, KF = estrutura.compute_Kf_Kt(LFnode,LRnode,RFnode,RRnode,LCnode,RCnode)
        t6 = time.perf_counter()

        penalty += penalidade_chassi(KT, KF, massa, von, frequencies)
        #penalty += penalidades_geometricas(nodes, elements, base_nodes)
        penalty += penalidades_tipo_tubo(nodes, elements, mapeamento_chassi, tubos_SAE)

        t7 = time.perf_counter()

        #print(f"""[EVALUATE]
        #- Instância:     {t1 - t0:.4f}s
        #- Montagem K/M:  {t2 - t1:.4f}s
        #- Modal:         {t3 - t2:.4f}s
        #- Estática:      {t4 - t3:.4f}s
        #- Tensão/Von:    {t5 - t4:.4f}s
        #- Massa/Rigidez: {t6 - t5:.4f}s
        #- Penalidade:    {t7 - t6:.4f}s
        #- Total:         {t7 - t0:.4f}s
        #""",end="\r")

        return (penalty), massa, KT, KF

    except Exception:
        traceback.print_exc()
        return 1e9, 1e9, 0.0, 0.0

def penalidade_chassi(KT, KF, massa, tensoes, frequencias,W_RIGIDEZ = 1.5,W_MASSA = 5.5,W_TENSAO = 5.0,W_FREQUENCIA = 10.0):
    """
        Calcula a função de penalidade composta (custo) baseada no desempenho físico da estrutura.

        Esta função implementa o método de "Soft Constraints" (Restrições Suaves), onde violações dos limites de projeto
        são convertidas em valores positivos adicionados à aptidão (fitness). Quanto maior a violação, maior a penalidade,
        pressionando o otimizador a corrigir o design.

        A função agrega quatro tipos de penalidades:
        1. Rigidez (KT/KF):
        - KT (Torcional): Penalidade Exponencial. Alta sensibilidade para garantir rigidez torcional mínima.
        - KF (Flexional): Penalidade Logarítmica. Crescimento mais suave, permitindo exploração próxima ao limite.
        
        

        2. Massa: Define uma "janela de alvo" (22kg a 25kg). Penaliza exponencialmente tanto chasssis
        muito leves (possivelmente frágeis/ilegais) quanto muito pesados.

        3. Tensão (Von Mises): Penaliza se a tensão máxima em qualquer elemento exceder o limite de escoamento/admissível
        (250 MPa).

        4. Vibração (Ressonância): Penaliza frequências naturais que caiam dentro de uma banda de 10% ao redor
        da frequência de excitação do motor (aprox. 80Hz). A penalidade cresce quadraticamente conforme a frequência
        se aproxima do pico de 80Hz.

        Args:
            KT (float): Rigidez Torcional calculada (Nm/deg).
            KF (float): Rigidez Flexional calculada (N/mm ou N/m, dependendo da escala).
            massa (float): Massa total da estrutura (kg).
            tensoes (np.ndarray): Array contendo as tensões máximas de Von Mises em cada elemento.
            frequencias (list | np.ndarray): As primeiras frequências naturais (modos de vibrar).
            W_... (float): Pesos (Weights) para calibrar a importância relativa de cada penalidade.

        Returns:
            float: O valor acumulado de penalidade (`fitness`). Se for 0.0, o chassi cumpre todos os requisitos físicos.
        """
        
    # Limites
    KT_min = 2000.0   # Rigidez torcional mínima (N·m/deg)
    KF_min = 1.0e6    # Rigidez flexional mínima (N/m)
    massa_min, massa_max = 22.0, 25.0  # faixa aceitável (kg)
    tensao_adm = 250.0e6 # Tensão admissível do material (Pa)
    
    # Fatores de penalidade
    alpha_exp = 4.0  # Fator de sensibilidade exponencial, ajustado para ser mais agressivo
    beta_log = 5.0   # Fator de escala logarítmica
      
    # Parâmetros de Frequência
    freq_motor_hz = 4800 / 60
    
    fitness = 0.0

    # 2. Penalidade por Rigidez (KT e KF)
    if KT < KT_min:
        deficit_kt = (KT_min - KT) / KT_min
        fitness += W_RIGIDEZ * (np.exp(alpha_exp * deficit_kt) - 1)

    if KF < KF_min:
        deficit_kf = (KF_min - KF) / KF_min
        fitness += W_RIGIDEZ * beta_log * np.log(1 + deficit_kf)

    # --- Massa ---
    if massa < massa_min:
        deficit = (massa_min - massa) / massa_min
        fitness += W_MASSA * (np.exp(alpha_exp * deficit) - 1)
    elif massa > massa_max:
        excesso = (massa - massa_max) / massa_max
        fitness += W_MASSA * (np.exp(alpha_exp * excesso) - 1)

    # 4. Penalidade por Tensões (Von Mises)
    if len(tensoes) > 0:
        tensao_max = np.max(tensoes)
        if tensao_max > tensao_adm:
            excesso_tensao = (tensao_max - tensao_adm) / tensao_adm
            fitness += W_TENSAO * (np.exp(alpha_exp * excesso_tensao) - 1)

    # 6. Penalidade por Frequências Naturais
    for f in frequencias:
        # A banda crítica é de 5% acima e abaixo da frequência do motor
        if abs(f - freq_motor_hz) / freq_motor_hz < 0.10:
            # Penalidade aumenta quadraticamente à medida que se aproxima do pico de ressonância
            distancia_relativa = abs(f - freq_motor_hz) / freq_motor_hz
            fitness += W_FREQUENCIA * (1 - distancia_relativa)**2

    return fitness

if __name__ == "__main__":
    nodes = np.array([
        [-0.160, -0.10,   0.350],   # 00 (Antigo 00)
        [-0.160, -0.10,   0.000],   # 01 (Antigo 01)
        [-0.245,  0.275,  0.220],   # 02 (Antigo 02)
        [-0.245,  0.605,  0.000],   # 03 (Antigo 03)
        [-0.245,  0.255,  0.000],   # 04 (Antigo 04)
        [-0.245,  0.585,  0.220],   # 05 (Antigo 05)
        [-0.210,  0.555,  0.500],   # 06 (Antigo 06)
        [-0.290,  1.370,  0.250],   # 07 (Antigo 07)
        [-0.268,  1.350,  0.000],   # 08 (Antigo 08)
        [-0.250,  1.670,  0.240],   # 09 (Antigo 09)
        [-0.250,  1.665,  0.030],   # 10 (Antigo 10)
        [-0.200,  1.950,  0.030],   # 11 (Antigo 11)
        [-0.200,  2.230,  0.250],   # 12 (Antigo 12)
        [-0.200,  2.230,  0.030],   # 13 (Antigo 13)
        [-0.170,  1.400,  0.965],   # 14 (Antigo 14)
        [-0.200,  1.950,  0.250],   # 15 (Antigo 15)
        [-0.200,  0.420,  0.469],   # 16 (Antigo 18)
        [-0.200,  2.090,  0.250],   # 17
        [ 0.000,  1.410,  1.105],   # 18 (Antigo 16) 
        [ 0.000,  0.555,  0.550],   # 19 (Antigo 17) 
    ])

    # ============================================================
    # CONEXÕES ORIGINAIS
    # ============================================================
    connections = [
    (0, 1), (0, 2), (0, 16),(6,16), (1, 2), (1, 4), (2, 3), (2, 4), (2, 5), (2, 6),
    (3, 4), (3, 5), (3, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 19), (7, 8),
    (7, 9), (8, 9), (7, 14), (8, 10), (9, 10), (9, 15), (10, 11), (9, 11),
    (11, 15), (11, 13), (12, 13), (12, 17),(15,17), (11, 12), (14, 18), (14, 15)]

    indices = [0,1,6,7,8,14]

    indices_suspensao = [2, 3, 4, 5, 11, 12, 13, 15 ,16 , 17]

    tubos_fixos = {
#    0: 'Tubo A',
#    1: 'Tubo A',
#    2: 'Tubo A',
#    'espelho_0': 'Tubo C',
#    'espelho_1': 'Tubo C',
#    'espelho_3': 'Tubo C',
#    'espelho_5': 'Tubo C',
}
    
    connect_with_mirror = [0,1,3,8,7,11,12,13,16,17]

    #Criando mapeamento do chassi completo base para análise (uma vez para ser aplicado em todos os indivíduos)
    mapeamento_chassi = {
        #Classificação: [elemento1, ..., 'Tubo X]
        "Front_Bulkhead": [(find_new_index(0, nodes)[0], find_new_index(1, nodes)[0]), (find_new_index(0, nodes)[0], find_new_index(0, nodes)[1]), (find_new_index(1, nodes)[0], find_new_index(1, nodes)[1]), 'Tubo B'],   
        "Front_Bulkhead_Support": [(find_new_index(0, nodes)[0], find_new_index(2, nodes)[0]), (find_new_index(1, nodes)[0], find_new_index(2, nodes)[0]), (find_new_index(1, nodes)[0], find_new_index(16, nodes)[0]), (find_new_index(2, nodes)[0], find_new_index(16, nodes)[0]), (find_new_index(3, nodes)[0], find_new_index(16, nodes)[0]), (find_new_index(2, nodes)[0], find_new_index(3, nodes)[0]), (find_new_index(2, nodes)[0], find_new_index(4, nodes)[0]), (find_new_index(2, nodes)[0], find_new_index(5, nodes)[0]), 'Tubo C'],    
        "Front_Hoop": [(find_new_index(4, nodes)[0], find_new_index(5, nodes)[0]), (find_new_index(3, nodes)[0], find_new_index(4, nodes)[0]), (find_new_index(5, nodes)[0], find_new_index(18, nodes)), 'Tubo A'],                 
        "Front_Hoop_Bracing": [(find_new_index(0, nodes)[0], find_new_index(5, nodes)[0]), 'Tubo B'],                   
        "Side_Impact_Structure": [(find_new_index(6, nodes)[0], find_new_index(5, nodes)[0]), (find_new_index(4, nodes)[0], find_new_index(6, nodes)[0]), (find_new_index(3, nodes)[0], find_new_index(7, nodes)[0]), (find_new_index(4, nodes)[0], find_new_index(7, nodes)[0]), 'Tubo B'],                   
        "Bent/Multi_Upper_Side_Impact_Member": ['Tubo D'],                   
        "Main_Hoop": [(find_new_index(6, nodes)[0], find_new_index(14, nodes)[0]), (find_new_index(14, nodes)[0], find_new_index(17, nodes)), 'Tubo A'],                  
        "Main_Hoop_Bracing": [(find_new_index(14, nodes)[0], find_new_index(15, nodes)[0]), 'Tubo B'],                  
        "Main_Hoop_Bracing_Supports": ['Tubo C'],
        "Driver_Restraint_Harness_Attachment": ['Tubo B'],
        "Shoulder_Harness_Mounting_Bar": ['Tubo A'],
        "Shoulder_Harness_Mounting_Bar_Bracing": ['Tubo C'],
        "Accumulator_Mounting_and_Protection": [(find_new_index(7, nodes)[0], find_new_index(6, nodes)[0]), (find_new_index(6, nodes)[0], find_new_index(8, nodes)[0]), (find_new_index(7, nodes)[0], find_new_index(8, nodes)[0]), (find_new_index(8, nodes)[0], find_new_index(10, nodes)[0]), (find_new_index(8, nodes)[0], find_new_index(9, nodes)[0]), (find_new_index(6, nodes)[0], find_new_index(9, nodes)[0]), (find_new_index(9, nodes)[0], find_new_index(10, nodes)[0]), 'Tubo B'],
        "Component_Protection": ['Tubo C'],
        "Structural_Tubing": [(find_new_index(10, nodes)[0], find_new_index(11, nodes)[0]), (find_new_index(9, nodes)[0], find_new_index(11, nodes)[0]), (find_new_index(9, nodes)[0], find_new_index(15, nodes)[0]), (find_new_index(11, nodes)[0], find_new_index(15, nodes)[0]), (find_new_index(15, nodes)[0], find_new_index(12, nodes)[0]), (find_new_index(12, nodes)[0], find_new_index(13, nodes)[0]), (find_new_index(11, nodes)[0], find_new_index(13, nodes)[0]), (find_new_index(13, nodes)[0], find_new_index(15, nodes)[0]), 'Tubo C']
    }

    #Dicionário de tubos com parâmetros adequados em relação às regras da SAE para serem aplicados nas análides de cada indivíduo
    tubos_SAE = {
        #Tubo X : [Diâmetro, Espessura, Área da secção transversal, Momento de Inércia de Área]
        'Tubo A': [0.025, 0.002, 0.000173, 1.13E-08],   
        'Tubo B': [0.025, 0.0012, 0.000114, 8.51E-09],
        'Tubo C': [0.025, 0.0012, 0.000091, 6.70E-09],
        'Tubo D': [0.035,0.0012, 0.000126, 1.80E-08]
    }

    # Criar diretório para resultados
    timestamp = datetime.now().strftime("%Y-%m-%d %H%M")
    max_gen = 300
    pop_size = 50
    otimizador = ChassisDEOptimizer(
        base_nodes=nodes,
        base_connections=connections,
        mandatory_indices=indices,
        fixed_nodes=indices_suspensao,
        pop_size=pop_size,
        max_generations=max_gen,
        use_parallel=False,                         # ou False, se quiser forçar modo serial
        n_workers=3,                                # defina quantos núcleos usar
        tubos_fixos=tubos_fixos,
        connect_with_mirror=connect_with_mirror
        )
    
    results_dir = f"Resultados_Otimizacao_LOCAL__{timestamp}_{max_gen}GEN_{pop_size}POP"
    os.makedirs(results_dir, exist_ok=True)

    best_indiv, best_cost, best_mass, best_KT, best_KF, history, duration = otimizador.optimize()
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