import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
from scipy.integrate import cumulative_trapezoid as cumtrapz

# Verifica se a biblioteca imageio está instalada
try:
    import imageio
except ImportError:
    print("A biblioteca imageio não está instalada. Instalando...")
    import subprocess
    subprocess.check_call(["pip", "install", "imageio"])
    import imageio

# Função para simular o movimento do carro
def simulate_car(F_pedal, pos_acel_max, pos_freio, pos_desejada):
    # Penalização se pos_acel_max > pos_freio ou pos_freio > pos_desejada
    if pos_acel_max >= pos_freio or pos_freio >= pos_desejada:
        # Retorna uma diferença grande para penalizar a solução inválida
        return float('inf'), pos_acel_max, pos_freio, 0, None, None, None, None

    # Parâmetros Físicos Fixos
    V_limite_carro = 26.7  # Velocidade máxima mecânica do carro (aprox 100 km/h)
    a_motor = 4.0          # Aceleração confortável/realista do motor (m/s²) - VALOR FIXO
    # Restante da função (cálculo da dinâmica do carro)
    V1 = 26.7  # Velocidade máxima (m/s)
    tr = 1.0   # Tempo de reação (s)
    ta = 0.25  # Tempo de aplicação do freio (s)
    tb = 0.3   # Tempo de transição da desaceleração (s)
    m = 230  # Massa do carro (kg)
    # 1. Calcula a velocidade atingida baseada na distância de aceleração disponível
    # Torricelli: V^2 = Vo^2 + 2*a*d  ->  V = sqrt(2 * a * d)
    V_atingida = np.sqrt(2 * a_motor * pos_acel_max)
    
    # 2. A velocidade V1 agora é dinâmica. Se a pista for curta, V1 será menor que 100km/h.
    # Se a pista for longa, limitamos a V1 ao máximo do carro (26.7 m/s)
    V1 = min(V_atingida, V_limite_carro)
    
    # 3. Recalculamos a aceleração efetiva e tempos
    # Se o carro foi limitado pela velocidade máxima, a aceleração para de atuar antes,
    # mas para simplificar a otimização, assumimos aceleração constante 'a_motor' até atingir V1.
    a_acel = a_motor 
    amax = F_pedal / m  # Desaceleração máxima (m/s²)
    V2 = V1 - (0.5 * amax * tb)

    # Tempo para atingir a velocidade máxima (aceleração constante)
    t_acel = V1 / a_acel  # Tempo para atingir V1

    # Tempo para percorrer de pos_acel_max até pos_freio (velocidade constante)
    dist_const = pos_freio - pos_acel_max  # Distância percorrida com velocidade constante
    t_const = dist_const / V1  # Tempo com velocidade constante

    # Tempo de desaceleração até parar
    ts = V2 / amax
    t_desacel = tr + ta + tb + ts  # Tempo para parar após começar a frear

    # Tempo total da simulação
    t_total = t_acel + t_const + t_desacel

    # Criando os vetores de tempo
    t_values = np.linspace(0, t_total, 1000)

    # Calculando a aceleração em cada instante
    acceleration = np.zeros_like(t_values)
    for i, t in enumerate(t_values):
        if t < t_acel:
            # Fase de aceleração (0 a pos_acel_max)
            acceleration[i] = a_acel
        elif t < t_acel + t_const:
            # Fase de velocidade constante (pos_acel_max a pos_freio)
            acceleration[i] = 0
        elif t < t_acel + t_const + tr + ta:
            # Fase de transição da desaceleração
            acceleration[i] = 0
        elif t < t_acel + t_const + tr + ta + tb:
            tb_elapsed = t - (t_acel + t_const + tr + ta)
            acceleration[i] = -amax * (tb_elapsed / tb)
        else:
            # Fase de desaceleração máxima
            acceleration[i] = -amax  # Desaceleração máxima constante

    # Calculando a velocidade pela integral da aceleração
    velocity = cumtrapz(acceleration, t_values, initial=0)

    # Calculando a posição pela integral da velocidade
    position = cumtrapz(velocity, t_values, initial=0)

    # Determinando a posição final onde o carro para
    pos_final = position[-1]

    # Retorna a diferença entre a posição final e a posição desejada
    return abs(pos_final - pos_desejada), pos_acel_max, pos_freio, pos_final, t_values, acceleration, velocity, position, V1

# Algoritmo de Evolução Diferencial Autoadaptativa
def differential_evolution(pos_desejada, F_min=445, F_max=823, pos_min=0, pos_max=1000, pop_size=4000, max_iter=6400, CR=0.85, F_scale=0.5, tol=10e-5):
    # Inicialização da população
    population_F = np.random.uniform(F_min, F_max, pop_size)
    population_pos_acel = np.random.uniform(pos_min, pos_desejada, pop_size)  # pos_acel_max < pos_desejada
    population_pos_freio = np.random.uniform(population_pos_acel, pos_desejada)  # pos_freio > pos_acel_max e < pos_desejada

    best_fitness = float('inf')
    best_F = None
    best_pos_acel = None
    best_pos_freio = None
    best_pos_final = None

    # Listas para armazenar os dados das iterações
    fitness_history = []
    F_history = []
    pos_acel_history = []
    pos_freio_history = []

    for iteration in range(max_iter):
        for i in range(pop_size):
            # Seleção de três indivíduos aleatórios
            idxs = [idx for idx in range(pop_size) if idx != i]
            a_F, b_F, c_F = population_F[np.random.choice(idxs, 3, replace=False)]
            a_pos_acel, b_pos_acel, c_pos_acel = population_pos_acel[np.random.choice(idxs, 3, replace=False)]
            a_pos_freio, b_pos_freio, c_pos_freio = population_pos_freio[np.random.choice(idxs, 3, replace=False)]

            # Mutação
            mutant_F = a_F + F_scale * (b_F - c_F)
            mutant_pos_acel = a_pos_acel + F_scale * (b_pos_acel - c_pos_acel)
            mutant_pos_freio = a_pos_freio + F_scale * (b_pos_freio - c_pos_freio)

            # Recombinação
            trial_F = np.where(np.random.rand() < CR, mutant_F, population_F[i])
            trial_pos_acel = np.where(np.random.rand() < CR, mutant_pos_acel, population_pos_acel[i])
            trial_pos_freio = np.where(np.random.rand() < CR, mutant_pos_freio, population_pos_freio[i])

            # Garantir que os valores estejam dentro dos limites
            trial_F = np.clip(trial_F, F_min, F_max)
            trial_pos_acel = np.clip(trial_pos_acel, pos_min, pos_desejada)
            trial_pos_freio = np.clip(trial_pos_freio, trial_pos_acel, pos_desejada)  # pos_freio > pos_acel_max

            # Seleção
            fitness_trial, pos_acel_trial, pos_freio_trial, pos_final_trial, _, _, _, _, _= simulate_car(trial_F, trial_pos_acel, trial_pos_freio, pos_desejada)
            fitness_current, pos_acel_current, pos_freio_current, pos_final_current, _, _, _, _, _ = simulate_car(population_F[i], population_pos_acel[i], population_pos_freio[i], pos_desejada)

            if fitness_trial < fitness_current:
                population_F[i] = trial_F
                population_pos_acel[i] = trial_pos_acel
                population_pos_freio[i] = trial_pos_freio
                if fitness_trial < best_fitness:
                    best_fitness = fitness_trial
                    best_F = trial_F
                    best_pos_acel = trial_pos_acel
                    best_pos_freio = trial_pos_freio
                    best_pos_final = pos_final_trial

        # Armazenando os dados da iteração
        fitness_history.append(best_fitness)
        F_history.append(best_F)
        pos_acel_history.append(best_pos_acel)
        pos_freio_history.append(best_pos_freio)

        # Critério de parada por convergência
        if best_fitness < tol:
            print(f"Convergência atingida na iteração {iteration + 1}!")
            break

    # Plotando os gráficos da evolução
    plt.figure(figsize=(12, 8))

    # Gráfico 1: Evolução da diferença para a posição desejada
    plt.subplot(3, 1, 1)
    plt.plot(fitness_history, label='Diferença para a posição desejada (m)')
    plt.xlabel('Iteração')
    plt.ylabel('Diferença (m)')
    plt.title('Evolução da Diferença para a Posição Desejada')
    plt.legend()
    plt.grid()

    # Gráfico 2: Evolução da força de frenagem (F_pedal)
    plt.subplot(3, 1, 2)
    plt.plot(F_history, label='Força de Frenagem (N)', color='orange')
    plt.xlabel('Iteração')
    plt.ylabel('F_pedal (N)')
    plt.title('Evolução da Força de Frenagem')
    plt.legend()
    plt.grid()

    # Gráfico 3: Evolução das posições de aceleração máxima e frenagem
    plt.subplot(3, 1, 3)
    plt.plot(pos_acel_history, label='Posição de Aceleração Máxima (m)', color='green')
    plt.plot(pos_freio_history, label='Posição de Frenagem (m)', color='red')
    plt.xlabel('Iteração')
    plt.ylabel('Posição (m)')
    plt.title('Evolução das Posições de Aceleração Máxima e Frenagem')
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.show()

    return best_F, best_pos_acel, best_pos_freio, best_fitness, best_pos_final

# Função para gerar a animação e os gráficos
def generate_animation_and_plots(F_pedal, pos_acel_max, pos_freio, pos_desejada, V1):
    # Parâmetros do sistema
    tr = 1.0   # Tempo de reação (s)
    ta = 0.25  # Tempo de aplicação do freio (s)
    tb = 0.3   # Tempo de transição da desaceleração (s)
    m = 230  # Massa do carro (kg)
    amax = F_pedal / m  # Desaceleração máxima (m/s²)
    V2 = V1 - (0.5 * amax * tb)  # Velocidade após a transição da desaceleração
    ts = V2 / amax  # Tempo para parar completamente após a desaceleração máxima

    # Simular o movimento do carro com os parâmetros otimizados
    _, _, _, _, t_values, acceleration, velocity, position, _ = simulate_car(F_pedal, pos_acel_max, pos_freio, pos_desejada)

    # Definindo os tempos correspondentes a S1, S2 e S3
    t_S0 = t_values[np.argmax(position >= pos_acel_max)]
    t_S1 = t_values[np.argmax(position >= pos_freio)]
    t_S2 = t_S1 + tr + ta
    t_S3 = t_S2 + tb + ts

    # Encontrando os índices correspondentes aos tempos S1, S2 e S3
    idx_S0 = np.argmin(np.abs(t_values - t_S0))
    idx_S1 = np.argmin(np.abs(t_values - t_S1))
    idx_S2 = np.argmin(np.abs(t_values - t_S2))
    idx_S3 = np.argmin(np.abs(t_values - t_S3))

    # Obtendo as posições correspondentes a S1, S2 e S3
    S0 = position[idx_S0]
    S1 = position[idx_S1]
    S2 = position[idx_S2]
    S3 = position[idx_S3]

    # Criando um diretório para salvar as imagens
    output_dir = "frames_carro"
    os.makedirs(output_dir, exist_ok=True)

    # Gerando as imagens
    for i, p in enumerate(position):
        if i % 10 == 0:  # Salvar a cada 5 pontos (ajustado conforme solicitado)
            plt.figure(figsize=(10, 2))
            plt.axhline(0, color='black', linewidth=2, linestyle='--', label='Centro da pista')  # Linha central da pista
            plt.fill_between([0, pos_desejada + 50], -1, 1, color='gray', alpha=0.3, label='Pista')  # Área da pista
            plt.plot(p, 0, 'bo', markersize=10, label='Carro')  # Carro
            plt.xlim(0, pos_desejada + 50)  # Limites do eixo x (pista)
            plt.ylim(-1, 1)  # Limites do eixo y (apenas para visualização)
            plt.xlabel('Posição (m)')
            plt.title(f'Posição do carro: {p:.2f} m')
            plt.legend()
            plt.grid(True)
            
            # Salvando a imagem
            plt.savefig(f"{output_dir}/frame_{i:04d}.png")
            plt.close()

    # Diretório onde os frames estão salvos
    input_dir = "frames_carro"

    # Lista de arquivos de frames
    frames = sorted([os.path.join(input_dir, img) for img in os.listdir(input_dir) if img.endswith(".png")])

    # Criando a animação (GIF ou MP4)
    output_file = "animacao_carro.gif"  # Ou "animacao_carro.mp4" para vídeo
    fps = 25  # Frames por segundo (ajuste conforme necessário)

    # Lendo os frames e criando a animação
    with imageio.get_writer(output_file, mode='I', fps=fps) as writer:
        for frame in frames:
            image = imageio.imread(frame)
            writer.append_data(image)

    print(f"Animação salva como: {output_file}")

    # Plotando os gráficos
    fig, ax = plt.subplots(3, 1, figsize=(8, 12))

    ax[0].plot(t_values, velocity, label='Velocidade (m/s)', color='b')
    ax[0].set_ylabel('Velocidade (m/s)')
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(t_values, acceleration, label='Aceleração (m/s²)', color='r')
    ax[1].set_ylabel('Aceleração (m/s²)')
    ax[1].grid()
    ax[1].legend()

    ax[2].plot(t_values, position, label='Posição (m)', color='g')
    ax[2].plot(t_S0, S0, 'yo', label=f'S0 = {S0:.2f} m')
    ax[2].plot(t_S1, S1, 'ro', label=f'S1 = {S1:.2f} m')  # Ponto S1
    ax[2].plot(t_S2, S2, 'bo', label=f'S2 = {S2:.2f} m')  # Ponto S2
    ax[2].plot(t_S3, S3, 'ko', label=f'S3 = {S3:.2f} m')  # Ponto S3
    ax[2].set_xlabel('Tempo (s)')
    ax[2].set_ylabel('Posição (m)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    plt.show()

    # Print da posição final
    print(f"Posição final onde o carro para: {position[-1]:.2f} m")
    print(f'S0 = {S0:.2f} m')
    print(f"S1 = {S1:.2f} m")
    print(f"S2 = {S2:.2f} m")
    print(f"S3 = {S3:.2f} m")
    
def prints():

    # 1. PRIMEIRO pegamos o input
    pos_desejada = int(input('Em que posição o carro deve parar?: '))

    # 2. SEGUNDO rodamos a otimização para descobrir os melhores parametros (best_pos_acel)
    best_F, best_pos_acel, best_pos_freio, best_fitness, best_pos_final = differential_evolution(pos_desejada)

    a_motor_fixo = 4.0 # TEM QUE SER O MESMO VALOR QUE ESTÁ DENTRO DA simulate_car
    v_limite = 26.7
    
    # Torricelli para descobrir a velocidade atingida na aceleração
    v_atingida = np.sqrt(2 * a_motor_fixo * best_pos_acel)
    
    # O V1 final é o menor entre o limite do carro e o que ele conseguiu atingir
    best_V1 = min(v_atingida, v_limite)
    # -------------------------------
    # Resultados finais
    print("\nResultados finais:")
    print(f"Melhor F_pedal encontrado: {best_F} N")
    print(f"Velocidade Máxima (V1) atingida: {best_V1:.2f} m/s ({best_V1*3.6:.0f} km/h)")
    print(f"Melhor pos_acel_max encontrada: {best_pos_acel} m")
    print(f"Melhor pos_freio encontrada: {best_pos_freio} m")
    print(f"Posição final do carro: {best_pos_final} m")
    print(f"Diferença para a posição desejada: {best_fitness} m")

    generate_animation_and_plots(best_F, best_pos_acel, best_pos_freio, pos_desejada, best_V1)

prints()