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
        # Retorna V1 como 0 no final para manter consistencia de 9 retornos
        return float('inf'), pos_acel_max, pos_freio, 0, None, None, None, None, 0

    # Parâmetros Físicos Fixos
    V_limite_carro = 26.7  # Velocidade máxima mecânica do carro (aprox 100 km/h)
    a_motor = 4.0          # Aceleração confortável/realista do motor (m/s²)
    
    # 1. Calcula a velocidade atingida baseada na distância de aceleração disponível
    V_atingida = np.sqrt(2 * a_motor * pos_acel_max)
    
    # 2. A velocidade V1 agora é dinâmica
    V1 = min(V_atingida, V_limite_carro)
    
    # 3. Recalculamos a aceleração efetiva e tempos
    a_acel = a_motor 
    t_acel = V1 / a_acel 

    # Outros parametros
    tr = 1.0   # Tempo de reação (s)
    ta = 0.25  # Tempo de aplicação do freio (s)
    tb = 0.3   # Tempo de transição da desaceleração (s)
    m = 230    # Massa do carro (kg)
    
    amax = F_pedal / m  # Desaceleração máxima (m/s²)
    V2 = V1 - (0.5 * amax * tb)

    # Tempo para percorrer de pos_acel_max até pos_freio (velocidade constante)
    dist_const = pos_freio - pos_acel_max 
    t_const = dist_const / V1 if V1 > 0 else 0

    # Tempo de desaceleração até parar
    if V2 < 0: V2 = 0 # Proteção
    ts = V2 / amax if amax > 0 else 0
    
    t_desacel = tr + ta + tb + ts 
    t_total = t_acel + t_const + t_desacel

    # Criando os vetores de tempo
    t_values = np.linspace(0, t_total, 1000)

    # Calculando a aceleração em cada instante
    acceleration = np.zeros_like(t_values)
    for i, t in enumerate(t_values):
        if t < t_acel:
            acceleration[i] = a_acel
        elif t < t_acel + t_const:
            acceleration[i] = 0
        elif t < t_acel + t_const + tr + ta:
            acceleration[i] = 0
        elif t < t_acel + t_const + tr + ta + tb:
            tb_elapsed = t - (t_acel + t_const + tr + ta)
            acceleration[i] = -amax * (tb_elapsed / tb)
        else:
            acceleration[i] = -amax 

    velocity = cumtrapz(acceleration, t_values, initial=0)
    velocity[velocity < 0] = 0 # Correção numérica
    position = cumtrapz(velocity, t_values, initial=0)
    pos_final = position[-1]

    # !!! IMPORTANTE: Retornando 9 valores (V1 é o último) !!!
    return abs(pos_final - pos_desejada), pos_acel_max, pos_freio, pos_final, t_values, acceleration, velocity, position, V1

# Algoritmo de Evolução Diferencial Autoadaptativa
def differential_evolution(pos_desejada, F_min=445, F_max=823, pos_min=0, pos_max=1000, pop_size=50, max_iter=200, CR=0.85, F_scale=0.5, tol=1e-3):
    # Inicialização da população
    population_F = np.random.uniform(F_min, F_max, pop_size)
    population_pos_acel = np.random.uniform(pos_min, pos_desejada, pop_size) 
    population_pos_freio = np.random.uniform(population_pos_acel, pos_desejada)

    best_fitness = float('inf')
    best_F = None
    best_pos_acel = None
    best_pos_freio = None
    best_pos_final = None

    fitness_history = []
    F_history = []
    pos_acel_history = []
    pos_freio_history = []

    for iteration in range(max_iter):
        for i in range(pop_size):
            idxs = [idx for idx in range(pop_size) if idx != i]
            a_F, b_F, c_F = population_F[np.random.choice(idxs, 3, replace=False)]
            a_pos_acel, b_pos_acel, c_pos_acel = population_pos_acel[np.random.choice(idxs, 3, replace=False)]
            a_pos_freio, b_pos_freio, c_pos_freio = population_pos_freio[np.random.choice(idxs, 3, replace=False)]

            mutant_F = a_F + F_scale * (b_F - c_F)
            mutant_pos_acel = a_pos_acel + F_scale * (b_pos_acel - c_pos_acel)
            mutant_pos_freio = a_pos_freio + F_scale * (b_pos_freio - c_pos_freio)

            trial_F = np.where(np.random.rand() < CR, mutant_F, population_F[i])
            trial_pos_acel = np.where(np.random.rand() < CR, mutant_pos_acel, population_pos_acel[i])
            trial_pos_freio = np.where(np.random.rand() < CR, mutant_pos_freio, population_pos_freio[i])

            trial_F = np.clip(trial_F, F_min, F_max)
            trial_pos_acel = np.clip(trial_pos_acel, pos_min, pos_desejada)
            trial_pos_freio = np.clip(trial_pos_freio, trial_pos_acel, pos_desejada)

            # !!! AQUI O UNPACKING TEM QUE TER UM "_" EXTRA PARA O V1 !!!
            fitness_trial, _, _, pos_final_trial, _, _, _, _, _ = simulate_car(trial_F, trial_pos_acel, trial_pos_freio, pos_desejada)
            fitness_current, _, _, _, _, _, _, _, _ = simulate_car(population_F[i], population_pos_acel[i], population_pos_freio[i], pos_desejada)

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

        fitness_history.append(best_fitness)
        F_history.append(best_F)
        pos_acel_history.append(best_pos_acel)
        pos_freio_history.append(best_pos_freio)

        if best_fitness < tol:
            print(f"Convergência atingida na iteração {iteration + 1}!")
            break

    # Plotting simplificado para nao ocupar muito espaço (pode usar o seu completo se preferir)
    plt.figure(figsize=(12, 8))
    plt.subplot(3, 1, 1); plt.plot(fitness_history); plt.title('Fitness')
    plt.subplot(3, 1, 2); plt.plot(F_history); plt.title('F_pedal')
    plt.subplot(3, 1, 3); plt.plot(pos_acel_history, label='Acel'); plt.plot(pos_freio_history, label='Freio'); plt.legend()
    plt.tight_layout(); plt.show()

    return best_F, best_pos_acel, best_pos_freio, best_fitness, best_pos_final

# Função para gerar a animação e os gráficos
def generate_animation_and_plots(F_pedal, pos_acel_max, pos_freio, pos_desejada, V1):
    # Parâmetros do sistema
    tr = 1.0   # Tempo de reação (s)
    ta = 0.25  # Tempo de aplicação do freio (s)
    tb = 0.3   # Tempo de transição da desaceleração (s)
    m = 230  # Massa do carro (kg)
    amax = F_pedal / m 
    V2 = V1 - (0.5 * amax * tb) 
    if V2 < 0: V2 = 0
    ts = V2 / amax if amax > 0 else 0

    # Note o "_" extra no final aqui também
    _, _, _, _, t_values, acceleration, velocity, position, _ = simulate_car(F_pedal, pos_acel_max, pos_freio, pos_desejada)

    t_S0 = t_values[np.argmax(position >= pos_acel_max)] if np.any(position >= pos_acel_max) else 0
    t_S1 = t_values[np.argmax(position >= pos_freio)] if np.any(position >= pos_freio) else 0
    t_S2 = t_S1 + tr + ta
    t_S3 = t_S2 + tb + ts

    idx_S0 = np.argmin(np.abs(t_values - t_S0))
    idx_S1 = np.argmin(np.abs(t_values - t_S1))
    idx_S2 = np.argmin(np.abs(t_values - t_S2))
    idx_S3 = np.argmin(np.abs(t_values - t_S3))

    S0 = position[idx_S0]
    S1 = position[idx_S1]
    S2 = position[idx_S2]
    S3 = position[idx_S3]

    output_dir = "frames_carro"
    os.makedirs(output_dir, exist_ok=True)

    for i, p in enumerate(position):
        if i % 5 == 0:  
            plt.figure(figsize=(10, 2))
            plt.axhline(0, color='black', linewidth=2, linestyle='--')  
            plt.fill_between([0, pos_desejada + 50], -1, 1, color='gray', alpha=0.3)  
            plt.plot(p, 0, 'bo', markersize=10)  
            plt.xlim(0, pos_desejada + 50)  
            plt.ylim(-1, 1)  
            plt.xlabel('Posição (m)')
            curr_vel = velocity[i] if i < len(velocity) else 0
            plt.title(f'Pos: {p:.2f} m | Vel: {curr_vel:.1f} m/s')
            plt.grid(True)
            plt.savefig(f"{output_dir}/frame_{i:04d}.png")
            plt.close()

    input_dir = "frames_carro"
    frames = sorted([os.path.join(input_dir, img) for img in os.listdir(input_dir) if img.endswith(".png")])
    output_file = "animacao_carro.gif" 
    
    with imageio.get_writer(output_file, mode='I', fps=25) as writer:
        for frame in frames:
            image = imageio.imread(frame)
            writer.append_data(image)

    print(f"Animação salva como: {output_file}")

    fig, ax = plt.subplots(3, 1, figsize=(8, 12))
    ax[0].plot(t_values, velocity, label=f'Velocidade (Max: {V1:.2f} m/s)', color='b')
    ax[0].grid(); ax[0].legend()
    ax[1].plot(t_values, acceleration, label='Aceleração', color='r')
    ax[1].grid(); ax[1].legend()
    ax[2].plot(t_values, position, label='Posição', color='g')
    ax[2].plot(t_S0, S0, 'yo', label='S0'); ax[2].plot(t_S1, S1, 'ro', label='S1')
    ax[2].plot(t_S2, S2, 'bo', label='S2'); ax[2].plot(t_S3, S3, 'ko', label='S3')
    ax[2].grid(); ax[2].legend()
    plt.tight_layout(); plt.show()

    print(f"Posição final onde o carro para: {position[-1]:.2f} m")

def prints():
    # 1. PRIMEIRO pegamos o input
    pos_desejada = int(input('Em que posição o carro deve parar?: '))

    # 2. SEGUNDO rodamos a otimização para descobrir os melhores parametros (best_pos_acel)
    best_F, best_pos_acel, best_pos_freio, best_fitness, best_pos_final = differential_evolution(pos_desejada)

    # 3. TERCEIRO calculamos a V1 baseado nos parametros descobertos acima
    a_motor_fixo = 12.0 
    v_limite = 26.7
    
    v_atingida = np.sqrt(2 * a_motor_fixo * best_pos_acel)
    best_V1 = min(v_atingida, v_limite)
    
    # 4. AGORA SIM podemos printar e gerar animação
    print("\nResultados finais:")
    print(f"Melhor F_pedal encontrado: {best_F:.2f} N")
    print(f"Melhor pos_acel_max encontrada: {best_pos_acel:.2f} m")
    print(f"Melhor pos_freio encontrada: {best_pos_freio:.2f} m")
    print(f"Velocidade Máxima (V1) atingida: {best_V1:.2f} m/s ({best_V1*3.6:.0f} km/h)")
    print(f"Posição final do carro: {best_pos_final:.2f} m")
    print(f"Diferença para a posição desejada: {best_fitness:.4f} m")

    generate_animation_and_plots(best_F, best_pos_acel, best_pos_freio, pos_desejada, best_V1)

prints()