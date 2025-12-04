import numpy as np
import sys
import os
import io
import contextlib # Import necessário para silenciar os prints

# Garante que o Python encontre os módulos (models, simulation)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from models import BatteryPack, Inversor, Tire, Transmission, Vehicle, Motor
from simulation.Simulation import Simulation

# -----------------------------------------------------------
# Função de simulação 
# -----------------------------------------------------------
def simulation(i_ratio, distancia_desejada):
    log_capture = io.StringIO() 
    
    try:
        # Tenta rodar escondendo os prints normais
        with contextlib.redirect_stdout(log_capture):
            
            # --- CONFIGURAÇÃO (Mantida igual ao seu código) ---
            transmission = Transmission.Transmission(final_drive_ratio=float(i_ratio), efficiency=0.95)
            vehicle = Vehicle.Vehicle(mass=230.0, wheel_radius=0.275, drag_coeff=0.7789,
                                      frontal_area=0.68, rolling_resistance=0.015, road_grade=0.0)
            battery = BatteryPack.BatteryPack(tipo_celula='Li-ion', n_serie=264, n_paralelo=1, soc_inicial=1.0)
            tire = Tire.Tire(pacejka_params=[0.333, 1.627, 1, 4.396, 931.4, 366.4], tire_friction_coef=1.45)
            inversor = Inversor.Inversor(eficiencia=0.95, freq_chaveamento=10000)
            
            motor = Motor.Motor(rs=0.00706, ld=0.0000965, lq=0.0000965, jm=0.02521, kf=0.1,
                                lambda_m=0.04748, p=10, valor_mu=0.99,
                                TL=False, torque=0.0,
                                speed_ref=1200.0)

            t_max_simulacao = 30.0
            sim = Simulation(motor=motor, vehicle=vehicle, transmission=transmission,
                             battery=battery, tire=tire, inversor=inversor, 
                             tmax=t_max_simulacao, steps=10000)
            
            results = sim.simulate(t0=0, tf=t_max_simulacao)

        # --- ANÁLISE ---
        posicoes = results.get('vehicle_position')
        tempos = results.get('t')

        # SE FALHAR: Imprime o log que estava escondido para debug
        if posicoes is None or len(posicoes) == 0:
            print(f"\n[DEBUG ERRO] Falha na simulação com i={i_ratio:.4f}")
            print("--- Log da Simulação ---")
            print(log_capture.getvalue()) # <--- ISSO VAI MOSTRAR O ERRO REAL
            print("------------------------")
            return 1e9, None
            
        if posicoes[-1] < 1.0:
            return 1e5, None # Carro não andou

        if posicoes[-1] < distancia_desejada:
            return 100.0 + (distancia_desejada - posicoes[-1]), None
        
        idx_cruzamento = np.argmax(posicoes >= distancia_desejada)
        tempo_final = tempos[idx_cruzamento]
        
        if np.isnan(tempo_final):
            return 1e9, None

        return tempo_final, tempo_final

    except Exception as e:
        print(f"\n🔴 [EXCEÇÃO PYTHON | i={i_ratio:.4f}]: {str(e)}")
        # Em caso de exceção real, também printamos o log interno
        print(log_capture.getvalue())
        return 1e9, None

# -----------------------------------------------------------
# Evolução Diferencial 
# -----------------------------------------------------------
def diferentialEvo(
                    distancia_desejada,
                    i_min=2.0,
                    i_max=20.0,
                    pop_size=20,    
                    max_iter=30,    
                    CR=0.9,
                    F_scale=0.8,
                    tol=1e-4):

    print(f"--- Iniciando Otimização para {distancia_desejada} metros ---")
    
    population_i = np.random.uniform(i_min, i_max, pop_size)

    best_fitness = float('inf')
    best_i = None

    fitness_history = []
    i_history = []

    fitness_pop = np.zeros(pop_size)
    
    print("Avaliando população inicial...")

    for k in range(pop_size):
       
        fitness_pop[k], _ = simulation(population_i[k], distancia_desejada)

        if fitness_pop[k] < best_fitness:
            best_fitness = fitness_pop[k]
            best_i = population_i[k]
           

    for iteration in range(max_iter):
        
        # Barra de progresso limpa
        progresso = (iteration + 1) / max_iter * 100
        sys.stdout.write(f"\rProgresso: {progresso:.1f}% | Melhor i: {best_i:.4f} ({best_fitness:.4f}s)")
        sys.stdout.flush()

        for i in range(pop_size):

            idxs = [idx for idx in range(pop_size) if idx != i]
            if len(idxs) < 3: continue
            
            a, b, c = np.random.choice(idxs, 3, replace=False)
            mutant_i = population_i[a] + F_scale * (population_i[b] - population_i[c])
            trial_i = np.where(np.random.rand() < CR, mutant_i, population_i[i])
            trial_i = np.clip(trial_i, i_min, i_max)

            fitness_trial, _ = simulation(trial_i, distancia_desejada)

            if fitness_trial < fitness_pop[i]:
                population_i[i] = trial_i
                fitness_pop[i] = fitness_trial

                if fitness_trial < best_fitness:
                    best_fitness = fitness_trial
                    best_i = trial_i

        fitness_history.append(best_fitness)
        i_history.append(best_i)

        if len(fitness_history) > 10:
            if np.std(fitness_history[-10:]) < tol:
                print(f"\n✅ Convergência atingida na iteração {iteration + 1}")
                break

    print("\n\n---------------- RESULTADO FINAL ----------------")
    print(f"Distância Alvo: {distancia_desejada} m")
    print(f"Melhor relação i encontrada = {best_i:.4f}")
    print(f"Tempo de percurso = {best_fitness:.4f} s\n")

    return best_i, best_fitness, fitness_history, i_history

if __name__ == "__main__":
    distancia = 75.0  
    best_i, best_fitness, fit_hist, i_hist = diferentialEvo(distancia_desejada=distancia)