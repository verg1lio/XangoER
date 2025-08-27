import numpy as np
import matplotlib.pyplot as plt
import pygad  # Importa a biblioteca do Algoritmo Genético
import sys    # Usado para criar a barra de progresso
from pwt import Motor               
from tire import Tire

class Drivetrain:
    def __init__(self, massa, massa_roda,
                 raio_pneu, i, tempo_i,
                 tire_friction_coef, tire_Fz):
        
        self.massa = massa
        self.massa_roda = massa_roda
        self.raio_pneu = raio_pneu
        self.i = i
        self.tempo_i = tempo_i
        self.tire_friction_coef = tire_friction_coef
        self.tire_Fz = tire_Fz

        self.motor = Motor(
            rs=0.04585, ld=0.00067, lq=0.00067, jm=0.05769,
            kf=0.1, lambda_m=0.13849, p=10, valor_mu=0.9
        )
        self.tire = Tire(
            tire_Fz=850, tire_Sa=0, tire_Ls=1,
            tire_friction_coef=self.tire_friction_coef, tire_Ca=0
        )

    def CarPerformanceDistancia(self, distancia_final=75.0):
        # --- parâmetros fixos / físicos ---
        peso = self.massa * 9.81
        coeficiente_arrasto = 0.54
        densidade_ar = 1.162
        area_frontal = 1.06
        raio = self.raio_pneu * 0.001
        c_r = 0.012
        b = 0.01
        eficiencia_transmissao = 0.95
        num_roda_motriz = 2
        num_rodas = 4
        Fz_por_roda = self.tire_Fz if self.tire_Fz is not None else (peso / num_rodas)
        pacejka_params = [0.333, 1.627, 1, 4.396, 931.4, 366.4]
        
        dt = 0.001
        tempo = self.tempo_i
        v_linear_ms = 0.0
        omega_roda = 0.0
        deslocamento = 0.0
        
        inercia_motor_refletida = self.motor.jm * (self.i**2)
        Jt = (0.5 * self.massa_roda * raio**2) + (inercia_motor_refletida / num_roda_motriz)

        eps_v = 1e-3
        performance = []

        # Loop de segurança para evitar simulações infinitas
        max_time = 20 # segundos
        
        while deslocamento < distancia_final:
            

            # forças resistivas
            forca_rolamento = c_r * peso
            forca_arrasto = 0.5 * densidade_ar * coeficiente_arrasto * area_frontal * v_linear_ms**2
            forca_resistiva_total = forca_rolamento + forca_arrasto

            # slip ratio
            denom_v = max(v_linear_ms, eps_v)
            slip_ratio_roda = (omega_roda * raio - v_linear_ms) / denom_v

            # modelo do pneu
            self.tire.tire_Fz = Fz_por_roda
            self.tire.tire_Ls = slip_ratio_roda
            _, _, Fx_roda = self.tire.Tire_forces(pacejka_params)

            # torque de carga no motor
            torque_rr = forca_rolamento * raio
            torque_fa = forca_arrasto * raio
            torque_atr_mec = b * (self.motor.wm) 
            torque_carga = torque_rr + torque_fa + torque_atr_mec
            torque_carga_motor = torque_carga / (self.i * eficiencia_transmissao)

            self.motor.set_load(torque_carga_motor)
            torque_motor_gerado = self.motor.step(dt)

            # torque na roda
            torque_trativo_rodas_total = torque_motor_gerado * self.i * eficiencia_transmissao
            torque_por_roda = torque_trativo_rodas_total / num_roda_motriz

            # EDOs
            forca_trativa_real = num_roda_motriz * Fx_roda
            dvdt = (forca_trativa_real - forca_resistiva_total) / self.massa
            domegadt = (torque_por_roda - raio * Fx_roda) / Jt

            # integração
            v_linear_ms += dvdt * dt
            omega_roda += domegadt * dt
            deslocamento += v_linear_ms * dt
            tempo += dt

            # salvar histórico
            self.motor.store_history(tempo)
            performance.append({
                "tempo": float(tempo),
                "s": float(deslocamento),           # posição acumulada [m]
                "vlk": float(v_linear_ms * 3.6),    # km/h
                "vlm": float(v_linear_ms),          # m/s
                "va": float(omega_roda),            # rad/s
                "fa": float(forca_arrasto),
                "rr": float(forca_rolamento),
                "ff": float(forca_trativa_real),
                "fp": float(torque_por_roda / raio * num_roda_motriz),
                "sr": float(slip_ratio_roda),
                "ft_roda": float(Fx_roda)
            })

        return performance, np.array([p["tempo"] for p in performance]), tempo

    def printCarPerformance(self, performance):
        """
        Gera gráficos de desempenho do veículo com base no histórico do motor
        e nos dados de simulação armazenados em performance.
        """

        if len(performance) == 0:
            print("Sem dados de performance para plotar.")
            return

        motor_history = self.motor.history
        perf_dict = {k: [dic[k] for dic in performance] for k in performance[0]}

        # --- Figura 1: Resumo geral ---
        plt.figure(figsize=(14, 10))

        # Velocidade linear
        plt.subplot(2, 2, 1)
        plt.plot(perf_dict['tempo'], perf_dict['vlk'])
        plt.title("Velocidade Linear do Veículo")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade [km/h]")
        plt.grid(True)

        # Forças
        plt.subplot(2, 2, 2)
        plt.plot(perf_dict['tempo'], perf_dict['ff'], label='Força Trativa')
        plt.plot(perf_dict['tempo'], perf_dict['fa'], label='Força de Arrasto')
        plt.plot(perf_dict['tempo'], perf_dict['rr'], label='Resist. Rolamento')
        plt.title("Forças Atuantes")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Força [N]")
        plt.legend()
        plt.grid(True)

        # RPM do motor
        plt.subplot(2, 2, 3)
        rpm_motor = np.array(motor_history['velocidade']) * 60 / (2 * np.pi)
        plt.plot(motor_history['tempo'], rpm_motor)
        plt.title("RPM do Motor")
        plt.xlabel("Tempo [s]")
        plt.ylabel("RPM")
        plt.grid(True)

        # Torques do motor
        plt.subplot(2, 2, 4)
        plt.plot(motor_history['tempo'], motor_history['torque_eletrico'], label='Torque Gerado')
        plt.plot(motor_history['tempo'], motor_history['torque_carga'], label='Torque de Carga', linestyle='--')
        plt.title("Torques no Motor")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Torque [Nm]")
        plt.legend()
        plt.grid(True)

        plt.tight_layout()
        plt.show()

         # Impressões resumidas
        print("Tempo [s]\tVelocidade Angular [rad/s]\tVelocidade Linear [km/h]")
        for i in range(0, len(performance), max(1, len(performance)//10)):
            t = performance[i]["tempo"]
            va = performance[i]["va"]
            vlk = performance[i]["vlk"]
            print(f"{t:.2f}\t\t{va:.2f}\t\t{vlk:.2f}")

        print("\nForça de Arrasto [N]\tForça Trativa Potencial [N]\tForça Trativa [N]")
        for i in range(0, len(performance), max(1, len(performance)//10)):
            fa = performance[i]["fa"]
            fp = performance[i]["fp"]
            ff = performance[i]["ff"]
            print(f"{fa:.2f}\t\t{fp:.2f}\t\t{ff:.2f}")

        # Plots secundários: velocidade angular e linear x tempo
        tempo = perf_dict["tempo"]
        velocidade_angular = perf_dict["va"]
        velocidade_linear = perf_dict["vlk"]

        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.plot(tempo, velocidade_angular)
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade Angular [rad/s]")
        plt.title("Velocidade Angular x Tempo")
        plt.grid(True)

        plt.subplot(1, 2, 2)
        plt.plot(tempo, velocidade_linear)
        plt.xlabel("Tempo [s]")
        plt.ylabel("Velocidade Linear [km/h]")
        plt.title("Velocidade Linear x Tempo")
        plt.grid(True)

        plt.tight_layout()
        plt.show()

# ==============================================================================
# INÍCIO DA IMPLEMENTAÇÃO DO ALGORITMO GENÉTICO
# ==============================================================================

# 1. FUNÇÃO DE APTIDÃO (FITNESS)
# Esta função recebe uma solução (relação de transmissão) e retorna um valor
# de aptidão. Queremos MINIMIZAR o tempo, então a aptidão será 1/tempo.
def fitness_func(ga_instance, solution, solution_idx):
    # Extrai a relação de transmissão da solução
    transmission_ratio = solution[0]
    
    # Cria uma instância do Drivetrain com a relação de transmissão atual
    dt_model = Drivetrain(
        massa=347, massa_roda=6, raio_pneu=220,
        i=transmission_ratio,  # Relação sendo testada pelo AG
        tempo_i=0, tire_friction_coef=1.45, tire_Fz=850
    )
    
    # Roda a simulação
    _, _, final_time = dt_model.CarPerformanceDistancia(distancia_final=75)
    
    # Calcula a aptidão. Adicionamos um valor pequeno ao denominador para evitar divisão por zero.
    fitness = 1.0 / (final_time + 1e-6)
    
    return fitness

# 2. FUNÇÃO DE CALLBACK PARA MOSTRAR O PROGRESSO
# Esta função será chamada ao final de cada geração.
def on_generation_callback(ga_instance):
    generation = ga_instance.generations_completed
    total_generations = ga_instance.num_generations
    progress = (generation / total_generations) * 100
    
    # Cria uma barra de progresso simples
    bar_length = 50
    filled_length = int(bar_length * generation // total_generations)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    
    # Escreve na mesma linha para atualizar a barra
    sys.stdout.write(f'\rProgresso da Otimização: |{bar}| {progress:.1f}% Concluído')
    sys.stdout.flush()

# 3. CONFIGURAÇÃO E EXECUÇÃO DO AG
if __name__ == "__main__":
    print("--- Otimização da Relação de Transmissão com Algoritmo Genético ---")

    # Parâmetros do Algoritmo Genético
    num_generations = 200       # Número de gerações (iterações)
    num_parents_mating = 5     # Número de soluções selecionadas como pais
    sol_per_pop = 40           # Número de soluções (cromossomos) em cada população
    num_genes = 1              # Estamos otimizando apenas um parâmetro: a relação de transmissão

    # Define o espaço de busca para a relação de transmissão. Ex: entre 3 e 15.
    gene_space = [{'low': 3.0, 'high': 15.0}]

    # Criação da instância do AG
    ga_instance = pygad.GA(
        num_generations=num_generations,
        num_parents_mating=num_parents_mating,
        fitness_func=fitness_func,
        sol_per_pop=sol_per_pop,
        num_genes=num_genes,
        gene_space=gene_space,
        gene_type=float,          # O gene é um número de ponto flutuante
        on_generation=on_generation_callback, # Função para mostrar o progresso
        # Parâmetros de mutação para explorar melhor o espaço de busca
        mutation_type="random",
        mutation_percent_genes=100, # Mutacionar o único gene
        random_mutation_min_val= -1.0, # Variação máxima da mutação
        random_mutation_max_val= 1.0
    )

    # Executa o AG
    ga_instance.run()
    
    print("\n\n--- Otimização Concluída! ---")

    # Extrai a melhor solução encontrada
    solution, solution_fitness, solution_idx = ga_instance.best_solution()
    best_ratio = solution[0]
    best_time = 1.0 / solution_fitness

    print(f"Melhor Relação de Transmissão Encontrada: {best_ratio:.3f}")
    print(f"Tempo Mínimo Estimado para 75m: {best_time:.4f} segundos")

    # Plota o gráfico de evolução da aptidão
    ga_instance.plot_fitness()

    # Roda e plota a simulação final com a melhor relação encontrada
    print("\n--- Gerando gráficos para a melhor solução encontrada... ---")
    
    best_drivetrain = Drivetrain(
        massa=347, massa_roda=6, raio_pneu=220,
        i=best_ratio,
        tempo_i=0, tire_friction_coef=1.45, tire_Fz=850
    )
    
    performance_final, tempo_vetor, _ = best_drivetrain.CarPerformanceDistancia(distancia_final=75)
    best_drivetrain.printCarPerformance(performance_final)
