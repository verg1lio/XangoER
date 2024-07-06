import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class DifferentialEvolutionOptimizer(BaseEstimator, TransformerMixin):
    def __init__(self, force=100, displacement=0.1234, NP=4, F=0.8, CR=0.9, max_iterations=500, patience=100):
        self.force = force
        #self.mass = mass
        self.displacement = displacement
        self.NP = NP
        self.F = F
        self.CR = CR
        self.max_iterations = max_iterations
        self.patience = patience
        self.best_solution_ = None
        self.best_fitness_ = np.inf
        self.best_solutions_history = []

    def objective_function(self, k):
        if k <= 0 or self.force == 0:
            return np.inf

        x = self.force / k
        fitness = np.abs(x - self.displacement)
        return fitness

    def fit(self):
        # Verificação inicial para força zero
        if self.force == 0:
            print("Força é zero. Não é possível realizar a otimização.")
            self.best_solution_ = None
            self.best_fitness_ = np.inf
            return self

        population = np.random.uniform(1, 1000, size=(self.NP,))
        fitness = np.array([self.objective_function(k) for k in population])
        best_index = np.argmin(fitness)
        self.best_solution_ = population[best_index]
        self.best_fitness_ = fitness[best_index]
        self.best_solutions_history.append((self.best_solution_, self.best_fitness_))

        no_improvement_counter = 0
        previous_best_fitness = self.best_fitness_

        for iteration in range(self.max_iterations):
            for i in range(self.NP):
                candidates = np.random.choice(self.NP, 3, replace=False)
                a, b, c = population[candidates]
                R = np.random.randint(0, 1)
                new_solution = np.copy(population[i])
                if np.random.rand() < self.CR or R == 0:
                    new_solution = a + self.F * (b - c)
                new_solution = np.clip(new_solution, 1, 1000000000000000000000)
                new_fitness = self.objective_function(new_solution)
                if new_fitness <= fitness[i]:
                    population[i] = new_solution
                    fitness[i] = new_fitness
                    if new_fitness < self.best_fitness_:
                        self.best_solution_ = new_solution
                        self.best_fitness_ = new_fitness
                        no_improvement_counter = 0

            if (iteration + 1) % 10 == 0:
                self.best_solutions_history.append((self.best_solution_, self.best_fitness_))

            if previous_best_fitness == self.best_fitness_:
                no_improvement_counter += 1
            else:
                no_improvement_counter = 0
            previous_best_fitness = self.best_fitness_

            if no_improvement_counter >= self.patience:
                print(f"Parando antecipadamente na iteração {iteration + 1} devido à convergência.")
                break

        if (iteration + 1) % 10 != 0:
            self.best_solutions_history.append((self.best_solution_, self.best_fitness_))

        return self

    def transform(self, X):
        return X

    def fit_transform(self, X):
        self.fit()
        return self.transform(X)


# Exemplo de uso
optimizer = DifferentialEvolutionOptimizer()  # A força pode ser definida aqui
optimizer.fit()

# Exibindo o melhor indivíduo encontrado ao final
best_solution_final = optimizer.best_solution_
best_fitness_final = optimizer.best_fitness_
if best_solution_final is not None:
    print(f"Melhor valor de k encontrado ao final: {best_solution_final}")
    print(f"Valor de aptidão da melhor solução ao final: {best_fitness_final}")
else:
    print("Nenhuma solução válida encontrada.")

# Exibindo os melhores indivíduos a cada 10 iterações
print("\nMelhores indivíduos a cada 10 iterações:")
for i, (solution, fitness) in enumerate(optimizer.best_solutions_history):
    print(f"Iteração {i * 10}: Melhor k = {solution}, Aptidão = {fitness}")
