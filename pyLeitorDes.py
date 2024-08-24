import numpy as np
import matplotlib.pyplot as plt

NPAR = 8  # Número de variáveis dependentes (número de colunas)
NPT = 100  # Número de pontos (número de linhas)
NTOT = NPAR * NPT

# Lendo o arquivo binário
with open("invtd.des", "rb") as file:
    labels = np.fromfile(file, dtype=np.int32, count=NPAR)
    abscisse = np.fromfile(file, dtype=np.float32, count=NPT)
    coordinate = np.fromfile(file, dtype=np.float32, count=NTOT).reshape(NPAR, NPT)

# Plotando os dados
plt.figure(figsize=(10, 6))

for i in range(NPAR):
    plt.plot(abscisse, coordinate[i], label=f'Var {i+1}')

plt.xlabel('Tempo (t)')
plt.ylabel('Valores')
plt.title('Plot dos Dados Gerados')
plt.legend()
plt.grid(True)
plt.show()
