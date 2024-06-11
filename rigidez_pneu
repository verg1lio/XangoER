import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

# Constantes de rigidez do pneu
k = 190.35  # Constante de rigidez do pneu (N/mm)

# Intervalo de deformação (em mm)
deformacao = np.linspace(0, 30, 1000) 

# Pontos fornecidos
pontos_forca = np.array([-29.41995, 323.28045, 558.96405, 921.0981, 1214.5456, 1579.17965, 1875.12015, 1893.33045, 2149.96135, 2435.2792, 2689.9101, 2817.92255, 3058.2738, 3312.9047, 3567.5356, 3778.49725, 4033.12815, 4287.75905, 4542.38995, 4785.8328, 5029.27565])
pontos_deformacao = np.array([0.0, 2.2, 3.6, 5.7, 7.2, 8.8, 9.7, 10.3, 11.3, 12.5, 13.6, 14.4, 15.3, 16.7, 17.7, 18.9, 20.1, 21.2, 22.7, 23.9, 25.2])

# Função de reta: y = m * x + b
def funcao_reta(x, m, b):
    return m * x + b

# Ajuste de curva usando scipy.optimize.curve_fit
params_otimizados, _ = opt.curve_fit(funcao_reta, pontos_deformacao, pontos_forca)
m_otimizado, b_otimizado = params_otimizados

print(f"Parâmetros Otimizados: m = {m_otimizado}, b = {b_otimizado}")

# Valores ajustados
pontos_forca_ajustados = funcao_reta(deformacao, m_otimizado, b_otimizado)

# Plotando gráfico deformação x carga
plt.figure(figsize=(18, 6))
plt.scatter(pontos_deformacao, pontos_forca, color='red', label='Pontos Fornecidos')
plt.plot(deformacao, pontos_forca_ajustados, color='green', label='Ajuste de Curva')
plt.xlabel('Deformação (mm)')
plt.ylabel('Carga Aplicada (N)')
plt.title('Deformação vs Carga Aplicada')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
