import numpy as np
import matplotlib.pyplot as plt

# Parâmetros do modelo
tempo_total = 100  # Tempo total da simulação
passo_tempo = 0.1  # Passo de tempo
coef_atrito = 0.3  # Coeficiente de atrito
pressao = 10  # Pressão aplicada na pastilha

# Inicialização das variáveis
tempo = np.arange(0, tempo_total, passo_tempo)
desgaste_pastilha = np.zeros_like(tempo)
forca_atrito = np.zeros_like(tempo)

# Simulação do desgaste da pastilha e força de atrito
for i in range(1, len(tempo)):
    # Modelo simplificado de desgaste da pastilha (apenas para demonstração)
    desgaste_pastilha[i] = desgaste_pastilha[i-1] + 0.01 * passo_tempo

    # Cálculo da força de atrito
    forca_atrito[i] = coef_atrito * pressao

# Plotagem dos resultados
plt.figure(figsize=(10, 6))

plt.subplot(2, 1, 1)
plt.plot(tempo, desgaste_pastilha, label='Desgaste da Pastilha')
plt.title('Desgaste da Pastilha de Freio ao Longo do Tempo')
plt.xlabel('Tempo')
plt.ylabel('Desgaste')
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(tempo, forca_atrito, label='Força de Atrito')
plt.title('Força de Atrito Gerada ao Longo do Tempo')
plt.xlabel('Tempo')
plt.ylabel('Força de Atrito')
plt.legend()

plt.tight_layout()
plt.show()
