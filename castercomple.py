import numpy as np
import matplotlib.pyplot as plt

def calcular_torque_autoalinhante(C, slip_angle, p, r, angulos_caster):
    """
   Calcula o torque auto-alinhante e filtra apenas os valores onde 
    o mechanical trail está dentro do intervalo permitido [0.0098, 0.02] m.
    """
    angulos_rad = np.radians(angulos_caster)
    mechanical_trail = r * np.tan(angulos_rad)

    # Filtrar apenas os valores válidos
    mascara_valida = (mechanical_trail >= 0.0098) & (mechanical_trail <= 0.02)
    angulos_validos = angulos_caster[mascara_valida]
    mechanical_trail_validos = mechanical_trail[mascara_valida]

    # Calcular força lateral
    F_y = C * np.radians(slip_angle)  # slip angle em radianos

    # Calcular torque auto-alinhante
    Mz_validos = F_y * (p + mechanical_trail_validos)

    return angulos_validos, Mz_validos

# Parâmetros fixos
C = 9170     # N/rad   gera aproximadadamente uma f_y = 800 N
p = 0.01       # m
r = 0.24       # m
angulos_caster = np.linspace(0, 5, 100)

# Slip angles a testar
slip_angles = [1, 3, 5]

# Plot
plt.figure(figsize=(8, 5))

for slip in slip_angles:
    angulos_filtrados, Mz = calcular_torque_autoalinhante(C, slip, p, r, angulos_caster)
    plt.plot(angulos_filtrados, Mz, label=f'Slip angle = {slip}°')

plt.xlabel('Ângulo de Caster (°)')
plt.ylabel('Torque Auto-Alinhante (Nm)')
plt.title('Torque Auto-Alinhante vs Ângulo de Caster para Vários Slip Angles')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
