import math
import matplotlib.pyplot as plt

# Funções para cálculos aerodinâmicos
def calcular_lift(CL, rho, v, A):
    return -CL * (0.5 * rho * v**2 * A)

def calcular_CL_alpha(AR):
    return 2 * math.pi / (1 + 2 / AR)

def calcular_CL(CL_alpha, alpha, alpha_L0):
    return CL_alpha * (alpha + alpha_L0)

def calcular_CD_induzido(CL, AR):
    return (CL ** 2) / (math.pi * AR)

def calcular_arrasto_viscoso(rho, v, Csf, Upper_Area, Cf, Front_Area):
    Skin_Friction_Drag = 0.5 * rho * v**2 * Csf * Upper_Area
    Form_Drag = 0.5 * rho * v**2 * Cf * Front_Area
    return Skin_Friction_Drag + Form_Drag

def calcular_arrasto_induzido(CL, A, AR, rho, v): 
    return calcular_CD_induzido(CL, AR) * (0.5 * rho * v**2 * A)

# Dados do exemplo
rho = 1.184                     # Densidade do ar em kg/m^3
A = 1.5                         # Área da asa em m^2
alpha = 8 * math.pi / 180       # Ângulo de ataque em radianos
alpha_L0 = 0                    # Valor de CL quando o ângulo de ataque é zero
AR = 4                          # Razão de aspecto das asas
e = 0.8                         # Eficiência da Asa (0.7 - 0.9)
Csf = 0.005                     # Coeficiente de Fricção da Superfície (Skin-Friction Coefficient)
Cf = 0.3                        # Coeficiente de Forma (Drag Form Coefficient)
Front_Area = 1.7                # Área de referência frontal em m²
Upper_Area = 2.2                # Área de superfície molhada em m²

# Cálculos
CL_alpha = calcular_CL_alpha(AR)
CL = calcular_CL(CL_alpha, alpha, alpha_L0)

# Listas para armazenar os resultados
velocidades = range(0, 61, 5)  # Velocidades de 0 m/s a 60 m/s com passos de 5 m/s
lifts = []
arrastos_viscosos = []
arrastos_induzidos = []
arrastos_totais = []

for v in velocidades:
    lift = calcular_lift(CL, rho, v, A)
    CDi = calcular_CD_induzido(CL, AR)
    Viscous_Drag = calcular_arrasto_viscoso(rho, v, Csf, Upper_Area, Cf, Front_Area)
    Induced_Drag = calcular_arrasto_induzido(CL, A, AR, rho, v)          
    Arrasto_Total = Induced_Drag + Viscous_Drag

    lifts.append(lift)
    arrastos_viscosos.append(Viscous_Drag)
    arrastos_induzidos.append(Induced_Drag)
    arrastos_totais.append(Arrasto_Total)

print("Arrasto Total:", Arrasto_Total)

# Plotando os gráficos
plt.figure(figsize=(14, 7))

# Gráfico Arrasto x Velocidade
plt.subplot(1, 3, 1)
plt.plot(velocidades, arrastos_totais, label='Arrasto Total')
plt.plot(velocidades, arrastos_viscosos, label='Arrasto Viscoso', linestyle='--')
plt.plot(velocidades, arrastos_induzidos, label='Arrasto Induzido', linestyle=':')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Arrasto (N)')
plt.title('Arrasto x Velocidade')
plt.legend()

# Gráfico Lift x Velocidade
plt.subplot(1, 3, 2)
plt.plot(velocidades, lifts, label='Lift', color='green')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Lift (N)')
plt.title('Lift x Velocidade')
plt.legend()

# Gráfico Arrasto e Lift x Velocidade
plt.subplot(1, 3, 3)
plt.plot(velocidades, arrastos_totais, label='Arrasto Total')
plt.plot(velocidades, lifts, label='Lift', color='green')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.title('Arrasto e Lift x Velocidade')
plt.legend()

plt.tight_layout()
plt.show()
