import math

# Funções para cálculos aerodinâmicos
def calcular_lift(CL, rho, v, A):
    return CL * (0.5 * rho * v**2 * A)

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
v = 60                          # Velocidade do ar em m/s
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
lift = calcular_lift(CL, rho, v, A)
CDi = calcular_CD_induzido(CL, AR)
Viscous_Drag = calcular_arrasto_viscoso(rho, v, Csf, Upper_Area, Cf, Front_Area)
Induced_Drag = calcular_arrasto_induzido(CL, A, AR, rho, v)
Arrasto_Total = Induced_Drag + Viscous_Drag

# Resultados
print("CL_alpha:", CL_alpha)
print("CL:", CL)
print("Lift:", lift)
print("CDi:", CDi)
print("Arrasto Viscoso:", Viscous_Drag)
print("Arrasto Induzido:", Induced_Drag)
print("Arrasto Total:", Arrasto_Total)
