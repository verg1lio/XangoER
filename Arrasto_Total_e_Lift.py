import math

# Primeiro vamos definir O Coef de Lift Induzido para obter o Lift da Asa (downforce = Lift negativo).
# Depois iremos 

# Função para calcular o coeficiente de sustentação de uma asa finita
def calcular_coeficiente_sustentacao(CL, rho, v, A):
    lift = CL * (0.5 * rho * v**2 * A)
    return lift

# Função para calcular a inclinação da sustentação em função da razão de aspecto
def calcular_inclinacao_sustentacao_AR(AR):
    CL_alpha = 2 * math.pi / (1 + 2 / AR)
    return CL_alpha

# Função para calcular o coeficiente de sustentação CL
def calcular_coeficiente_sustentacao_CLalpha(CL_alpha, alpha, alpha_L0):
    CL = CL_alpha * (alpha + alpha_L0)
    return CL

# Dados do exemplo
rho = 1.225                     # Densidade do ar em kg/m^3
v = 28                 # Velocidade do ar em m/s
A = 1.5                         # Área da asa em m^2
alpha =  8 * math.pi / 180      # Ângulo de ataque em graus
alpha_L0 = 0 * math.pi / 180    # Valor de CL quando o ângulo de ataque é zero
AR = 4                          # Razão de aspecto das duas asas

# Calcular CL_alpha
CL_alpha = calcular_inclinacao_sustentacao_AR(AR)

# Calcular CL
CL = calcular_coeficiente_sustentacao_CLalpha(CL_alpha, alpha, alpha_L0)

# Calcular lift
lift = calcular_coeficiente_sustentacao(CL, rho, v, A)

print("CL_alpha:", CL_alpha)
print("CL:", CL)
# print("lift:", lift)  

# Agora o CL obtido será fundamental, assim como o AR determinado, para o calculo do Arrasto Induzido.
# Vale Lembrar que o Arassto Total será o Arassto Induzido + o Arrasto Viscoso (Fricção + Forma). Ou seja: CD = CDi + CD_Viscous

# Calculo do Coeficiente de Arrasto Induzido
def Calculo_Induced_Coef (CL,AR):
    CDi = ((CL ** 2) / (math.pi * AR))
    return CDi

CDi = Calculo_Induced_Coef (CL,AR)

## ETAPA Validação: Arrasto Viscoso + Arrasto Induzido
# Arrasto Viscoso
def Calculo_Coef_Viscous_Drag (rho, v, Csf, Upper_Aerea, Cf, Front_Area): 

    # Calculo do arrasto por fricção
    Skin_Friction_Drag = 0.5 * rho * v**2 * Csf * Upper_Aerea

    # Calculo do arrasto por pressão
    Form_Drag = 0.5 * rho * v**2 * Cf * Front_Area

    # Arrasto Viscoso
    Viscous_Drag = Skin_Friction_Drag + Form_Drag
    return Viscous_Drag

# Arrasto induzido
def calc_arrasto_induzido(lift, e, A, AR):
    Induced_Drag = (lift ** 2) / (math.pi * e * A * AR)
    return Induced_Drag
print("CDi = ", CDi)
print(" ")

e = 0.8             # Eficiência da Asa (0.7 - 0.9)
Csf = 0.005         # coeficiente de Fricção da Superfície (Skin-Friction Coefficient)
Cf = 0.3            # coeficiente de Fomra (Drag Form Coefficient)
Front_Area = 1.8    # área de referência frontal em m²(exemplo)
Upper_Aerea = 2.2   # área de superfície molhada em m² (exemplo)
CD = Csf + Cf + CDi

Viscous_Drag = Calculo_Coef_Viscous_Drag (rho, v, Csf, Upper_Aerea, Cf, Front_Area)
Induced_Drag = calc_arrasto_induzido (lift, e, A, AR)

# Arrasto Total
Arrasto_Total = Induced_Drag + Viscous_Drag
print("Arrasto: ", Arrasto_Total)
print("Lift = ", lift)

#                                       Calculo do Drag por CDi, Csf, Cf

#def Calculo_Drag_por_Coefs (CD, rho, v, A):
#    Drag = CD * (0.5 * rho * v**2 * A)
#    return Drag

#Drag = Calculo_Drag_por_Coefs (CD, rho, v, A)
#print("O Coef de Arrasto Total é: ", CD)
#print("Arrato total 1: ", Drag)
