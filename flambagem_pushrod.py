
import math

D = 0.022225  # diâmetro externo [m]
d = 0.019177    # diâmetro interno [m]
L = 0.5    # comprimento entre apoios [m]

E = 205e9     # módulo de elasticidade [Pa]
sigma_y = 435e6  # tensão de escoamento [Pa]

# Condição de apoio (defina K)
# K = 1.0 -> ambas articuladas
# K = 0.7 -> engastada-articulada
# K = 0.5 -> ambas engastadas
# K = 2.0 -> engastada-livre
K = 1

# Carga aplicada (N)
P_aplicada = -7000  # defina a carga real aqui [N]


# =======================

# Área da seção
A = (math.pi / 4) * (D**2 - d**2)

# Momento de inércia
I = (math.pi / 64) * (D**4 - d**4)

# Raio de giração
r = math.sqrt(I / A)

# Esbeltez
lambda_ = (K * L) / r

# Limite de esbeltez para transição (Johnson x Euler)
# lambda_lim = math.pi * math.sqrt(E / sigma_y)       
lambda_lim = math.sqrt(2*(math.pi**2)*E / sigma_y)

# --- Flambagem de Euler ---
Pcr_Euler = (math.pi**2 * E * I) / ((K * L)**2)
sigma_cr_Euler = Pcr_Euler / A

# --- Flambagem de Johnson (coluna intermediária) ---
sigma_cr_Johnson = sigma_y * (1 - (sigma_y / (4 * math.pi**2 * E)) * (lambda_**2))
Pcr_Johnson = sigma_cr_Johnson * A

# --- Carga de escoamento simples (esmagamento) ---
P_esmag = sigma_y * A


# =======================
if lambda_ > lambda_lim:
    Pcrit = Pcr_Euler
    modo = "Flambagem elástica (Euler)"
else:
    Pcrit = Pcr_Johnson
    modo = "Flambagem inelástica (Johnson)"


# =======================
if abs(P_aplicada) >= Pcrit:
    resultado = f"⚠️ FALHA por {modo}! (P_aplicada = {P_aplicada:.1f} N > Pcr = {Pcrit:.1f} N)"
elif abs(P_aplicada) >= P_esmag:
    resultado = f"⚠️ FALHA por ESMAGAMENTO! (P_aplicada = {P_aplicada:.1f} N > P_esmag = {P_esmag:.1f} N)"
else:
    resultado = f"✅ Seguro: não ocorre flambagem nem esmagamento (P_aplicada = {P_aplicada:.1f} N)."


# =======================
print("=== ANÁLISE DE FLAMBAGEM E ESMAGAMENTO ===")
print(f"Diâmetro externo D = {D*1000:.2f} mm")
print(f"Diâmetro interno d = {d*1000:.2f} mm")
print(f"Comprimento L = {L*1000:.3f} mm")
print(f"Módulo de elasticidade E = {E/1e9:.1f} GPa")
print(f"Tensão de escoamento σy = {sigma_y/1e6:.1f} MPa")
print(f"Fator de comprimento efetivo K = {K}")
print("--------------------------------------------")
print(f"Área da seção A = {A*1e6:.3f} mm²")
print(f"Momento de inércia I = {I*1e12:.3f} mm⁴")
print(f"Raio de giração r = {r*1000:.3f} mm")
print(f"Esbeltez λ = {lambda_:.1f}")
print(f"Limite de esbeltez λ_lim = {lambda_lim:.1f}")
print("--------------------------------------------")
print(f"Carga crítica de Euler = {Pcr_Euler:.2f} N")
print(f"Carga crítica de Johnson = {Pcr_Johnson:.2f} N")
print(f"Carga de escoamento (esmagamento) = {P_esmag:.2f} N")
print(f"Carga aplicada = {P_aplicada:.2f} N")
print("--------------------------------------------")
print(f"Modo governante: {modo}")
print(resultado)
