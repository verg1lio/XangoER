import math

def CalculateParams(RedP, a, Mt, psi, μl, HCG, μ, Rdp, Npast, red, Dwc, Dcm, L, pi, atrito_coeficiente, FzF, FzR):
    BF = 2 * μl
    χ = HCG / L
    W = Mt * 9.81 
    FzF_dyn = (1 - psi + a * χ) * W
    FzR_dyn = (psi - a * χ) * W
    τF = FzF_dyn * μ * Rdp
    τR = FzR_dyn * μ * Rdp
    FnF = τF / Npast * RedP * red
    FnR = τR / Npast * RedP * red
    Awc = (pi * (Dwc ** 2)) / 4
    Acm = (pi * (Dcm ** 2)) / 4
    PF = FnF / Awc
    PR = FnR / Awc
    FaCM = PF * Acm
    lF = psi * L
    lR = (1 - psi) * L
    return BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR

# Definindo os parâmetros
params = {
    'RedP': 4,  # Redução do pedal 
    'a': 0.8, # Desaceleração [g]
    'psi': 0.40, # Distribuição de carga estática por eixo [adm]
    'μl': 0.45, # Coeficiente de atrito do contato pastilha/disco
    'pi': 3.14,
    'HCG': 0.5, # Altura do centro de gravidade [m]
    'μ': 0.60, # Coeficiente de atrito do pneu/solo
    'FzF': 1471.5, # Força de reação estática na dianteira [N]
    'FzR': 981.0, # Força de reação estática na traseira [N]
    'Rdp': 0.30, # Raio do pneu [m]
    'Dcm': 0.02,  # Diâmetro do cilindro mestre em metros
    'Dwc': 0.032,  # Diâmetro do cilindro da roda em metros
    'Npast': 2,  # Número de pastilhas por disco
    'atrito_coeficiente': 0.35,  # Coeficiente de atrito da pastilha
    'red': 0.75,  # Raio efetivo do disco [m]
    'Mt': 250,  # Massa total do veículo [kg]
    'L': 1.5,  # Distância entre eixos [m]
}

# Chamando a função e atribuindo os resultados
BF, χ, W, FzF_dyn, FzR_dyn, τF, τR, FnF, FnR, Awc, Acm, PF, PR, FaCM, lF, lR = CalculateParams(**params)

# Imprimindo os resultados
print("BF:", BF)
print("χ:", χ)
print("W:", W, "N")
print("FzF_dyn:", FzF_dyn, "N")
print("FzR_dyn:", FzR_dyn, "N")
print("τF:", τF, "Nm")
print("τR:", τR, "Nm")
print("FnF:", FnF, "N")
print("FnR:", FnR, "N")
print("Awc:", Awc, "mm²")
print("Acm:", Acm, "mm²")
print("PF:", PF, "MPa")
print("PR:", PR, "Mpa")
print("FaCM:", FaCM, "N")
print("lF:", lF, "m")
print("lR:", lR, "m")


    
