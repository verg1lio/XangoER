import math

def titulo(texto):
    print(f"\n{'='*60}\n {texto}\n{'='*60}")

def subtitulo(texto):
    print(f"\n--- {texto} ---")

class CalculadoraSuspensaoMetrica:
    def __init__(self):
        # Fator de conversão APENAS para consulta de catálogo
        self.NMM_TO_LBS = 5.71 

    def calcular_rigidez(self, massa_kg, freq_hz, k_pneu_nmm):
        # 1. Rigidez de Rodagem (K_ride)
        # K = M * w^2
        k_ride = massa_kg * (2 * math.pi * freq_hz) ** 2 / 1000.0 # /1000 para N/mm
        
        # 2. Rigidez na Roda (K_roda) - Correção do Pneu em série
        k_roda = (k_ride * k_pneu_nmm) / (k_pneu_nmm - k_ride)
        
        return k_ride, k_roda

    def selecionar_mola_comercial(self, k_roda_tras_nmm, mr_tras_padrao):
        # Calcula a mola teórica necessária na traseira (N/mm)
        k_mola_teorica = k_roda_tras_nmm / (mr_tras_padrao ** 2)
        
        # Converte para libras APENAS para achar a peça no catálogo
        lbs_teorica = k_mola_teorica * self.NMM_TO_LBS
        
        # Arredonda para o valor comercial comum mais próximo (ex: 300, 325, 350)
        # Vamos assumir incrementos de 25 lbs
        lbs_comercial = round(lbs_teorica / 25) * 25
        
        # Converte de volta para N/mm para usar nas contas
        k_mola_real_nmm = lbs_comercial / self.NMM_TO_LBS
        
        print(f"   > Mola Teórica Traseira: {k_mola_teorica:.2f} N/mm ({lbs_teorica:.0f} lb/in)")
        print(f"   > Mola Comercial Escolhida: {lbs_comercial} lb/in ({k_mola_real_nmm:.2f} N/mm)")
        print(f"   (Usaremos esta mola de {k_mola_real_nmm:.2f} N/mm nas 4 rodas)")
        
        return k_mola_real_nmm

    def calcular_proporcao_balancim(self, k_roda_alvo, k_mola_real, ir_real):
        # 1. Descobre o Motion Ratio (MR) que precisamos na frente para essa mola funcionar
        # MR = raiz(K_roda / K_mola)
        mr_necessario = math.sqrt(k_roda_alvo / k_mola_real)
        
        # 2. Calcula a Razão do Balancim (Rocker Ratio - RR)
        # RR = MR / IR
        rocker_ratio = mr_necessario / ir_real
        
        return mr_necessario, rocker_ratio

# ============================================================================
# DADOS DE ENTRADA (MÉTRICOS)
# ============================================================================

# Massas e Frequências
M_FRENTE = 67.5  # kg
M_TRAS = 82.5    # kg
FREQ_FRENTE = 2.5 # Hz
FREQ_TRAS = 2.8   # Hz

# Geometria e Pneu
K_PNEU = 120.0       # N/mm
MR_TRAS_PADRAO = 0.80 # MR Traseiro (Design eficiente)
IR_REAL = 0.86        # Installation Ratio (Pushrod a 280mm / Braço 325mm)

# Instancia a calculadora
calc = CalculadoraSuspensaoMetrica()

titulo("DEFINIÇÃO DA MOLA PADRÃO (BASEADA NA TRASEIRA)")
# Passo 1: Calcular o que a traseira precisa
_, k_roda_tras = calc.calcular_rigidez(M_TRAS, FREQ_TRAS, K_PNEU)

# Passo 2: Escolher a mola que atende a traseira
k_mola_padrao = calc.selecionar_mola_comercial(k_roda_tras, MR_TRAS_PADRAO)


titulo("PROJETO DOS BALANCINS (PROPORÇÕES)")

# --- ANÁLISE TRASEIRA ---
subtitulo("1. Balancim Traseiro (Mola Correta)")
# O MR é o padrão 0.80
rr_tras = MR_TRAS_PADRAO / IR_REAL

print(f"Objetivo: Manter MR de {MR_TRAS_PADRAO}")
print(f"PROPORÇÃO DO BALANCIM (Saída : Entrada):")
print(f"   1 : {1/rr_tras:.2f}")
print(f"   (Ou seja: O braço do Amortecedor deve ser {rr_tras*100:.1f}% do tamanho do braço do Pushrod)")


# --- ANÁLISE DIANTEIRA ---
subtitulo("2. Balancim Dianteiro (Adaptação para Mola Dura)")
# Passo 1: O que a frente precisa?
_, k_roda_frente = calc.calcular_rigidez(M_FRENTE, FREQ_FRENTE, K_PNEU)

# Passo 2: Qual a proporção do balancim para transformar a mola dura nisso?
mr_frente, rr_frente = calc.calcular_proporcao_balancim(k_roda_frente, k_mola_padrao, IR_REAL)

print(f"Objetivo: Atingir rigidez de roda de {k_roda_frente:.2f} N/mm usando a mola dura.")
print(f"Motion Ratio Ajustado Necessário: {mr_frente:.3f}")
print(f"PROPORÇÃO DO BALANCIM (Saída : Entrada):")
print(f"   1 : {1/rr_frente:.2f}")
print(f"   (Ou seja: O braço do Amortecedor deve ser apenas {rr_frente*100:.1f}% do tamanho do braço do Pushrod)")

titulo("RESUMO DAS PROPORÇÕES PARA O CAD")
print(f"TRASEIRA: Proporção {rr_tras:.2f} (Quase 1 para 1)")
print(f"   -> Exemplo visual: Se a entrada tiver 100mm, a saída tem {100*rr_tras:.1f}mm")
print(f"\nDIANTEIRA: Proporção {rr_frente:.2f} (Desmultiplicado)")
print(f"   -> Exemplo visual: Se a entrada tiver 100mm, a saída tem {100*rr_frente:.1f}mm")
