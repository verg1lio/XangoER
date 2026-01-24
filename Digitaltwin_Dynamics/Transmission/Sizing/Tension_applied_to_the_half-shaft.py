import math

def calcular_momento_inercia_polar(diametro):
    raio = diametro / 2
    return (math.pi * (raio**4)) / 2

def calcular_momento_inercia_area(diametro):
    return (math.pi * (diametro**4)) / 64

def calcular_tensao_torcao(torque, diametro, kt_torcao=1.0):
    raio = diametro / 2
    j = calcular_momento_inercia_polar(diametro)
    if j == 0:
        return float('inf')  # Evita divisão por zero
    tensao_torcao = (torque * raio) / j
    return tensao_torcao * kt_torcao

def calcular_tensao_flexao(momento_fletor, diametro, kt_flexao=1.0):
    c = diametro / 2
    i = calcular_momento_inercia_area(diametro)
    if i == 0:
        return float('inf')  # Evita divisão por zero
    tensao_flexao = (momento_fletor * c) / i
    return tensao_flexao * kt_flexao

def calcular_tensao_von_mises(tensao_normal, tensao_cisalhamento):
    return math.sqrt(tensao_normal**2 + 3 * tensao_cisalhamento**2)

def analisar_semieixo_fsae(secoes_semieixo, torque_max, momento_fletor_max, material_propriedades):
    resultados_analise = {}

    print("Analise de Tensoes do Semieixo")
    print(f"Torque Max: {torque_max:.2f} Nm")
    print(f"Momento Fletor m_ap: {momento_fletor_max:.2f} Nm")
    print(f"Limite de Escoamento: {material_propriedades['limite_escoamento'] / 1e6:.2f} MPa")
    print(f"Limite de Resistencia: {material_propriedades['limite_resistencia'] / 1e6:.2f} MPa\n")

    for i, secao in enumerate(secoes_semieixo):
        diametro = secao['diametro'] / 1000  # Converter mm para metros
        kt_flexao = secao.get('kt_flexao', 1.0)  # Padrão 1.0 se não especificado
        kt_torcao = secao.get('kt_torcao', 1.0)  # Padrão 1.0 se não especificado

        print(f"Secao {i+1} (Diametro: {secao['diametro']} mm)")

        # Calcular tensões
        tensao_torcao = calcular_tensao_torcao(torque_max, diametro, kt_torcao)
        tensao_flexao = calcular_tensao_flexao(momento_fletor_max, diametro, kt_flexao)
        tensao_von_mises = calcular_tensao_von_mises(tensao_flexao, tensao_torcao)

        # Converter para MPa para melhor visualização
        tensao_torcao_mpa = tensao_torcao / 1e6
        tensao_flexao_mpa = tensao_flexao / 1e6
        tensao_von_mises_mpa = tensao_von_mises / 1e6

        # Calcular fator de segurança
        fator_seguranca_escoamento = material_propriedades['limite_escoamento'] / tensao_von_mises if tensao_von_mises > 0 else float('inf')
        fator_seguranca_resistencia = material_propriedades['limite_resistencia'] / tensao_von_mises if tensao_von_mises > 0 else float('inf')

        print(f"  T tau: {tensao_torcao_mpa:.2f} MPa (Kt={kt_torcao})")
        print(f"  T sigma: {tensao_flexao_mpa:.2f} MPa (Kt={kt_flexao})")
        print(f"  T equivalente: {tensao_von_mises_mpa:.2f} MPa")
        print(f"  FS Escoamento: {fator_seguranca_escoamento:.2f}")
        print(f"  FS Resistencia: {fator_seguranca_resistencia:.2f}")

        # Armazenar resultados
        resultados_analise[f"secao_{i+1}"] = {
            'diametro_mm': secao['diametro'],
            'tensao_torcao_mpa': tensao_torcao_mpa,
            'tensao_flexao_mpa': tensao_flexao_mpa,
            'tensao_von_mises_mpa': tensao_von_mises_mpa,
            'fator_seguranca_escoamento': fator_seguranca_escoamento,
            'fator_seguranca_resistencia': fator_seguranca_resistencia
        }
        print("-" * 40)

    return resultados_analise

# Teste
if __name__ == "__main__":
    material_propriedades = {  # Aço 4340
        'limite_escoamento': 1000e6,  # 1000 MPa (aproximado)
        'limite_resistencia': 1200e6  # 1200 MPa (aproximado)
    }

    torque_max_semieixo = 500  # Nm 
    momento_fletor_max_semieixo = 300  # Nm 

    secoes_semieixo = [
        {'diametro': 30, 'kt_flexao': 1.2, 'kt_torcao': 1.1},  # Seção maior próxima ao diferencial
        {'diametro': 20, 'kt_flexao': 2.0, 'kt_torcao': 1.8},  # Seção menor próxima à roda
    ]

    # Análise
    resultados = analisar_semieixo_fsae(secoes_semieixo, torque_max_semieixo, momento_fletor_max_semieixo, material_propriedades)

    print("\n--- Resumo dos Resultados ---")
    for secao_nome, dados in resultados.items():
        print(f"{secao_nome.replace('_', ' ').capitalize()}:")
        print(f"  Diametro: {dados['diametro_mm']} mm")
        print(f"  Tensao: {dados['tensao_von_mises_mpa']:.2f} MPa")
        print(f"  FS Escoamento: {dados['fator_seguranca_escoamento']:.2f}")
        print(f"  FS Resistencia: {dados['fator_seguranca_resistencia']:.2f}")
        print("-" * 20)
