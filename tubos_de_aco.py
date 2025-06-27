import math

def calcular_propriedades_tubo(d_ext, espessura):
    """
    Calcula propriedades mecânicas de um tubo de aço 1020.

    Parâmetros:
    - d_ext (float): Diâmetro externo do tubo em metros.
    - espessura (float): Espessura da parede do tubo em metros.

    Retorna:
    - dict com propriedades calculadas 
    modo de usar: digite no terminal  "python tubosaco.py"  e depois coloque o diâmetro e a espessura desejada 
    """

    # Propriedades do aço 1020
    E = 2.05e11  # Pa
    nu = 0.29    # coeficiente de Poisson
    G = E / (2 * (1 + nu))  # Módulo de cisalhamento (Pa)

    # Geometria
    d_int = d_ext - 2 * espessura
    if d_int <= 0:
        raise ValueError("Espessura grande demais: o diâmetro interno é zero ou negativo.")

    # Área da seção transversal
    A = (math.pi / 4) * (d_ext**2 - d_int**2)

    # Momento de inércia (Ix)
    Ix = (math.pi / 64) * (d_ext**4 - d_int**4)

    # Momento polar de inércia (J)
    J = (math.pi / 32) * (d_ext**4 - d_int**4)

    return {
        "Módulo de Young (Pa)": E,
        "Módulo de Cisalhamento (Pa)": G,
        "Área da seção (m²)": A,
        "Momento de Inércia (m⁴)": Ix,
        "Momento Polar de Inércia (m⁴)": J
    }

def main():
    print("=== Cálculo de propriedades de tubo de aço 1020 ===")
    try:
        d_ext = float(input("Digite o diâmetro externo do tubo (em metros): "))
        espessura = float(input("Digite a espessura da parede do tubo (em metros): "))

        propriedades = calcular_propriedades_tubo(d_ext, espessura)

        print("\n--- Resultados ---")
        for nome, valor in propriedades.items():
            print(f"{nome}: {valor:.4e}")

    except ValueError as e:
        print(f"Erro: {e}")
    except Exception as e:
        print(f"Erro inesperado: {e}")

if __name__ == "__main__":
    main()
