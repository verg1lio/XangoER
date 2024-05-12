opcao = int(input("Matrizes rigidez flexional (1), massa (2), axial (3) ou torcional (4) :"))
L = float(input("Digite o tamanho do elemento L: "))
num_linhas = int(input("Quantas linhas tem a matriz? "))
num_colunas = int(input("Quantas colunas tem a matriz? "))

def construir_matriz():
    matriz = []

    # Loop para solicitar os elementos da matriz
    for i in range(num_linhas):
        linha = []
        for j in range(num_colunas):
            elemento = input(f"Qual é o elemento da linha {i+1} e coluna {j+1}? ")
            linha.append(calculadora_rebuscada(elemento))
        matriz.append(linha)

    return matriz

def calculadora_rebuscada(expressao):
    try:
        resultado = eval(expressao)
        # L=2
        # 2*L**2
        # 8
        return resultado
    except Exception as e:
        return "Erro: " + str(e)
    
# Construção da matriz

matriz = construir_matriz()

if opcao == 1:
    modulo_de_young = float(input("Digite o módulo de Young do material:"))
    momento_de_inercia = float(input("Digite o momento de inércia do elemento:"))
    Kf = 2*modulo_de_young*momento_de_inercia/(L**3)
    print(f"A constante fora da matriz é igual a {Kf}")
    C=Kf
    for i in range(num_linhas):
        for j in range(num_colunas):
            matriz[i][j]=matriz[i][j]*C
elif opcao == 2:
    massa_especifica = float(input("Digite a massa específica do material:"))
    area_do_elemento = float(input("Digite a área do elemento:"))
    Km = massa_especifica*area_do_elemento*L/420
    print(f"A constante fora da matriz é igual a {Km}")
    C=Km
    for i in range(num_linhas):
        for j in range(num_colunas):
            matriz[i][j]=matriz[i][j]*C
elif opcao == 3:
    modulo_de_young = float(input("Digite o módulo de Young do material:"))
    area_do_elemento = float(input("Digite a área do elemento:"))
    Ka = modulo_de_young*area_do_elemento/(3*L)
    print(f"A constante fora da matriz é igual a {Ka}")
    C=Ka
    for i in range(num_linhas):
        for j in range(num_colunas):
            matriz[i][j]=matriz[i][j]*C
elif opcao == 4:
    constante_torcao = float(input("Digite a constante de torcao do material:"))
    elasticidade_transversal = float(input("Digite o módulo de elasticidade transversal do material:"))
    momento_polar = float(input("Digite o momento polar de inércia do elemento:"))
    Kt = constante_torcao*elasticidade_transversal*momento_polar/L
    print(f"A constante fora da matriz é igual a {Kt}")
    C=Kt
    for i in range(num_linhas):
        for j in range(num_colunas):
            matriz[i][j]=matriz[i][j]*C
else:
    print("Erro")
print(matriz)