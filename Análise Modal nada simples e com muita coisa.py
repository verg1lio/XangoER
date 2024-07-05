import numpy as np
from scipy.linalg import eigh
import gmsh
import sys

#Criando uma matriz de massa de massa de um elemento de barra apenas
def bar_elem_mass_matrix(n_elem, ρ, A, L):
    #Matriz de massa elementar
    m = (ρ*A*L/6*n_elem) * np.array([[2,1], [1,2]])

    return m            #Retornando apenas a matriz elementar pois é isso que é necessário agora

#Parâmetros da barra
E = 210e9               #Módulo de Elasticidade lontigudinal, em Pa
A = 0.225               #Área da seção transversal, em m²
ρ = 7850                #Massa específica do aço, em kg/m³

#MALHA DE ELEMENTO QUALQUER DE TESTE + MATRIZES DE LOCALIZAÇÃO E CONECTIVIDADE
#Parâmetros iniciais
n_div = 2

L = 10                      #Comprimento total da barra, m
H = 1                       #Altura total da barra, m
W = 2                       #Espessura total da barra, m

l_elem = L/n_div            #Comprimento do elemento
h_elem = H/n_div            #Altura do elemento
w_elem = W/n_div            #Espessura do elemento

#CRIANDO UMA MALHA QUALQUER

def generate_mesh():
    #Inicializando o gmsh
    gmsh.initialize()

    #Criando o modelo
    gmsh.model.add("isso aqui é só pra teste por enquanto, deve sair uma barra 3D agora")

    #Definindo os pontos do elemento e linhas do elemento
    pontos = [(0, 0, 0),
              (L, 0, 0),
              (L, H, 0),
              (0, H, 0),
              (0, 0, W),
              (L, 0, W),
              (L, H, W),
              (0, H, W)]

    point_tags = []
    for x, y, z in pontos:
        point_tags.append(gmsh.model.geo.addPoint(x, y, z))
    
    linhas = [(point_tags[0], point_tags[1]),
              (point_tags[1], point_tags[2]),
              (point_tags[2], point_tags[3]),
              (point_tags[3], point_tags[0]),
              (point_tags[4], point_tags[5]),
              (point_tags[5], point_tags[6]),
              (point_tags[6], point_tags[7]),
              (point_tags[7], point_tags[4]),
              (point_tags[0], point_tags[4]),
              (point_tags[1], point_tags[5]),
              (point_tags[2], point_tags[6]),
              (point_tags[3], point_tags[7])]

    for p1, p2 in linhas:
        gmsh.model.geo.addLine(p1, p2)
    #Ideia para melhorar essa parte: montar numa planilha de excel os pontos dos nós do chassi com coordenadas e suas conectividades
    #e importar esses dados e montar o código definindo os pontos e linhas acima a partir dos dados dessa planilha

    #Sincronizando a gemoetria
    gmsh.model.geo.synchronize()

    #Comprimento característico
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", l_elem)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", l_elem)

    #Gerar a malha
    gmsh.model.mesh.generate(1)         #Elementos de 1 dimensão (barras simples)

    #Obtendo informações da malha
    node_tags, node_coords, node_parametric_coords = gmsh.model.mesh.getNodes()     #Índices dos nós, coordenadas e coordenadas paramétricas
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()           #Tipos de elemento, índicies dos elementos e conecctividades

    #Tratando os dados para extrair o necessário
    node_coords = node_coords.reshape((-1,3))           #Mostra as coordenadas x, y e z de cada nó ordenadamente
    elements = elem_node_tags[0].reshape((-1, 2))       #Conectividade dos elementos de barra 1D
    
    #Visualizar a malha
    gmsh.fltk.run()

    #Encerrar isso daqui
    gmsh.finalize()

    #Coisas para retornar
    return node_tags, node_coords, elements

#MATRIZ DE LOCALIZAÇÃO DOS NÓS
def node_loc_matrix(node_tags, node_coord):
    #Número de Nós
    num_nodes = len(node_tags)

    #Gerando uma matrix número de nos x 4:                          (Para cada linha, teremos: [índice do nó, x, y, z])
    node_loc_matrix = np.zeros((num_nodes,4), dtype = float)        #Na primeira coluna, teremos os índices, nas seguintes, x y e z.
    
    #Preenchendo a matriz de zeros
    for i, (x, y, z) in enumerate(node_coord, start=0):
        node_loc_matrix[i][0] = node_tags[i]                        #Número do nó na primeira coluna     
        node_loc_matrix[i][1] = x                                   #Coordenada x na segunda coluna
        node_loc_matrix[i][2] = y                                   #Coordenada y na terceira coluna
        node_loc_matrix[i][3] = z                                   #Coordenada z na quarta coluna
    
    print("   Nó   x   y   z")
    print(node_loc_matrix)

#MATRIZ DE CONECTIVIDADE DOS NÓS
def connect_matrix(elements):
    #Número de conexões
    num_connect = len(elements)

    #Gerando a matriz de conectividade:                             (Para cada linha, teremos: [índice a coneção, 1º nó, 2º nó])
    CM = np.zeros((num_connect, 3), dtype = int)

    #Preenchendo a matriz:
    for i, (no1, no2) in enumerate(elements, start=0):
        CM[i][0] = i+1
        CM[i][1] = no1
        CM[i][2] = no2
    
    print("Conexão   1º Nó   2º Nó")
    print(CM)

# Matrizes de rigidez qualquer para exemplo, assumir no momento de forma extremamente errada que todas as matrizes de rigidez eelementar são iguais
K = np.array([[12, -4], 
              [-4, 12]], dtype=float)

#Função de análise modal elemento a elemento
def modal_analysis(K, node_coord, elements, num_modes = 20):
    #Número de elementos a serem analizados
    num_elements = len(elements)

    #Inicializando vetores de autovalores, autovetores e frequências
    unsorted_eigenvalues = []
    unsorted_eigenvectors = []
    unsorted_frequencies = []

    #Inicializando o loop de análise dos elementos
    for i in range(num_elements):
        no_a , no_b = elements[i] - 1               #Estabelecendo nós a serem utilizados

        xa, ya, za = node_coord[no_a]               #Coordenadas do 1º nó que forma o elemento
        xb, yb, zb = node_coord[no_b]               #Coordenadas do 2º nó que forma o elemento

        #Calculando o comprimento do elemento (distância entre os pontos)
        L = np.sqrt((xb - xa)**2 + (yb - ya)**2 + (zb - za)**2)

        m = bar_elem_mass_matrix(1, ρ, A, L)        #Número de elementos = 1 elemento por barra

        # Resolver o problema de autovalores e autovetores generalizado
        temp_eigenvalues, temp_eigenvectors = eigh(K, m)

        # Frequências naturais (raiz quadrada dos autovalores)
        temp_frequencies = np.sqrt(temp_eigenvalues) / (2*np.pi)            #Divisão por 2*pi para converter para hertz

        #Adicionando os resultados de cada análise no vetor original
        unsorted_eigenvalues.extend(temp_eigenvalues)
        unsorted_eigenvectors.extend(temp_eigenvectors)                       
        unsorted_frequencies.extend(temp_frequencies)

    #Tratando os dados (tomando apenas as 20 primeiras frequências naturais)
    sorted_indices = np.argsort(unsorted_frequencies)                           #Ordena as frequências em ordem crescente
    top_indices = sorted_indices[:num_modes]                                    #Seleciona os índices dos primeiros n modos

    eigenvalues = np.array(unsorted_eigenvalues)[top_indices]                   #Filtra os primeiros n autovalores
    eigenvectors = np.array(unsorted_eigenvectors)[top_indices]                 #Filtra os primeiros n autovetores
    frequencies = np.array(unsorted_frequencies)[top_indices]                   #Filtra as primeiras n frequências

    return eigenvalues, eigenvectors, frequencies

#Parâmetros da Malha
node_tags, node_coord, elements = generate_mesh()

#Análise modal dos elementos da malha
autovalores, autovetores, frequências = modal_analysis(K, node_coord, elements)

# Exibir resultados
print("Matriz de localização dos nós:")
node_loc_matrix(node_tags, node_coord)
print("")

print("Matriz de conectividade dos elementos")
connect_matrix(elements)
print("")

#print("Autovalores (λ):")
#print(autovalores)
#print("")

print("Frequências Naturais (ω):")
print(frequências)
print("")

print("Autovetores (Modos de Vibração):")
print(autovetores)