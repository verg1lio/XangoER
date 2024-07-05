#Inicializando as bibliotecas
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import sys

#Parâmetros iniciais
n_div = 5

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


#Fazendo a coisa funcionar
node_tags, node_coord, elements = generate_mesh()

print("Matriz de localização dos nós:")
node_loc_matrix(node_tags, node_coord)
print("")

print("Matriz de conectividade dos elementos")
connect_matrix(elements)