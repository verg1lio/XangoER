import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=200, suppress=True)

class Estrutura:
    def __init__(self, elements, nodes, m, Id, Ip):
        self.elements = elements                                              #Matriz de elementos conectados
        self.num_elements = len(elements)                                     #Número de elementos
        self.nodes = nodes                                                    #Matriz de nós com suas posiççoes
        self.num_nodes = len(nodes)                                           #Número total de nós
        self.massa = m                                                        #Massa do carro (Kg)
        self.momento_inercia_direcao = Id                                     #Momento de inércia em relação à direção (kg.m^2)
        self.momento_inercia_plano = Ip                                       #Momento de inércia em relação ao plano (kg.m^2)
        self.num_dofs_per_node = 6                                            #Graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node               #Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((len(nodes) * 6, len(nodes) * 6))  # Tamanho adequado para a matriz global
        self.M_global = np.zeros((len(nodes) * 6, len(nodes) * 6))
        self.num_modes = 12                                                   #Número de modos de vibração a serem retornados


    def node_loc_matrix(self, node_tags, node_coord):
        # Gerando uma matriz número de nos x 4: (Para cada linha, teremos: [índice do nó, x, y, z])
        node_loc_matrix = np.zeros((self.num_nodes, 4), dtype=float)

        # Preenchendo a matriz de zeros
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i] + 1  # Número do nó na primeira coluna
            node_loc_matrix[i][1] = x  # Coordenada x na segunda coluna
            node_loc_matrix[i][2] = y  # Coordenada y na terceira coluna
            node_loc_matrix[i][3] = z  # Coordenada z na quarta coluna

        print("\n   Nó   x   y   z")
        print(node_loc_matrix)


    def connect_matrix(self):
        # Número de conexões
        num_connect = len(self.elements)

        # Gerando a matriz de conectividade: (Para cada linha, teremos: [índice a coneção, 1º nó, 2º nó])
        CM = np.zeros((num_connect, 3), dtype=int)

        # Preenchendo a matriz:
        for i, (no1, no2) in enumerate(self.elements, start=0):
            CM[i][0] = i + 1
            CM[i][1] = no1 +1
            CM[i][2] = no2 +1

        print("\n Conexão   1º Nó   2º Nó")
        print(CM)


    def calcular_comprimento(self, element):
        node1, node2 = element
        x1, y1, z1 = self.nodes[node1] 
        x2, y2, z2 = self.nodes[node2] 
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)


    def element(self,elemento):
        # Variáveis e constantes físicas do modelo
        E = 210e9   	#Modulo de Young (Pa)
        I = 1.6667e-5 	#Momento de inercia (m^4)
        G = 81.2e9  	#Modulo de Cisalhamento(Pa)
        A= 0.0125	    #Área da seção do elemento (m^2)	
        J = I/2     	#Momento polar de inércia (m^4) 
        kappa=0.9       #Fator de correção para cisalhamento 
        L_e = self.calcular_comprimento(elemento)
        Phi = (12 * E * I) / (kappa * G * A * L_e**2)
        rho = 7850  # kg/m^3

        c1 = E * A / L_e
        c2 = G * J / L_e
        c3 = E * I / L_e**3             #Euler-Bernoulli
        c4 = (E*I)/(L_e**3*(1+Phi))      #Timoshenko
        t1=(4+Phi)
        t2=(2-Phi)
        d1 = rho*A*L_e
        d2 = (I*L_e)/6
        d3 = (rho*A*L_e)/420

        # Matriz de Rigidez Elementar (Euler-Bernoulli)
        # Para converter para timoshenko basta trocar c3 por c4,onde tem (4 * L_e**2 * c3) substitui por (t1* L_e**2 * c3) e onde tiver (2 * L_e**2 * c3) por (t2* L_e**2 * c3) )
        k_e= np.array([
                [12 * c3, 0, 0, 6 * L_e * c3, 0, 0, -12 * c3, 0, 0, 6 * L_e * c3, 0, 0],
                [0, c1, 0, 0, 0, 0, 0, -c1, 0, 0, 0, 0],
                [0, 0, 12 * c3, 0, 0, 6 * L_e* c3, 0, 0, -12 * c3, 0, 0, 6 * L_e * c3],
                [6 * L_e * c3, 0, 0, 4 * L_e**2 * c3, 0, 0, -6 * L_e * c3, 0, 0, 2 * L_e**2 * c3, 0, 0],
                [0, 0, 0, 0, c2, 0, 0, 0, 0, 0, -c2, 0],
                [0, 0, 6 * L_e * c3, 0, 0, 4 * L_e**2 * c3, 0, 0, -6 * L_e * c3, 0, 0, 2 * L_e**2 * c3],
                [-12 * c3, 0, 0, -6 * L_e * c3, 0, 0, 12 * c3, 0, 0, -6 * L_e * c3, 0, 0],
                [0, -c1, 0, 0, 0, 0, 0, c1, 0, 0, 0, 0],
                [0, 0, -12 * c3, 0, 0, -6 * L_e * c3, 0, 0, 12 * c3, 0, 0, -6 * L_e * c3],
                [6 * L_e * c3, 0, 0, 2 * L_e**2 * c3, 0, 0, -6 * L_e * c3, 0, 0, 4 * L_e**2 * c3, 0, 0],
                [0, 0, 0, 0, -c2, 0, 0, 0, 0, 0, c2, 0],
                [0, 0, 6 * L_e * c3, 0, 0, 2 * L_e**2 * c3, 0, 0, -6 * L_e * c3, 0, 0, 4 * L_e**2 * c3]
            ])
        
        # Matriz de Massa Elementar
        m_e= np.array([
                [156 * d3, 0, 0, 22 * L_e * d3, 0, 0, 54 * d3, 0, 0, -13 * L_e * d3, 0, 0],
                [0,2*d1, 0, 0, 0, 0, 0, d1, 0, 0, 0, 0],
                [0, 0, 156 * d3, 0, 0, 22 * L_e* d3, 0, 0, 54 * d3, 0, 0, -13 * L_e * d3],
                [22 * L_e * d3, 0, 0, 4 * L_e**2 * d3, 0, 0, 13 * L_e * d3, 0, 0, -3 * L_e**2 * d3, 0, 0],
                [0, 0, 0, 0, 2*d2, 0, 0, 0, 0, 0, d2, 0],
                [0, 0, 22 * L_e * d3, 0, 0, 4 * L_e**2 * d3, 0, 0, 13 * L_e * d3, 0, 0, -3 * L_e**2 * d3],
                [54 * d3, 0, 0, 13 * L_e * d3, 0, 0, 156* d3, 0, 0, -22 * L_e * d3, 0, 0],
                [0, d1, 0, 0, 0, 0, 0, 2*d1, 0, 0, 0, 0],
                [0, 0, 54 * d3, 0, 0, 13 * L_e * d3, 0, 0, 156 * d3, 0, 0, -22 * L_e * d3],
                [-13 * L_e * d3, 0, 0, -3 * L_e**2 * d3, 0, 0, -22 * L_e * d3, 0, 0, 4 * L_e**2 * d3, 0, 0],
                [0, 0, 0, 0, d2, 0, 0, 0, 0, 0, 2*d2, 0],
                [0, 0, -13 * L_e * d3, 0, 0,-3 * L_e**2 * d3, 0, 0, -22 * L_e * d3, 0, 0, 4 * L_e**2 * d3]
            ])
        return k_e,m_e
    
    def matrizes_global(self):
        #Calculando as matrizes de rigidez e massa de cada elemento
        for element in self.elements:
            node1,node2=element
            k_e,m_e= self.element(element)
            #Montando as matrizes globais
            #           x1      y1          z1        rx1       ry1         rz1       x2        y2          z2          rx2     ry2         rz2        
            dofs = [6*node1, 6*node1+1, 6*node1+2, 6*node1+3, 6*node1+4, 6*node1+5, 6*node2, 6*node2+1, 6*node2+2, 6*node2+3, 6*node2+4, 6*node2+5]
            for i in range(len(dofs)):
                for j in range(len(dofs)):
                    self.K_global[dofs[i], dofs[j]] += k_e[i, j]
                    self.M_global[dofs[i], dofs[j]] += m_e[i, j]
        
        #Aplicando engastes nas extremidades FAZER AS CONDIÇÕES DE CONTORNO NUMA DEF À PARTE
        self.K_global[0, 0] = 10**7
        self.K_global[1, 1] = 10**7                            
        self.K_global[-1, -1] = 10**7
        self.K_global[-2, -2] = 10**7                           
        
        return self.K_global,self.M_global
    
    def shape_fun(self, F_flexao, F_axial, F_torcao):
        E = 210e9   	#Modulo de Young (Pa)
        I = 1.6667e-5 	#Momento de inercia (m^4)
        G = 81.2e9  	#Modulo de Cisalhamento(Pa)
        A= 0.0125	    #Área da seção do elemento (m^2)	
        J = I/2     	#Momento polar de inércia (m^4) 
        torcao, deformacao, flexao = [], [], []
        for element in self.elements:
            L_e = self.calcular_comprimento(element)
            # Equação de torsão
            torcao_val = (F_torcao * L_e) / (G * J)         #Fonte[1]
            torcao.append(torcao_val)
            # Equação  para deformação axial
            deformacao_val = (F_axial* L_e / (A * E))       #Fonte[2]
            deformacao.append(deformacao_val)
            # Equação para flexão
            flexao_val = (F_flexao*L_e**3)/(3 * E * I)      #Fonte[3]
            flexao.append(flexao_val)
        return np.array(torcao), np.array(deformacao), np.array(flexao)
    
    def modal_analysis(self):
       # Análise modal por resolução do problema de autovalor e autovetor
        unsorted_eigenvalues, unsorted_eigenvectors = eigh(self.K_global, self.M_global)

        # Frequências naturais (raiz quadrada dos autovalores)
        unsorted_frequencies = np.sqrt(unsorted_eigenvalues) / (2 * np.pi)  # Divisão por 2*pi para converter para hertz

        # Tratando os dados (tomando apenas as 20 primeiras frequências naturais)
        sorted_indices = np.argsort(unsorted_frequencies)  # Ordena as frequências em ordem crescente
        top_indices = sorted_indices[:self.num_modes]  # Seleciona os índices dos primeiros n modos

        self.eigenvalues = np.array(unsorted_eigenvalues)[top_indices]  # Filtra os primeiros n autovalores
        self.eigenvectors = np.array(unsorted_eigenvectors)[:, top_indices]  # Filtra os primeiros n autovetores
        self.frequencies = np.array(unsorted_frequencies)[top_indices]  # Filtra as primeiras n frequências

        return self.eigenvalues, self.eigenvectors, self.frequencies  

    def structure_plot(self):
        # Plotando o gráfico 3D da estrutura
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Adicionando os pontos
        for i, (x, y, z) in enumerate(self.nodes):
            ax.scatter(x, y, z, color='b', s=100)
            ax.text(x, y, z, f'  {i+1}', color='black', fontsize=12)

        # Adicionando as linhas de ligação entre os nós
        for node1, node2 in elements:
            x = [nodes[node1][0], nodes[node2][0]]
            y = [nodes[node1][1], nodes[node2][1]]
            z = [nodes[node1][2], nodes[node2][2]]
            ax.plot(x, y, z, marker='o')

        # Configurações adicionais
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Estrutura 3D')

        plt.show()        
    
    def modal_analysis_plot(self):
        #Criando plot para os modos de vibração separadamente
        for mode_idx in range(len(self.eigenvalues)):
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            ax.set_title(f'Modo {mode_idx + 1} Viga')

            #Colocando a legenda dos nós no gráfico
            for i, (x, y, z) in enumerate(self.nodes):
                ax.scatter(x, y, z, color='b', s=100)
                ax.text(x, y, z, f'  {i+1}', color='black', fontsize=8)

            #Colocando os nós no gráfico
            for node1, node2 in elements:
                x = [nodes[node1][0], nodes[node2][0]]
                y = [nodes[node1][1], nodes[node2][1]]
                z = [nodes[node1][2], nodes[node2][2]]
                ax.plot(x, y, z, 'b--')

            #Inicializando os modos de vibração e vetor de deslocamentos
            mode_shape = self.eigenvectors[:, mode_idx]
            displacements = np.zeros((len(self.nodes), 3))

            #Aplicando os deslocamentos na estrutura
            for j, (x, y, z) in enumerate(self.nodes):
                if j > 0 and j < len(self.nodes) - 1:
                    displacements[j, 0] = mode_shape[2 * j]
                    displacements[j, 1] = mode_shape[2 * j + 1]
                    displacements[j, 2] = 0

            deformed_nodes = np.array(nodes) + displacements

            #Plotando a estrutura deformada
            for node1, node2 in elements:
                x = [deformed_nodes[node1][0], deformed_nodes[node2][0]]
                y = [deformed_nodes[node1][1], deformed_nodes[node2][1]]
                z = [deformed_nodes[node1][2], deformed_nodes[node2][2]]
                ax.plot(x, y, z, 'r-')

            #Configurações adicionais
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')

            plt.tight_layout()
            plt.show()
    
    def shape_fun_plot(self, F_flexao, F_axial, F_torcao, flexao, axial, torcao):
        # Plotando os resultados das deformações
        fig, axs = plt.subplots(3, 1, figsize=(12, 18))

        # Plot da Torção
        axs[0].plot(torcao, 'o-', label=[f'Força {F}N' for F in F_torcao])
        axs[0].set_title('Deformação por Torção de cada Elemento')
        axs[0].set_xlabel('Elemento')
        axs[0].set_ylabel('Torção (rad)')
        axs[0].legend()

        # Plot da Deformação Axial
        axs[1].plot(axial, 's-', label=[f'Força {F}N' for F in F_axial])
        axs[1].set_title('Deformação Axial de cada Elemento')
        axs[1].set_xlabel('Elemento')
        axs[1].set_ylabel('Deformação (m)')
        axs[1].legend()

        # Plot da Flexão
        axs[2].plot(flexao,'o-', label=[f'Força {F}N' for F in F_flexao])
        axs[2].set_title('Deformação por Flexão de cada Elemento')
        axs[2].set_xlabel('Posição ao longo do elemento')
        axs[2].set_ylabel('Flexão(m)')
        axs[2].legend()

        plt.tight_layout()
        plt.show()

#Coordenadas dos nós (x, y, z)
nodes = [
    (0, 0, 0),
    (0, 0.375, 0),
    (0, 0.700, 0),
    (1.500, 0.375, 0),
    (1.500, 0, 0),
    (1.500, 0.700, 0)
]

#Conectividade dos elementos (índices dos nós)
elements = [
    (0, 1),  # Elemento entre os nós 0 e 1
    (1, 2),  # Elemento entre os nós 1 e 2
    (4, 3),  # Elemento entre os nós 4 e 3
    (3, 5),  # Elemento entre os nós 3 e 5
    (1, 3)  # Elemento entre os nós 1 e 3
]

#Criar a estrutura e montar as matrizes de rigidez e massa globais, atribuir forças
F_flexao = np.array([1000, 2000, 3000, 4000, 5000])
F_axial = np.array([1000, 2000, 3000, 4000, 5000])
F_torcao = np.array([1000, 2000, 3000, 4000, 5000])
estrutura = Estrutura(elements, nodes, 1500, 8.33e-6, 8.33e-6)
K_global, M_global = estrutura.matrizes_global()

#Gerar as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
estrutura.node_loc_matrix(node_tags, nodes)
estrutura.connect_matrix()

print("\n Matriz de Rigidez Global")
print(K_global)

print("\n Matriz de Massa Global")
print(M_global)

#Gerar autovalores, autovetores e frequências naturais
autovalores, autovetores, frequencias = estrutura.modal_analysis()

#Gerar resultados de torção, axial e flexão
torcao, axial, flexao = estrutura.shape_fun(F_flexao, F_axial, F_torcao)

#Plotando os resultados das deformações
estrutura.shape_fun_plot(F_flexao, F_axial, F_torcao, flexao, axial, torcao)

# Plotando o gráfico 3D da estrutura
estrutura.structure_plot()

#Exibindo as frequências naturais e modos de vibração da estrutura
print("\n Frequências Naturais (ω) da estrutura:")
print(frequencias)

#Plotagem dos modos de vibração para a estrutura de vigas
estrutura.modal_analysis_plot()