#find_new_index + class Estrutura

import numpy as np
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback
from scipy.spatial.distance import pdist
#from rascunho_KT_new import Estrutura #Substituir pelo nome correto da main
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import pandas as pd
import os
from datetime import datetime
import time

class Estrutura:

    def __init__(self, elements, nodes):
        """
        Initializes the structure with elements, nodes, and physical properties.
        
        Inputs:
            - elements: connectivity matrix between nodes (tuples of node indices).
            - nodes: node coordinates (Nx3 array, where N is the number of nodes).
            - m: total mass of the system (float).
            - Id: directional moment of inertia (float).
            - Ip: planar moment of inertia (float).

        Outputs: None.

        Code authors: 
            - Patrícia Nascimento Vaccarezza; 
            - Eduardo Almeida Menezes; 
            - Cayque Lemos Souza; 
            - Antônio Marcos Lopes Brito Junior; 
            - Larissa Pereira Leanor;
            - Alexandre Duque;
            - Vergílio Torezan Silingardi Del Claro.
        """
        self.elements = elements                                             #Matriz de elementos conectados
        self.num_elements = len(elements)                                    #Número de elementos
        self.nodes = nodes                                                   #Matriz de nós com suas posições
        self.num_nodes = len(nodes)                                          #Número total de nós
        self.num_dofs_per_node = 6                                           #6 graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node              #Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))             #Matriz de rigidez global
        self.M_global = np.zeros((self.num_dofs, self.num_dofs))             #Matriz de massa global
        self.num_modes = 1                                                  #Número de modos de vibração a serem retornados
        self.car_mass = 0

    def calcular_comprimento(self, element):    
        """
        Calculates the length of an element based on node coordinates.
        Inputs:
            - element: tuple (start node index, end node index).
        Outputs:
            - element length (float).
        """                    
        node1, node2, tube_type = element
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    def momento_inercia_area_e_polar(self, diametro, espessura):                            #Função para calculo dos momentos de inércia de área e polar
            outer_radius = (diametro / 2)
            inner_radius = (outer_radius - espessura)
            I = (np.pi * 0.25) * (outer_radius ** 4 - inner_radius ** 4)
            J = (np.pi * 0.5) * (outer_radius ** 4 - inner_radius ** 4)
            return I, J

    def area_seccao_transversal(self, diameter, espessura):                                 #Função para calcular a área da secção transversal do tubo (diâmetro externo)
        outer_radius = diameter / 2
        inner_radius = (outer_radius - espessura) 
        A = (outer_radius ** 2 - inner_radius ** 2) * np.pi
        return A 
    
    def mass(self):
        """
        Calculate the mass of the entire structure.

        Parameters:
        None (uses class attributes):


        outputs:
            - self.car_mass (float): Mass of the entire structure
        """
        for element in self.elements:
            L_e = self.calcular_comprimento(element)
            d = self.obter_propriedades(element[2])[2]         #Diâmetro do Tubo (m)
            e = self.obter_propriedades(element[2])[3]         #Espessura do Tubo (m)
            A = self.area_seccao_transversal(d,e)                #Área da secção transversal (m^2)
            rho = self.obter_propriedades(element[2])[4]       #Densidade do material (kg/m^3)
            raio_externo = d / 2
            raio_interno = raio_externo - e
            volume = np.pi*L_e* (raio_externo**2 - raio_interno**2)
            self.car_mass+= volume*rho

        return self.car_mass
    
    def obter_propriedades(self, tube_nome):                                                #Função que lê a planilha 'tubos.csv' e extrai as propriedades de cada tipo de tubo lá presentes
        df = pd.read_csv('tubos.csv')  
        tubo = df[df['Tube'] == tube_nome]
    
        if tubo.empty:
            raise ValueError(f"Tubo '{tube_nome}' não encontrado.")
    
       
        propriedades = tubo.iloc[0]
        E = propriedades['E']
        G = propriedades['G']
        diametro = propriedades['Diametro(m)']
        espessura = propriedades['Espessura(m)']
        densidade = propriedades['Densidade(kg/m^3)']
        poisson = propriedades['Poisson']

        return float(E), float(G), float(diametro), float(espessura), float(densidade), float(poisson)

    def node_loc_matrix(self, node_tags, node_coord): 
        """
        Creates a matrix with node locations for visualization.
        Inputs:
            - node_tags: list of node identifiers.
            - node_coord: matrix of node coordinates.
        Outputs: None."""

        num_nodes = len(node_tags)
        node_loc_matrix = np.zeros((num_nodes, 4), dtype=float)
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i] + 1
            node_loc_matrix[i][1] = x
            node_loc_matrix[i][2] = y
            node_loc_matrix[i][3] = z
        
        # print("\n   Nó   x   y   z")
        # print(node_loc_matrix)

    def connect_matrix(self):
        """
        Generates and prints the connectivity matrix of elements.
        Inputs: None (uses class attributes).
        Outputs: None.
        """
        # Inicializar uma lista para armazenar as conexões
        connections = []

        # Criando a lista a partir de Connections para monstar a matriz connect
        for i, element in enumerate(self.elements):
            node_start, node_end, tube_type = element
            connections.append([i + 1, node_start, node_end])
            
        # Converter a lista em um array numpy
        connections_matrix = np.array(connections)

        # print("Matriz de conectividade:")
        # print(connections_matrix)

    def element(self, element):
        """
        Computes the element stiffness and mass matrices.
        Inputs:
            - element: tuple (start node index, end node index).
        Outputs:
            - k_e: element stiffness matrix.
            - m_e: element mass matrix.
        """
        # Variáveis e constantes físicas do modelo
        d = self.obter_propriedades(element[2])[2]         #Diâmetro do Tubo (m)
        e = self.obter_propriedades(element[2])[3]         #Espessura do Tubo (m)
        E = self.obter_propriedades(element[2])[0]         #Modulo de Young (Pa)      
        G = self.obter_propriedades(element[2])[1]         #Modulo de Cisalhamento (Pa)
        A = self.area_seccao_transversal(d, e)             #Área da seção do elemento (m^2)
        I = self.momento_inercia_area_e_polar(d,e)[0]      #Momento de inercia (m^4)
        J = self.momento_inercia_area_e_polar(d,e)[1]      #Momento polar de inércia (m^4)
        rho = self.obter_propriedades(element[2])[4]       #Densidade do material (kg/m^3)
        kappa=0.9                                          #Fator de correção para cisalhamento 
        L_e = self.calcular_comprimento(element)
        Phi = (12 * E * I) / (kappa * G * A * L_e**2)

        c1 = E * A / L_e
        c2 = G * J / L_e
        c3 = E * I / L_e**3             #Euler-Bernoulli
        c4 = (E*I)/(L_e**3*(1+Phi))     #Timoshenko
        t1 = (4+Phi)
        t2 = (2-Phi)
        d1 = rho*A*L_e
        d2 = (I*L_e)/6
        d3 = (rho*A*L_e)/420
        #print(A, I, J, d, e)
        # Matriz de Rigidez Elementar (Euler-Bernoulli)
        # Para converter para timoshenko basta trocar c3 por c4,onde tem (4 * L_e**2 * c3) substitui por (t1* L_e**2 * c4) e onde tiver (2 * L_e**2 * c3) por (t2* L_e**2 * c4))
        k_e= np.array([
                [12 * c4, 0, 0, 6 * L_e * c4, 0, 0, -12 * c4, 0, 0, 6 * L_e * c4, 0, 0],
                [0, c1, 0, 0, 0, 0, 0, -c1, 0, 0, 0, 0],
                [0, 0, 12 * c4, 0, 0, 6 * L_e* c4, 0, 0, -12 * c4, 0, 0, 6 * L_e * c4],
                [6 * L_e * c4, 0, 0, t1* L_e**2 * c4, 0, 0, -6 * L_e * c4, 0, 0, t2 * L_e**2 * c4, 0, 0],
                [0, 0, 0, 0, c2, 0, 0, 0, 0, 0, -c2, 0],
                [0, 0, 6 * L_e * c4, 0, 0, t1* L_e**2 * c4, 0, 0, -6 * L_e * c4, 0, 0, t2 * L_e**2 * c4],
                [-12 * c4, 0, 0, -6 * L_e * c4, 0, 0, 12 * c4, 0, 0, -6 * L_e * c4, 0, 0],
                [0, -c1, 0, 0, 0, 0, 0, c1, 0, 0, 0, 0],
                [0, 0, -12 * c4, 0, 0, -6 * L_e * c4, 0, 0, 12 * c4, 0, 0, -6 * L_e * c4],
                [6 * L_e * c4, 0, 0, t2 * L_e**2 * c4, 0, 0, -6 * L_e * c4, 0, 0, t1* L_e**2 * c4, 0, 0],
                [0, 0, 0, 0, -c2, 0, 0, 0, 0, 0, c2, 0],
                [0, 0, 6 * L_e * c4, 0, 0, t2 * L_e**2 * c4, 0, 0, -6 * L_e * c4, 0, 0, t1* L_e**2 * c4]
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

    def aplicar_engastes(self, nodes, dofs):
        """
        Applies constraints (fixed DOFs) on specific nodes.
        Inputs:
            - nodes: list of node indices to be constrained.
            - dofs: list of degrees of freedom to be fixed.
        Outputs: None.
        """
        for node in nodes:                                          # Laço para selecionar cada nó que será engastado
            for dof in dofs:                                        # Laço para selecionar quais graus de liberdade serão fixados
                index = node * self.num_dofs_per_node + dof         # Identificação da entrada da matriz que precisa ser restringida pelo engaste        
                self.K_global[index, index] = 10**9                # Um valor suficientemente grande para simular um engaste 

    def matrizes_global(self):
        """
        Assembles the global stiffness and mass matrices.
        Inputs: None (uses class attributes).
        Outputs:
            - K_global: global stiffness matrix.
            - M_global: global mass matrix.
        """
        for element in self.elements:
            node1, node2, tube_type = element
            k_e, m_e = self.element(element)
            # DOFs associados ao elemento            
            dofs = [6 * node1, 6 * node1 + 1, 6 * node1 + 2, 6 * node1 + 3, 6 * node1 + 4, 6 * node1 + 5,
                    6 * node2, 6 * node2 + 1, 6 * node2 + 2, 6 * node2 + 3, 6 * node2 + 4, 6 * node2 + 5]

            # Atualizando as matrizes globais
            self.K_global[np.ix_(dofs, dofs)] += k_e
            self.M_global[np.ix_(dofs, dofs)] += m_e

        #self.aplicar_engastes([30,31,32,33], [0, 1, 2, 3, 4, 5])                             #Por enquanto não estaremos considerando engastes
        pd.DataFrame(self.K_global).to_csv('Matriz_Global_Rigidez.csv', index=True, header=True)
        pd.DataFrame(self.M_global).to_csv('Matriz_Global_Massa.csv', index=True, header=True)        

        # print (self.K_global)
        # print (self.M_global)

        #plt.figure(figsize=(6, 6))
        #plt.spy(self.K_global, markersize=10)  # Adjust markersize for visibility
        #plt.title("Spy Plot of the Kg")
        #plt.xlabel("Columns")
        #plt.ylabel("Rows")
        #plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        #plt.show()

        #plt.figure(figsize=(6, 6))
        #plt.spy(self.M_global, markersize=10)  # Adjust markersize for visibility
        #plt.title("Spy Plot of the Mg")
        #plt.xlabel("Columns")
        #plt.ylabel("Rows")
        #plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        #plt.show()

        return self.K_global,self.M_global

    def shape_fun(self, F_flexao1=1000, F_flexao2=1000, F_axial=1000,F_torcao=1000): 
        """
        Calculates deformations and stiffness of elements under loads.
        Inputs:
            - F_flex1: array of point bending forces.
            - F_flex2: array of distributed bending forces.
            - F_axial: array of axial forces.
            - F_torsion: array of torsion forces.
        Outputs:
            - Arrays of torsion, deformations, and stiffness (bending and torsional).
        """
        #E = 2.1e11  	#Modulo de Young (Pa)
        #I = 1.6667e-5 	#Momento de inercia (m^4)
        #G = 81.2e9  	#Modulo de Cisalhamento (Pa)
        #A= 0.0125	    #Área da seção do elemento (m^2)	
        #J = I/2     	#Momento polar de inércia (m^4)
        KF_total = 0
        KT_total = 0
        KF_elements = []
        KT_elements = [] 
        torcao, deformacao, flexao1, flexao2, flexao3 = [], [], [], [], []
        for element in self.elements:
            d = self.obter_propriedades(element[2])[2]         #Diâmetro do Tubo (m)
            e = self.obter_propriedades(element[2])[3]         #Espessura do Tubo (m)
            E = self.obter_propriedades(element[2])[0]         #Modulo de Young (Pa)      
            G = self.obter_propriedades(element[2])[1]         #Modulo de Cisalhamento (Pa)
            A = self.area_seccao_transversal(d, e)             #Área da seção do elemento (m^2)
            I = self.momento_inercia_area_e_polar(d,e)[0]      #Momento de inercia (m^4)
            J = self.momento_inercia_area_e_polar(d,e)[1]      #Momento polar de inércia (m^4)
            L_e = self.calcular_comprimento(element)
            # Equação de torsão
            torcao_val = (F_torcao * L_e) / (G * J)         #Fonte[1]
            torcao.append(torcao_val)
            # Equação  para deformação axial
            deformacao_val = (F_axial* L_e / (A * E))       #Fonte[2]
            deformacao.append(deformacao_val)
            # Equação para flexão
            flexao_val1 = (F_flexao1*L_e**3)/(48 * E * I)      #Fonte[3.1] (carga pontual no meio do elemento biapoiado)
            flexao_val2 = (5*F_flexao2*L_e**4)/(384 * E * I)   #Fonte[3.2] (carga distribuída ao longo de todo o elemento biapoiado)
            flexao_val3 = flexao_val1 + flexao_val2            #Fonte[3.3] (tentativa de carregamento misto)
            flexao1.append(flexao_val1)
            flexao2.append(flexao_val2)
            flexao3.append(flexao_val3)

            # Rigidez flexional
            KF = E * I / L_e
            # Rigidez torsional
            KT = G * J / L_e

            KF_total += KF
            KT_total += KT

            KF_elements.append(KF)
            KT_elements.append(KT)
            
        return (np.array(torcao), np.array(deformacao), np.array(flexao1), 
            np.array(flexao2), np.array(flexao3), KF_total, KT_total, KF_elements, KT_elements)

    def modal_analysis(self):
        """
        Performs modal analysis to compute natural frequencies and mode shapes.
        Inputs: None.
        Outputs:
            - eigenvalues: eigenvalues (squared natural frequencies).
            - eigenvectors: eigenvectors (mode shapes).
            - frequencies: natural frequencies (Hz).
        """
        # Análise modal por resolução do problema de autovalor e autovetor
        unsorted_eigenvalues, unsorted_eigenvectors = eigh(self.K_global, self.M_global)

        # Frequências naturais (raiz quadrada dos autovalores)
        unsorted_frequencies = np.sqrt(unsorted_eigenvalues) / (2 * np.pi)  # Divisão por 2*pi para converter para hertz

        # Tratando os dados (tomando apenas as 20 primeiras frequências naturais)
        sorted_indices = np.argsort(unsorted_frequencies)  # Ordena as frequências em ordem crescente
        top_indices = sorted_indices[:self.num_modes]  # Seleciona os índices dos primeiros n modos

        eigenvalues = np.array(unsorted_eigenvalues)[top_indices]  # Filtra os primeiros n autovalores
        eigenvectors = np.array(unsorted_eigenvectors)[:, top_indices]  # Filtra os primeiros n autovetores
        frequencies = np.array(unsorted_frequencies)[top_indices]  # Filtra as primeiras n frequências

        return eigenvalues, eigenvectors, frequencies
    
    def static_analysis(self, F_global, fixed_dofs):
        """
        Perform static analysis by solving Ku = F with boundary conditions.

        Parameters:
            F_global (ndarray): Global force vector (N).
            fixed_dofs (list): List of DOF indices to be fixed.

        Returns:
            displacements (ndarray): Displacement vector (N).

        Resolve uma análise estática para deslocamentos em DOFs livres.
        Entradas:
            - K_global: matriz de rigidez global.
            - F_global: vetor de forças globais.
            - fixed_dofs: índices de graus de liberdade fixos.
        Saídas:
            - displacements: vetor de deslocamentos nos DOFs.
        """
        # Total number of DOFs
        n_dofs = self.K_global.shape[0]

        # Create a mask for free DOFs (DOFs not constrained)
        free_dofs = np.array([i for i in range(n_dofs) if i not in fixed_dofs])

        # Reduce the stiffness matrix and force vector
        K_reduced = self.K_global[np.ix_(free_dofs, free_dofs)]
        F_reduced = F_global[free_dofs]

        # Solve for displacements at free DOFs
        # USE "linalg.lstsq" FOR NEAR SINGULAR MATRICES (ALL OF THEM)
        u_reduced = np.linalg.lstsq(K_reduced, F_reduced, rcond=None)[0]  # Get only the solution vector

        # Construct full displacement vector
        displacements = np.zeros(n_dofs)
        displacements[free_dofs] = u_reduced
        
        return displacements

#    def calcular_B_Elementar(self, displacements, element):
#         """
#         Calcula a matriz B de um elemento individual.
#         """
#         L_e = self.calcular_comprimento(element)
#         node1, node2, tube_type = element
#         # Extract displacements for the nodes of the element
#         dofs_node1 = displacements[node1 * self.num_dofs_per_node: (node1 + 1) * self.num_dofs_per_node]
#         dofs_node2 = displacements[node2 * self.num_dofs_per_node: (node2 + 1) * self.num_dofs_per_node]
#         # Construct the B matrix (simplified for axial and bending)
#         B = np.array([
#             [-(dofs_node2[1] - dofs_node1[1]) / L_e**2, 0, 0, 0, 0, 0, (dofs_node2[1] - dofs_node1[1]) / L_e**2, 0, 0, 0, 0, 0],       # Y
#             [0, -(dofs_node2[2] - dofs_node1[2]) / L_e**2, 0, 0, 0, 0, 0, (dofs_node2[2] - dofs_node1[2]) / L_e**2, 0, 0, 0, 0],       # Z
#             [0, 0, -(dofs_node2[3] - dofs_node1[3]) / L_e**2, 0, 0, 0, 0, 0, (dofs_node2[3] - dofs_node1[3]) / L_e**2, 0, 0, 0],       # RX
#             [0, 0, 0, -(dofs_node2[4] - dofs_node1[4]) / L_e**2, 0, 0, 0, 0, 0, (dofs_node2[4] - dofs_node1[4]) / L_e**2, 0, 0],       # RY
#             [-(dofs_node2[0] - dofs_node1[0]) / L_e, 0, 0, 0, 0, 0, (dofs_node2[0] - dofs_node1[0]) / L_e, 0, 0, 0, 0, 0],             # X
#             [0, 0, 0, 0, -(dofs_node2[5] - dofs_node1[5]) / L_e, 0, 0, 0, 0, 0, (dofs_node2[5] - dofs_node1[5]) / L_e, 0]              # RZ
#         ])
#         return B

#    def calcular_B_Elementar(self, displacements, element):
#        """
#        Calcula a matriz B de um elemento individual.
#        """
#        L_e = self.calcular_comprimento(element)
#        node1,node2,_= element
#        # Extract displacements for the nodes of the element
#        dofs_node1 = displacements[node1 * self.num_dofs_per_node: (node1 + 1) * self.num_dofs_per_node]
#        dofs_node2 = displacements[node2 * self.num_dofs_per_node: (node2 + 1) * self.num_dofs_per_node]
#
#        #Construct the B matrix
#        B = np.array([
#        [  0   ,-1/L_e,   0  ,    0 ,    0 ,    0 ,  0  , 1/L_e,   0 ,    0 ,    0 ,   0    ],
#
#        [  0   ,   0  ,   0  ,    0 ,-1/L_e,   0  ,  0  ,  0   ,  0  ,   0  ,1/L_e ,   0    ],
#
#        [  0   ,   0  ,   0  ,-1/L_e,   0  ,   0  ,  0  ,  0   ,  0  ,1/L_e ,  0   ,   0    ],
#
#        [  0   ,   0  ,   0  ,    0 ,    0 ,-1/L_e,  0  ,  0   ,  0  ,   0  ,   0  , 1/L_e  ],
#
#        [-1/L_e,   0  ,   0  ,    0 ,    0 ,-1+dofs_node1[5]/L_e,1/L_e,  0   ,  0  ,   0  ,   0  , -dofs_node2[5]/L_e],
#            
#        [  0   ,   0  ,-1/L_e, 1-dofs_node1[3]/L_e,   0  ,   0  ,  0  ,  0   ,1/L_e, dofs_node2[3]/L_e,  0   ,   0    ],
#        ])
#        return B

    def calcular_B_Elementar(self, element):
        """
        Calcula a matriz B (6x12) para um elemento de viga de Timoshenko
        com o eixo longitudinal no eixo Y.
        """
        L = self.calcular_comprimento(element)
        B =  [  [0,       -1/L,     0,       0,     0,       0,       0,       1/L,     0,       0,     0,       0     ],  # ε_yy (axial)
                [0,        0,       0,    -1/L,     0,       0,       0,        0,      0,     1/L,     0,       0     ],  # κ_x (curvatura sobre X)
                [0,        0,       0,       0,     0,   -1/L,        0,        0,      0,       0,     0,     1/L     ],  # κ_z (curvatura sobre Z)
                [-1/L,     0,       0,       0,     0,   -0.5,     1/L,        0,      0,       0,     0,    -0.5      ],  # γ_xy (cisalhamento XY)
                [0,        0,   -1/L,     0.5,      0,     0,        0,        0,    1/L,     0.5,     0,       0      ],  # γ_zy (cisalhamento ZY)
                [0,        0,       0,       0,   -1/L,     0,        0,        0,      0,       0,   1/L,       0     ]]   # φ_y (torção)

        return B

    def compute_strain(self, displacements):
        """
        Compute strains for all elements. The B_matrices (Strain-displacement 
        matrices for each element are computed here, locally, by the will of
        Newton.

        Parameters:
            displacements (ndarray): Displacement vector for all nodes.

        Returns:
            strains (list of ndarray): Strain tensors for all elements.
        """

        strains = []
        for element in self.elements:
            B = self.calcular_B_Elementar(element)
            node1, node2, tube_type = element
            element_dofs = []
            for node in [node1, node2]:
                for dof in range(self.num_dofs_per_node):
                    element_dofs.append(node * self.num_dofs_per_node + dof)
            element_displacements = displacements[element_dofs]
            strain = np.dot(B, element_displacements)  # B-matrix times element's displacement vector
            strains.append(strain)
        return strains
    
#    def compute_stress(self, strains):
#        """
#        Compute stresses for all elements using Hooke's law.
#
#        Parameters:
#            strains (list of ndarray): Strain tensors for all elements.
#            E (float): Young's modulus.
#            nu (float): Poisson's ratio.
#
#        Returns:
#            stresses (list of ndarray): Stress tensors for all elements.
#        """
#        stresses = [] 
#        # Construct constitutive matrix (isotropic 3D elasticity)
#        for i in range(len(self.elements)):
#            element=self.elements[i]
#            E = self.obter_propriedades(element[2])[0]         #Modulo de Young (Pa)      
#            G = self.obter_propriedades(element[2])[1]         #Modulo de Cisalhamento (Pa)
#            nu = self.obter_propriedades(element[2])[5]
#            lambda_ = (E * nu) / ((1 + nu) * (1 - 2 * nu))
#            C = np.array([
#                [lambda_ + 2*G  , lambda_       , lambda_       ,   0,  0,  0],
#                [lambda_        , lambda_ + 2*G , lambda_       ,   0,  0,  0],
#                [lambda_        , lambda_       , lambda_ + 2*G ,   0,  0,  0],
#                [              0,              0,              0,   G,  0,  0],
#                [              0,              0,              0,   0,  G,  0],
#                [              0,              0,              0,   0,  0,  G]
#            ])
#
#            strain=strains[i]
#            stress = np.dot(C, strain)  # Hooke's law: C times strain
#            stresses.append(stress)
#        return stresses

    def compute_stress(self, strains):
        """
        Compute stresses for all elements using a Timoshenko-beam diagonal constitutive matrix.
        """
        stresses = []
        for i, element in enumerate(self.elements):
            # extrai propriedades
            E, G, d, e, rho, nu = self.obter_propriedades(element[2])

            # matriz constitutiva diagonal (6×6)
            C = np.array([
                [E,    0.0,  0.0,  0.0,  0.0,  0.0],  # axial ε_yy → σ_yy = E * ε_yy
                [0.0,  E,    0.0,  0.0,  0.0,  0.0],  # curvatura κ_x → M_x/E (não usado para tensão direta)
                [0.0,  0.0,  E,    0.0,  0.0,  0.0],  # curvatura κ_z → M_z/E
                [0.0,  0.0,  0.0,  G,    0.0,  0.0],  # cisalhamento γ_xy → τ_xy = G * γ_xy
                [0.0,  0.0,  0.0,  0.0,  G,    0.0],  # cisalhamento γ_zy → τ_zy = G * γ_zy
                [0.0,  0.0,  0.0,  0.0,  0.0,  G   ]   # torção φ_y → τ_y = G * φ_y
            ])

            stress = C @ strains[i]
            stresses.append(stress)
        return stresses
        
    def compute_Kf_Kt(self, K_global, lfNode, lrNode, rfNode, rrNode, lcNode, rcNode):
        """
        Calculates Torsional and Beaming stiffness of chassis
        Inputs:
        - K_global: Global stiffness matrix
        - lfNode: Left front suspension node index
        - lrNode: Left rear suspension node index
        - rfNode: Right front suspension node index
        - rrNode: Right rear suspension node index
        - lcNode: Left central node index
        - rcNode: Right central node index
        Outputs:
        - Kt: Torsional stiffness of the chassis
        - Kf: Beaming stiffness of the chassis
        """
        #Start Kt simulation
        F_global = np.zeros(K_global.size)
        F_global[2+lfNode*6] = 250  #Forças aplicadas nos nós onde estaria a suspensão dianteira
        F_global[2+rfNode*6] = -250  #Mesmo módulo para gerar um torque no eixo longitudinal do chassi
        
        fixed_nodes=[lrNode, rrNode]  #Fixação dos nós onde estaria a suspensão traseira
        fixed_dofs=[(node*6+i) for node in fixed_nodes for i in range(6)]  #Lista com os dofs fixados
        
        displacements = self.static_analysis(F_global, fixed_dofs)  #Calcula os displacements com as condições de contorno aplicadas acima
        mi1=np.abs(displacements[lfNode*6+2])  
        mi2=np.abs(displacements[rfNode*6+2])  
        L = np.abs(self.nodes[lfNode][0] - self.nodes[rfNode][0])  #0 - pega a coordenada no eixo x
        alpha= np.degrees(np.atan((mi1+mi2)/(L)))  #Ângulo de torção do chassi após aplicação do torque
        tau = (np.abs(F_global[2+rfNode*6]))*(L)  #Cálculo do torque aplicado
        
        Kt = tau/alpha

        #Start Kf simulation
        F_global = np.zeros(K_global.size)
        F = 5000  #Módulo da força que vai gerar a flexão
        F_global[2+lcNode*6] = -F/2  #Força distribuída nos nós centrais do chassi (nodes 22 e 23)
        F_global[2+rcNode*6] = -F/2  #Sinal negativo por conta da direção de aplicação da força
        
        fixed_nodes=[rfNode, lfNode, rrNode, lrNode]  #Fixação dos nós onde estaria a suspensão dianteira e traseira
        fixed_dofs=[(node*6+i) for node in fixed_nodes for i in range(6)]  #Lista com os dofs fixados
        
        displacements = self.static_analysis(F_global, fixed_dofs)  #Calcula os displacements com as condições de contorno aplicadas acima
        dY=np.abs(displacements[2+lcNode*6])  #Deslocamento em Y de um dos nós onde foi aplicado a força

        Kf=F/dY

        return Kt, Kf

    def compute_von_mises(self, stresses):
        """
        Compute von Mises stress for all elements.

        Parameters:
            stresses (list of ndarray): Stress tensors for all elements.

        Returns:
            von_mises_stresses (list of float): Von Mises stress for each element.
        """
        von_mises_stresses = []
        for stress in stresses:
            sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_zx = stress
            von_mises = np.sqrt(
                0.5 * (
                    (sigma_xx - sigma_yy)**2 +
                    (sigma_yy - sigma_zz)**2 +
                    (sigma_zz - sigma_xx)**2 +
                    6 * (tau_xy**2 + tau_yz**2 + tau_zx**2)
                )
            )
            von_mises_stresses.append(von_mises/10**6)
        return np.array(von_mises_stresses)

    def Mesh(self):
        """
        Generates a `.geo` file for the structure mesh in GMSH.
        Inputs: None (uses class attributes and user-provided file name).
        Outputs: None.
        """
        
        filename = input("Insira o nome do arquivo: ") + ".geo"
        diretorio = input("Insira o diretorio onde o arquivo .geo deve ser salvo: ")

        if not os.path.exists(diretorio):
            os.makedirs(diretorio)

        filepath = os.path.join(diretorio, filename)

        with open(filepath, 'w') as geo_file:
            for i, (x, y, z) in enumerate(self.nodes):
                geo_file.write(f'Point({i + 1}) = {{{x}, {y}, {z}, 1.0}};\n')

            for i, (start, end) in enumerate(self.elements):
                geo_file.write(f'Line({i + 1}) = {{{start + 1}, {end + 1}}};\n')

            if len(self.elements) > 2:
                line_loop_indices = ', '.join(str(i + 1) for i in range(len(self.elements)))
                geo_file.write(f'Line Loop(1) = {{{line_loop_indices}}};\n')
                geo_file.write('Plane Surface(1) = {1};\n')

            geo_file.write('Mesh.Algorithm = 6;\n')
            geo_file.write('Mesh.ElementOrder = 1;\n')
            geo_file.write('Mesh.Format = 1;\n')

        print(f'O arquivo foi salvo em: {filepath}, basta abrir o GMSH, e abrir o arquivo')

    def plot_colored_wireframe(nodes, elements, scalar_values, graphtitle='Wireframe plot', scalelabel='Your variable here', colormap='jet'):
        """
        Plots a 3D wireframe of the structure with color mapping based on scalar values.
        
        Parameters:
            nodes (array): Array of node coordinates (N x 3).
            elements (list): List of tuples defining connections between nodes.
            scalar_values (array): 1D array of scalar values (e.g., strain) at each node.
            graphtitle (str): Graph title.
            scalelabel (str): LColormap scale label description.
            colormap (str): Colormap name for visualization.
        """
        # Normalize scalar values to [0, 1] for colormap
        norm = plt.Normalize(vmin=np.min(scalar_values), vmax=np.max(scalar_values))
        cmap = plt.get_cmap(colormap)
        
        # Create the plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot each element with color based on scalar values
        for node1, node2, tube_type in elements:
            # Get coordinates for the two nodes
            x = [nodes[node1][0], nodes[node2][0]]
            y = [nodes[node1][1], nodes[node2][1]]
            z = [nodes[node1][2], nodes[node2][2]]
            
            # Get the scalar value for the midpoint of the element
            scalar_midpoint = (scalar_values[node1] + scalar_values[node2]) / 2
            
            # Map scalar value to color
            color = cmap(norm(scalar_midpoint))
            
            # Plot the line segment with the corresponding color
            ax.plot(x, y, z, color=color, linewidth=2)

        # Add a colorbar
        mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        mappable.set_array(scalar_values)
        cbar = plt.colorbar(mappable, ax=ax, orientation='vertical', shrink=0.8, pad=0.1)
        cbar.set_label(scalelabel, fontsize=12)

        # Set axis labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(graphtitle)
        ax.set_box_aspect([1,3,2])
        #plt.xlim([min(nodes)*1.1,max(nodes)*1.1])
        #plt.ylim([min(nodes)*1.1,max(nodes)*1.1])
        plt.tight_layout()
        plt.show()

    def modal_analysis_plot(self):
        """
        Plots the modal shapes of the structure in 3D, showing the deformed and original configurations.

        This function performs the following:
        - Computes eigenvalues and eigenvectors using the modal analysis method.
        - Extracts the translational displacements (x, y, z) for each node from the modal shapes.
        - Scales and visualizes the deformed structure for each mode, overlaying it on the original structure.
        - Annotates the nodes of both the original and deformed structures.

        Parameters:
            None (uses class attributes):
                - coordinates (array): Array of node coordinates (N x 3).
                - connections (list): List of tuples defining element connections between nodes.
                - modal_analysis (method): Computes eigenvalues and eigenvectors.
        """     
        autovalores, autovetores, _ = self.modal_analysis()

        for mode_idx in range(len(autovalores)):
            mode_shape = autovetores[:, mode_idx]
            displacements = np.zeros((len(self.nodes), 3))  # Assuming we want to visualize x, y, z displacements only

            # Loop through nodes to extract the translations
            for j, (x, y, z) in enumerate(self.nodes):
                # 6 DOFs per node: [u_x, u_y, u_z, theta_x, theta_y, theta_z]
                dof_start = 6 * j  # Start index of DOFs for node j
                displacements[j, 0] = mode_shape[dof_start]     # u_x
                displacements[j, 1] = mode_shape[dof_start+1] # u_y
                displacements[j, 2] = mode_shape[dof_start+1] # u_z

            # Scale displacements for plots
            scale_factor = 1  # Adjust as needed
            deformed_nodes = np.array(self.nodes) + displacements * scale_factor

            # Plot deformed
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')

            # Plot deformed structure
            for i, (node1, node2,_) in enumerate(self.elements):
                x = [deformed_nodes[node1][0], deformed_nodes[node2][0]]
                y = [deformed_nodes[node1][1], deformed_nodes[node2][1]]
                z = [deformed_nodes[node1][2], deformed_nodes[node2][2]]
                ax.plot(x, y, z, 'r-', label="Deformed" if i == 0 else "")  # Add label only once

            #Colocando a legenda dos nós no gráfico
            for i, (x, y, z) in enumerate(self.nodes):
                ax.scatter(x, y, z, color='b', s=20)
                ax.text(x, y, z, f'  {i+1}', color='black', fontsize=8)

            #Colocando a legenda dos nós após a deformação no gráfico
            #for i, (x, y, z) in enumerate(deformed_nodes):
            #    ax.scatter(x, y, z, color='r', s=25)
            #    ax.text(x, y, z, f'  {i+1}', color='black', fontsize=8)

            # Plot original structure
            for i, (node1, node2,_) in enumerate(self.elements):
                x = [self.nodes[node1][0], self.nodes[node2][0]]
                y = [self.nodes[node1][1], self.nodes[node2][1]]
                z = [self.nodes[node1][2], self.nodes[node2][2]]
                ax.plot(x, y, z, 'k--', label="Original" if i == 0 else "")  # Add label only once

            # Add labels and title
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_title(f'Forma modal nº: {mode_idx + 1}')
            ax.legend()  # Ensure the legend is displayed
            ax.set_zlim([-0.5,1.5])
            ax.set_ylim([0,2])
            plt.tight_layout()
            plt.show()

    def shape_fun_plot(self, F_flexao1, F_flexao2, F_axial, F_torcao):
        """
            Generates plots of deformations and stiffness values for each element based on the given forces.

            This function calls `shape_fun` to calculate torsional, axial, and flexural deformations, as well as stiffness values.
            It then plots these results across subplots to visualize the behavior of each element under the applied forces.

            Parameters:
                F_flexao1 (float): Force applied for flexural deformation (point load at mid-span).
                F_flexao2 (float): Force applied for distributed flexural deformation.
                F_axial (float): Axial force applied to the elements.
                F_torcao (float): Torsional force applied to the elements.

            Plots:
                1. Torsional deformation for each element.
                2. Axial deformation for each element.
                3. Flexural deformation due to point load for each element.
                4. Flexural deformation due to distributed load for each element.
                5. Combined flexural deformation for each element.
                6. Flexural and torsional stiffness for each element.

            Notes:
                - Total flexural stiffness (KF_total) and torsional stiffness (KT_total) are displayed in the overall plot title.
                - The function uses subplots to organize the visuals, and the layout is adjusted for clarity.
            """
        torcao,deformacao_axial,flexao1,flexao2,flexao3,KF_total,KT_total,KF_elements,KT_elements= self.shape_fun(F_flexao1, F_flexao2, F_axial, F_torcao)
        
        # Configuração dos subplots
        fig, axs = plt.subplots(6, 1, figsize=(12, 22))

        # Plot da Torção
        axs[0].plot(torcao, 'o-', label=[f'Força {F}N' for F in F_torcao])
        axs[0].set_title('Deformação por Torção de cada Elemento')
        axs[0].set_xlabel('Elemento')
        axs[0].set_ylabel('Torção (rad)')
        axs[0].legend()

        # Plot da Deformação Axial
        axs[1].plot(deformacao_axial, 's-', label=[f'Força {F}N' for F in F_axial])
        axs[1].set_title('Deformação Axial de cada Elemento')
        axs[1].set_xlabel('Elemento')
        axs[1].set_ylabel('Deformação (m)')
        axs[1].legend()

        # Plot da Flexão por Carga Pontual
        axs[2].plot(flexao1, 'o-', label=[f'Força {F}N' for F in F_flexao1])
        axs[2].set_title('Deformação por Carga Pontual de cada Elemento')
        axs[2].set_xlabel('Elemento')
        axs[2].set_ylabel('Deflexão (m)')
        axs[2].legend()

        # Plot da Flexão por Carga Distribuída
        axs[3].plot(flexao2, 'o-', label=[f'Força {F}N' for F in F_flexao2])
        axs[3].set_title('Deformação por Carga Distribuída de cada Elemento')
        axs[3].set_xlabel('Elemento')
        axs[3].set_ylabel('Deflexão (m)')
        axs[3].legend()

        # Plot da Flexão Mista
        axs[4].plot(flexao3, 'o-', label='Carregamento misto')
        axs[4].set_title('Deformação por Flexão Mista de cada Elemento')
        axs[4].set_xlabel('Elemento')
        axs[4].set_ylabel('Deflexão (m)')
        axs[4].legend()

        # Plot da Rigidez Flexional e Torsional por Elemento
        axs[5].plot(KF_elements, 'o-', label='Rigidez Flexional (KF)')
        axs[5].plot(KT_elements, 's-', label='Rigidez Torsional (KT)')
        axs[5].set_title('Rigidez Flexional e Torsional de cada Elemento')
        axs[5].set_xlabel('Elemento')
        axs[5].set_ylabel('Rigidez (N/m)')
        axs[5].legend()

        # Mostrando os totais no título geral
        plt.suptitle(f'KF Total: {KF_total:.2e} N/m, KT Total: {KT_total:.2e} N/m', fontsize=16)

        # Ajustes de layout
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()
        print(f'KF Total: {KF_total:.2e} N/m \nKT Total: {KT_total:.2e} N/m')



def init_worker(optimizer_instance):
    """
    Inicializador de cada worker: atribui a instância global OPTIMIZER
    """
    global OPTIMIZER
    OPTIMIZER = optimizer_instance

class ChassisDEOptimizer:   
    """
    Otimizador de geometria de chassi tubular usando Differential Evolution (DE).

    Atributos:
    - base_nodes: np.ndarray, coordenadas iniciais dos nós de um lado (x>=0).
    - base_connections: list of tuple, pares de índices representando arestas base.
    - mandatory_indices: índices de nós que têm deslocamento limitado a radius_mand.
    - pop_size, F, CR, max_gens: parâmetros do DE (tamanho pop., taxa de mutação, crossover, gerações).
    - radius_mand, radius_opt: raios máximos de deslocamento para nós mandatórios e opcionais.
    - tipos_tubos: lista de strings com perfis de tubo permitidos.
    """

    def __init__(
        self,
        base_nodes: np.ndarray,
        base_connections: list,
        mandatory_indices: list,
        pop_size: int = 10,
        F: float = 0.5,
        CR: float = 0.9,
        max_generations: int = 20,
        radius_mand: float = 0.025,
        radius_opt: float = 0.05,
        use_parallel: bool = True,
        n_workers: int = None,
    ):
        """
        Inicializa o otimizador.

        Entradas:
        - base_nodes: array (n,3) de floats.
        - base_connections: lista de tuplas (i,j).
        - mandatory_indices: lista de inteiros.
        - pop_size, F, CR, max_generations: parâmetros DE.
        - radius_mand, radius_opt: floats definindo limites de deslocamento.

        Retorno:
        - Nenhum (configura atributos internos).
        """
        self.base_nodes = base_nodes.copy()
        self.n = base_nodes.shape[0]
        self.n_tubes = len(base_connections)
        self.base_connections = base_connections
        self.mandatory = set(mandatory_indices)
        self.radius_mand = radius_mand
        self.radius_opt = radius_opt

        self.tipos_tubos = ['Tubo A', 'Tubo B', 'Tubo C', 'Tubo D']

        self.pop_size = pop_size
        self.F = F
        self.CR = CR
        self.max_gens = max_generations

        self.dim_coords = 3 * self.n
        self.dim_tubes = self.n_tubes
        self.dim = self.dim_coords + self.dim_tubes

        self.use_parallel = use_parallel
        self.n_workers = max(1, os.cpu_count() if n_workers is None else n_workers)

        self.is_central = np.isclose(self.base_nodes[:, 0], 0.0)
        self.central_alignment_map = self._create_alignment_map()

        self.global_decoded_cache = {}
        
    def _create_alignment_map(self,exceptions=[15,20]):
        """
        Cria um dicionário mapeando cada nó central ao nó lateral que ele deve espelhar.
        exceptions não são mapeados.
        """
        alignment_map = {}
        central_nodes = [i for i in range(self.n) if self.is_central[i]]
        for i, j in self.base_connections:
            if i in central_nodes and not self.is_central[j] and i not in exceptions:
                alignment_map.setdefault(i, []).append(j)
            elif j in central_nodes and not self.is_central[i] and j not in exceptions:
                alignment_map.setdefault(j, []).append(i)
        return alignment_map

    def apply_central_alignment(self, coords: np.ndarray) -> np.ndarray:
        """
        Força nós centrais a terem o mesmo Y e Z dos nós laterais conectados.
        Para múltiplas conexões, usa a média das coordenadas.
        """
        aligned_coords = coords.copy()
        
        for central_idx, lateral_indices in self.central_alignment_map.items():
            aligned_coords[central_idx, 1] = np.mean(coords[lateral_indices, 1])
            aligned_coords[central_idx, 2] = np.mean(coords[lateral_indices, 2])
            aligned_coords[central_idx, 0] = 0.0
            
        return aligned_coords
    
    def enforce_bounds(self, coords: np.ndarray) -> np.ndarray:
        """
        Aplica limites de deslocamento e arredonda as coordenadas.

        Entrada:
        - coords: array (n,3) de floats (proposta de deslocamento).
        Saída:
        - adjusted: array (n,3) de floats, cada deslocamento limitado a radius_mand ou radius_opt
          e arredondado para 3 casas decimais.
        """
        adjusted = coords.copy()
        for i in range(self.n):
            orig = self.base_nodes[i]

            if self.is_central[i]:
                adjusted[i, 0] = 0.0
                delta = coords[i, 1:] - orig[1:]
                dist = np.linalg.norm(delta)
                r = self.radius_mand if i in self.mandatory else self.radius_opt
                if dist > r:
                    adjusted[i, 1:] = orig[1:] + (delta / dist) * r
            else:
                delta = coords[i] - orig
                dist = np.linalg.norm(delta)
                r = self.radius_mand if i in self.mandatory else self.radius_opt
                if dist > r:
                    adjusted[i] = orig + (delta / dist) * r
        adjusted[self.is_central, 0] = 0.0
        return np.round(adjusted, 3)
    
    def validate_min_distance(self, coords: np.ndarray, min_dist: float = 0.05) -> bool:
        """
        Verifica se todas as distâncias entre pares de nós são maiores ou iguais a uma distância mínima.
        Vetorizado usando pdist.

        Parâmetros:
        - coords: np.ndarray (N, 3), coordenadas dos nós.
        - min_dist: float, distância mínima permitida.

        Retorno:
        - bool: True se válido, False se algum par violar a distância mínima.
        """
        return np.all(pdist(coords) >= min_dist)
 
    def decode_individual(self, x: np.ndarray):
        """
        Converte um vetor genotípico em nós completos e lista de elementos.

        Entrada:
        - x: array (dim_coords+dim_tubes,) de floats.
        Saída:
        - nodes_full: array (N,3) de floats com nós de ambos os lados.
        - elements: lista de tuplas (i, j, perfil), com índices nos nodes_full.

        Processos:
        1. Extrai coords e tube_vars do vetor x.
        2. Aplica enforce_bounds às coords.
        3. Para nós com base x≈0, inclui apenas um nó (central).
           Para outros, inclui coord e seu espelho.
        4. Usa mapeamento para gerar conexões simétricas com tipo de tubo.
        """
        # Separa coords e tube_vars
        coords = x[:self.dim_coords].reshape((self.n, 3))
        coords = self.enforce_bounds(coords)
        coords = self.apply_central_alignment(coords)
        tube_vars = x[self.dim_coords:]

        full_nodes = []
        mapping = {}
        for i, coord in enumerate(coords):
            if self.is_central[i]:
                idx = len(full_nodes)
                full_nodes.append(coord)
                mapping[i] = [idx]
            else:
                idx1 = len(full_nodes)
                full_nodes.append(coord)
                mirrored = coord.copy()
                mirrored[0] *= -1
                idx2 = len(full_nodes)
                full_nodes.append(mirrored)
                mapping[i] = [idx1, idx2]

        nodes_full = np.array(full_nodes)
        elements = []
        for idx_conn, (i, j) in enumerate(self.base_connections):
            perfil = self.tipos_tubos[int(np.clip(np.floor(tube_vars[idx_conn]), 0, len(self.tipos_tubos) - 1))]
            ids_i, ids_j = mapping[i], mapping[j]
            if len(ids_i) == 1 and len(ids_j) == 1:
                elements.append((ids_i[0], ids_j[0], perfil))
            elif len(ids_i) == 2 and len(ids_j) == 2:
                elements.append((ids_i[0], ids_j[0], perfil))
                elements.append((ids_i[1], ids_j[1], perfil))
            else:
                cent, lats = (ids_i[0], ids_j) if len(ids_i) == 1 else (ids_j[0], ids_i)
                for lat in lats:
                    elements.append((cent, lat, perfil))

        return nodes_full, elements

    def decode_population(self, pop):
        decoded = []
        for x in pop:
            key = hash(x.tobytes())
            if key not in self.global_decoded_cache:
                nodes, elems = self.decode_individual(x)
                if not self.validate_min_distance(nodes):
                    continue
                self.global_decoded_cache[key] = (nodes, elems)
            decoded.append(self.global_decoded_cache[key])
        return decoded 
   
    def initialize_individual(self) -> np.ndarray:
        """
        Gera um indivíduo válido.

        Saída:
        - x: array (dim,) de floats, contendo coords arredondadas e tube_vars iniciais.
        """
        while True:
            deltas = np.random.normal(size=(self.n, 3))
            deltas[self.is_central, 0] = 0.0
            norms = np.linalg.norm(deltas, axis=1, keepdims=True)
            radii = (np.random.rand(self.n, 1) ** (1/3)) * self.radius_opt
            coords = self.base_nodes + (deltas / norms) * radii
            coords = self.enforce_bounds(coords)
            coords = self.apply_central_alignment(coords)
            tube_vars = np.random.uniform(0, len(self.tipos_tubos), size=(self.n_tubes,))
            x = np.concatenate([coords.reshape(-1), tube_vars])
            nodes, _ = self.decode_individual(x)
            if self.validate_min_distance(nodes):
                return x

    def initialize_population(self):
        pop = []
        seen = set()
        while len(pop) < self.pop_size:
            x = self.initialize_individual()
            key = hash(x.tobytes())
            if key in seen:
                continue
            seen.add(key)
            pop.append(x)
        pop = np.array(pop)
        decoded = self.decode_population(pop)
        nodes_list, elems_list = zip(*decoded)
        results = self.parallel_evaluate(nodes_list, elems_list)
        fitness = np.array([res[0] for res in results])
        return pop, fitness

    def mutate_and_crossover(self, pop):
        new_pop = pop.copy()
        for i in range(self.pop_size):
            a, b, c = np.random.choice(np.delete(np.arange(self.pop_size), i), 3, replace=False)
            v = pop[a] + self.F * (pop[b] - pop[c])
            v[0:self.n * 3].reshape(self.n, 3)[self.is_central, 0] = 0.0
            j_rand = np.random.randint(self.dim)
            u = np.where((np.random.rand(self.dim) < self.CR) | (np.arange(self.dim) == j_rand), v, pop[i])
            if np.allclose(u, pop[i], atol=1e-4): continue
            key = hash(u.tobytes())
            if key not in self.global_decoded_cache:
                try:
                    nodes, _ = self.decode_individual(u)
                    if not self.validate_min_distance(nodes): continue
                    self.global_decoded_cache[key] = (nodes, _)
                except: continue
            new_pop[i] = u
        return new_pop
    
    def parallel_evaluate(self, nodes_list, elements_list):
        if self.use_parallel:
            with ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=init_worker,
                initargs=(self,)
            ) as executor:
                results = list(executor.map(evaluate, nodes_list, elements_list))
        else:
            results = [evaluate(nodes_list[i], elements_list[i]) for i in range(len(nodes_list))]
        return results  
      
    def optimize(self):
        """
        Executa o loop principal de otimização, com avaliação paralela.

        Retorna:
        - best_solution: tupla (nodes, elements)
        - best_cost: float
        - best_mass, best_KT, best_KF: métricas físicas
        - history: dicionário com histórico da evolução
        - duration: float (tempo total em segundos)
        """
        print("Otimização Iniciada")
        pop, fit = self.initialize_population()

        history = {
            'best_fit': [], 'avg_fit': [], 'std_dev': [],
            'best_mass': [], 'avg_mass': [],
            'best_KT': [], 'avg_KT': [],
            'best_KF': [], 'avg_KF': []
        }

        convergence_count = 0
        convergence_threshold = 3
        start_time = time.time()
        # Abre pool de processos apenas uma vez
        executor = None
        if self.use_parallel:
            executor = ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=init_worker,
                initargs=(self,)
            )

        # Loop de gerações
        for gen in range(1, self.max_gens + 1):
            t_gen = time.time()
            
            t1 = time.time()
            cand_pop = self.mutate_and_crossover(pop)
            t2 = time.time()

            decoded = self.decode_population(cand_pop)
            t3 = time.time()

            nodes_list, elems_list = zip(*decoded)
            t4 = time.time()
            # Avaliação (paralela ou sequencial)
            if executor:
                results = list(executor.map(evaluate, nodes_list, elems_list))
            else:
                results = [evaluate(n, e) for n, e in zip(nodes_list, elems_list)]
            t5 = time.time()

            t6 = time.time()
            new_pop = pop.copy()
            new_fit = fit.copy()

            new_mass = np.full(self.pop_size, np.nan)
            new_KT   = np.full(self.pop_size, np.nan)
            new_KF   = np.full(self.pop_size, np.nan)

            for i, (f, m, kt, kf) in enumerate(results):
                # Armazena métricas para análise estatística
                new_mass[i] = m
                new_KT[i]   = kt
                new_KF[i]   = kf

                # Substituição por seleção de DE
                if f <= fit[i]:
                    new_pop[i] = cand_pop[i]
                    new_fit[i] = f

            pop, fit = new_pop, new_fit

            best_idx = np.argmin(fit)
            best_fit = fit[best_idx]
            avg_fit  = np.mean(fit)
            std_dev  = float(np.std(pop))

            history['best_fit'].append(best_fit)
            history['avg_fit'].append(avg_fit)
            history['std_dev'].append(std_dev)
            history['best_mass'].append(new_mass[best_idx])
            history['avg_mass'].append(np.nanmean(new_mass))
            history['best_KT'].append(new_KT[best_idx])
            history['avg_KT'].append(np.nanmean(new_KT))
            history['best_KF'].append(new_KF[best_idx])
            history['avg_KF'].append(np.nanmean(new_KF))

            if std_dev < 0.02:
                convergence_count += 1
            else:
                convergence_count = 0

            t7 = time.time()

            print(
                f"[Gen {gen}] mutate: {t2 - t1:.2f}s | decode: {t3 - t2:.2f}s | eval: {t5 - t4:.2f}s | stats: {t7 - t6:.2f}s | total: {t7 - t_gen:.2f}s"
            )

            status = (f"Gen {gen}/{self.max_gens} — Best={best_fit:.4e} Std={std_dev:.4f}"
                    if convergence_count < convergence_threshold
                    else f"Convergência após {gen} gerações")
            print(status, end='\r')

            if convergence_count >= convergence_threshold:
                print()
                break

        # Fecha pool de processos
        if executor:
            executor.shutdown()

        duration = time.time() - start_time

        best_idx       = np.argmin(fit)
        best_solution  = self.decode_individual(pop[best_idx])
        best_cost      = fit[best_idx]
        best_mass      = history['best_mass'][-1]
        best_KT        = history['best_KT'][-1]
        best_KF        = history['best_KF'][-1]

        return best_solution, best_cost, best_mass, best_KT, best_KF, history, duration

    def plotar(self, individuo, save_path=None):
        """
        Plota a geometria do chassi evoluído em 3D.

        Entrada:
        - individuo: tupla (nodes, elements)
        - save_path: caminho opcional para salvar a imagem

        Saída:
        - exibe o plot e/ou salva a imagem se o caminho for fornecido.
        """
        nodes, elements = individuo
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')

        num_nodes = len(nodes)
        elements_valid = [(i, j, t) for i, j, t in elements if 0 <= i < num_nodes and 0 <= j < num_nodes]

        xs, ys, zs = zip(*nodes)
        ax.scatter(ys, xs, zs, s=25, c='black')
        for i, j, t in elements_valid:
            ni, nj = nodes[i], nodes[j]
            ax.plot([ni[1], nj[1]], [ni[0], nj[0]], [ni[2], nj[2]], 'b-')

        ax.set_xlabel('Y')
        ax.set_ylabel('X')
        ax.set_zlabel('Z')
        ax.set_box_aspect([3, 1, 2])
        plt.title("Chassi Evoluído")

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Visualização 3D salva em: {save_path}")

        plt.show()

    def plot_metrics(self, history, save_path=None, show=True):
        """
        Gera gráfico de evolução das métricas: massa, rigidez torcional e flexional.

        Entradas:
        - history: dict com as listas 'best_mass', 'best_KT', 'best_KF'.
        - save_path: caminho opcional para salvar o PNG.
        - show: bool, se True exibe o plot.
        """
        plt.figure(figsize=(10, 6))
        gens = range(len(history['best_mass']))
        
        # Configura eixo primário (massa)
        ax1 = plt.gca()
        ax1.plot(gens, history['best_mass'], 'b-', linewidth=2, label='Massa (kg)')
        ax1.set_xlabel('Geração')
        ax1.set_ylabel('Massa (kg)', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        
        # Configura eixo secundário (rigidezes)
        ax2 = ax1.twinx()
        ax2.plot(gens, history['best_KT'], 'r-', linewidth=2, label='Rigidez Torcional (KT)')
        ax2.plot(gens, history['best_KF'], 'g-', linewidth=2, label='Rigidez Flexional (KF)')
        ax2.set_ylabel('Rigidez (N·m/rad ou N/m)', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        
        plt.title('Evolução das Métricas do Melhor Indivíduo')
        
        # Unificar legendas
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
        
        plt.grid(True)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            print(f"Gráfico de métricas salvo em: {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()

    def plot_convergence(self, history, save_path=None, show=True):
        """
        Gera gráfico de convergência: fitness vs geração e std_dev.

        Entradas:
        - history: dict com 'best_fit', 'avg_fit', 'std_dev'.
        - save_path: caminho opcional para salvar o PNG.
        - show: bool, se True exibe o plot.
        """
        plt.figure(figsize=(10, 6))
        
        # Gráfico duplo eixo Y
        ax1 = plt.gca()
        ax2 = ax1.twinx()
        
        # Curva de fitness
        ax1.plot(history['best_fit'], 'b-', linewidth=2, label='Melhor Fitness')
        ax1.plot(history['avg_fit'], 'b--', alpha=0.7, label='Fitness Médio')
        ax1.set_ylabel('Fitness', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        ax1.set_yscale('log')
        
        # Curva de desvio padrão
        ax2.plot(history['std_dev'], 'r-', linewidth=2, label='Desvio Padrão')
        ax2.axhline(y=0.2, color='g', linestyle='--', label='Limite Convergência')
        ax2.set_ylabel('Desvio Padrão', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        
        plt.title('Progresso da Otimização')
        plt.xlabel('Geração')
        plt.grid(True)
        
        # Unificar legendas
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            print(f"Gráfico de convergência salvo em: {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()

    def plot_tradeoff(self, history, axises=('mass','KT'), save_path=None, show=True):
        x_key = {'mass':'best_mass','KT':'best_KT','KF':'best_KF'}[axises[0]]
        y_key = {'mass':'best_mass','KT':'best_KT','KF':'best_KF'}[axises[1]]
        x, y = history[x_key], history[y_key]
        plt.figure(figsize=(8,6))
        sc = plt.scatter(x, y, c=range(len(x)), cmap='viridis', s=50)
        plt.colorbar(sc, label='Geração')
        plt.xlabel(axises[0].upper())
        plt.ylabel(axises[1].upper())
        plt.title(f"Trade-off: {axises[0]} × {axises[1]}")
        plt.grid(True)
        if save_path: plt.savefig(save_path, dpi=300)
        if show: plt.show()
        else: plt.close()

    def plot_fitness_vs_metric(self, history, metric='mass', save_path=None, show=True):
        mkey = f"best_{metric}"
        x, y = history[mkey], history['best_fit']
        plt.figure(figsize=(8,6))
        sc = plt.scatter(x, y, c=range(len(x)), cmap='plasma', s=40)
        plt.colorbar(sc, label='Geração')
        plt.xlabel(metric.upper())
        plt.ylabel('Fitness')
        plt.title(f"Fitness × {metric.upper()}")
        plt.grid(True)
        if save_path: plt.savefig(save_path, dpi=300)
        if show: plt.show()
        else: plt.close()

    def plot_all(self, history, results_dir, show=False):
        os.makedirs(results_dir, exist_ok=True)
        self.plot_convergence(history, save_path=os.path.join(results_dir,'convergencia.png'), show=show)
        self.plot_metrics(history, save_path=os.path.join(results_dir,'evolucao_metricas.png'), show=show)
        #self.plot_tradeoff(history, axises=('mass','KT'), save_path=os.path.join(results_dir,'tradeoff_mass_KT.png'), show=show)
        #self.plot_tradeoff(history, axises=('mass','KF'), save_path=os.path.join(results_dir,'tradeoff_mass_KF.png'), show=show)
        #self.plot_fitness_vs_metric(history,'mass', save_path=os.path.join(results_dir,'fitness_vs_mass.png'), show=show)
        #self.plot_fitness_vs_metric(history,'KT',   save_path=os.path.join(results_dir,'fitness_vs_KT.png'), show=show)
        #self.plot_fitness_vs_metric(history,'KF',   save_path=os.path.join(results_dir,'fitness_vs_KF.png'), show=show)
        print(f"Todos os gráficos salvos em: {results_dir}")

    def save_solution(self, nodes, elements, file_path, fitness=None, mass=None, KT=None, KF=None, duration=None):
        """
        Salva relatório completo em TXT, incluindo duração e valores do melhor indivíduo.
        """
        with open(file_path, 'w') as f:
            f.write("SOLUÇÃO FINAL DO CHASSI\n")
            f.write("="*60 + "\n\n")
            f.write(f"{self.max_gens} Gerações, {self.pop_size} individuos por população \n")
            if duration is not None:
                hrs, rem = divmod(duration, 3600)
                mins, secs = divmod(rem, 60)
                f.write(f"Duração da Otimização: {int(hrs)}h {int(mins)}m {secs:.1f}s\n")
            if fitness is not None:
                f.write(f"Melhor Fitness: {fitness:.6e}\n")
            if mass is not None:
                f.write(f"Massa: {mass:.3f} kg\n")
            if KT is not None and KF is not None:
                f.write(f"Rigidez Torcional (KT): {KT:.3e} N·m/rad\n")
                f.write(f"Rigidez Flexional (KF): {KF:.3e} N/m\n")
            f.write("\nNÓS (X, Y, Z):\n")
            for i, node in enumerate(nodes):
                f.write(f"{i:3d}: {node[0]:6.3f}, {node[1]:6.3f}, {node[2]:6.3f}\n")
            f.write("\nELEMENTOS (Conexões):\n")
            for i, (ni, nj, tp) in enumerate(elements):
                f.write(f"{i:3d}: Nó {ni} - Nó {nj} | Perfil: {tp}\n")
        print(f"Relatório completo salvo em: {file_path}")

    def export_solution(self, nodes: np.ndarray, elements: list, directory: str, filename: str = "melhor_solucao"):
        """
        Exporta a solução para arquivos CSV dentro do diretório especificado.
        
        Args:
            nodes: array (N,3) com coordenadas dos nós
            elements: lista de tuplas (i, j, perfil)
            directory: diretório onde os arquivos serão salvos
            filename: nome base para os arquivos (sem extensão)
        """
        # Garante que o diretório existe
        os.makedirs(directory, exist_ok=True)
        
        # Cria caminhos completos
        nodes_path = os.path.join(directory, f"{filename}_nodes.csv")
        elements_path = os.path.join(directory, f"{filename}_elements.csv")
        
        # Exportar nós
        with open(nodes_path, 'w') as f_nodes:
            f_nodes.write("index,x,y,z\n")
            for i, node in enumerate(nodes):
                f_nodes.write(f"{i},{node[0]:.6f},{node[1]:.6f},{node[2]:.6f}\n")
        
        # Exportar elementos
        with open(elements_path, 'w') as f_elements:
            f_elements.write("node_i,node_j,perfil\n")
            for elem in elements:
                f_elements.write(f"{elem[0]},{elem[1]},{elem[2]}\n")
        
        print(f"Solucao exportada para:")
        print(f" - Nodes: {nodes_path}")
        print(f" - Elements: {elements_path}")

def find_new_index(old_index, nodes):
    """
    Calcula o indice do nó após passar pela otimização
    Retorna o novo indice do nó e do seu espelhado correspondente
    Caso seja um nó central retorna apenas o novo indice do nó
    
    Entradas:
    - old_index: Índice do nó no base nodes
    - nodes: lista de nós pós otimização
    Retorno:
    - new_index: Índice do nó após passar pela otimização
    - mirrored_index: índice do nó espelhado correspondente
    - new_central_index: novo index se for um nó central
    """
    is_central = np.isclose(nodes[:, 0], 0.0)
    
    first_central = np.sum(~is_central)

    if (old_index*2)>=first_central:
        new_central_index=first_central+(old_index-(first_central//2))
        return new_central_index
    else:
        new_index=old_index*2
        mirrored_index=new_index+1
        return new_index,mirrored_index
    
def penalidades_geometricas(nodes, elements):

    penalidade = 0

    # Penalidade em relação ao FH e FHB
    def penalidade_fh_fhb(nodes):
        pen = 0
        fronthoop_node = nodes[find_new_index(20, nodes)]                          #declara o nó do fronthoop com indice novo
        fhb_node = nodes[find_new_index(5, nodes)[0]]                              #declara o nó de um front hoop bracing com indice novo
        dist_fh_fhb = fronthoop_node[2] - fhb_node[2]                              #declara a distância no eixo z entre esses dois nós 
        if dist_fh_fhb > 0.05:                                                     #condição retirada do regulamento
            pen += ((dist_fh_fhb - 0.05) ** 2 ) * 1e6                              #aplicação de penalidade

        return pen

    # Penalidade em relação ao MH e o MHB
    def penalidade_mh_mhb(nodes):
        pen = 0
        mainhoop_node = nodes[find_new_index(16, nodes)]                                        #declara o nó do main hoop com indice novo
        mhb_node = nodes[find_new_index(14, nodes)[0]]                                          #declara o nó do main hoop bracing com indice novo não espelhado
        deltax_mh_mhb = mainhoop_node[0] - mhb_node[0]                                          #diferença das coordenadas "x" em ambos os nós
        deltay_mh_mhb = mainhoop_node[1] - mhb_node[1]                                          #diferença das coordenadas "y" em ambos os nós
        deltaz_mh_mhb = mainhoop_node[2] - mhb_node[2]                                          #diferença das coordenadas "z" em ambos os nós
        dist_mh_mhb = np.sqrt(deltax_mh_mhb ** 2 + deltay_mh_mhb ** 2 + deltaz_mh_mhb ** 2)     #declara a distância entre os dois nós      
        if dist_mh_mhb > 0.16:                                                                  #condição retirada do regulamento
            pen += ((dist_mh_mhb - 0.16) ** 2) * 1e6                                            #aplicação de penalidade
        
        return pen
    
    # Penalidade ângulo entre o Main Hoop e o Main Hoop Bracing
    def penalidade_angulo_mh_mhb(nodes):
        pen = 0
        x_porcao_mh = nodes[find_new_index(14, nodes)[0]][0] - nodes[find_new_index(6, nodes)[0]][0]                                                                  #coordenada x do vetor formado pelos nós do elemento da porção do mainhoop analisada
        y_porcao_mh = nodes[find_new_index(14, nodes)[0]][1] - nodes[find_new_index(6, nodes)[0]][1]                                                                  #coordenada y do vetor formado pelos nós do elemento da porção do mainhoop analisada
        z_porcao_mh = nodes[find_new_index(14, nodes)[0]][2] - nodes[find_new_index(6, nodes)[0]][2]                                                                  #coordenada z do vetor formado pelos nós do elemento da porção do mainhoop analisada
        x_mhb = nodes[find_new_index(14, nodes)[0]][0] - nodes[find_new_index(15, nodes)[0]][0]                                                                       #coordenada x do vetor formado pelos nós do elemento do Main Hoop Bracing
        y_mhb = nodes[find_new_index(14, nodes)[0]][1] - nodes[find_new_index(15, nodes)[0]][1]                                                                       #coordenada y do vetor formado pelos nós do elemento do Main Hoop Bracing
        z_mhb = nodes[find_new_index(14, nodes)[0]][2] - nodes[find_new_index(15, nodes)[0]][2]                                                                       #coordenada z do vetor formado pelos nós do elemento do Main Hoop Bracing
        vetor_porcao_mh = (x_porcao_mh, y_porcao_mh, z_porcao_mh)                                                                                               #vetor formado pelos nós do elemento da porção do mainhoop analisada
        vetor_mhb = (x_mhb, y_mhb, z_mhb)                                                                                                                       #vetor formado pelos nós do elemento do Main Hoop Bracing
        modulo_vetor_porcao_mh = np.sqrt(vetor_porcao_mh[0] ** 2 + vetor_porcao_mh[1] ** 2 + vetor_porcao_mh[2] ** 2 )                                          #módulo do vetor formado pelos nós do elemento da porção do mainhoop analisada
        modulo_vetor_mhb = np.sqrt(vetor_mhb[0] ** 2 + vetor_mhb[1] ** 2 + vetor_mhb[2] ** 2 )                                                                  #módulo do vetor formado pelos nós do elemento do Main Hoop Bracing
        produto_escalar_mh_porcao_e_mhb = (vetor_porcao_mh[0] * vetor_mhb[0]) + (vetor_porcao_mh[1] * vetor_mhb[1]) + (vetor_porcao_mh[2] * vetor_mhb[2])       #produto escalar entre os dois vetores criados
        cos_theta_mh_mhb = produto_escalar_mh_porcao_e_mhb / (modulo_vetor_porcao_mh * modulo_vetor_mhb)                                                        #valor do cosseno do ângulo formado pelos dois vetores
        theta_mh_mhb = np.degrees(np.acos(cos_theta_mh_mhb))                                                                                                    #valor do ângulo formado pelos dois vetores
        if theta_mh_mhb < 30:                                                                                                                                   #condição retirada do regulamento
            pen += ((theta_mh_mhb - 30) ** 2) * 1e6                                                                                                             #aplicação da penalidade
    
        return pen

    # Penalidade ângulo com a vertical da parte do Front Hoop que fica acima da Upper Side Impact Structure
    def penalidade_angulo_vetical_fh(nodes):
        pen = 0
        x_porcao_fh = nodes[find_new_index(5, nodes)[0]][0] - nodes[find_new_index(4, nodes)[0]][0]                                          #coordenada x do vetor formado pelos nós do elemento da porção do mainhoop analisada
        y_porcao_fh = nodes[find_new_index(5, nodes)[0]][1] - nodes[find_new_index(4, nodes)[0]][1]                                          #coordenada y do vetor formado pelos nós do elemento da porção do mainhoop analisada
        z_porcao_fh = nodes[find_new_index(5, nodes)[0]][2] - nodes[find_new_index(4, nodes)[0]][2]                                          #coordenada z do vetor formado pelos nós do elemento da porção do mainhoop analisada
        vetor_porcao_fh = (x_porcao_fh, y_porcao_fh, z_porcao_fh)                                                                      #vetor formado pelos nós que formam a porção do front hoop analisada
        modulo_vetor_porcao_fh = np.sqrt(vetor_porcao_fh[0] ** 2 + vetor_porcao_fh[1] ** 2 + vetor_porcao_fh[2] **2)                   #módulo do vetor formado pelos nós que formam a porção do front hoop analisada
        produto_escalar_porcao_fh_e_vertical = vetor_porcao_fh[2]                                                                      #produto escalar entre o vetor formado pelos nós que formam a porção do front hoop analisada e o versor da vertical
        cos_theta_fh_porcao_vertical = produto_escalar_porcao_fh_e_vertical / modulo_vetor_porcao_fh                                   #cosseno ângulo formado entre a porção do front hoop analisada e a vertical
        theta_fh_porcao_vertical = np.degrees(np.acos(cos_theta_fh_porcao_vertical))                                                   #cálculo do ângulo através do cosseno
        if theta_fh_porcao_vertical > 20:                                                                                              #condição retirada do regulamento
            pen += ((theta_fh_porcao_vertical - 20) ** 2) * 1e6                                                                        #aplicação da penalidade
    
        return pen

    # Penalidade ângulo com a vertical da parte do Main Hoop que fica acima do ponto que o conecta ao Upper Side Impact Tube 
    def penalidade_angulo_vertical_mh(nodes):
        pen = 0
        x_porcao_mh = nodes[find_new_index(14, nodes)[0]][0] - nodes[find_new_index(6, nodes)[0]][0]                                          #coordenada x do vetor formado pelos nós do elemento da porção do mainhoop analisada
        y_porcao_mh = nodes[find_new_index(14, nodes)[0]][1] - nodes[find_new_index(6, nodes)[0]][1]                                          #coordenada y do vetor formado pelos nós do elemento da porção do mainhoop analisada
        z_porcao_mh = nodes[find_new_index(14, nodes)[0]][2] - nodes[find_new_index(6, nodes)[0]][2]                                          #coordenada z do vetor formado pelos nós do elemento da porção do mainhoop analisada
        vetor_porcao_mh = (x_porcao_mh, y_porcao_mh, z_porcao_mh)                                                                       #vetor formado pelos nós do elemento da porção do mainhoop analisada
        modulo_vetor_porcao_mh = np.sqrt(vetor_porcao_mh[0] ** 2 + vetor_porcao_mh[1] ** 2 + vetor_porcao_mh[2] ** 2 )                  #módulo do vetor formado pelos nós do elemento da porção do mainhoop analisada
        produto_escalar_porcao_mh_e_vertical = vetor_porcao_mh[2]                                                                       #produto escalar do vetor formado pelo elemento da porção do Main Hoop trabalhada com o versor da vertical
        cos_theta_mh_porcao_vertical = produto_escalar_porcao_mh_e_vertical / modulo_vetor_porcao_mh                                    #cosseno do ângulo formado entre este vetor mencionado e a vertical
        theta_mh_porcao_vertical = np.degrees(np.acos(cos_theta_mh_porcao_vertical))                                                    #ângulo formado entre este vetor mencionado e a vertical
        if theta_mh_porcao_vertical > 10:                                                                                               #condição retirada do regulamento
            pen += ((theta_mh_porcao_vertical -10) ** 2) * 1e6                                                                          #aplicação da penalidade

        return pen
    
     # Penalidade ângulo com a vertical da parte do Front Bulkhead gaarantindo a retidão
    def penalidade_angulo_vertical_fb(nodes):
        pen = 0
        x_porcao_fb = nodes[find_new_index(0, nodes)[0]][0] - nodes[find_new_index(1, nodes)[0]][0]                                          #coordenada x do vetor formado pelos nós do elemento da porção do front bulkhead analisada
        y_porcao_fb = nodes[find_new_index(0, nodes)[0]][1] - nodes[find_new_index(1, nodes)[0]][1]                                          #coordenada y do vetor formado pelos nós do elemento da porção do front bulkhead analisada
        z_porcao_fb = nodes[find_new_index(0, nodes)[0]][2] - nodes[find_new_index(1, nodes)[0]][2]                                          #coordenada z do vetor formado pelos nós do elemento da porção do front bulkhead analisada
        vetor_porcao_fb = (x_porcao_fb, y_porcao_fb, z_porcao_fb)                                                                       #vetor formado pelos nós do elemento da porção do front bulkhead analisada
        modulo_vetor_porcao_fb = np.sqrt(vetor_porcao_fb[0] ** 2 + vetor_porcao_fb[1] ** 2 + vetor_porcao_fb[2] ** 2 )                  #módulo do vetor formado pelos nós do elemento da porção do front bulkhead analisada
        produto_escalar_porcao_fb_e_vertical = vetor_porcao_fb[2]                                                                       #produto escalar do vetor formado pelo elemento da porção do front bulkhead trabalhada com o versor da vertical
        cos_theta_fb_porcao_vertical = produto_escalar_porcao_fb_e_vertical / modulo_vetor_porcao_fb                                    #cosseno do ângulo formado entre este vetor mencionado e a vertical
        theta_fb_porcao_vertical = np.degrees(np.acos(cos_theta_fb_porcao_vertical))                                                    #ângulo formado entre este vetor mencionado e a vertical
        if theta_fb_porcao_vertical != 0:                                                                                               #condição retirada do regulamento
            pen += ((theta_fb_porcao_vertical -0) ** 2) * 1e6                                                                          #aplicação da penalidade

        return pen   

     # Penalidade ângulo com a vertical da parte do Rear Bulkhead gaarantindo a retidão
    def penalidade_angulo_vertical_rb(nodes):
        pen = 0
        x_porcao_rb = nodes[find_new_index(12, nodes)[0]][0] - nodes[find_new_index(13, nodes)[0]][0]                                          #coordenada x do vetor formado pelos nós do elemento da porção do rear bulkhead analisada
        y_porcao_rb = nodes[find_new_index(12, nodes)[0]][1] - nodes[find_new_index(13, nodes)[0]][1]                                          #coordenada y do vetor formado pelos nós do elemento da porção do rear bulkhead analisada
        z_porcao_rb = nodes[find_new_index(12, nodes)[0]][2] - nodes[find_new_index(13, nodes)[0]][2]                                          #coordenada z do vetor formado pelos nós do elemento da porção do rear bulkhead analisada
        vetor_porcao_rb = (x_porcao_rb, y_porcao_rb, z_porcao_rb)                                                                       #vetor formado pelos nós do elemento da porção do rear bulkhead analisada
        modulo_vetor_porcao_rb = np.sqrt(vetor_porcao_rb[0] ** 2 + vetor_porcao_rb[1] ** 2 + vetor_porcao_rb[2] ** 2 )                  #módulo do vetor formado pelos nós do elemento da porção do rear bulkhead analisada
        produto_escalar_porcao_rb_e_vertical = vetor_porcao_rb[2]                                                                       #produto escalar do vetor formado pelo elemento da porção do rear bulkhead trabalhada com o versor da vertical
        cos_theta_rb_porcao_vertical = produto_escalar_porcao_rb_e_vertical / modulo_vetor_porcao_rb                                    #cosseno do ângulo formado entre este vetor mencionado e a vertical
        theta_rb_porcao_vertical = np.degrees(np.acos(cos_theta_rb_porcao_vertical))                                                    #ângulo formado entre este vetor mencionado e a vertical
        if theta_rb_porcao_vertical != 0:                                                                                               #condição retirada do regulamento
            pen += ((theta_rb_porcao_vertical -0) ** 2) * 1e6                                                                          #aplicação da penalidade

        return pen  

    # Penalidade altura mínima da Side Impact Structure (aparentemente, essa regra é para monocoque)
    '''def penalidade_altura_min_SIS(nodes):
        pen = 0
        z_zone_impact_bottom_back = nodes[find_new_index(7, nodes)[0]][2]                                              #coordenada vertical do ponto mais baixo da parte posterior da side impact structure
        z_zone_impact_top_back = nodes[find_new_index(6, nodes)[0]][2]                                                 #coordenada vertical do ponto mais alto da parte posterior da side impact structure
        z_zone_impact_bottom_front = nodes[find_new_index(3, nodes)[0]][2]                                             #coordenada vertical do ponto mais baixo da parte frontal da side impact structure
        z_zone_impact_top_front = nodes[find_new_index(4, nodes)[0]][2]                                                #coordenada vertical do ponto mais alto da parte frontal da side impact structure
        dist_bottom_top_back = z_zone_impact_top_back - z_zone_impact_bottom_back                                   #distância vertical entre os extremos da parte posterior desta estrutura
        dist_bottom_top_front = z_zone_impact_top_front - z_zone_impact_bottom_front                                #distância vertical entre os extremos da parte frontal desta estrutura
        if dist_bottom_top_back < 0.29 and dist_bottom_top_front < 0.29:                                            #condição retirada do regulamento
            pen += ((dist_bottom_top_front - 0.29) ** 2 + (dist_bottom_top_back - 0.29) ** 2) * 1e6                 #aplicação da penalidade

        return pen'''
    
    # Penalidade para posição do Upper Side Impact Member
    def penalidade_posicao_upper_SI_member(nodes):
        pen = 0
        upper_SI_node_front_z = nodes[find_new_index(4, nodes)[0]][2]
        upper_SI_node_back_z = nodes[find_new_index(6, nodes)[0]][2]
        lowest_point_SI_z =  nodes[find_new_index(7, nodes)[0]][2]
        dist_lowestpoint_front_usim = upper_SI_node_front_z - lowest_point_SI_z
        dist_lowestpoint_back_usim = upper_SI_node_back_z - lowest_point_SI_z

        if dist_lowestpoint_front_usim * 1000 not in range(240, 321) and dist_lowestpoint_back_usim * 1000 not in range(240, 321):
            pen += ((dist_lowestpoint_front_usim) ** 2 + (dist_lowestpoint_back_usim) ** 2) * 1e6


        return pen
    
    # Controle das penalidades
    penalidade += penalidade_fh_fhb(nodes)
    penalidade += penalidade_mh_mhb(nodes)
    penalidade += penalidade_angulo_vetical_fh(nodes)
    penalidade += penalidade_angulo_vertical_mh(nodes)
    penalidade += penalidade_angulo_mh_mhb(nodes)
    #penalidade += penalidade_altura_min_SIS(nodes)
    penalidade += penalidade_posicao_upper_SI_member(nodes)
    penalidade += penalidade_angulo_vertical_fb(nodes)
    penalidade += penalidade_angulo_vertical_rb(nodes)

    return penalidade



def evaluate(nodes,elements) -> float:
    """
    Avalia o custo de um indivíduo.

    Entrada:
    - nodes
    - elements
    Saída:
    - float, valor de penalidade (fitness).

    Processo:
    1. Decodifica para nodes e elements.
    2. Checa distância mínima.
    3. Monta e analisa pela classe Estrutura (FEA).
    4. Calcula penalidade via penalidade_chassi.
    """
    try:
        t0 = time.perf_counter()

        penalty = 0
        penalty += penalidades_geometricas(nodes, elements)

        # Instanciamento da estrutura
        estrutura = Estrutura(elements, nodes)
        t1 = time.perf_counter()

        estrutura.matrizes_global()
        t2 = time.perf_counter()

        fixed = list(range(6))
        Fg = np.zeros(estrutura.num_dofs)
        for i in [7, 8, 21, 22]:
            Fg[i * 6 + 2] = (65 * 9.81) / 4

        _, _, frequencies = estrutura.modal_analysis()
        t3 = time.perf_counter()

        disp = estrutura.static_analysis(Fg, fixed)
        t4 = time.perf_counter()

        strain = estrutura.compute_strain(disp)
        stresses = estrutura.compute_stress(strain)
        von = estrutura.compute_von_mises(stresses)
        t5 = time.perf_counter()

        massa = estrutura.mass()
        LFnode,RFnode = find_new_index(3,nodes) #Left front suspension node index (3)
        LRnode,RRnode = find_new_index(11,nodes) #Left rear suspension node index (11)
        LCnode,RCnode = find_new_index(7,nodes) #Left central node index (7)
        KT, KF = estrutura.compute_Kf_Kt(estrutura.K_global,LFnode,LRnode,RFnode,RRnode,LCnode,RCnode)
        print(KT,KF)
        t6 = time.perf_counter()

        penalty += penalidade_chassi(KT, KF, massa, von, frequencies)
        t7 = time.perf_counter()

        print(f"""[EVALUATE]
        - Instância:     {t1 - t0:.4f}s
        - Montagem K/M:  {t2 - t1:.4f}s
        - Modal:         {t3 - t2:.4f}s
        - Estática:      {t4 - t3:.4f}s
        - Tensão/Von:    {t5 - t4:.4f}s
        - Massa/Rigidez: {t6 - t5:.4f}s
        - Penalidade:    {t7 - t6:.4f}s
        - Total:         {t7 - t0:.4f}s
        """)

        return penalty, massa, KT, KF

    except Exception:
        traceback.print_exc()
        return 1e9, 1e9, 0.0, 0.0

def penalidade_chassi(KT, KF, massa, tensoes, frequencias):
    """
    Calcula penalidade total do chassi.

    Entradas:
    - KT, KF: floats de rigidezes.
    - massa: float.
    - tensoes: array de floats.
    - frequencias: array de frequencias naturais.
    Saída:
    - penalidade_total: float.
    """

    # Limites e parâmetros
    KT_min = 1e7                # Rigidez torcional mínima (N·m/rad)
    KF_min = 1e6                # Rigidez flexão mínima (N/m)
    massa_ideal = 23            # Massa alvo (kg)
    K_mola = 5e5                # Constante da mola do amortecedor (N/m)
    tensao_adm = 250e6          # Tensão admissível do material (Pa)
    alpha = 0.5                 # Fator de sensibilidade exponencial
    beta = 10                   # Fator de escala logarítmica
    freq_motor = 4800           # Rotação do motor dada em (RPM) e convertida para (Hz)
    penalidade_total = 0

    # 1. Rigidez Torcional (Função Exponencial)
    if KT < KT_min:
        deficit = (KT_min - KT) / KT_min
        # Penalidade cresce exponencialmente com o déficit
        penalidade_total += np.exp(alpha * deficit) - 1

    # 2. Rigidez em Flexão (Função Logarítmica)
    if KF < KF_min:
        deficit = (KF_min - KF) / KF_min
        # Penalidade logarítmica: suave para pequenas violações, forte para grandes
        penalidade_total += beta * np.log(1 + deficit)

    # 3. Massa (Função Híbrida)
    if massa > massa_ideal:
        excesso = (massa - massa_ideal) / massa_ideal
        # Combina resposta linear inicial com crescimento exponencial
        penalidade_total += excesso + np.exp(alpha * excesso) - 1

    # 4. Compatibilidade com Mola (Lógica Aprimorada)
    ratio_KT = K_mola / KT if KT > 0 else float('inf')
    ratio_KF = K_mola / KF if KF > 0 else float('inf')
    
    if ratio_KT > 25 or ratio_KF > 25:
        # Penalidade proporcional ao nível de incompatibilidade
        violacao = max(ratio_KT/25, ratio_KF/25) - 1
        penalidade_total += 100 * violacao**2

    # 5. Tensões (Abordagem Baseada em Risco)
    tensao_max = max(tensoes)
    if tensao_max > tensao_adm:
        # Penalidade exponencial para tensões acima do admissível
        excesso = (tensao_max - tensao_adm) / tensao_adm
        penalidade_total += np.exp(5 * excesso) - 1
    
    # Penalidade por distribuição desigual de tensões (logarítmica)
    razao_tensoes = np.ptp(tensoes) / np.mean(tensoes) if np.mean(tensoes) > 0 else 0
    penalidade_total += np.log(1 + razao_tensoes)

    # 7. Frequências Naturais em Zona de Ressonância com o Motor (Função exponencial)
    freq_crit_min = (0.95*freq_motor)/60        #Frequência crítica mínima convertida para Hz
    freq_crit_max = (1.05*freq_motor)/60        #Frequência crítica máxima convertida para Hz

    for f in frequencias:
        # Penalidade exponencial mais severa para frequências próximas a do motor
        if freq_crit_min <= f <= freq_crit_max:
            severidade = np.exp(alpha * (1 - abs((f - freq_motor) / freq_motor)))
            penalidade_total += 100 * severidade 
    return penalidade_total * 100  # Fator de escala global

if __name__ == "__main__":                                                      #centrais começando depois dos nós passíveis de espelhamento
    nodes = np.array([[-0.181,  0.000,  0.360],             #00
    [-0.181,  0.000,  0.050],                               #01
    [-0.280,  0.275,  0.240],                               #02
    [-0.285,  0.495,  0.045],                               #03
    [-0.285,  0.555,  0.270],                               #04
    [-0.212,  0.555,  0.550],                               #05
    [-0.293,  1.370,  0.250],                               #06
    [-0.268,  1.350,  0.000],                               #07
    [-0.268,  1.495,  0.015],                               #08
    [-0.302,  1.670,  0.240],                               #09
    [-0.271,  1.665,  0.030],                               #10
    [-0.271,  1.835,  0.070],                               #11
    [-0.183,  2.015,  0.285],                               #12
    [-0.183,  2.015,  0.060],                               #13
    [-0.170,  1.400,  0.965],                               #14
    [-0.293,  1.950,  0.250],                               #15
    [ 0.000,  1.410,  1.105],                               #16
    [0.000,  0.000,  0.360],                                #conection between the 2 sides 
    [0.000,  0.000,  0.050],                                #conection between the 2 sides
    [0.000,  0.495,  0.045],                                #conection between the 2 sides
    [0.000,  0.555,  0.550],                                #conection between the 2 sides
    [0.000,  1.350,  0.000],                                #conection between the 2 sides
    [0.000,  1.495,  0.015],                                #conection between the 2 sides
    [0.000,  1.665,  0.030],                                #conection between the 2 sides
    [0.000,  1.835,  0.070],                                #conection between the 2 sides
    [0.000,  2.015,  0.285],                                #conection between the 2 sides
    [0.000,  2.015,  0.060],                                #conection between the 2 sides
    [0.000,  1.370,  0.250]])                               #conection between the 2 sides 
                                
    connections = [(0,1)  ,(0,2)  ,(1,2)  ,(0,5)  ,(2,5)  ,(2,4)  ,(2,3)  ,(1,3)  ,
                (3,7)  ,(3,4)  ,(4,7)  ,(4,6)  ,(5,6)  ,(7,6)  ,(7,8)  ,(6,8)  ,
                (6,9)  ,(8,9)  ,(8,10) ,(9,10) ,(9,11) ,(10,11),(11,15),(15,13),
                (11,12),(11,13),(12,13),(12,15),(15,14),(14,16),(14,6) ,(9,15) ,
                (4,5)  ,(0,17) ,(1,18) ,(3,19) ,(5,20) ,(7,21) ,(8,22) ,(10,23),
                (11,24),(12,25),(6,27),(13,26)]

    indices = [0,1,3,4,5,6,7,8,11,16,14,15,20]

    # Criar diretório para resultados
    timestamp = datetime.now().strftime("%Y-%m-%d %H%M")
    max_gen = 40
    pop_size = 10
    otimizador = ChassisDEOptimizer(
        base_nodes=nodes,
        base_connections=connections,
        mandatory_indices=indices,
        pop_size=pop_size,
        max_generations=max_gen,
        use_parallel=True,     # ou False, se quiser forçar modo serial
        n_workers=4            # defina quantos núcleos usar
    )
    
    results_dir = f"Resultados_Otimizacao_LOCAL__{timestamp}_{max_gen}GEN_{pop_size}POP"
    os.makedirs(results_dir, exist_ok=True)

    best_indiv, best_cost, best_mass, best_KT, best_KF, history,duration = otimizador.optimize()
    print(f"\nRESULTADOS FINAIS:")
    print(f"Melhor custo: {best_cost:.6e}")
    print(f"Massa: {best_mass:.2f} kg")
    print(f"Rigidez Torcional (KT): {best_KT:.2f} N·m/rad")
    print(f"Rigidez Flexional (KF): {best_KF:.2e} N/m")
    nodes_final, elements_final = best_indiv

    # Salvar solução em arquivo TXT
    solution_path = os.path.join(results_dir, "solucao_final.txt")
    otimizador.save_solution(
        nodes_final, elements_final, solution_path,
        fitness=best_cost, mass=best_mass, KT=best_KT, KF=best_KF, duration=duration)
    
    # Geração de todos os gráficos de uma única vez
    otimizador.plot_all(history, results_dir, show=False)
    # visualização 3D final
    otimizador.plotar(best_indiv, save_path=os.path.join(results_dir,'visualizacao_3d.png'))

    # Exportar para arquivos
    otimizador.export_solution(
        nodes_final, 
        elements_final, 
        directory=results_dir,
        filename="solucao_otimizada")
    