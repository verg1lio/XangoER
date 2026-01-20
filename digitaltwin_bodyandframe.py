import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import pandas as pd
import sys
import os
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=200, suppress=True)

class Estrutura:


    def __init__(self, elements, nodes, m, Id, Ip):
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
            - Alexandre Duque Gondim Pires;
            - Arthur Miguel Costa Pimentel Brandão;
            - Enzo De Andrade Magalhães;
            - João Felipe Sanders Pinto;
            - Vergílio Torezan Silingardi Del Claro.
        """
        self.elements = elements                                             #Matriz de elementos conectados
        self.num_elements = len(elements)                                    #Número de elementos
        self.nodes = nodes                                                   #Matriz de nós com suas posições
        self.num_nodes = len(nodes)                                          #Número total de nós
        self.massa = 30                                                      #Massa do carro (Kg)
        self.momento_inercia_direcao = Id                                    #Momento de inércia em relação à direção (kg.m^2)
        self.momento_inercia_plano = Ip                                      #Momento de inércia em relação ao plano (kg.m^2)
        self.num_dofs_per_node = 6                                           #6 graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node              #Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))             #Matriz de rigidez global
        self.M_global = np.zeros((self.num_dofs, self.num_dofs))             #Matriz de massa global
        self.num_modes = 6                                                  #Número de modos de vibração a serem retornados


    def calcular_comprimento(self, element):    
        """
        Calculates the length of an element based on node coordinates.
        Inputs:
            - element: tuple (start node index, end node index).
        Outputs:
            - element length (float).
        """                    
        node1, node2 = element
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)


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
            node_start, node_end = element
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
        E = 210e9   	#Modulo de Young (Pa)
        I = 1.6667e-5 	#Momento de inercia (m^4)
        G = 81.2e9  	#Modulo de Cisalhamento(Pa)
        A = 0.0125	    #Área da seção do elemento (m^2)	
        J = I/2     	#Momento polar de inércia (m^4) 
        kappa=0.9       #Fator de correção para cisalhamento 
        L_e = self.calcular_comprimento(element)
        Phi = (12 * E * I) / (kappa * G * A * L_e**2)
        rho = 7850  # kg/m^3

        c1 = E * A / L_e
        c2 = G * J / L_e
        c3 = E * I / L_e**3             #Euler-Bernoulli
        c4 = (E*I)/(L_e**3*(1+Phi))     #Timoshenko
        t1 = (4+Phi)
        t2 = (2-Phi)
        d1 = rho*A*L_e
        d2 = (I*L_e)/6
        d3 = (rho*A*L_e)/420

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
                self.K_global[index, index] = 10**10                # Um valor suficientemente grande para simular um engaste 


    def matrizes_global(self):
        """
        Assembles the global stiffness and mass matrices.
        Inputs: None (uses class attributes).
        Outputs:
            - K_global: global stiffness matrix.
            - M_global: global mass matrix.
        """
        for element in self.elements:
            node1, node2 = element
            k_e, m_e = self.element(element)
            # DOFs associados ao elemento            
            dofs = [6 * node1, 6 * node1 + 1, 6 * node1 + 2, 6 * node1 + 3, 6 * node1 + 4, 6 * node1 + 5,
                    6 * node2, 6 * node2 + 1, 6 * node2 + 2, 6 * node2 + 3, 6 * node2 + 4, 6 * node2 + 5]

            # Atualizando as matrizes globais
            self.K_global[np.ix_(dofs, dofs)] += k_e
            self.M_global[np.ix_(dofs, dofs)] += m_e

        # self.aplicar_engastes([0, 2, 4, 5], [0, 1, 2, 3, 4, 5])                             #Por enquanto não estaremos considerando engastes
        pd.DataFrame(self.K_global).to_csv('Matriz_Global_Rigidez.csv', index=True, header=True)
        pd.DataFrame(self.M_global).to_csv('Matriz_Global_Massa.csv', index=True, header=True)        

        # print (self.K_global)
        # print (self.M_global)

        plt.figure(figsize=(6, 6))
        plt.spy(self.K_global, markersize=10)  # Adjust markersize for visibility
        plt.title("Spy Plot of the Kg")
        plt.xlabel("Columns")
        plt.ylabel("Rows")
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.show()

        plt.figure(figsize=(6, 6))
        plt.spy(self.M_global, markersize=10)  # Adjust markersize for visibility
        plt.title("Spy Plot of the Mg")
        plt.xlabel("Columns")
        plt.ylabel("Rows")
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.show()

        return self.K_global,self.M_global


    def shape_fun(self, F_flexao1, F_flexao2, F_axial,F_torcao): 
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
        E = 2.1e11  	#Modulo de Young (Pa)
        I = 1.6667e-5 	#Momento de inercia (m^4)
        G = 81.2e9  	#Modulo de Cisalhamento (Pa)
        A= 0.0125	    #Área da seção do elemento (m^2)	
        J = I/2     	#Momento polar de inércia (m^4)
        KF_total = 0
        KT_total = 0
        KF_elements = []
        KT_elements = [] 
        torcao, deformacao, flexao1, flexao2, flexao3 = [], [], [], [], []
        for element in self.elements:
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
    
    def static_analysis(self, K_global, F_global, fixed_dofs):
        """
        Perform static analysis by solving Ku = F with boundary conditions.

        Parameters:
            K_global (ndarray): Global stiffness matrix (N x N).
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
        n_dofs = K_global.shape[0]

        # Create a mask for free DOFs (DOFs not constrained)
        free_dofs = np.array([i for i in range(n_dofs) if i not in fixed_dofs])

        # Reduce the stiffness matrix and force vector
        K_reduced = K_global[np.ix_(free_dofs, free_dofs)]
        F_reduced = F_global[free_dofs]

        # Solve for displacements at free DOFs
        # USE "linalg.lstsq" FOR NEAR SINGULAR MATRICES (ALL OF THEM)
        u_reduced = np.linalg.lstsq(K_reduced, F_reduced, rcond=None)[0]  # Get only the solution vector

        # Construct full displacement vector
        displacements = np.zeros(n_dofs)
        displacements[free_dofs] = u_reduced
        
        return displacements

    def calcular_B_Elementar(self, displacements, element):
        """
        Calcula a matriz B de um elemento individual.
        """
        L_e = self.calcular_comprimento(element)
        node1, node2 = element

        # Extract displacements for the nodes of the element
        dofs_node1 = displacements[node1 * self.num_dofs_per_node: (node1 + 1) * self.num_dofs_per_node]
        dofs_node2 = displacements[node2 * self.num_dofs_per_node: (node2 + 1) * self.num_dofs_per_node]

        # Construct the B matrix (simplified for axial and bending)
        B = np.array([
            [-(dofs_node2[1] - dofs_node1[1]) / L_e**2, 0, 0, 0, 0, 0, (dofs_node2[1] - dofs_node1[1]) / L_e**2, 0, 0, 0, 0, 0],
            [0, -(dofs_node2[2] - dofs_node1[2]) / L_e**2, 0, 0, 0, 0, 0, (dofs_node2[2] - dofs_node1[2]) / L_e**2, 0, 0, 0, 0],
            [0, 0, -(dofs_node2[3] - dofs_node1[3]) / L_e**2, 0, 0, 0, 0, 0, (dofs_node2[3] - dofs_node1[3]) / L_e**2, 0, 0, 0],
            [0, 0, 0, -(dofs_node2[4] - dofs_node1[4]) / L_e**2, 0, 0, 0, 0, 0, (dofs_node2[4] - dofs_node1[4]) / L_e**2, 0, 0],
            [-(dofs_node2[0] - dofs_node1[0]) / L_e, 0, 0, 0, 0, 0, (dofs_node2[0] - dofs_node1[0]) / L_e, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -(dofs_node2[5] - dofs_node1[5]) / L_e, 0, 0, 0, 0, 0, (dofs_node2[5] - dofs_node1[5]) / L_e, 0]
        ])

        return B


#    def calcular_B_Elementar(self, displacements, element):
#        """
#        Calcula a matriz B de um elemento individual.
#        """
#        L_e = self.calcular_comprimento(element)
#        N1 = 1-displacements[element]/L_e
#        N2 = displacements[element]/L_e
#
#        B = np.array([
#            [  0   ,-1/L_e,   0  ,    0 ,    0 ,    0 ,  0  , 1/L_e,   0 ,    0 ,    0 ,   0    ],
#            [  0   ,   0  ,   0  ,    0 ,-1/L_e,   0  ,  0  ,  0   ,  0  ,   0  ,1/L_e ,   0    ],
#            [  0   ,   0  ,   0  ,-1/L_e,   0  ,   0  ,  0  ,  0   ,  0  ,1/L_e ,  0   ,   0    ],
#            [  0   ,   0  ,   0  ,    0 ,    0 ,-1/L_e,  0  ,  0   ,  0  ,   0  ,   0  , 1/L_e  ],
#            [-1/L_e,   0  ,   0  ,    0 ,    0 ,-N1[5],1/L_e,  0   ,  0  ,   0  ,   0  , -N2[11]],
#            [  0   ,   0  ,-1/L_e, N1[3],   0  ,   0  ,  0  ,  0   ,1/L_e, N2[9],  0   ,   0    ],
#        ])
#
#        return B

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
            B = self.calcular_B_Elementar(displacements, element)
            node1, node2 = element
            element_dofs = []
            for node in [node1, node2]:
                for dof in range(self.num_dofs_per_node):
                    element_dofs.append(node * self.num_dofs_per_node + dof)
            element_displacements = displacements[element_dofs]
            strain = np.dot(B, element_displacements)  # B-matrix times element's displacement vector
            strains.append(strain)
        return strains

#        strains = []
#        for element in self.elements:
#            B = self.calcular_B_Elementar(displacements, element)
#            strain = np.dot(B, displacements)  # B-matrix times displacement vector
#            strains.append(strain)
#        return strains
    

    def compute_stress(self, strains, E, nu):
        """
        Compute stresses for all elements using Hooke's law.

        Parameters:
            strains (list of ndarray): Strain tensors for all elements.
            E (float): Young's modulus.
            nu (float): Poisson's ratio.

        Returns:
            stresses (list of ndarray): Stress tensors for all elements.
        """
        # Construct constitutive matrix (isotropic 3D elasticity)
        lambda_ = (E * nu) / ((1 + nu) * (1 - 2 * nu))
        G = E / (2 * (1 + nu))
        C = np.array([
            [lambda_ + 2*G  , lambda_       , lambda_       ,   0,  0,  0],
            [lambda_        , lambda_ + 2*G , lambda_       ,   0,  0,  0],
            [lambda_        , lambda_       , lambda_ + 2*G ,   0,  0,  0],
            [              0,              0,              0,   G,  0,  0],
            [              0,              0,              0,   0,  G,  0],
            [              0,              0,              0,   0,  0,  G]
        ])
        
        stresses = []
        for strain in strains:
            stress = np.dot(C, strain)  # Hooke's law: C times strain
            stresses.append(stress)
        return stresses


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
            von_mises_stresses.append(von_mises)
        return von_mises_stresses


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
                line_loop_indices = ', '.join(str(i + 1) for i in range(len(elements)))
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
        for node1, node2 in elements:
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
        ax.set_box_aspect([3,1,2])
        #plt.xlim([min(nodes)*1.1,max(nodes)*1.1])
        #plt.ylim([min(nodes)*1.1,max(nodes)*1.1])
        plt.tight_layout()
        plt.show()



# THE EXAMPLE FUN STARTS HERE, JUST RUN IN INTERACTIVE WINDOW

#Coordenadas dos nós (x, y, z)
#i = 1.7
#j = 1.5
#k = 1.8

i = 0.3
j = 0.4
k = 0.3

nodes = np.array ([[64*i, 0*j, 0*k] , [64*i, 16*j, 0*k] ,[64*i, 0*j, 16*k] , [64*i, 16*j, 16*k] ,[59*i, 0*j, 7*k] , [59*i, 16*j, 7*k] , [64*i, 0*j, 3*k] , [64*i, 16*j, 3*k] , [50*i, 0*j, 1*k] , [50*i, 16*j, 1*k] , [38*i, 2*j, 1*k] , [38*i, 14*j, 1*k] , [38*i, 0*j, 3*k] , [38*i, 16*j, 3*k] , [38*i, 0*j, 12*k] , [41*i, 16*j, 12*k] , [38*i, 1*j, 24*k] , [38*i, 15*j, 24*k] , [21*i, 0*j, 18*k] , [21*i, 16*j, 18*k] , [23*i, 0*j, 8*k] , [23*i, 16*j, 8*k] , [23*i, 0*j, 0*k] , [23*i, 16*j, 0*k] , [15*i, 0*j, 7*k] , [15*i, 16*j, 7*k] , [8*i, 0*j, 3*k] , [8*i, 16*j, 3*k] , [0*i, 4*j, 7*k] , [0*i, 12*j, 7*k] , [0*i, 4*j, 3*k] , [0*i, 12*j, 3*k] , [0*i, 4*j, 14*k],[0*i, 12*j, 14*k] , [11*i, 1*j, 22*k] , [11*i, 15*j, 22*k] , [19*i, 1*j, 40*k] , [19*i, 15*j, 40*k] , [18*i, 8*j, 45*k] , [38*i, 8*j, 26*k]])  
elements = [(0,1),(0,2),(1,3),(2,3),(4,0),(4,2),(5,1),(5,3),(4,5),(6,7),(0,8),(1,9),(4,8),(5,9),(8,9),(10,8),(10,4),(11,9),(11,5),(10,11),(12,10),(12,4),(13,11),(13,5),(14,12),(14,4),(15,13),(15,5),(16,14),(16,4),(17,15),(17,5),(2,16),(3,17),(16,18),(17,19),(20,18),(20,16),(20,14),(20,10),(21,19),(21,17),(21,15),(21,11),(22,10),(22,20),(23,11),(23,21),(22,23),(24,18),(24,20),(24,22),(25,19),(25,21),(25,23),(26,22),(26,24),(27,23),(27,25),(26,27),(28,30),(28,32),(29,31),(29,33),(30,26),(31,27),(30,31),(28,24),(29,25),(32,24),(32,18),(33,25),(33,19),(32,33),(34,18),(34,32),(35,19),(35,33),(34,35),(36,34),(36,18),(37,35),(37,19),(36,38),(37,38),(16,39),(17,39)]

F_flexao1 = np.array([1000, 2000, 3000, 4000, 5000])
F_flexao2 = np.array([1000, 1000, 1000, 1000, 1000])
F_axial = np.array([1000, 2000, 3000, 4000, 5000])
F_torcao = np.array([1000, 2000, 3000, 4000, 5000])

estrutura = Estrutura(elements, nodes, 180, 4.18e-6, 8.33e-6)

K_global, M_global = estrutura.matrizes_global()

#Gera as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
estrutura.node_loc_matrix(node_tags, nodes)
estrutura.connect_matrix()

#Gerar autovalores, autovetores e frequências naturais
autovalores, autovetores, frequencias = estrutura.modal_analysis()

# Chamando a função shape_fun
torcao, deformacao_axial, flexao1, flexao2, flexao3, KF_total, KT_total, KF_elements, KT_elements = estrutura.shape_fun(F_flexao1, F_flexao2, F_axial, F_torcao)

# # Plotando os resultados das deformações
# fig, axs = plt.subplots(6, 1, figsize=(12, 22))
# 
# # Plot da Torção
# axs[0].plot(torcao, 'o-', label=[f'Força {F}N' for F in F_torcao])
# axs[0].set_title('Deformação por Torção de cada Elemento')
# axs[0].set_xlabel('Elemento')
# axs[0].set_ylabel('Torção (rad)')
# axs[0].legend()
# 
# # Plot da Deformação Axial
# axs[1].plot(deformacao_axial, 's-', label=[f'Força {F}N' for F in F_axial])
# axs[1].set_title('Deformação Axial de cada Elemento')
# axs[1].set_xlabel('Elemento')
# axs[1].set_ylabel('Deformação (m)')
# axs[1].legend()
# 
# # Plot da Flexão por Carga Pontual
# axs[2].plot(flexao1,'o-', label=[f'Força {F}N' for F in F_flexao1])
# axs[2].set_title('Deformação por Carga Pontual de cada Elemento')
# axs[2].set_xlabel('Elemento')
# axs[2].set_ylabel('Deflexão(m)')
# axs[2].legend()
# 
# # Plot da Flexão por Carga Distribuída
# axs[3].plot(flexao2,'o-', label=[f'Força {F}N' for F in F_flexao2])
# axs[3].set_title('Deformação por Carga Distribuída de cada Elemento')
# axs[3].set_xlabel('Elemento')
# axs[3].set_ylabel('Deflexão(m)')
# axs[3].legend()
# 
# # Plot da Flexão Mista
# axs[4].plot(flexao3, 'o-', label='Carregamento misto')
# axs[4].set_title('Deformação por Flexão Mista de cada Elemento')
# axs[4].set_xlabel('Elemento')
# axs[4].set_ylabel('Deflexão (m)')
# axs[4].legend()
# 
# # Plot da Rigidez Flexional e Torsional por Elemento
# axs[5].plot(KF_elements, 'o-', label='Rigidez Flexional (KF)')
# axs[5].plot(KT_elements, 's-', label='Rigidez Torsional (KT)')
# axs[5].set_title('Rigidez Flexional e Torsional de cada Elemento')
# axs[5].set_xlabel('Elemento')
# axs[5].set_ylabel('Rigidez (N/m)')
# axs[5].legend()
# 
# # Mostrando os totais no título geral
# plt.suptitle(f'KF Total: {KF_total:.2e} N/m, KT Total: {KT_total:.2e} N/m', fontsize=16)
# 
# plt.tight_layout(rect=[0, 0, 1, 0.96])
# plt.show()

# Plotando o gráfico 3D da estrutura
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plotar os nós
ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], c='b')   

# Numerar os nós
for i, node in enumerate(nodes):
    ax.text(node[0], node[1], node[2], str(i), color='black')

# Conectando os nós
for element in elements:
    node_start, node_end = element
    ax.plot([nodes[node_start, 0], nodes[node_end, 0]],  # X
            [nodes[node_start, 1], nodes[node_end, 1]],  # Y
            [nodes[node_start, 2], nodes[node_end, 2]])  # Z

# Configurações adicionais
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Estrutura 3D')
ax.legend()
ax.set_box_aspect([3,1,2])
plt.tight_layout()
plt.show()
#Exibindo as frequências naturais e modos de vibração da estrutura
print("\n Frequências Naturais (ω) da estrutura montada por vigas:")
print(frequencias)

#Plotagem dos modos de vibração para a estrutura de vigas
for mode_idx in range(len(autovalores)):
    mode_shape = autovetores[:, mode_idx]
    displacements = np.zeros((len(nodes), 3))  # Assuming we want to visualize x, y, z displacements only

    # Loop through nodes to extract the translations
    for j, (x, y, z) in enumerate(nodes):
        # 6 DOFs per node: [u_x, u_y, u_z, theta_x, theta_y, theta_z]
        dof_start = 6 * j  # Start index of DOFs for node j
        displacements[j, 0] = mode_shape[dof_start]     # u_x
        displacements[j, 1] = mode_shape[dof_start + 1] # u_y
        displacements[j, 2] = mode_shape[dof_start + 2] # u_z

    # Scale displacements for plots
    scale_factor = 1000  # Adjust as needed
    deformed_nodes = np.array(nodes) + displacements * scale_factor

    # Plot deformed
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for node1, node2 in elements:
        x = [deformed_nodes[node1][0], deformed_nodes[node2][0]]
        y = [deformed_nodes[node1][1], deformed_nodes[node2][1]]
        z = [deformed_nodes[node1][2], deformed_nodes[node2][2]]
        ax.plot(x, y, z, 'r-', label="Deformed" if node1 == 0 else "")  # Add label for the first line

    # Add original for comparison
    for node1, node2 in elements:
        x = [nodes[node1][0], nodes[node2][0]]
        y = [nodes[node1][1], nodes[node2][1]]
        z = [nodes[node1][2], nodes[node2][2]]
        ax.plot(x, y, z, 'k--', label="Original" if node1 == 0 else "")  # Add label for the first line

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Forma modal nº: {mode_idx}')
    ax.legend()
    ax.set_box_aspect([3,1,2])
    plt.tight_layout()
    plt.show()

# estrutura.Mesh()
F_global = np.zeros(K_global.size)  # Force vector
F_global[2+5*6] = 100
F_global[2+5*9] = -50
fixed_dofs = [0, 1, 2, 3, 4, 5]

# Perform deformation analysis
displacements = estrutura.static_analysis(K_global, F_global, fixed_dofs)
# print("Displacement Vector:", displacements)
Estrutura.plot_colored_wireframe(nodes, elements, displacements, 'Displacements', 'Displacements [m]')

# Perform equivalent von mises stress determination
strains = estrutura.compute_strain(displacements)
stresses = estrutura.compute_stress(strains, 2.1e11, 0.27)
eq_von_mises = estrutura.compute_von_mises(stresses)
# print("Equivalent Von-Mises Stress:", eq_von_mises)
Estrutura.plot_colored_wireframe(nodes, elements, eq_von_mises, 'Stress', 'Equivalent Von-Mises Stress [Pa]')


"""
Referências
[1] FERREIRA, L. C. M.; SILVA, M. C. F.; GOMES, L. A. Análise da Rigidez Torsional do Chassi de um Protótipo de Competição Formula SAE. Anais do Congresso Brasileiro de Engenharia Mecânica (COBEM), 2021. Disponível em: https://www.scielo.br. Acesso em: 7 jul. 2024.

[2] ZAVATTI, A.; JARDIM, A.; BALTHAZAR, J. M. Estudo Experimental e Numérico da Rigidez Torsional de um Chassi Tubular de Competição. Revista Brasileira de Engenharia Mecânica, v. 22, n. 2, p. 89-101, 2022. Disponível em: https://doi.org. Acesso em: 7 jul. 2024.

[3] KRZIKALLA, F.; et al. Analysis of Torsional Stiffness of the Frame of a Formula Student Vehicle. Journal of Applied Mechanical Engineering, v. 8, n. 1, p. 1-6, 2019. Disponível em: https://www.walshmedicalmedia.com. Acesso em: 7 jul. 2024.

[4] MONTEIRO, R. B.; CORREIA, M. D.; SANTOS, T. A. Determinação Experimental da Rigidez Torsional em Estruturas Tubulares. Anais do Simpósio de Engenharia Automotiva, 2020. Disponível em: https://www.researchgate.net. Acesso em: 7 jul. 2024.

[5] AZEVEDO CANUT, Felipe; MALCHER, Lucival; HENRIQUES, Antonio Manuel Dias. Structural Analysis of a Formula SAE Chassis Under Rollover Loads. In: Proceedings of the 23rd ABCM International Congress of Mechanical Engineering. ABCM, 2015. Disponível em: ABCM Proceedings. Acesso em: 7 jul. 2024.

[6] SETHUPATHI, P. Baskara et al. Design and Optimization of FSAE Chassis Using FEA. IOP Conference Series: Materials Science and Engineering, vol. 402, 2018. Disponível em: IOPscience. Acesso em: 7 jul. 2024.

[7] PATEL, Vijaykumar V.; PATEL, Vibhu. Finite Element Analysis of Truck Chassis Frame. World Journal of Science and Technology, vol. 2, n. 4, p. 5-8, 2012. Disponível em: ResearchGate. Acesso em: 7 jul. 2024.

[8] CALLISTER, William D.; RETHWISCH, David G. Ciência e Engenharia de Materiais: Uma Introdução. 9. ed. Rio de Janeiro: LTC, 2016.

[9] ZIENKIEWICZ, O. C.; TAYLOR, R. L.; ZHU, J. Z. The Finite Element Method for Solid and Structural Mechanics. 6. ed. Oxford: Butterworth-Heinemann, 2005.

[10] EER, Ferdinand P.; JOHNSTON Jr., E. Russell; DEWOLF, John T.; MAZUREK, David F. Resistência dos Materiais. 7. ed. São Paulo: McGraw-Hill Brasil, 2016.
"""
