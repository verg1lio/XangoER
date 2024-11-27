import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import pandas as pd
import sys
import os
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=200, suppress=True)

class Estrutura:
    def __init__(self, elements,element_properties, nodes, m, Id, Ip):
        self.elements = elements                                              #Matriz de elementos conectados
        self.num_elements = len(elements)                                     #Número de elementos
        self.element_properties = element_properties                          #Matriz de propriedades dos elementos
        self.nodes = nodes                                                    #Matriz de nós com suas posições
        self.coordinates = np.array(nodes[['x', 'y', 'z']])
        self.connections = np.array(elements[['Node a', 'Node b']])
        self.num_nodes = len(nodes)                                           #Número total de nós
        self.massa = m                                                        #Massa do carro (Kg)
        self.momento_inercia_direcao = Id                                     #Momento de inércia em relação à direção (kg.m^2)
        self.momento_inercia_plano = Ip                                       #Momento de inércia em relação ao plano (kg.m^2)
        self.num_dofs_per_node = 6                                            #Graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node               #Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((len(nodes) * 6, len(nodes) * 6))            #Tamanho adequado para a matriz global
        self.M_global = np.zeros((len(nodes) * 6, len(nodes) * 6))
        self.num_modes = 12                                                   #Número de modos de vibração a serem retornados

    def node_loc_matrix(self):
        print(self.nodes)

    def connect_matrix(self):
        print(self.elements)

    def calcular_comprimento(self, index):                                 #Função auxiliar para cálculo de comprimento dos elementos
        
        node1 , node2 = int(self.elements["Node a"][index]-1) , int(self.elements["Node b"][index]-1)
        x1, y1, z1 = self.nodes['x'][node1], self.nodes['y'][node1], self.nodes['z'][node1]
        x2, y2, z2 = self.nodes['x'][node2], self.nodes['y'][node2], self.nodes['z'][node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    
    def element(self, index):
        #Tomando as propriedades de cada elemento
        E = self.element_properties['E'][index]
        A = self.element_properties['A'][index]
        I = self.element_properties['I'][index]
        J = self.element_properties['J'][index]
        G = self.element_properties['G'][index]

        kappa=0.9       #Fator de correção para cisalhamento 

        L_e = self.calcular_comprimento(index)
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
        for node in nodes:                                          #Laço para selecionar cada nó que será engastado
            for dof in dofs:                                        #Laço para selecionar quais graus de liberdade serão fixados
                index = node * self.num_dofs_per_node + dof         #Identificação da entrada da matriz que precisa ser restringida pelo engaste        
                self.K_global[index, index] = 10**10                # Um valor suficientemente grande para simular um engaste 
                   
    def matrizes_global(self):
        #Calculando as matrizes de rigidez e massa de cada elemento
        for index in range(self.num_elements):

            k_e, m_e= self.element(index)
            node1 , node2 = int(self.elements["Node a"][index]-1) , int(self.elements["Node b"][index]-1)
            # DOFs associados ao elemento            
            dofs = [6 * node1, 6 * node1 + 1, 6 * node1 + 2, 6 * node1 + 3, 6 * node1 + 4, 6 * node1 + 5,
                    6 * node2, 6 * node2 + 1, 6 * node2 + 2, 6 * node2 + 3, 6 * node2 + 4, 6 * node2 + 5]

            # Atualizando as matrizes globais
            self.K_global[np.ix_(dofs, dofs)] += k_e
            self.M_global[np.ix_(dofs, dofs)] += m_e

        #self.aplicar_engastes([0, 2, 4, 5], [0, 1, 2, 3, 4, 5])                             #Por enquanto não estaremos considerando engastes
        pd.DataFrame(self.K_global).to_csv('Matriz_Global_Rigidez.csv', index=True, header=True)
        pd.DataFrame(self.M_global).to_csv('Matriz_Global_Massa.csv', index=True, header=True)        

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
        KF_total = 0
        KT_total = 0
        KF_elements = []
        KT_elements = [] 
        torcao, deformacao, flexao1, flexao2, flexao3 = [], [], [], [], []
        for index in range(len(self.elements)):
            E = self.element_properties['E'][index]
            A = self.element_properties['A'][index]
            I = self.element_properties['I'][index]
            J = self.element_properties['J'][index]
            G = self.element_properties['G'][index]
            L_e = self.calcular_comprimento(index)
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

    def static_analysis(self,K_global, F_global, fixed_dofs):
        """
        Perform static analysis by solving Ku = F with boundary conditions.

        Parameters:
            K_global (ndarray): Global stiffness matrix (N x N).
            F_global (ndarray): Global force vector (N).
            fixed_dofs (list): List of DOF indices to be fixed.

        Returns:
            displacements (ndarray): Displacement vector (N).
        """
        # Total number of DOFs
        n_dofs = K_global.shape[0]

        # Create a mask for free DOFs (DOFs not constrained)
        free_dofs = np.array([i for i in range(n_dofs) if i not in fixed_dofs])

        # Reduce the stiffness matrix and force vector
        K_reduced = K_global[np.ix_(free_dofs, free_dofs)]
        F_reduced = F_global[free_dofs]

        # Solve for displacements at free DOFs
        u_reduced = np.linalg.solve(K_reduced, F_reduced)

        # Construct full displacement vector
        displacements = np.zeros(n_dofs)
        displacements[free_dofs] = u_reduced

        return displacements

    def compute_strain(self, B_matrices,F_global, fixed_dofs):
        """
        Compute strains for all elements.

        Parameters:
            displacements (ndarray): Displacement vector for all nodes.
            B_matrices (list of ndarray): Strain-displacement matrices for each element.

        Returns:
            strains (list of ndarray): Strain tensors for all elements.
        """
        #### AQUI A MATRIZ B PRECISA SER CALCULADA DIRETO DAS FUNÇÕES DE FORMA
        #### NÃO É TRIVIAL E VAI LEVAR TEMPO, MAS O RESTO DAS CONTAS ESTÃO OK
        displacements = self.static_analysis(F_global, fixed_dofs)
        strains = []
        for B in B_matrices:
            strain = np.dot(B, displacements)  # B-matrix times displacement vector
            strains.append(strain)
        return strains
    
    def compute_stress(self,B_matrices,F_global, fixed_dofs, E, nu):
        """
        Compute stresses for all elements using Hooke's law.

        Parameters:
            strains (list of ndarray): Strain tensors for all elements.
            E (float): Young's modulus.
            nu (float): Poisson's ratio.

        Returns:
            stresses (list of ndarray): Stress tensors for all elements.
        """
        strains = self.compute_strain(B_matrices,F_global, fixed_dofs)
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

    def compute_von_mises(self,B_matrices,F_global, fixed_dofs, E, nu):
        """
        Compute von Mises stress for all elements.

        Parameters:
            stresses (list of ndarray): Stress tensors for all elements.

        Returns:
            von_mises_stresses (list of float): Von Mises stress for each element.
        """
        stresses = self.compute_stress(B_matrices,F_global, fixed_dofs, E, nu)
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

    def structure_plot(self):
        # Plotando o gráfico 3D da estrutura
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Adicionando os pontos (nós)
        for i, (x, y, z) in enumerate(self.coordinates):
            ax.scatter(x, y, z, color='b', s=50)
            ax.text(x, y, z, f' {i+1}', color='black', fontsize=8)

        # Adicionando as linhas de ligação entre os nós
        for node1, node2 in self.connections:
            x_coords = [self.coordinates[node1 - 1][0], self.coordinates[node2 - 1][0]]
            y_coords = [self.coordinates[node1 - 1][1], self.coordinates[node2 - 1][1]]
            z_coords = [self.coordinates[node1 - 1][2], self.coordinates[node2 - 1][2]]
            ax.plot(x_coords, y_coords, z_coords, 'r-', marker='o')

        # Configurações adicionais do gráfico
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Estrutura 3D do Chassi')

        plt.show()

    def plot_colored_wireframe(self,scalar_values, colormap='jet'):
        """
        Plots a 3D wireframe of the structure with color mapping based on scalar values.
        
        Parameters:
            nodes (array): Array of node coordinates (N x 3).
            elements (list): List of tuples defining connections between nodes.
            scalar_values (array): 1D array of scalar values (e.g., strain) at each node.
            colormap (str): Colormap name for visualization.
        """
        # Normalize scalar values to [0, 1] for colormap
        norm = plt.Normalize(vmin=np.min(scalar_values), vmax=np.max(scalar_values))
        cmap = plt.get_cmap(colormap)
        
        # Create the plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Adicionando os pontos (nós)
        for i, (x, y, z) in enumerate(self.coordinates):
            ax.scatter(x, y, z, color='b', s=50)
            ax.text(x, y, z, f' {i+1}', color='black', fontsize=8)

        # Adicionando as linhas de ligação entre os nós
        for node1, node2 in self.connections:
            x_coords = [self.coordinates[node1 - 1][0], self.coordinates[node2 - 1][0]]
            y_coords = [self.coordinates[node1 - 1][1], self.coordinates[node2 - 1][1]]
            z_coords = [self.coordinates[node1 - 1][2], self.coordinates[node2 - 1][2]]
            ax.plot(x_coords, y_coords, z_coords, 'r-', marker='o')
            # Get the scalar value for the midpoint of the element
            scalar_midpoint = (scalar_values[node1-1] + scalar_values[node2-1]) / 2
            
            # Map scalar value to color
            color = cmap(norm(scalar_midpoint))
            
            # Plot the line segment with the corresponding color
            ax.plot(x_coords, y_coords, z_coords, color=color, linewidth=2)

        # Add a colorbar
        mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        mappable.set_array(scalar_values)
        cbar = plt.colorbar(mappable, ax=ax, orientation='vertical', shrink=0.8, pad=0.1)
        cbar.set_label("Strain (or other variable)", fontsize=12)

        # Set axis labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Wireframe with Scalar Color Mapping')
        plt.tight_layout()
        plt.show()
           
    def modal_analysis_plot(self):
        autovalores, autovetores, _ = self.modal_analysis()

        for mode_idx in range(len(autovalores)):
            mode_shape = autovetores[:, mode_idx]
            displacements = np.zeros((len(self.coordinates), 3))  # Assuming we want to visualize x, y, z displacements only

            # Loop through nodes to extract the translations
            for j, (x, y, z) in enumerate(self.coordinates):
                # 6 DOFs per node: [u_x, u_y, u_z, theta_x, theta_y, theta_z]
                dof_start = 6 * j  # Start index of DOFs for node j
                displacements[j, 0] = mode_shape[dof_start]     # u_x
                displacements[j, 1] = mode_shape[dof_start + 1] # u_y
                displacements[j, 2] = mode_shape[dof_start + 2] # u_z

            # Scale displacements for plots
            scale_factor = 1  # Adjust as needed
            deformed_nodes = np.array(self.coordinates) + displacements * scale_factor

            # Plot deformed
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')

            # Plot deformed structure
            for i, (node1, node2) in enumerate(self.connections):
                x = [deformed_nodes[node1-1][0], deformed_nodes[node2-1][0]]
                y = [deformed_nodes[node1-1][1], deformed_nodes[node2-1][1]]
                z = [deformed_nodes[node1-1][2], deformed_nodes[node2-1][2]]
                ax.plot(x, y, z, 'r-', label="Deformed" if i == 0 else "")  # Add label only once

            # Plot original structure
            for i, (node1, node2) in enumerate(self.connections):
                x = [self.coordinates[node1-1][0], self.coordinates[node2-1][0]]
                y = [self.coordinates[node1-1][1], self.coordinates[node2-1][1]]
                z = [self.coordinates[node1-1][2], self.coordinates[node2-1][2]]
                ax.plot(x, y, z, 'k--', label="Original" if i == 0 else "")  # Add label only once

            # Add labels and title
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_title(f'Forma modal nº: {mode_idx + 1}')
            ax.legend()  # Ensure the legend is displayed
            ax.set_zlim([-0.5,2])
            plt.tight_layout()
            plt.show()

    def shape_fun_plot(self, F_flexao1, F_flexao2, F_axial, F_torcao):
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


nodes_file_path = "C:\\Users\\dudua\\OneDrive\\Documentos\\GitHub\\EduardoChassi\\Nós e Elementos modelo de chassi básico - Nodes.csv"
elements_file_path = "C:\\Users\\dudua\\OneDrive\\Documentos\\GitHub\\EduardoChassi\\Nós e Elementos modelo de chassi básico - Elements.csv"

# Carregar o arquivo para inspecionar seu conteúdo
nodes = pd.read_csv(nodes_file_path)
element_data = pd.read_csv(elements_file_path)

# Selecionando as colunas de conectividades e propriedades dos elementos

elements = element_data[['Element ID','Node a', 'Node b']]
element_properties = element_data[['Element ID','A', 'I', 'J', 'E', 'G']]

#Criar a estrutura e montar as matrizes de rigidez e massa globais
#Dados: n = len(nodes), 
#       m = 1500 kg, 
#       rho = 7850 kg/m^3
#       A = 0.225 m^2
#       E = 210e9  # Módulo de elasticidade em Pa
#       I = 8.33e-6  # Momento de inércia em m^4
#       Ip = Id = 8.33e-6 kg.m^2

#Criar a estrutura e montar as matrizes de rigidez e massa globais, atribuir forças
F_flexao1 = np.array([1000, 2000, 3000, 4000, 5000])
F_flexao2 = np.array([1000, 1000, 1000, 1000, 1000])
F_axial   = np.array([1000, 2000, 3000, 4000, 5000])
F_torcao  = np.array([1000, 2000, 3000, 4000, 5000])

#Inicializando a Estrutura
estrutura = Estrutura(elements, element_properties, nodes, 1500, 8.33e-6, 8.33e-6)

#Gerar as matrizes de localização dos nós e de conectividade
estrutura.node_loc_matrix()
estrutura.connect_matrix()

# Plotando o gráfico 3D da estrutura
estrutura.structure_plot()

K_global, M_global = estrutura.matrizes_global()

#Plotando os resultados das deformações
estrutura.shape_fun_plot(F_flexao1, F_flexao2, F_axial,F_torcao)

#Gerar autovalores, autovetores e frequências naturais
autovalores, autovetores, frequencias = estrutura.modal_analysis()

#Exibindo as frequências naturais e modos de vibração da estrutura
print("\\n Frequências Naturais (ω) da estrutura:")
print(frequencias)

#Plotagem dos modos de vibração para a estrutura de vigas
estrutura.modal_analysis_plot()

F_global = np.zeros(K_global.size)  # Force vector
F_global[2+5*6] = 100
F_global[2+5*9] = -50
fixed_dofs = [0, 1, 2, 3, 4, 5]

# Perform analysis
displacements = estrutura.static_analysis(K_global,F_global, fixed_dofs)
print("Displacement Vector:", displacements)

estrutura.plot_colored_wireframe(displacements)
print(nodes.size)
print(displacements.size)
print(F_torcao.size)
print(F_flexao1.size)

"""
Estrutura.plot_colored_wireframe(nodes, elements, torcao/(np.max(np.max(torcao))))
Estrutura.plot_colored_wireframe(nodes, elements, deformacao_axial)
Estrutura.plot_colored_wireframe(nodes, elements, flexao1)
Estrutura.plot_colored_wireframe(nodes, elements, flexao2)
Estrutura.plot_colored_wireframe(nodes, elements, flexao3)

"""