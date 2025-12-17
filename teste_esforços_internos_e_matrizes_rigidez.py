import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import pandas as pd
import seaborn as sns
import sys
import os

np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=200, suppress=True)

class Estrutura:

    def __init__(self, elements, nodes, planilha = 'tubos_atualizado.csv'):
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
        self.elements = sorted(elements)                                                        #Matriz de elementos conectados
        self.num_elements = len(elements)                                                       #Número de elementos
        self.planilha = planilha
        self.nodes = nodes                                                                      #Matriz de nós com suas posições
        self.num_nodes = len(nodes)                                                             #Número total de nós
        self.num_dofs_per_node = 6                                                              #6 graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node                                 #Total de Graus de liberdade (gdls)
        self.K_global = lil_matrix((self.num_dofs, self.num_dofs), dtype=np.float64)             #Matriz de rigidez global
        self.M_global = lil_matrix((self.num_dofs, self.num_dofs), dtype=np.float64)             #Matriz de massa global
        self.num_modes = 12                                                                     #Número de modos de vibração a serem retornados
        self.car_mass = 0

        # Cache for tube properties to minimize CSV access
        self.tube_properties_cache = {}                                                         #Inicializa armazenamento em cache das propriedades do tubo
        self.load_tube_properties()                                                             #Carrega as propriedades dos tudos antes de realizar as demais funções

# ---------------------
# FUNÇÕES AUXILIARES
# ---------------------

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

    def momento_inercia_area_e_polar(self, dimension, espessura, shape):
        """
        Calculate the momentum of inertia polar and momentum of area of a tube.

        Parameters:
            - dimension:float (Diameter for circular tube or side for square tube)
            - thickness: float (Thickness of the wall)
            - shafe: string (specify the type of tube) 


        outputs:
            - momentum of inertia polar and momentum of area
        """
        
        if shape == 'Circular':
            outer_radius = (dimension / 2)
            inner_radius = (outer_radius - espessura)
            I = (np.pi * 0.25) * (outer_radius ** 4 - inner_radius ** 4)
            J = (np.pi * 0.5) * (outer_radius ** 4 - inner_radius ** 4)
        elif shape == 'Square':
            outer_side = dimension
            inner_side = dimension - 2 * espessura
            if inner_side <= 0:
                raise ValueError("Inner side of square tube is zero or negative. Check dimensions and thickness.")
            I = (outer_side**4 - inner_side**4) / 12
            J = 2 * I # Approximation for polar moment of inertia of a square section
        else:
            raise ValueError(f"Shape '{shape}' not supported for inertia calculation.")
        return I, J

    def area_seccao_transversal(self, dimension, thickness, shape):
        """
        Calculate the area of the tube.

        Parameters:
            - dimension:float (Diameter for circular tube or side for square tube)
            - thickness: float (Thickness of the wall)
            - shafe: string (specify the type of tube) 


        outputs:
            - Section area of the element
        """
        if shape == 'Circular':
            outer_radius = dimension / 2
            inner_radius = (outer_radius - thickness)
            A = (outer_radius ** 2 - inner_radius ** 2) * np.pi
        elif shape == 'Square':
            outer_side = dimension
            inner_side = dimension - 2 * thickness
            if inner_side <= 0:
                raise ValueError("Inner side of square tube is zero or negative. Check dimensions and thickness.")
            A = outer_side**2 - inner_side**2
        else:
            raise ValueError(f"Shape '{shape}' not supported for area calculation.")
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
            # Unpack all properties including 'shape'
            E, G, dimension, e, rho, poisson, shape = self.obter_propriedades(element[2])
            
            A = self.area_seccao_transversal(dimension, e, shape)
            
            volume = L_e * A
            self.car_mass += volume * rho
        return self.car_mass

    def load_tube_properties(self):
        """Load all tube properties from CSV into memory once"""
        df = pd.read_csv(self.planilha)
        for _, row in df.iterrows():
            tube_name = row['Tube']
            self.tube_properties_cache[tube_name] = {
                'E': float(row['E']),
                'G': float(row['G']),
                'Densidade': float(row['Densidade(kg/m^3)']),
                'Poisson': float(row['Poisson']),
                'Espessura': float(row['Espessura(m)']),
                'Shape': row['Shape'],
                'Diametro': float(row['Diametro(m)']) if 'Diametro(m)' in row else None,
                'Side': float(row['Side(m)']) if 'Side(m)' in row else None
            }

    def obter_propriedades(self, tube_name):
        """
        Gets tube properties from csv file

        Parameters:
            - tube_name: string (Contains the type of tubes)

        outputs:
            - Diverses properties of the tube
        """
        if tube_name not in self.tube_properties_cache:
            raise ValueError(f"Tubo '{tube_name}' não encontrado.")
            
        props = self.tube_properties_cache[tube_name]
        E = props['E']
        G = props['G']
        densidade = props['Densidade']
        poisson = props['Poisson']
        espessura = props['Espessura']
        shape = props['Shape']
        
        # Get the relevant dimension based on shape
        if shape == 'Circular':
            dimension = props['Diametro']
        elif shape == 'Square':
            dimension = props['Side']
        else:
            raise ValueError(f"Shape '{shape}' not supported for tube '{tube_name}'.")

        return E, G, dimension, espessura, densidade, poisson, shape
   
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

    def calcular_comprimentos_e_quantidades_por_tipo(self):
        """
        Calculate the amount of tubing of each type

        Output:
            - Total length of each type of tube
            - Quantity of tubes of each type

        """
        
        resultados = {}
        for element in self.elements:
            comprimento = self.calcular_comprimento(element)
            _,_,tipo = element
            if tipo not in resultados:
                resultados[tipo] = {'comprimento_total': 0.0, 'quantidade': 0}

            resultados[tipo]['comprimento_total'] += comprimento
            resultados[tipo]['quantidade'] += 1

        return resultados

    def calcular_K_Elementar(self,element):
        """
        Calcula a matriz de rigidez elementar (12x12) usando a matriz B.
        K_e = Integral(B^T * D * B * dV)

        Para elementos de Timoshenko com B constante ao longo do comprimento,
        a integração se simplifica para: K_e = B^T * D * B 
        (Assumindo que D e B são constantes sobre a seção e ao longo do comprimento,
         ou que as integrais foram pré-calculadas nos termos de B).
        """
        # Variáveis e constantes físicas do modelo
        E, G, dimension, e, rho, poisson, shape = self.obter_propriedades(element[2])

        A = self.area_seccao_transversal(dimension, e, shape)
        I, J = self.momento_inercia_area_e_polar(dimension, e, shape)
        k=0.9                                          #Fator de correção para shear 
        L_e = self.calcular_comprimento(element)
        D = np.array([
            [E * A, 0, 0, 0, 0, 0],                 # Força axial (N)
            [0, E * I, 0, 0, 0, 0],                 # Momento fletor (M_x)
            [0, 0, E * I, 0, 0, 0],                 # Momento fletor (M_z)
            [0, 0, 0, G * A * k , 0, 0],            # Força cortante (V_x)
            [0, 0, 0, 0, G * A * k , 0],            # Força cortante (V_z)
            [0, 0, 0, 0, 0, G * J]])                # Momento torsor (T)

        # 3-point Gauss on xi in [0,1]
        xi_pts = np.array([0.5 - np.sqrt(3.0/5.0)/2.0, 0.5, 0.5 + np.sqrt(3.0/5.0)/2.0])
        w_pts  = np.array([5.0/18.0, 8.0/18.0, 5.0/18.0])

        Ke = np.zeros((12,12), dtype=float)

        for xi, w in zip(xi_pts, w_pts):
            B = self.calcular_B_Elementar(element, xi=xi)  
            # integrando em xi no intervalo [0,1]; para converter para y: integral dy = L * integral dxi
            Ke += w * (B.T @ D @ B)

        Ke *= L_e  # fator Jacobiano (dy = L dxi)
        node1, node2, tube_type = element
        # Coordenadas dos nós
        coord1 = self.nodes[node1]
        coord2 = self.nodes[node2]
        t=self.construir_T(coord1,coord2)
        K_e= t.T @ Ke @ t
        return K_e

# ---------------------
# FUNÇÕES PRNICIPAIS
# ---------------------

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
        E, G, dimension, e, rho, poisson, shape = self.obter_propriedades(element[2])

        A = self.area_seccao_transversal(dimension, e, shape)
        I, J = self.momento_inercia_area_e_polar(dimension, e, shape)
        kappa=0.9                                          #Fator de correção para shear 
        L_e = self.calcular_comprimento(element)
        Phi = (12 * E * I) / (kappa * G * A * L_e**2)

        c1 = (E * A) / L_e
        c2 = (G * J) / L_e
        c3 = (E * I) / L_e**3               #Euler-Bernoulli
        c4 = (E*I)/(L_e**3*(1+Phi))         #Timoshenko
        t1 = (4+Phi)
        t2 = (2-Phi)
        d1 = rho*A*L_e
        d2 = (I*L_e)/6
        d3 = (rho*A*L_e)/420
        
        k_e = np.array([
            # x1          y1      z1       rx1      ry1       rz1       x2          y2        z2       rx2       ry2      rz2
            [12*c4,      0,      0,       0,        0,     6*L_e*c4, -12*c4,       0,        0,        0,       0,     6*L_e*c4],           # x1     0
            [0,          c1,     0,       0,        0,        0,       0,          -c1,       0,        0,       0,       0],               # y1     1
            [0,          0,      12*c4,   6*L_e*c4, 0,        0,       0,           0,      -12*c4,   6*L_e*c4, 0,       0],                # z1     2
            [0,          0,      6*L_e*c4,t1*L_e**2*c4, 0,    0,       0,           0,       -6*L_e*c4, t2*L_e**2*c4, 0,      0],           # rx1    3
            [0,          0,      0,       0,        c2,       0,       0,           0,        0,        0,       -c2,     0],               # ry1    4
            [6*L_e*c4,   0,      0,       0,        0,     t1*L_e**2*c4,-6*L_e*c4,   0,        0,        0,       0,    t2*L_e**2*c4],      # rz1    5
            [-12*c4,     0,      0,       0,        0,     -6*L_e*c4, 12*c4,       0,        0,        0,       0,     -6*L_e*c4],          # x2     6
            [0,          -c1,    0,       0,        0,        0,       0,           c1,       0,        0,       0,       0],               # y2     7   
            [0,          0,      -12*c4,  -6*L_e*c4,  0,       0,       0,           0,      12*c4,    -6*L_e*c4,  0,       0],             # z2     8
            [0,          0,      6*L_e*c4,t2*L_e**2*c4, 0,    0,       0,           0,       -6*L_e*c4, t1*L_e**2*c4, 0,      0],           # rx2    9
            [0,          0,      0,       0,        -c2,      0,       0,           0,        0,        0,        c2,      0],              # ry2    10
            [6*L_e*c4,   0,      0,       0,        0,     t2*L_e**2*c4,-6*L_e*c4,   0,        0,        0,       0,    t1*L_e**2*c4]       # rz2    11
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

    def construir_T(self,coord1, coord2):
        """
        Constrói a matriz de transformação T 12x12 do sistema local para global
        para um elemento 3D com 2 nós e 6 GDL por nó (3 translação + 3 rotação).
        
        Args:
            coord1, coord2: arrays [x, y, z] das coordenadas dos nós
        
        Returns:
            T: matriz de transformação 12x12
        """
        # Vetor axial do elemento
        v = np.array(coord2) - np.array(coord1)
        L = np.linalg.norm(v)
        if L < 1e-12:
            raise ValueError("Comprimento do elemento muito pequeno!")

        y_local = v / L  # eixo axial do elemento
        up = np.array([0, 1, 0], dtype=float)

        if np.abs(np.dot(y_local, up)) > 0.99: 
            up = np.array([1, 0, 0], dtype=float)
        if np.abs(np.dot(y_local, np.array([0, 0, 1]))) > 0.99:
            up = np.array([1, 0, 0], dtype=float)


        z_local = np.cross(y_local, up)
        z_norm = np.linalg.norm(z_local)
        if z_norm < 1e-12:

            raise ValueError("Eixo z_local mal definido!")
        z_local /= z_norm

        x_local = np.cross(y_local, z_local)

        # Monta a matriz de rotação 3x3
        R = np.column_stack((x_local, y_local, z_local))

        # Monta a matriz de transformação 12x12
        T = np.zeros((12, 12))
        for i in range(4):
            T[3 * i:3 * (i + 1), 3 * i:3 * (i + 1)] = R

        return T

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
            ke, me = self.element(element)
            # Coordenadas dos nós
            coord1 = self.nodes[node1]
            coord2 = self.nodes[node2]

            # Matriz de rotação (12x12)
            T = self.construir_T(coord1, coord2)

            # Transformação para o sistema global
            k_e = T.T @ ke @ T
            m_e = T.T @ me @ T

            # DOFs associados ao elemento            
            dofs = [6 * node1, 6 * node1 + 1, 6 * node1 + 2, 6 * node1 + 3, 6 * node1 + 4, 6 * node1 + 5,
                    6 * node2, 6 * node2 + 1, 6 * node2 + 2, 6 * node2 + 3, 6 * node2 + 4, 6 * node2 + 5]
            
            # Update sparse matrices
            for i, dof_i in enumerate(dofs):
                for j, dof_j in enumerate(dofs):
                    self.K_global[dof_i, dof_j] += k_e[i, j]
                    self.M_global[dof_i, dof_j] += m_e[i, j]

        #pd.DataFrame(self.K_global).to_csv('Matriz_Global_Rigidez.csv', index=True, header=True)
        #pd.DataFrame(self.M_global).to_csv('Matriz_Global_Massa.csv', index=True, header=True)        
        #self.K_global = self.K_global.tocsr()
        #self.M_global = self.M_global.tocsr()
        
        # print (self.K_global)
        # print (self.M_global)

        #lt.figure(figsize=(6, 6))
        #lt.spy(self.K_global, markersize=10)  # Adjust markersize for visibility
        #lt.title("Spy Plot of the Kg")
        #lt.xlabel("Columns")
        #lt.ylabel("Rows")
        #lt.grid(True, which="both", linestyle="--", linewidth=0.5)
        #lt.show()

        #lt.figure(figsize=(6, 6))
        #lt.spy(self.M_global, markersize=10)  # Adjust markersize for visibility
        #lt.title("Spy Plot of the Mg")
        #lt.xlabel("Columns")
        #lt.ylabel("Rows")
        #lt.grid(True, which="both", linestyle="--", linewidth=0.5)
        #lt.show()

        return self.K_global,self.M_global

         
    def static_analysis(self, F_global, fixed_dofs):
        """
        Perform static analysis by solving Ku = F with boundary conditions.

        Parameters:
            F_global (ndarray): Global force vector (N).
            fixed_dofs (list): List of DOF indices to be fixed.

        Returns:
            displacements (ndarray): Displacement vector (N).
        """
        # Total number of DOFs
        n_dofs = self.K_global.shape[0]

        # Create a mask for free DOFs (DOFs not constrained)
        free_dofs = np.array([i for i in range(n_dofs) if i not in fixed_dofs])

        # Reduce the stiffness matrix and force vector
        K_reduced = self.K_global[free_dofs, :][:, free_dofs].tocsc()   # CSC format for efficient solve
        F_reduced = F_global[free_dofs]

        try:
            # First try sparse solver
            u_reduced = spsolve(K_reduced, F_reduced)
        except RuntimeError as e:
            print("Sparse solver failed, falling back to dense solver")
            # If sparse solver fails (common with constraints), convert to dense
            u_reduced = np.linalg.solve(K_reduced.toarray(), F_reduced)

        # Construct full displacement vector
        displacements = np.zeros(n_dofs)
        displacements[free_dofs] = u_reduced
        
        return displacements
    
    def funcoes_de_forma(self, xi, L, alpha, beta):
        """
        Funções de forma consistentes para pares de flexão Timoshenko 3D
        (u,ψ) e (w,φ), válidas quando Ix = Iz.
        """
        # H (as in paper) -- note H3,H4 include factor L in definition
        H1 = beta*(2*xi**3 - 3*xi**2 + alpha*xi + 1.0 - alpha)                  #Correto
        H2 = beta*(-2*xi**3 + 3*xi**2 - alpha*xi)                               #Correto
        H3 = L*beta*(xi**3 + (alpha/2.0 - 2.0)*xi**2 + (1.0 - alpha/2.0)*xi)    #Correto
        H4 = L*beta*(xi**3 - (1.0 + alpha/2.0)*xi**2 + (alpha/2.0)*xi)          #Correto
        # G (paper) - note G1,G2 have 1/L factor
        G1 = (6.0*beta/L)*(xi**2 - xi)                                          #correto
        G2 = (6.0*beta/L)*(-xi**2 + xi)                                         #correto
        G3 = beta*(3*xi**2 + (alpha - 4.0)*xi + 1.0 - alpha)                    #correto
        G4 = beta*(3*xi**2 - (alpha + 2.0)*xi)                                  #correto

        # Derivadas
        dH1_dxi = beta*(6*xi**2 - 6*xi + alpha)                                 #Correto
        dH2_dxi = beta*(-6*xi**2 + 6*xi - alpha)                                #Correto
        dH3_dxi = L*beta*(3*xi**2 + (alpha - 4.0)*xi + 1.0 - alpha/2.0)         #Correto
        dH4_dxi = L*beta*(3*xi**2 - (alpha + 2.0)*xi + alpha/2.0)               #Correto
        dG1_dxi = (6.0*beta/L)*(2*xi - 1.0)
        dG2_dxi = (6.0*beta/L)*(1.0 - 2*xi)
        dG3_dxi = beta*(6*xi + (alpha - 4.0))
        dG4_dxi = beta*(6*xi - (alpha + 2.0))

        return {
            'H1':H1,'H2':H2,'H3':H3,'H4':H4,
            'G1':G1,'G2':G2,'G3':G3,'G4':G4,
            'dH1_dxi':dH1_dxi,'dH2_dxi':dH2_dxi,'dH3_dxi':dH3_dxi,'dH4_dxi':dH4_dxi,
            'dG1_dxi':dG1_dxi,'dG2_dxi':dG2_dxi,'dG3_dxi':dG3_dxi,'dG4_dxi':dG4_dxi
        }
    
    def calcular_B_Elementar(self, element,xi=0.5):
        """
        Retorna a matriz B (matriz deformação-deslocamento) para um elemento de viga de Timoshenko 3D
        com 6 graus de liberdade por nó (2 nós por elemento),
        onde a direção axial é o eixo Y. A torção é a última linha.

        Args:
            xi (float): Coordenada normalizada ao longo da viga (y/L), variando de 0 a 1.

        Returns:
            numpy.ndarray: A matriz B 6x12.
        """  

        # Parâmetros do elemento

        E, G, dimension, e, rho, poisson, shape = self.obter_propriedades(element[2])
        L_e = self.calcular_comprimento(element)
        A = self.area_seccao_transversal(dimension, e, shape)
        I, J = self.momento_inercia_area_e_polar(dimension, e, shape)
        k = 0.9  
        alfa = 12.0 * E * I / (k * G * A * L_e**2)
        beta = 1.0 / (1.0 + alfa)        
        # Funções consistentes (uma vez só!)
        p = self.funcoes_de_forma(xi, L_e,alfa,beta)
        invL = 1.0 / L_e
        # convert derivatives to d/dy exactly once:
        dH1_dy = p['dH1_dxi'] * invL
        dH2_dy = p['dH2_dxi'] * invL
        dH3_dy = p['dH3_dxi'] * invL
        dH4_dy = p['dH4_dxi'] * invL
        dG1_dy = p['dG1_dxi'] * invL
        dG2_dy = p['dG2_dxi'] * invL
        dG3_dy = p['dG3_dxi'] * invL
        dG4_dy = p['dG4_dxi'] * invL
        G1 = p['G1']; G2 = p['G2']; G3 = p['G3']; G4 = p['G4']

        dN1_dy = -1.0 / L_e
        dN2_dy = 1.0 / L_e

        B = np.zeros((6,12))
        # order [u1,v1,w1,theta_x1,theta_y1,theta_z1, u2,v2,w2,theta_x2,theta_y2,theta_z2]
        B[0,1] = dN1_dy
        B[0,7] = dN2_dy

        # kappa_x (theta_x' from w/phi family)
        B[1,2] = dG1_dy
        B[1,3] = dG3_dy
        B[1,8] = dG2_dy
        B[1,9] = dG4_dy

        # kappa_z = - theta_z'
        B[2,0] =  (dG1_dy)
        B[2,5] =  (dG3_dy)
        B[2,6] =  (dG2_dy)
        B[2,11] = ( dG4_dy)

        # gamma_xy = u' - theta_z
        B[3,0] = dH1_dy - G1
        B[3,5] = dH3_dy - G3
        B[3,6] = dH2_dy - G2
        B[3,11] = dH4_dy - G4

    # gamma_yz = w' + theta_x
        B[4,2] = dH1_dy - G1
        B[4,3] = dH3_dy - G3
        B[4,8] = dH2_dy - G2
        B[4,9] = dH4_dy - G4

        # torsion theta_y'
        B[5,4] = dN1_dy
        B[5,10] = dN2_dy
            
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
            # B é calculada em coordenadas LOCAIS do elemento
            B = self.calcular_B_Elementar(element)
            
            node1, node2, tube_type = element
            
            # Precisamos das coordenadas dos nós para construir a matriz T
            coord1 = self.nodes[node1]
            coord2 = self.nodes[node2]
            T = self.construir_T(coord1, coord2)

            # Define os DOFs GLOBAIS associados a este elemento
            dofs = [6 * node1 + i for i in range(6)] + [6 * node2 + i for i in range(6)]
            
            element_displacements_global = displacements[dofs]
            element_displacements_local = T @ element_displacements_global

            # Agora sim, multiplica a matriz B (LOCAL) pelos deslocamentos (LOCAIS)
            strain = np.dot(B, element_displacements_local) 

            strains.append(strain)
        return strains

    def compute_stress(self, strains):
        """
        Compute stresses for all elements using Hooke's law.

        Parameters:
            strains (list of ndarray): Strain tensors for all elements.
            E (float): Young's modulus.
            nu (float): Poisson's ratio.

        Returns:
            stresses (list of ndarray): Stress tensors for all elements.
        """
        esforcos = [] 
        # Construct constitutive matrix (isotropic 3D elasticity)
        for i in range(len(self.elements)):
            element=self.elements[i]
            E, G, dimension, e, rho, nu, shape = self.obter_propriedades(element[2])
            L_e = self.calcular_comprimento(element)
            A = self.area_seccao_transversal(dimension, e, shape)
            I, J = self.momento_inercia_area_e_polar(dimension, e, shape)
            k=0.9
            """
            lambda_ = (E * nu) / ((1 + nu) * (1 - 2 * nu))
            C = np.array([
                [lambda_ + 2*G  , lambda_       , lambda_       ,   0,  0,  0],
                [lambda_        , lambda_ + 2*G , lambda_       ,   0,  0,  0],
                [lambda_        , lambda_       , lambda_ + 2*G ,   0,  0,  0],
                [              0,              0,              0,   G,  0,  0],
                [              0,              0,              0,   0,  G,  0],
                [              0,              0,              0,   0,  0,  G]])
            """
            D = np.array([
            [E * A, 0, 0, 0, 0, 0],                 # Força axial (N)
            [0, E * I, 0, 0, 0, 0],                 # Momento fletor (M_x)
            [0, 0, E * I, 0, 0, 0],                 # Momento fletor (M_z)
            [0, 0, 0, G * A * k , 0, 0],            # Força cortante (V_x)
            [0, 0, 0, 0, G * A * k , 0],            # Força cortante (V_z)
            [0, 0, 0, 0, 0, G * J]])                # Momento torsor (T)
            
            strain=strains[i]
        
            esforco = np.dot(D, strain)
            esforcos.append(esforco)
        return esforcos

    def calcular_tensoes_reais(self,esforcos):
        """
        Calcula tensões reais máximas a partir dos esforços internos em um tubo da estrutura.

        Parâmetros:
            esforcos (ndarray): Vetor com 6 esforços internos [N, Mx, Mz, Vx, Vz, T]
            element (tuple): Tupla (nó1, nó2, tipo_tubo)

        Retorno:
            sigma (ndarray): Tensões reais [σ_axial, σ_flex_x, σ_flex_z, τ_cis_x, τ_cis_z, τ_torção]
        """
        tensoes=[]
        componentes=[]
        # Propriedades do tubo
        for i in range(len(self.elements)):
            element=self.elements[i]
            E, G, dimension, e, rho, poisson, shape = self.obter_propriedades(element[2])
            A = self.area_seccao_transversal(dimension, e, shape)
            I, J = self.momento_inercia_area_e_polar(dimension, e, shape)

            D_ext = dimension  # Diâmetro externo
            r = D_ext / 2
            c = r  # Fibra extrema para flexão
            k = 0.9  # Fator de shear/cisalhamento/

            N, Mx, Mz, Vx, Vz, T = esforcos[i]

            # Tensões reais
            sigma_axial = N / A             # Normal axial em y
            sigma_flex_x = Mx * c / I       # Flexão em torno de x
            sigma_flex_z = Mz * c / I       # Flexão em torno de z
            tau_shear_x = Vx / (k * A)      # Cisalhamento em x
            tau_shear_z = Vz / (k * A)      # Cisalhamento em z
            tau_torcao_max = T * r / J          # Torção axial

            # Calculando a tensão normal devido à flexão máxima no eixo y
            M_comb_max = np.sqrt(Mx**2 + Mz**2)
            sigma_flex_max = M_comb_max * c / I
            
            # Calculando a tensão normal máxima no eixo y
            sigma_yy_max = sigma_axial + sigma_flex_max
            
            # Calculando as tensões transversais por efeito poisson (seção simétrica nos dois eixos)
            sigma_xx = -poisson * sigma_yy_max
            sigma_zz = -poisson * sigma_yy_max

            # Calculando a tensão de cisalhamento cortante combinado máxima
            V_res = np.sqrt(Vx**2 + Vz**2)
            tau_V_res_max = V_res / (k * A)

            # Análise do cisalhamento resultante
            # Tensão cisalhante cortante máxima p/ tubo oco ocorre:
                # - no linha média da espessura (centro da parede)
                # - para espessuras pequenas, é aproximadamente igual ao longo da espessura da parede
            
            # Tensão cisalhante torcional é máxima quanto r = r_ext, e mínima quando r = r_int

            # Cálculo da tensão de cisalhamento total máxima
            tau_res_max = np.sqrt(tau_V_res_max**2 + tau_torcao_max**2)

            tensao=np.array([sigma_xx,sigma_yy_max,sigma_zz,tau_res_max])
            tensoes.append(tensao)

            componente = np.array([sigma_axial, sigma_flex_x, sigma_flex_z, tau_shear_x, tau_shear_z, tau_torcao_max])
            componentes.append(componente)

        return tensoes, componentes

# --- 2. Definição da Geometria  ---
# Eixo Y é o axial. Z é vertical.
'''
nodes = np.array([
    [0.0, 0.0, 0.0],  # Nó 0 (Engaste)
    [0.0, 0.5, 0.0],  # Nó 1 (Meio)
    [0.0, 1.0, 0.0]   # Nó 2 (Ponta com Carga)
])

# Conectividade: (Nó Inicial, Nó Final, Nome do Tubo no CSV)
elements = [
    (0, 1, 'Tubo A'), # Elemento 0
    (1, 2, 'Tubo A')  # Elemento 1
]
'''


nodes = np.array([
    [0.0, 0.0, 0.0],  # Nó 0 (Engaste)
    [0.0, 0.1, 0.0],  # Nó 1 
    [0.0, 0.2, 0.0],
    [0.0, 0.3, 0.0],
    [0.0, 0.4, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.6, 0.0],
    [0.0, 0.7, 0.0],
    [0.0, 0.8, 0.0],
    [0.0, 0.9, 0.0],
    [0.0, 1.0, 0.0]   # Nó 10 (Ponta com Carga)
])

# Conectividade: (Nó Inicial, Nó Final, Nome do Tubo no CSV)
elements = [
    (0, 1, 'Tubo A'), # Elemento 0
    (1, 2, 'Tubo A'),
    (2, 3, 'Tubo A'),
    (3, 4, 'Tubo A'),
    (4, 5, 'Tubo A'),
    (5, 6, 'Tubo A'),
    (6, 7, 'Tubo A'),
    (7, 8, 'Tubo A'),
    (8, 9, 'Tubo A'),
    (9, 10, 'Tubo A')  # Elemento 9
]


# --- 3. Inicialização da Estrutura ---
estrutura = Estrutura(elements, nodes)

# Montagem das matrizes globais
K, M = estrutura.matrizes_global()

# --- 4. Condições de Contorno e Carga ---
num_dofs = estrutura.num_dofs
F = np.zeros(num_dofs)

# A) Engaste no Nó 0 (Todos os 6 GDLS travados)
fixed_dofs_0 = [0, 1, 2, 3, 4, 5] 
estrutura.aplicar_engastes([0], fixed_dofs_0) 

# B) Aplicação de Cargas no Nó 2
# Nó 2 inicia no índice de DOF: 2 * 6 = 12
# Fórmula genérica: Nó*GDLs + [0, 1, 2, 3, 4, 5] para [x, y, z, rx, ry, rz]
# GDLs do Nó 2: [12: x, 13: y, 14: z, 15: rx, 16: ry, 17: rz]

# Carga Axial (Tração) de 1000 N no eixo Y
F[10*6 + 1] = 1000.0

# Força de 500 N para baixo no eixo Z
F[10*6 + 2] = -500.0

# Força de 500 N para a direita no eixo X
#F[10*6 + 0] = 500

# Momento de 200 N.m em torno de X aplicado no meio da barra (nó 5 para 10 elementos, nó 1 para 2 elementos)
#F[5*6 + 3] = 200

# Momento de 200 N.m em torno de Z aplicado no meio da barra (nó 5 para 10 elementos, nó 1 para 2 elementos)
#F[5*6 + 5] = 200

# Torque de 100 N.m aplicado no meio da barra (nó 5 para 10 elementos, nó 1 para 2 elementos)
#F[5*6 + 4] = 100

print("--- Cargas Aplicadas ---")
print(f"Força Axial (Y) no Nó 2: {F[13]} N")
print(f"Força Cortante (Z) no Nó 2: {F[14]} N\n")

# --- 5. Análise Estática ---
displacements = estrutura.static_analysis(F, fixed_dofs_0)

print("--- Deslocamentos no Nó 2 ---")
print(f"Deslocamento Y (Axial): {displacements[13]:.6e} m")
print(f"Deslocamento Z (Deflexão): {displacements[14]:.6e} m")
print(f"Rotação em X (Flexão): {displacements[15]:.6e} rad\n")

# --- 6. Cálculo de Esforços e Tensões ---
strains = estrutura.compute_strain(displacements)
internal_forces = estrutura.compute_stress(strains) # O nome é compute_stress mas retorna os strains >:(
real_stresses, stress_components = estrutura.calcular_tensoes_reais(internal_forces)

# --- 7. Exibição dos Resultados ---
print("="*60)
print(f"{'RESULTADOS FINAIS':^60}")
print("="*60)

titulos_esforcos = ['Normal (N)', 'Momento X (Nm)', 'Momento Z (Nm)', 'Cortante X (N)', 'Cortante Z (N)', 'Torção (Nm)']
titulos_tensoes = ['Sigma_XX', 'Sigma_YY_max (Axial+Flex)', 'Sigma_ZZ', 'Tau_res_max']
titulos_componentes = ['Sigma_y', 'Sigma_x', 'Sigma_z', 'Tau_x', 'Tau_z', 'Tau_T']

for i, elem in enumerate(elements):
    n1, n2, tipo = elem
    print(f"\nELEMENTO {i}: Nó {n1} -> Nó {n2}")
    
    # Print Esforços Internos
    print("-" * 30)
    print("ESFORÇOS INTERNOS (Vetor de Forças):")
    forces = internal_forces[i]
    for label, val in zip(titulos_esforcos, forces):
        print(f"  {label:<15}: {val:10.4f}")
        
    # Print Tensões Reais
    print("-" * 30)
    print("TENSÕES REAIS (Pascal):")
    stress = real_stresses[i]
    for label, val in zip(titulos_tensoes, stress):
        print(f"  {label:<15}: {val:10.4e} Pa")
    
    # Print Componentes Normais
    print("-" * 30)
    print("COMPONENTES DE TENSÃO (Pascal):")
    component = stress_components[i]
    for label, val in zip(titulos_componentes, component):
        print(f"  {label:<15}: {val:10.4e} Pa")

# Validando o código comparando as matrizes elementares calculadas pelas diferentes funções
def calcular_diferencas_matrizes_K(estrutura):
    """
    Itera sobre todos os elementos da estrutura, calculando a matriz de rigidez
    pelo método analítico e pelo numérico, e armazena a diferença.
    """
    resultados = []
    
    print(f"Calculando diferenças para {estrutura.num_elements} elementos...")
    
    for i, element in enumerate(estrutura.elements):
        # 1. Método Analítico (Local -> Global)
        # element() retorna k_e no sistema LOCAL
        k_local, _ = estrutura.element(element) 
        
        # Construir matriz de transformação T
        n1, n2, _ = element
        coord1 = estrutura.nodes[n1]
        coord2 = estrutura.nodes[n2]
        T = estrutura.construir_T(coord1, coord2)
        
        # Transformar para GLOBAL: K_global = T.T * K_local * T
        k_analitica = T.T @ k_local @ T
        
        # 2. Método Numérico (Já Global)
        k_numerica = estrutura.calcular_K_Elementar(element)
        
        # 3. Diferença
        # Importante: k_numerica pode vir como matriz esparsa ou densa dependendo da implementação
        # Garantir que ambos sejam numpy arrays densos para a subtração
        if hasattr(k_analitica, 'toarray'): k_analitica = k_analitica.toarray()
        if hasattr(k_numerica, 'toarray'): k_numerica = k_numerica.toarray()
            
        diff = k_numerica - k_analitica
        max_erro = np.max(np.abs(diff))
        
        resultados.append({
            'id': i,
            'analitica': k_analitica,
            'numerica': k_numerica,
            'diff': diff,
            'max_erro': max_erro,
            'nos': (n1, n2)
        })
        
    return resultados

def plot_interativo_matrizes_K(resultados):
    """
    Cria uma janela interativa com Slider para navegar entre as matrizes dos elementos.
    """
    # Labels dos Graus de Liberdade
    labels = ['u1','v1','w1','rx1','ry1','rz1','u2','v2','w2','rx2','ry2','rz2']
    num_elementos = len(resultados)
    
    # Configuração da Figura
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    plt.subplots_adjust(bottom=0.25) # Espaço para o slider
    
    # Títulos iniciais
    ax1, ax2, ax3 = axes
    ax1.set_title("Analítica (Global)")
    ax2.set_title("Numérica (Gauss)")
    ax3.set_title("Diferença")

    # Função auxiliar para plotar heatmap usando imshow (mais rápido para update que seaborn)
    def plot_matrix(ax, data, cmap):
        im = ax.imshow(data, cmap=cmap, aspect='auto', interpolation='nearest')
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=45)
        ax.set_yticklabels(labels)
        return im

    # Plot inicial (Elemento 0)
    data0 = resultados[0]
    im1 = plot_matrix(ax1, data0['analitica'], 'viridis')
    im2 = plot_matrix(ax2, data0['numerica'], 'viridis')
    im3 = plot_matrix(ax3, data0['diff'], 'coolwarm')
    
    # Barras de cores (Colorbars)
    cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
    cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)
    
    # Texto de informação
    info_text = fig.text(0.5, 0.95, 
                         f"Elemento 0 (Nós {data0['nos']}) | Max Erro: {data0['max_erro']:.4e}", 
                         ha='center', fontsize=12, fontweight='bold')

    # --- Configuração do Slider ---
    ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(
        ax=ax_slider,
        label='Elemento ID ',
        valmin=0,
        valmax=num_elementos - 1,
        valinit=0,
        valstep=1
    )

    # Função de atualização
    def update(val):
        idx = int(slider.val)
        data = resultados[idx]
        
        # Atualizar dados das imagens
        im1.set_data(data['analitica'])
        im2.set_data(data['numerica'])
        im3.set_data(data['diff'])
        
        # Atualizar escalas de cores (necessário pois a magnitude muda entre elementos)
        # Para K (Analítica e Numérica), usamos o mesmo range
        vmin_k = min(data['analitica'].min(), data['numerica'].min())
        vmax_k = max(data['analitica'].max(), data['numerica'].max())
        im1.set_clim(vmin_k, vmax_k)
        im2.set_clim(vmin_k, vmax_k)
        
        # Para a Diferença, centralizamos em 0
        max_abs_diff = np.max(np.abs(data['diff']))
        if max_abs_diff == 0: max_abs_diff = 1e-9 # Evitar erro se tudo for zero
        im3.set_clim(-max_abs_diff, max_abs_diff)
        
        # Atualizar texto
        info_text.set_text(f"Elemento {idx} (Nós {data['nos']}) | Max Erro: {data['max_erro']:.4e}")
        
        # Redesenhar
        fig.canvas.draw_idle()

    # Conectar slider à função update
    slider.on_changed(update)
    
    plt.show()

# Calculando as diferenças entre os métodos
k_diffs = calcular_diferencas_matrizes_K(estrutura)

# Exibindo os resultados
plot_interativo_matrizes_K(k_diffs)