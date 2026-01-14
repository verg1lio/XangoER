import numpy as np
from matplotlib.lines import Line2D
from scipy.linalg import eigh
import cadquery as cq
from cadquery import exporters
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix
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

    def calcular_centro_de_massa(self):
    
        soma_massa = 0.0
        soma_massa_pos = np.zeros(3)

        for element in self.elements:
            node1, node2, tube_type = element

            # Coordenadas dos nós
            p1 = np.array(self.nodes[node1])
            p2 = np.array(self.nodes[node2])

            # Comprimento do elemento
            L = np.linalg.norm(p2 - p1)

            # Propriedades do tubo (vindas do CSV)
            E, G, dimension, espessura, densidade, poisson, shape = self.obter_propriedades(tube_type)

            # Área da seção transversal
            A = self.area_seccao_transversal(dimension, espessura, shape)

            # Massa do elemento
            massa_elemento = densidade * A * L

            # Centro geométrico do tubo (meio do elemento)
            centro_elemento = 0.5 * (p1 + p2)

            # Acumulação
            soma_massa += massa_elemento
            soma_massa_pos += massa_elemento * centro_elemento

        if soma_massa == 0:
            raise ValueError("Massa total da estrutura é zero. Verifique os dados.")

        centro_massa = soma_massa_pos / soma_massa

        return centro_massa, soma_massa
    
    def adicionar_massas_concentradas(self, centro_massa_estrutura, massa_estrutura, massas_concentradas):

        soma_massa = massa_estrutura
        soma_massa_pos = massa_estrutura * np.array(centro_massa_estrutura)

        for item in massas_concentradas:
            m = item['massa']
            pos = np.array(item['posicao'])

            soma_massa += m
            soma_massa_pos += m * pos

        if soma_massa == 0:
            raise ValueError("Massa total é zero. Verifique os dados.")

        centro_massa_total = soma_massa_pos / soma_massa

        return centro_massa_total, soma_massa

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

        Para elementos de Timoshenko com B constante ao longo do comprimento (como a sua),
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

        Ke = lil_matrix((self.num_dofs, self.num_dofs), dtype=np.float64)

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

    def shear_correction_factor(self, dimension, e, poisson, shape):
        """
        Calcula o fator de correção de cisalhamento (kappa) para vigas de Timoshenko.
        
        Parâmetros:
            dimension : float
                Dimensão externa (diâmetro externo ou lado externo)
            e : float
                Espessura da parede
            poisson : float
                Coeficiente de Poisson do material
            shape : str
                'circular' ou 'quadrado'
        
        Retorno:
            kappa : float
                Fator de correção de cisalhamento
        """

        if shape.lower() == 'circular':
            # Tubo circular oco
            kappa = (6 * (1 + poisson)) / (7 + 6 * poisson)

        elif shape.lower() == 'square':
            # Tubo quadrado oco (Cowper, 1966)
            b = dimension
            t = e
            beta = (b - 2 * t) / b

            kappa = (
                10 * (1 + poisson) * (1 - beta**2)**2
            ) / (
                12 + 11 * beta**2 + 10 * poisson * (1 - beta**2)**2
            )

        else:
            raise ValueError("Forma de seção inválida. Use 'circular' ou 'quadrado'.")

        return kappa

    def compute_Kf_Kt(self, lfNode, lrNode, rfNode, rrNode, lcNode, rcNode):
            """
            Calculates
            Inputs:
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
            F_global_1 = np.zeros(self.K_global.shape[0])
            F_global_1[2+lfNode*6] = 250                                                                  #Forças aplicadas nos nós onde estaria a suspensão dianteira (nodes 8 e 9)
            F_global_1[2+rfNode*6] = -250                                                                #Mesmo módulo para gerar um torque no eixo longitudinal do chassi
            
            fixed_nodes_1=[lrNode, rrNode]                                                                #Fixação dos nós onde estaria a suspensão traseira
            fixed_dofs_1=[(node*6+i) for node in fixed_nodes_1 for i in range(6)]                           #Lista com os dofs fixados
            
            displacements = self.static_analysis( F_global_1, fixed_dofs_1)                                 #Calcula os displacements com as condições de contorno aplicadas acima
            mi1=np.abs(displacements[lfNode*6+2])                                                       #displacement do nó 9 em módulo
            mi2=np.abs(displacements[rfNode*6+2])                                                            #displacement do nó 8 em módulo
            L = np.abs(self.nodes[lfNode][0] - self.nodes[rfNode][0])                                   #Distancia entre o nó 8 e 9
            alpha= np.degrees(np.atan((mi1+mi2)/(L)))                                                   #Ângulo de torção do chassi após aplicação do torque
            tau = (np.abs(F_global_1[2+rfNode*6]))*(L)                                                    #Cálculo do torque aplicado
            
            Kt = tau/alpha

            #Start Kf simulation
            F_global = np.zeros(self.K_global.shape[0])
            F = 5000                                                                                #Módulo da força que vai gerar a flexão
            F_global[2+lcNode*6] = -F/2                                                             #Força distribuída nos nós centrais do chassi (nodes 22 e 23)
            F_global[2+rcNode*6] = -F/2                                                             #Sinal negativo por conta da direção de aplicação da força
            
            fixed_nodes=[rfNode, lfNode, rrNode, lrNode]                                            #Fixação dos nós onde estaria a suspensão dianteira e traseira
            fixed_dofs=[(node*6+i) for node in fixed_nodes for i in range(6)]                       #Lista com os dofs fixados
            
            displacements = self.static_analysis(F_global, fixed_dofs)                              #Calcula os displacements com as condições de contorno aplicadas acima
            dY=np.abs(displacements[2+lcNode*6])                                                    #Deslocamento em Y de um dos nós onde foi aplicado a força

            Kf=F/dY

            return Kt, Kf
   
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
        kappa= self.shear_correction_factor(dimension,e,poisson,shape)                                         #Fator de correção para shear 
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
                self.K_global[index, index] = 10**20                # Um valor suficientemente grande para simular um engaste 

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
        Beta deve ser beta = 1 / (1 + alfa), e as funções de forma devem ser consistentes com essa formulação
        """
        # Funções H
        # Condições de contorno: H1(0) = 1, H2(1) = 1, H3(0) = 1 e H4(1) = 1.
        H1 = beta*(2*xi**3 - 3*xi**2 + alpha*xi + (1.0 + alpha))                # Deslocamento no nó 1
        H2 = beta*(-2*xi**3 + 3*xi**2 + alpha*xi)                               # Deslocamento no nó 2
        H3 = L*beta*(xi**3 - (2.0 + alpha/2.0)*xi**2 + (1.0 + alpha/2.0)*xi)    # Rotação no nó 1
        H4 = L*beta*(xi**3 - (1.0 - alpha/2.0)*xi**2 - (alpha/2.0)*xi)          # Rotação no nó 2
        
        # Funções G
        G1 = (6.0*beta/L)*(xi**2 - xi)                                          # Influência do deslocamento H1 na rotação
        G2 = (6.0*beta/L)*(-xi**2 + xi)                                         # Influência do deslocamento H2 na rotação
        G3 = beta*(3*xi**2 - (4.0 + alpha)*xi + (1.0 + alpha))                  # Influência da rotação H3 na rotação
        G4 = beta*(3*xi**2 - (2.0 - alpha)*xi)                                  # Influência da rotação H4 na rotação

        # Derivadas
        dH1_dxi = beta*(6*xi**2 - 6*xi - alpha)                                 
        dH2_dxi = beta*(-6*xi**2 + 6*xi + alpha)                                
        dH3_dxi = L*beta*(3*xi**2 - (4.0 + alpha)*xi + (1.0 + alpha/2.0))         
        dH4_dxi = L*beta*(3*xi**2 - (2.0 - alpha)*xi - (alpha/2.0))               

        dG1_dxi = (6.0*beta/L)*(2*xi - 1.0)
        dG2_dxi = (6.0*beta/L)*(1.0 - 2*xi)
        dG3_dxi = beta*(6*xi - (4.0 + alpha))
        dG4_dxi = beta*(6*xi - (2.0 - alpha))

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
        k = self.shear_correction_factor(dimension,e,poisson,shape)                                         #Fator de correção para shear  
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
            E, G, dimension, e, rho, poisson, shape = self.obter_propriedades(element[2])
            L_e = self.calcular_comprimento(element)
            A = self.area_seccao_transversal(dimension, e, shape)
            I, J = self.momento_inercia_area_e_polar(dimension, e, shape)
            k= self.shear_correction_factor(dimension,e,poisson,shape)                                         #Fator de correção para shear
            """
            lambda_ = (E * poisson) / ((1 + poisson) * (1 - 2 * poisson))
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
            k = self.shear_correction_factor(dimension,e,poisson,shape)

            N, Mx, Mz, Vx, Vz, T = esforcos[i]

            # Tensões reais
            sigma_axial = N / A             # Normal axial em y
            sigma_flex_x = Mx * c / I       # Flexão em torno de x
            sigma_flex_z = Mz * c / I       # Flexão em torno de z
            tau_shear_x = Vx / (k * A)      # Cisalhamento em x
            tau_shear_z = Vz / (k * A)      # Cisalhamento em z
            tau_torcao_max = T * r / J      # Torção axial

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
            sigma_xx, sigma_yy, sigma_zz, tau_res = stress
            von_mises = np.sqrt(
                0.5 * (
                    (sigma_xx - sigma_yy)**2 +
                    (sigma_yy - sigma_zz)**2 +
                    (sigma_zz - sigma_xx)**2 +
                    6 * (tau_res**2)
                )
            )
            von_mises_stresses.append(von_mises/10**6)
        return np.array(von_mises_stresses)

# ---------------------
# DEFS DE PLOTAGEM DE GRAFICOS
# ---------------------

    def plotar_chassi_3d(self):
        """
        Plots a 3D wireframe of the structure with color mapping based on the type of tubes
        """
        
        cores_tubos = {
            'Tubo A': 'red',
            'Tubo B': 'blue',
            'Tubo C': 'green',
            'Tubo D': 'orange',
            'Tubo E': 'black'
        }

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Plotar os nós
        for i, (x, y, z) in enumerate(self.nodes):
            ax.scatter(x, y, z, color='black', s=30)
            ax.text(x, y, z, f'{i}', color='black', fontsize=9)

        # Para ajustar escala dos eixos depois
        xs, ys, zs = [], [], []

        # Plotar os tubos
        for idx1, idx2, tipo in self.elements:
            cor = cores_tubos.get(tipo, 'gray')
            x = [self.nodes[idx1][0], self.nodes[idx2][0]]
            y = [self.nodes[idx1][1], self.nodes[idx2][1]]
            z = [self.nodes[idx1][2], self.nodes[idx2][2]]
            ax.plot(x, y, z, color=cor, linewidth=2)
            # Correto
            xs.extend(x)
            ys.extend(y)
            zs.extend(z)

        # Ajustar a escala para que os eixos não fiquem distorcidos
        max_range = max(
            max(xs) - min(xs),
            max(ys) - min(ys),
            max(zs) - min(zs)
        ) / 2

        mid_x = (max(xs) + min(xs)) / 2
        mid_y = (max(ys) + min(ys)) / 2
        mid_z = (max(zs) + min(zs)) / 2

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

        # Legenda com base nos tipos usados no chassi
        tipos_usados = set(tipo for _, _, tipo in self.elements)
        legendas = [
            Line2D([0], [0], color=cores_tubos[tipo], lw=4, label=tipo)
            for tipo in sorted(tipos_usados)
        ]
        ax.legend(handles=legendas, loc='upper left')

        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('Visualização 3D da Estrutura do Chassi')
        plt.tight_layout()
        plt.show()

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

        # --- CÓDIGO PARA AJUSTAR A PROPORÇÃO DOS EIXOS ---
        x_coords, y_coords, z_coords = nodes[:, 0], nodes[:, 1], nodes[:, 2]

        max_range = np.array([np.ptp(x_coords), np.ptp(y_coords), np.ptp(z_coords)]).max()
        
        mid_x = (x_coords.min() + x_coords.max()) / 2.0
        mid_y = (y_coords.min() + y_coords.max()) / 2.0
        mid_z = (z_coords.min() + z_coords.max()) / 2.0

        ax.set_xlim(mid_x - max_range / 2.0, mid_x + max_range / 2.0)
        ax.set_ylim(mid_y - max_range / 2.0, mid_y + max_range / 2.0)
        ax.set_zlim(mid_z - max_range / 2.0, mid_z + max_range / 2.0)
        # --- FIM DO CÓDIGO DE AJUSTE ---

        # Set axis labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(graphtitle)
        
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
        num_modes = len(autovalores)

        # --- Pré-calcular deformações para um ajuste de eixo estável ---
        # (Esta parte do código permanece a mesma)
        all_nodes_for_scaling = [self.nodes]
        for mode_idx in range(num_modes):
            mode_shape = autovetores[:, mode_idx]
            displacements = np.zeros((len(self.nodes), 3))
            for j in range(len(self.nodes)):
                dof_start = 6 * j
                displacements[j, 0] = mode_shape[dof_start]
                displacements[j, 1] = mode_shape[dof_start + 1]
                displacements[j, 2] = mode_shape[dof_start + 2]

            max_dim = np.ptp(self.nodes, axis=0).max()
            max_displacement = np.linalg.norm(displacements, axis=1).max()
            scale_factor = max_dim / (max_displacement * 10) if max_displacement > 1e-9 else 0
            all_nodes_for_scaling.append(np.array(self.nodes) + displacements * scale_factor)
        
        global_nodes = np.vstack(all_nodes_for_scaling)

        # --- Configuração do Gráfico ---
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d', position=[0.05, 0.15, 0.9, 0.8])
        
        class ModeIndex:
            def __init__(self, start_index=0, max_index=num_modes):
                self.index = start_index
                self.max_index = max_index

        current_mode = ModeIndex()

        # --- Função de Atualização do Gráfico ---
        def update_plot():
            mode_idx = current_mode.index
            ax.clear()

            # Calcula e plota a deformação
            # ... (código de cálculo e plotagem das linhas) ...
            mode_shape = autovetores[:, mode_idx]
            displacements = np.zeros((len(self.nodes), 3))
            for j in range(len(self.nodes)):
                dof_start = 6 * j
                displacements[j, 0] = mode_shape[dof_start]
                displacements[j, 1] = mode_shape[dof_start + 1]
                displacements[j, 2] = mode_shape[dof_start + 2]
            
            max_dim = np.ptp(self.nodes, axis=0).max()
            max_displacement = np.linalg.norm(displacements, axis=1).max()
            scale_factor = max_dim / (max_displacement * 10) if max_displacement > 1e-9 else 0
            deformed_nodes = np.array(self.nodes) + displacements * scale_factor

            for i, (node1, node2, _) in enumerate(self.elements):
                x_o = [self.nodes[node1][0], self.nodes[node2][0]]
                y_o = [self.nodes[node1][1], self.nodes[node2][1]]
                z_o = [self.nodes[node1][2], self.nodes[node2][2]]
                ax.plot(x_o, y_o, z_o, 'k--', label="Original" if i == 0 else "", linewidth=1, alpha=0.7)
                x_d = [deformed_nodes[node1][0], deformed_nodes[node2][0]]
                y_d = [deformed_nodes[node1][1], deformed_nodes[node2][1]]
                z_d = [deformed_nodes[node1][2], deformed_nodes[node2][2]]
                ax.plot(x_d, y_d, z_d, 'r-', label="Deformada" if i == 0 else "", linewidth=2)
                
            # Ajuste de proporção
            x_coords, y_coords, z_coords = global_nodes[:, 0], global_nodes[:, 1], global_nodes[:, 2]
            max_range = np.array([np.ptp(x_coords), np.ptp(y_coords), np.ptp(z_coords)]).max()
            mid_x, mid_y, mid_z = np.mean([x_coords.min(), x_coords.max()]), np.mean([y_coords.min(), y_coords.max()]), np.mean([z_coords.min(), z_coords.max()])
            ax.set_xlim(mid_x - max_range / 2.0, mid_x + max_range / 2.0)
            ax.set_ylim(mid_y - max_range / 2.0, mid_y + max_range / 2.0)
            ax.set_zlim(mid_z - max_range / 2.0, mid_z + max_range / 2.0)

            # Labels e Legenda
            ax.set_xlabel('X (m)'); ax.set_ylabel('Y (m)'); ax.set_zlabel('Z (m)')
            ax.legend(loc='upper left')
            
            # =================================================================
            # MUDANÇA PRINCIPAL: Usando o Título da Figura
            # =================================================================
            freq = np.sqrt(autovalores[mode_idx]) / (2 * np.pi)
            fig.suptitle(f'Análise Modal - Modo {mode_idx + 1} (Frequência: {freq:.2f} Hz)', fontsize=16)
            
            # Redesenha a tela
            fig.canvas.draw_idle()

        # --- Funções e criação dos Botões ---
        def next_mode(event):
            current_mode.index = (current_mode.index + 1) % current_mode.max_index
            update_plot()

        def prev_mode(event):
            current_mode.index = (current_mode.index - 1 + current_mode.max_index) % current_mode.max_index
            update_plot()

        ax_prev = fig.add_axes([0.7, 0.05, 0.1, 0.05])
        ax_next = fig.add_axes([0.81, 0.05, 0.1, 0.05])
        button_prev = Button(ax_prev, 'Voltar')
        button_next = Button(ax_next, 'Avançar')
        button_prev.on_clicked(prev_mode)
        button_next.on_clicked(next_mode)

        # --- Exibição Inicial ---
        update_plot()
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

# ---------------------
# DEFS DE EXPORTAÇÃO
# ---------------------
    def tubo_oco(self, p1, p2, tube_type):
            """
            Create a hollow tube using cadquerry commands with intersection fix.

            Parameters:
                -p1: int (start node of the element)
                -p2: int (end node of the element)
                -tube_type: string (Type of the tube)
            """

            # Scale input nodes (which are in meters) to millimeters
            p1_mm = p1 * 1000
            p2_mm = p2 * 1000

            vec = p2_mm - p1_mm # Vector in millimeters
            original_length = np.linalg.norm(vec)

            if original_length < 1e-6:
                print(f"AVISO: Tubo entre {p1_mm} e {p2_mm} tem comprimento zero. Será ignorado.")
                return None

            overlap = 5.0 # 5mm de penetração em cada ponta
            
            # O comprimento final considera a extensão em ambas as extremidades
            final_length = original_length + (2 * overlap)

            # Vetor normalizado para cálculos de direção
            vec_normalized = vec / original_length 

            # Obter propriedades
            try:
                E, G, dimension_m, thickness_m, rho, poisson, shape = self.obter_propriedades(tube_type)
            except Exception as e:
                print(f"❌ Erro ao obter propriedades para {tube_type}: {e}")
                return None

            d = dimension_m * 1000  # Convert outer dimension to mm
            e = thickness_m * 1000  # Convert thickness to mm

            print(f"\n--- DEBUG TUBO_OCO ---")
            print(f"    Tipo: {tube_type}, Forma: {shape}")
            print(f"    Comprimento Original: {original_length:.2f}mm -> Comprimento C/ Overlap: {final_length:.2f}mm")

            outer_shape_workplane = None
            inner_shape_workplane = None

            # Definição da Forma 2D
            if shape == "Circular":
                outer_radius = d / 2
                inner_radius = outer_radius - e
                if inner_radius <= 0:
                    print(f"❌ ERRO GEOMÉTRICO: Raio interno ({inner_radius:.4f}mm) <= 0.")
                    return None
                outer_shape_workplane = cq.Workplane("XY").circle(outer_radius)
                inner_shape_workplane = cq.Workplane("XY").circle(inner_radius)
            
            elif shape == 'Square':
                outer_side = d
                inner_side = d - 2*e
                if inner_side <= 0:
                    print(f"❌ ERRO GEOMÉTRICO: Lado interno ({inner_side:.4f}mm) <= 0.")
                    return None
                outer_shape_workplane = cq.Workplane("XY").rect(outer_side, outer_side)
                inner_shape_workplane = cq.Workplane("XY").rect(inner_side, inner_side)
            
            else:
                print(f"Erro: Forma de tubo '{shape}' não suportada.")
                return None
                
            # --- CORREÇÃO 2: Extrusão com comprimento estendido ---
            # Usamos final_length em vez de length
            outer_solid = outer_shape_workplane.extrude(final_length)
            inner_solid = inner_shape_workplane.extrude(final_length)
            
            # Operação de Corte (Cut)
            tubo_oco_result = outer_solid.cut(inner_solid)

            if not tubo_oco_result.val().isValid(): 
                print(f"❌ ERRO: Objeto inválido após operação de corte.")
                return None

            # --- CORREÇÃO 3: Posicionamento Ajustado ---
            # Calculamos onde o tubo deve começar. 
            # Como aumentamos o tamanho, recuamos o ponto de início na direção oposta do vetor.
            start_point_adjusted = p1_mm - (vec_normalized * overlap)

            # Rotação e Translação
            z_axis = np.array([0, 0, 1])
            axis = np.cross(z_axis, vec_normalized)
            axis_length = np.linalg.norm(axis)

            tubo_final = None

            if axis_length < 1e-6: 
                # Se o tubo já está alinhado com Z (vertical)
                # Se o vetor aponta para baixo (0,0,-1), o produto escalar será -1
                dot_product = np.dot(z_axis, vec_normalized)
                if dot_product > 0:
                    # Aponta para cima, apenas translada
                    tubo_final = tubo_oco_result.translate(tuple(start_point_adjusted))
                else:
                    # Aponta para baixo, rotaciona 180 graus (ou espelha) e translada
                    # O jeito mais seguro no CadQuery para 180 puro é rotacionar em torno de X ou Y
                    tubo_final = tubo_oco_result.rotate((0,0,0), (1,0,0), 180).translate(tuple(start_point_adjusted))
            else:
                angle = np.degrees(np.arccos(np.dot(z_axis, vec_normalized)))
                tubo_final = tubo_oco_result.rotate((0, 0, 0), tuple(axis / axis_length), angle).translate(tuple(start_point_adjusted))
            
            if not tubo_final.val().isValid(): 
                print(f"❌ ERRO: Objeto inválido após rotação/translação.")
                return None

            return tubo_final
    
    def create_step_complete(self, nome="chassi"):
        """
        Create a .STEP file with the entire tubular structure    

        Parameters:

            - nome: string (Name of the file)      
        Outputs:

            - .STEP file
        """

        todos_tubos_validos = [] 
        
        print(f"\n--- INICIANDO GERAÇÃO DA ESTRUTURA ({len(self.elements)} elementos) ---")

        for i, element in enumerate(self.elements):
            node_a, node_b, tube_type = element
            p1 = self.nodes[node_a]
            p2 = self.nodes[node_b]

            # Cria o tubo já com a sobreposição correta
            tubo = self.tubo_oco(p1, p2, tube_type)

            if tubo and tubo.val() and tubo.val().isValid():
                # --- ATENÇÃO: Rotação Individual Removida ---
                # O tubo_oco já alinha o tubo do ponto A ao ponto B.
                # Rotacionar aqui individualmente costuma "desmontar" a estrutura.
                # Se precisar rotacionar o chassi todo, faça isso no final (na variável 'estrutura').
                
                todos_tubos_validos.append(tubo.val())
            else:
                print(f"❌ AVISO: Falha ao criar geometria para o Elemento {i+1} ({tube_type}).")

        if not todos_tubos_validos:
            print("\n❌ ERRO FINAL: Nenhum tubo válido gerado.")
            return

        print(f"\nTentando unir {len(todos_tubos_validos)} tubos...")

        estrutura = None
        
        # 1. Tenta a união booleana (Sólido Único)
        try:
            # Cria um Workplane com todos os objetos
            temp_workplane = cq.Workplane("XY")
            
            # Adiciona todos os objetos ao workplane. 
            # Dica: Adicionar um por um pode ser lento, mas o .add() aceita lista.
            for t in todos_tubos_validos:
                temp_workplane = temp_workplane.add(t)

            # --- MUDANÇA CRÍTICA AQUI ---
            # Removemos o tol=5e-1. Com o overlap físico, a união deve ser exata.
            # Se falhar, usamos combine(glue=True) que é mais rápido para peças que apenas se tocam,
            # mas como temos interseção, o combine() padrão é o correto.
            estrutura = temp_workplane.combine()
            
            # Validação pós-união
            if not estrutura.val() or not estrutura.val().isValid():
                raise Exception("Resultado do combine() é inválido.")
            
            print("\n✅ SUCESSO: Todos os tubos fundidos em um único sólido perfeito!")

        except Exception as ex:
            print(f"\n⚠️ FALHA NA UNIÃO BOOLEANA: {ex}")
            print("ℹ️  Motivo provável: Geometria complexa ou auto-interseção excessiva.")
            print("👉  Alternando para modo 'Compound' (vários sólidos agrupados, visualmente idêntico).")
            
            # 2. Fallback para Compound (Saco de peças)
            try:
                estrutura = cq.Compound.makeCompound(todos_tubos_validos)
                if not estrutura.isValid():
                    print("❌ ERRO GRAVE: Nem o Compound pôde ser criado.")
                    return
            except Exception as e_comp:
                print(f"❌ ERRO CRÍTICO no Fallback: {e_comp}")
                return

        # 3. Exportação
        try:
            # Se você precisar rotacionar O CHASSI INTEIRO (ex: Z-up para Y-up), faça aqui:
            # estrutura = estrutura.rotate((0,0,0), (1,0,0), -90)
            
            nome_arquivo = f"{nome}.step"
            exporters.export(estrutura, nome_arquivo)
            print(f"\n💾 ARQUIVO SALVO: {nome_arquivo}")
        except Exception as ex:
            print(f"\n❌ ERRO DE EXPORTAÇÃO: {ex}")

        # (Opcional) Visualização com Matplotlib mantém-se igual...
        # ...
        # Visualização 3D (This part plots lines from your nodes/elements data, not the 3D solids from CadQuery)
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        # Mostrar nós
        for idx, (x, y, z) in enumerate(self.nodes):
            ax.scatter(x, y, z, c='red', s=60)
            ax.text(x, y, z, f'{idx}', color='black', fontsize=10,
                            verticalalignment='bottom', horizontalalignment='right')

        # Mostrar vigas (as lines)
        for start_idx, end_idx, tube_type in self.elements:
            p1 = self.nodes[start_idx]
            p2 = self.nodes[end_idx]
            xs, ys, zs = zip(p1, p2)
            ax.plot(xs, ys, zs, color='blue', linewidth=5, alpha=0.7)

        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")
        ax.set_zlabel("Z [m]")
        ax.set_title("Frame")
        ax.set_box_aspect([1, 2, 1])

        plt.tight_layout()
        plt.show()

    def create_step_esboco(self, nome_arquivo="chassi_esboço.step", visualizar=True):
            """
            Gera o arquivo STEP da estrutura para CAD e plota a visualização 3D interativa.

            Parâmetros:
                nome_arquivo (str): Nome do arquivo .step a ser salvo.
                visualizar (bool): Se True, exibe o gráfico 3D ao final.
            """
            
            # ------------------------------------------------------------
            # 1. GERAÇÃO DO ARQUIVO STEP (CADQUERY)
            # ------------------------------------------------------------
            # Converter m -> mm para o CAD
            nodes_mm = self.nodes * 1000.0  # m → mm

            wp = cq.Workplane("XY")

            for element in self.elements:
                i = int(element[0])
                j = int(element[1])

                p1 = tuple(nodes_mm[i])
                p2 = tuple(nodes_mm[j])

                wp = wp.add(cq.Edge.makeLine(p1, p2))

            # ⚠️ FORÇAR UM SHAPE ÚNICO
            wireframe = wp.combine()

            exporters.export(
                wireframe,
                nome_arquivo,
                exportType="STEP",
                tolerance=1e-6
            )

            print(f"✅ STEP gerado corretamente: {nome_arquivo}")
            # ------------------------------------------------------------
            # 2. PLOTAGEM (MATPLOTLIB)
            # ------------------------------------------------------------
            if visualizar:
                fig = plt.figure(figsize=(12, 9))
                ax = fig.add_subplot(111, projection='3d')

                # Plotagem dos Nós
                # self.nodes[:, 0] = X, [:, 1] = Y, [:, 2] = Z
                ax.scatter(self.nodes[:, 0], self.nodes[:, 1], self.nodes[:, 2], c='red', s=45)

                # Adicionar rótulos (IDs) aos nós
                for i, (x, y, z) in enumerate(self.nodes):
                    ax.text(x, y, z, f'{i}', fontsize=9)

                # Plotagem das Conexões (Elementos)
                for element in self.elements:
                    start = int(element[0])
                    end = int(element[1])
                    
                    ax.plot(
                        [self.nodes[start, 0], self.nodes[end, 0]],
                        [self.nodes[start, 1], self.nodes[end, 1]],
                        [self.nodes[start, 2], self.nodes[end, 2]],
                        color='black', linewidth=2
                    )

                # Ajuste de aspecto e exibição
                # Tenta manter a proporção real dos eixos
                try:
                    ax.set_box_aspect([1, 1, 1]) 
                    # Ou use limites automáticos para "equal" se a versão do matplotlib suportar
                    max_range = np.array([self.nodes[:,0].max()-self.nodes[:,0].min(), 
                                        self.nodes[:,1].max()-self.nodes[:,1].min(), 
                                        self.nodes[:,2].max()-self.nodes[:,2].min()]).max() / 2.0
                    mid_x = (self.nodes[:,0].max()+self.nodes[:,0].min()) * 0.5
                    mid_y = (self.nodes[:,1].max()+self.nodes[:,1].min()) * 0.5
                    mid_z = (self.nodes[:,2].max()+self.nodes[:,2].min()) * 0.5
                    ax.set_xlim(mid_x - max_range, mid_x + max_range)
                    ax.set_ylim(mid_y - max_range, mid_y + max_range)
                    ax.set_zlim(mid_z - max_range, mid_z + max_range)
                except:
                    # Fallback para o aspecto fixo do seu código original caso o cálculo acima falhe
                    ax.set_box_aspect([1, 3, 2])

                plt.show() 

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

"""
# THE EXAMPLE FUN STARTS HERE, JUST RUN IN INTERACTIVE WINDOW

nodes = np.array([
    [1.92, 0.0, 0.0],    # 64*
    [1.92, 0.64, 0.0],   # 64*i
    [1.92, 0.0, 0.48],   # 64*i
    [1.92, 0.64, 0.48],  # 64*i,
    [1.77, 0.0, 0.21],   # 59*
    [1.77, 0.64, 0.21],  # 59*i
    [1.92, 0.0, 0.09],   # 64*
    [1.92, 0.64, 0.09],  # 64*i
    [1.5, 0.0, 0.03],    # 50*
    [1.5, 0.64, 0.03],   # 50*i
    [1.14, 0.08, 0.03],  # 38*
    [1.14, 0.56, 0.03],  # 38*i
    [1.14, 0.0, 0.09],   # 38*
    [1.14, 0.64, 0.09],  # 38*i
    [1.23, 0.0, 0.36],   # 41*i
    [1.23, 0.64, 0.36],  # 41*i,
    [1.14, 0.03, 0.72],  # 38*i
    [1.14, 0.6, 0.72],   # 38*i,
    [0.63, 0.0, 0.54],   # 21*i
    [0.63, 0.64, 0.54],  # 21*i,
    [0.69, 0.0, 0.24],   # 23*
    [0.69, 0.64, 0.24],  # 23*i
    [0.69, 0.0, 0.0],    # 23*
    [0.69, 0.64, 0.0],   # 23*i
    [0.45, 0.0, 0.21],   # 15*
    [0.45, 0.64, 0.21],  # 15*i
    [0.24, 0.0, 0.09],   # 8*
    [0.24, 0.64, 0.09],  # 8*i
    [0.0, 0.16, 0.21],   # 0*
    [0.0, 0.48, 0.21],   # 0*i
    [0.0, 0.16, 0.09],   # 0*
    [0.0, 0.48, 0.09],   # 0*i
    [0.0, 0.16, 0.42],    # 0*i
    [0.0, 0.48, 0.42],    # 0*i,
    [0.33, 0.03, 0.66],   # 11*i
    [0.33, 0.6, 0.66],    # 11*i
    [0.57, 0.03, 1.2],    # 19*i
    [0.57, 0.6, 1.2],     # 19*i
    [0.54, 0.32, 1.35],   # 18*i
    [1.14, 0.32, 0.78]    # 38*i
])
 


elements = [(0, 1, 'Tubo A'), (0, 6, 'Tubo A'), (6, 2, 'Tubo A'), (1, 7, 'Tubo A'), (7, 3, 'Tubo A'),
 (2, 3, 'Tubo A'), (4, 0, 'Tubo A'), (4, 2, 'Tubo A'), (5, 1, 'Tubo A'), (5, 3, 'Tubo A'),
 (4, 5, 'Tubo A'), (6, 7, 'Tubo A'), (0, 8, 'Tubo A'), (1, 9, 'Tubo A'), (4, 8, 'Tubo A'),
 (5, 9, 'Tubo A'), (8, 9, 'Tubo A'), (10, 8, 'Tubo A'), (10, 4, 'Tubo A'), (11, 9, 'Tubo A'),
 (11, 5, 'Tubo A'), (10, 11, 'Tubo A'), (12, 10, 'Tubo A'), (12, 4, 'Tubo A'), (13, 11, 'Tubo A'),
 (13, 5, 'Tubo A'), (14, 12, 'Tubo A'), (14, 4, 'Tubo A'), (15, 13, 'Tubo A'), (15, 5, 'Tubo A'),
 (16, 14, 'Tubo A'), (16, 4, 'Tubo A'), (17, 15, 'Tubo A'), (17, 5, 'Tubo A'), (2, 16, 'Tubo A'),
 (3, 17, 'Tubo A'), (16, 18, 'Tubo A'), (17, 19, 'Tubo A'), (20, 18, 'Tubo A'), (20, 16, 'Tubo A'),
 (20, 14, 'Tubo A'), (20, 10, 'Tubo A'), (21, 19, 'Tubo A'), (21, 17, 'Tubo A'), (21, 15, 'Tubo A'),
 (21, 11, 'Tubo A'), (22, 10, 'Tubo A'), (22, 20, 'Tubo A'), (23, 11, 'Tubo A'), (23, 21, 'Tubo A'),
 (22, 23, 'Tubo A'), (24, 18, 'Tubo A'), (24, 20, 'Tubo A'), (24, 22, 'Tubo A'), (25, 19, 'Tubo A'),
 (25, 21, 'Tubo A'), (25, 23, 'Tubo A'), (26, 22, 'Tubo A'), (26, 24, 'Tubo A'), (27, 23, 'Tubo A'),
 (27, 25, 'Tubo A'), (26, 27, 'Tubo A'), (28, 30, 'Tubo A'), (28, 32, 'Tubo A'), (29, 31, 'Tubo A'),
 (29, 33, 'Tubo A'), (30, 26, 'Tubo A'), (31, 27, 'Tubo A'), (30, 31, 'Tubo A'), (28, 24, 'Tubo A'),
 (29, 25, 'Tubo A'), (32, 24, 'Tubo A'), (32, 18, 'Tubo A'), (33, 25, 'Tubo A'), (33, 19, 'Tubo A'),
 (32, 33, 'Tubo A'), (34, 18, 'Tubo A'), (34, 32, 'Tubo A'), (35, 19, 'Tubo A'), (35, 33, 'Tubo A'),
 (34, 35, 'Tubo A'), (36, 34, 'Tubo A'), (36, 18, 'Tubo A'), (37, 35, 'Tubo A'), (37, 19, 'Tubo A'),
 (36, 38, 'Tubo A'), (37, 38, 'Tubo A'), (16, 39, 'Tubo A'), (17, 39, 'Tubo A')]

 
F_flexao1 = np.array([1000, 2000, 3000, 4000, 5000])
F_flexao2 = np.array([1000, 1000, 1000, 1000, 1000])
F_axial = np.array([1000, 2000, 3000, 4000, 5000])
F_torcao = np.array([1000, 2000, 3000, 4000, 5000])

estrutura = Estrutura(elements, nodes)
K_global, M_global = estrutura.matrizes_global()
#estrutura.aplicar_engastes([30,31,32,33],[0,1,2,3,4,5])
#Gera as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
estrutura.node_loc_matrix(node_tags, nodes)
estrutura.connect_matrix()
print(f"Massa do chassi: {estrutura.mass()} kg")
#Plotando os resultados das deformações
#estrutura.shape_fun_plot(F_flexao1, F_flexao2, F_axial,F_torcao)

#Gerar autovalores, autovetores e frequências naturais
autovalores, autovetores, frequencias = estrutura.modal_analysis()
#Exibindo as frequências naturais e modos de vibração da estrutura
print("\\n Frequências Naturais (ω) da estrutura:")
print(frequencias)
estrutura.modal_analysis_plot()

# estrutura.Mesh()
F_global = np.zeros(K_global.size)  # Force vector
#F_global[2+8*6] = -100
#F_global[2+9*6] = 100
#fixed_dofs = []
F_global[15*6+2] = 100
fixed_dofs = [0,1,2,3,4,5]
# Perform deformation analysis
displacements = estrutura.static_analysis(F_global, fixed_dofs)
print("Displacement Vector:", displacements)


# Perform torsional and flexional stiffness analysis
#Kt,Kf= estrutura.compute_Kf_Kt(8,27,9,26)
#print(f'Rigidez Torcional do chassi atual: {Kt}')
#print(f'Rigidez Flexional do chassi atual: {Kf:.2e}')

Estrutura.plot_colored_wireframe(nodes, elements, displacements, 'Displacements', 'Displacements [m]')

# Perform equivalent von mises stress determination
strains = estrutura.compute_strain(displacements)
#print("strain:", strains)
stresses = estrutura.compute_stress(strains)
#print("Stress:", stresses)
tensoes = estrutura.calcular_tensoes_reais(stresses)
von_mises = estrutura.compute_von_mises(tensoes)

print("Equivalent Von-Mises Stress:", von_mises)
Estrutura.plot_colored_wireframe(nodes, elements, von_mises, 'Stress', 'Equivalent Von-Mises Stress [MPa]')

estrutura.plotar_chassi_3d()

print(estrutura.calcular_comprimentos_e_quantidades_por_tipo())
"""

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

[11] LUO, Y. An efficient 3D Timoshenko beam element with consistent shape functions. Advances in Theoretical and Applied Mechanics, v. 1, n. 3, p. 95-106, out. 2008.
"""
