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
from matplotlib.lines import Line2D # Importe esta linha para a função plotar_chassi_3d

np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=200, suppress=True)

# --- Funções principais (copiadas da sua primeira mensagem) ---

def calcular_comprimento(node1, node2):
    return np.linalg.norm(node2 - node1)

def calcular_comprimentos_e_quantidades_por_tipo(nodes, elements):
    resultados = {}
    for idx1, idx2, tipo in elements:
        node1 = nodes[idx1]
        node2 = nodes[idx2]
        comprimento = calcular_comprimento(node1, node2)

        if tipo not in resultados:
            resultados[tipo] = {'comprimento_total': 0.0, 'quantidade': 0}

        resultados[tipo]['comprimento_total'] += comprimento
        resultados[tipo]['quantidade'] += 1

    return resultados

def set_axes_equal(ax):
    # Deixa a escala dos eixos x, y e z iguais
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max(x_range, y_range, z_range)

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)

    ax.set_xlim3d([x_middle - max_range / 2, x_middle + max_range / 2])
    ax.set_ylim3d([y_middle - max_range / 2, y_middle + max_range / 2])
    ax.set_zlim3d([z_middle - max_range / 2, z_middle + max_range / 2]) 


"""class Chassi3D:
    def __init__(self, nodes, elements):
        self.nodes = nodes
        self.elements = elements 

def plotar_chassi_3d(self):

        # Definindo cores para tipos diferentes de tubo
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

        # Plotar os tubos
        for idx1, idx2, tipo in self.elements:
            cor = cores_tubos.get(tipo, 'gray')
            x = [self.nodes[idx1][0], self.nodes[idx2][0]]
            y = [self.nodes[idx1][1], self.nodes[idx2][1]]
            z = [self.nodes[idx1][2], self.nodes[idx2][2]]
            ax.plot(x, y, z, color=cor, linewidth=2)

        ax.set_box_aspect([1, 2, 1])
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('Visualização 3D da Estrutura do Chassi')
        plt.tight_layout()
        plt.show() """ 


"""class Chassi3D:
    def __init__(self, nodes, elements):
        self.nodes = nodes
        self.elements = elements

    def plotar_chassi_3d(self):
        cores_tubos = {
            'Tubo A': 'red',
            'Tubo B': 'blue',
            'Tubo C': 'green',
            'Tubo D': 'orange',
            'Tubo E': 'black'
        }

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        for i, (x, y, z) in enumerate(self.nodes):
            ax.scatter(x, y, z, color='black', s=30)
            ax.text(x, y, z, f'{i}', color='black', fontsize=9)

        for idx1, idx2, tipo in self.elements:
            cor = cores_tubos.get(tipo, 'gray')
            x = [self.nodes[idx1][0], self.nodes[idx2][0]]
            y = [self.nodes[idx1][1], self.nodes[idx2][1]]
            z = [self.nodes[idx1][2], self.nodes[idx2][2]]
            ax.plot(x, y, z, color=cor, linewidth=2)

        ax.set_box_aspect([1, 2, 1])
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('Visualização 3D da Estrutura do Chassi')
        plt.tight_layout()
        plt.show()"""


class Chassi3D:
    def __init__(self, nodes, elements):
        self.nodes = nodes
        self.elements = elements

    def plotar_chassi_3d(self):
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


# --- Fim das funções principais ---


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
        self.elements = elements                                        #Matriz de elementos conectados
        self.num_elements = len(elements)                               #Número de elementos
        self.nodes = nodes                                              #Matriz de nós com suas posições
        self.num_nodes = len(nodes)                                     #Número total de nós
        self.num_dofs_per_node = 6                                      #6 graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node         #Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))        #Matriz de rigidez global
        self.M_global = np.zeros((self.num_dofs, self.num_dofs))        #Matriz de massa global
        self.num_modes = 1                                              #Número de modos de vibração a serem retornados
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
            d = self.obter_propriedades(element[2])[2]      #Diâmetro do Tubo (m)
            e = self.obter_propriedades(element[2])[3]      #Espessura do Tubo (m)
            A = self.area_seccao_transversal(d,e)           #Área da secção transversal (m^2)
            rho = self.obter_propriedades(element[2])[4]    #Densidade do material (kg/m^3)
            raio_externo = d / 2
            raio_interno = raio_externo - e
            volume = np.pi*L_e* (raio_externo**2 - raio_interno**2)
            self.car_mass+= volume*rho

        return self.car_mass
    
    def obter_propriedades(self, tube_nome):                                         #Função que lê a planilha 'tubos.csv' e extrai as propriedades de cada tipo de tubo lá presentes
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
            
        # print("\n    Nó    x    y    z")
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
        d = self.obter_propriedades(element[2])[2]      #Diâmetro do Tubo (m)
        e = self.obter_propriedades(element[2])[3]      #Espessura do Tubo (m)
        E = self.obter_propriedades(element[2])[0]      #Modulo de Young (Pa)      
        G = self.obter_propriedades(element[2])[1]      #Modulo de Cisalhamento (Pa)
        A = self.area_seccao_transversal(d, e)          #Área da seção do elemento (m^2)
        I = self.momento_inercia_area_e_polar(d,e)[0]   #Momento de inercia (m^4)
        J = self.momento_inercia_area_e_polar(d,e)[1]   #Momento polar de inércia (m^4)
        rho = self.obter_propriedades(element[2])[4]    #Densidade do material (kg/m^3)
        kappa=0.9                                       #Fator de correção para cisalhamento 
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
        for node in nodes:                                   # Laço para selecionar cada nó que será engastado
            for dof in dofs:                                 # Laço para selecionar quais graus de liberdade serão fixados
                index = node * self.num_dofs_per_node + dof  # Identificação da entrada da matriz que precisa ser restringida pelo engaste        
                self.K_global[index, index] = 10**9          # Um valor suficientemente grande para simular um engaste    

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

        #self.aplicar_engastes([30,31,32,33], [0, 1, 2, 3, 4, 5])                          #Por enquanto não estaremos considerando engastes
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
        #E = 2.1e11      #Modulo de Young (Pa)
        #I = 1.6667e-5   #Momento de inercia (m^4)
        #G = 81.2e9      #Modulo de Cisalhamento (Pa)
        #A= 0.0125       #Área da seção do elemento (m^2)    
        #J = I/2         #Momento polar de inércia (m^4)
        KF_total = 0
        KT_total = 0
        KF_elements = []
        KT_elements = []    
        torcao, deformacao, flexao1, flexao2, flexao3 = [], [], [], [], []
        for element in self.elements:
            d = self.obter_propriedades(element[2])[2]      #Diâmetro do Tubo (m)
            e = self.obter_propriedades(element[2])[3]      #Espessura do Tubo (m)
            E = self.obter_propriedades(element[2])[0]      #Modulo de Young (Pa)      
            G = self.obter_propriedades(element[2])[1]      #Modulo de Cisalhamento (Pa)
            A = self.area_seccao_transversal(d, e)          #Área da seção do elemento (m^2)
            I = self.momento_inercia_area_e_polar(d,e)[0]   #Momento de inercia (m^4)
            J = self.momento_inercia_area_e_polar(d,e)[1]   #Momento polar de inércia (m^4)
            L_e = self.calcular_comprimento(element)
            # Equação de torsão
            torcao_val = (F_torcao * L_e) / (G * J)        #Fonte[1]
            torcao.append(torcao_val)
            # Equação  para deformação axial
            deformacao_val = (F_axial* L_e / (A * E))      #Fonte[2]
            deformacao.append(deformacao_val)
            # Equação para flexão
            flexao_val1 = (F_flexao1*L_e**3)/(48 * E * I)       #Fonte[3.1] (carga pontual no meio do elemento biapoiado)
            flexao_val2 = (5*F_flexao2*L_e**4)/(384 * E * I)    #Fonte[3.2] (carga distribuída ao longo de todo o elemento biapoiado)
            flexao_val3 = flexao_val1 + flexao_val2             #Fonte[3.3] (tentativa de carregamento misto)
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

    def calcular_B_Elementar(self, element):
        """
        Calcula a matriz B (6x12) para um elemento de viga de Timoshenko
        com o eixo longitudinal no eixo Y.
        """
        L = self.calcular_comprimento(element)
        B = [   [0,         -1/L,        0,         0,         0,         0,         0,         1/L,         0,         0,         0,         0      ], # ε_yy (axial)
                [0,         0,         0,      -1/L,         0,         0,         0,         0,         0,       1/L,         0,         0      ], # κ_x (curvatura sobre X)
                [0,         0,         0,         0,         0,      -1/L,         0,         0,         0,         0,         0,       1/L      ], # κ_z (curvatura sobre Z)
                [-1/L,         0,         0,         0,         0,      -0.5,       1/L,         0,         0,         0,         0,       -0.5       ], # γ_xy (cisalhamento XY)
                [0,         0,      -1/L,       0.5,         0,         0,         0,         0,       1/L,       0.5,         0,         0      ], # γ_zy (cisalhamento ZY)
                [0,         0,         0,         0,      -1/L,         0,         0,         0,         0,         0,       1/L,         0      ]]    # φ_y (torção)

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
                [0.0,  0.0,  0.0,  0.0,  0.0,  G    ]    # torção φ_y → τ_y = G * φ_y
            ])

            stress = C @ strains[i]
            stresses.append(stress)
        return stresses
        
    def compute_Kf_Kt(self, K_global):
        """
        Calculates
        Inputs:
        - K_global: Global stiffness matrix
        Outputs:
        - Kt: Torsional stiffness of the chassis
        - Kf: Beaming stiffness of the chassis
        """
        #Start Kt simulation
        F_global = np.zeros(K_global.size)
        F_global[2+9*6] = 250   #Forças aplicadas nos nós onde estaria a suspensão dianteira (nodes 8 e 9)
        F_global[2+8*6] = -250  #Mesmo módulo para gerar um torque no eixo longitudinal do chassi
        
        fixed_nodes=[26, 27]    #Fixação dos nós onde estaria a suspensão traseira
        fixed_dofs=[(node*6+i) for node in fixed_nodes for i in range(6)]   #Lista com os dofs fixados
        
        displacements = self.static_analysis(F_global, fixed_dofs)  #Calcula os displacements com as condições de contorno aplicadas acima
        mi1=np.abs(displacements[9*6+2])                             #displacement do nó 9 em módulo
        mi2=np.abs(displacements[8*6+2])                             #displacement do nó 8 em módulo
        L = np.abs(self.nodes[9][1] - self.nodes[8][1])              #Distancia entre o nó 8 e 9
        alpha= np.degrees(np.atan((mi1+mi2)/(L)))                   #Ângulo de torção do chassi após aplicação do torque
        tau = (np.abs(F_global[2+8*6]))*(L)                          #Cálculo do torque aplicado
        
        Kt = tau/alpha

        #Start Kf simulation
        F_global_2 = np.zeros(K_global.size)
        F = 5000    #Módulo da força que vai gerar a flexão
        F_global_2[2+22*6] = -F/2   #Força distribuída nos nós centrais do chassi (nodes 22 e 23)
        F_global_2[2+23*6] = -F/2   #Sinal negativo por conta da direção de aplicação da força
        
        fixed_nodes=[8, 9, 26, 27]  #Fixação dos nós onde estaria a suspensão dianteira e traseira
        fixed_dofs=[(node*6+i) for node in fixed_nodes for i in range(6)]   #Lista com os dofs fixados
        
        displacements_f = self.static_analysis(F_global_2, fixed_dofs)  #Calcula os displacements com as condições de contorno aplicadas acima
        dY=np.abs(displacements_f[2+22*6])  #Deslocamento em Y de um dos nós onde foi aplicado a força

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

            for i, (start, end, _) in enumerate(self.elements): # Alterado para ignorar o 'tipo'
                geo_file.write(f'Line({i + 1}) = {{{start + 1}, {end + 1}}};\n')

            if len(self.elements) > 2: # Esta parte pode precisar de ajustes dependendo da sua topologia
                line_loop_indices = ', '.join(str(i + 1) for i in range(len(self.elements))) # Usar self.elements
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
        for node1_idx, node2_idx, _ in elements: # Desempacotar corretamente
            node1 = nodes[node1_idx]
            node2 = nodes[node2_idx]
            
            # Use average scalar value for element color for simplicity
            # Or you might want to assign scalar values to elements directly
            # For now, let's assume scalar_values corresponds to elements
            # You might need to adjust this logic based on how scalar_values are structured
            try:
                element_scalar_value = scalar_values[elements.index((node1_idx, node2_idx, _))] # This assumes scalar_values matches elements order
            except ValueError:
                # If element tuple is not directly in scalar_values index,
                # you might need a different mapping or calculation.
                # For this example, let's just pick a default or a simplistic average
                element_scalar_value = np.mean([scalar_values[node1_idx], scalar_values[node2_idx]])


            color = cmap(norm(element_scalar_value))
            
            ax.plot([node1[0], node2[0]], 
                    [node1[1], node2[1]], 
                    [node1[2], node2[2]], 
                    color=color, linewidth=2)

        # Add color bar
        m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        m.set_array(scalar_values)
        cbar = fig.colorbar(m, ax=ax, pad=0.1)
        cbar.set_label(scalelabel)

        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title(graphtitle)
        plt.tight_layout()
        plt.show()

# --- Carregar os dados de nós e elementos ---
df_nodes = pd.read_csv("solucao_otimizada_nodes.csv")
nodes = df_nodes[['x', 'y', 'z']].to_numpy()

# Carregar elementos (certifique-se de que 'perfil' esteja no CSV)
df_elements = pd.read_csv("solucao_otimizada_elements.csv")
elements = [(row.node_i, row.node_j, row.perfil) for row in df_elements.itertuples()]

# --- Executando as novas funções ---

# Calcular e exibir comprimentos e quantidades por tipo de tubo
print("--- Análise de Comprimentos e Quantidades por Tipo de Tubo ---")
resultados = calcular_comprimentos_e_quantidades_por_tipo(nodes, elements)

for tipo, dados in resultados.items():
    print(f"{tipo}:")
    print(f"  - Comprimento total: {dados['comprimento_total']:.3f} metros")
    print(f"  - Quantidade de tubos: {dados['quantidade']}")
print("------------------------------------------------------------\n")

# Plotar o chassi em 3D com cores por tipo de tubo
#plotar_chassi_3d(nodes,elements)  
chassi = Chassi3D(nodes, elements)
chassi.plotar_chassi_3d()

# Exemplo de uso da classe Estrutura (mantido do seu código original)
estrutura = Estrutura(elements, nodes)
K_global, M_global = estrutura.matrizes_global()

# Exemplo de uso da função plot_colored_wireframe (se você tiver dados de tensão, por exemplo)
# Para usar plot_colored_wireframe, você precisa de um array de valores escalares
# Vamos criar um exemplo de valores de tensão von Mises, assumindo que você já calculou isso.
# Por exemplo, após uma análise estática:
# F_exemplo = np.zeros(estrutura.num_dofs)
# F_exemplo[2] = 1000 # Aplica uma força de exemplo
# fixed_dofs_exemplo = [0, 1, 2, 3, 4, 5] # Fixa o nó 0
# displacements_exemplo = estrutura.static_analysis(F_exemplo, fixed_dofs_exemplo)
# strains_exemplo = estrutura.compute_strain(displacements_exemplo)
# stresses_exemplo = estrutura.compute_stress(strains_exemplo)
# von_mises_exemplo = estrutura.compute_von_mises(stresses_exemplo)

# Se você tiver os resultados de von Mises para cada elemento, pode plotar assim:
# plot_colored_wireframe(nodes, elements, von_mises_exemplo, 
#                        graphtitle='Tensões de Von Mises no Chassi', 
#                        scalelabel='Tensão Von Mises (MPa)', 
#                        colormap='viridis') # Escolha um colormap adequado

# Continue com o restante do seu código principal...
# ...