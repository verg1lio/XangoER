class BodyAndFrame:
    """A Body And Frame object.

    This class will create a Body And Frame digital twin with the chassis elements provided.

    Parameters AINDA A ALTERAR, ISSO É UM EXEMP
    ----------
    shaft_elements : list
        List with the shaft elements.
    disk_elements : list
        List with the disk elements.
    bearing_elements : list
        List with the bearing elements.
    automeshing : boolean
        Set it True to use the automeshing method. Default is False.
        If automeshing is True, the previous shaft_elements parameter is now defined by:
            shaft_elements : list
                List with the length, inner and outter diameter and material of each element, as follows:
                    [[length, inner, outter, material],[length, inner, outter, material],...]
        For the other parameters please check the respective bearing and disk classes for more information.
    **kwargs : dict, optional
        If automeshing is True, these parameters needs to be informed.
        The possible arguments are:
            alpha : float
                Proportional damping coefficient, associated to the element Mass matrix
            beta : float
                Proportional damping coefficient, associated to the element Stiffness matrix

    Returns
    -------
    A rotor object.

    Attributes
    ----------
    MM : array
        Global mass matrix.
    KK : array
        Global stiffness matrix.
    CCgyros: array
        Global gyroscopic matrix.
    CCtotal: array
        Global damping matrix

    Examples
    --------
    >>> import lmest_rotor as lm
    >>> rotor = lm.rotor_example()
    >>> rotor.MM
    array(30x30)
    """




    def __init__(self, n, m, Id, Ip):
        self.n = int(n)
        self.n_l = n
        self.n_r = n

        self.m = float(m)
        self.Id = float(Id)
        self.Ip = float(Ip)
        

    import numpy as np
import matplotlib.pyplot as plt

def matriz_conectividade(nos, elementos):
    num_nos = len(nos)
    matriz = np.zeros((num_nos, num_nos), dtype=int)

    for no1, no2 in elementos:
        if no1 in nos and no2 in nos:
            j1 = nos.index(no1)
            j2 = nos.index(no2)
            matriz[j1, j2] += 1
            matriz[j2, j1] += 1

    return matriz

# Definindo os nós
nos = [1, 2, 3, 4, 5]

# Definindo os elementos (arestas)
elementos = [(1, 2), (1, 4), (2, 3), (2, 5), (3, 5), (4, 5), (1, 2), (2, 4)]

# Definindo as coordenadas dos nós
coordenadas = {
    1: (0, 0),  # No 1 na coordenada (0, 0)
    2: (1, 0),  # No 2 na coordenada (1, 0)
    3: (2, 0),  # No 3 na coordenada (2, 0)
    4: (0, 1),  # No 4 na coordenada (0, 1)
    5: (1, 1)   # No 5 na coordenada (1, 1)
}

# Gerando a matriz de conectividade
matriz = matriz_conectividade(nos, elementos)

# Plotando o grafo
plt.figure(figsize=(8, 6))

# Plotando os nós
for no, (x, y) in coordenadas.items():
    plt.scatter(x, y, s=100)  # s é o tamanho do ponto
    plt.text(x, y, str(no), fontsize=12, ha='right')

# Plotando as arestas
for no1, no2 in elementos:
    x1, y1 = coordenadas[no1]
    x2, y2 = coordenadas[no2]
    plt.plot([x1, x2], [y1, y2], 'k-')  # 'k-' é a cor preta para as arestas

plt.title('Grafo dos Nós e Elementos')
plt.xlabel('Coordenada X')
plt.ylabel('Coordenada Y')
plt.grid(True)
plt.show()



    def ElementMat(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return Ele



    def GlobalMat(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return Glo

    

    def ShapeFun(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return SF



    def WeakForm(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return WF



    def Mesh(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return Mesh



    def EIGSolver(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return EIG



    def TIMESolver(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return TIME




