class BodyAndFrame:
    """A Body And Frame object.

    This class will create a Body And Frame digital twin with the chassis elements provided.

    Parameters AINDA A ALTERAR, ISSO É UM EXEMPLO
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




    def __init__(self, nos, massa, momento_inercia_direção, momento_inercia_plano):
        self.I = I #Tirar I e E e colocar em elemento de Barra e viga numa nova def
        self.E = E
        self.num_nodes = n  # Número total de nós
        self.massa = m  # Massa do carro (Kg)
        self.momento_inercia_direcao = iD  # Momento de inércia em relação à direção (Kg.m^2)
        self.momento_inercia_plano = IP  # Momento de inércia em relação ao plano (Kg.m^2)
        self.area_transversal = area_transversal  # Área da seção transversal da Estrutura
        self.num_dofs_per_node = 6  # Graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node  # Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))  # Matriz de rigidez global
        
        

    def Connect(self):
         
         def node_loc_matrix(self, node_tags, node_coord):
        # Número de Nós
        num_nodes = len(node_tags)

        # Gerando uma matrix número de nos x 4: (Para cada linha, teremos: [índice do nó, x, y, z])
        node_loc_matrix = np.zeros((num_nodes, 4),
                                   dtype=float)  # Na primeira coluna, teremos os índices, nas seguintes, x y e z.

        # Preenchendo a matriz de zeros
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i]  # Número do nó na primeira coluna
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
            CM[i][1] = no1
            CM[i][2] = no2

        print("\n Conexão   1º Nó   2º Nó")
        print(CM)




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




