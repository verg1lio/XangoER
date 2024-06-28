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
        self.nos = nos #Number of chassi's nodes
        self.massa = massa #Car mass (Kg)
        self.momento_inercia_direção = momento_inercia_direção #Moment of inertia relative to direction (Kg.m^2)
        self.momento_inercia_plano = momento_inercia_plano #Moment of inertia relative to the plane (Kg.m^2)
        

    def Connect(self):
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
        
        return Con



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




