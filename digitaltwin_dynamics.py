class Dynamics:
    """A Vehicle Dynamics object.

    This class will create a Vehicle Dynamics digital twin with the springs, suspension bars, damper cylinder, brakes, tires, and transmission elements provided.

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



     def __init__(self, spring_type, spring_k, spring_F, spring_non_lin_coef, tire_Fz, tire_Sa, tire_Ls, tire_type):
        # Modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # rigidez da mola [N/m]
        self.spring_F = spring_F  # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # coeficiente de ganho não-linear
        # Modelo de pneu
        self.tire_Fz = tire_Fz  # carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # slip angle do pneu [rad]
        self.tire_Ls = tire_Ls  # longitudinal slip do pneu [Admensional]
        self.tire_type = tire_type  # Default, Admensional
        

    def Spring(self):
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
        if self.spring_type == 'Hooke':
            spring_x = self.spring_F/self.spring_k
        if self.spring_type == 'Softening'
            spring_x = self.spring_F/(self.spring_non_lin_coef*(self.spring_k)**2)

        return spring_x



    def Suspension(self):
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
        
        return Sus



    def Damper(self):
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
        
        return Dam

    

    def Brake(self):
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
        
        return Bra



    def Tire(self):
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
        
        return Tir



    def Transmission(self):
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
        
        return Tra


