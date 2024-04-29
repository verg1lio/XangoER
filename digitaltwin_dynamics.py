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



    def __init__(self, spring_type, spring_k, spring_F, spring_non_lin_coef, tire_Fz, tire_Sa, tire_Ls):
        # Modelo de mola
        self.spring_type = spring_type # Hooke, 
        self.spring_k = spring_k # rigidez da mola [N/m]
        self.spring_F = spring_F # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef # coeficiente de ganho não-linear
        # Modelo de pneu
        self.tire_Fz = tire_Fz  # carga vertical no pneu [N]
        self.tire_Sa = tire_Sa  # slip angle do pneu [rad]
        self.tire_Ls = tire_Ls  # longitudinal slip do pneu [Admensional]
        self.tire_type = 'Default'
        

    def Spring(self):

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
        # Pacejka parâmetros
        E = -1
        Cy = 1.4  # C para força lateral
        Cx = 1.65  # C para força longitudinal
        Cz = 1  # C para momento de torque auto-alinhante
        c1 = 54000
        c2 = 6600
        Cs = c1 * np.sin(2 * np.arctan(self.tire_Fz / c2))

        if self.tire_type == 'Default':
            D = 1.45 * self.tire_Fz
            Bz = Cs / (Cz * D)
            Bx = Cs / (Cx * D)
            By = Cs / (Cy * D)
            tire_lateral_force = D * np.sin(Cy * np.arctan(By * self.tire_Sa - E * (By * self.tire_Sa - np.arctan(By * self.tire_Sa))))
            tire_longitudinal_force = D * np.sin(Cx * np.arctan(Bx * self.tire_Ls - E * (Bx * self.tire_Ls - np.arctan(Bx * self.tire_Ls))))
            tire_auto_align_moment = D * np.sin(Cz * np.arctan(Bz * self.tire_Sa - E * (Bz * self.tire_Sa - np.arctan(Bz * self.tire_Sa))))

        return tire_lateral_force, tire_longitudinal_force, (tire_auto_align_moment/55)



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


