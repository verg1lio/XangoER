import numpy as np
import matplotlib.pyplot as plt
import math as m

class BodyAndFrame:

    """parameters
    __________________
    
    A: cross section area
    E: elastic modulos
    L: length
    I: second moment area"""
    


     def __init__(self,A, E, L, I):
        self.A = float(A)
        self.E = float(E)
        self.L = float(L)
        self.I = float(I)
        
 def bar(self, E, A, L)
#matriz nodal de elemento de barra
#Ainda precisa integrar o angulo de torção na matriz
        self.result= np.array([[E*A/L, -E*A/L],
                               [-E*A/L, E*A/L]])
        
        return bar 
        
        
    def beam(self, E, I, L)
#matriz nodal elemento de viga
        self.result= np.array([[12*E*I/L**3, 6*E*I/L**2, -12*E*I/L**3, 6*E*L/L**2],
                               [6*E*L/L**2, 4*E*I/L, -6*E*L/L**2, 2*E*I/L],
                              [-12*E*I/L**3, -6*E*L/L**2, 12*E*I/L**3, -6*E*L/L**2],
                              [6*E*I/L**2, 2*E*I/L, -6*E*L/L**2, 4*E*I/L ]])
        
        return beam
    



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




