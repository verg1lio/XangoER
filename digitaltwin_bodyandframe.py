import numpy as np
import matplotlib.pyplot as plt


class BodyAndFrame:
    
    """
    parameters:
    ________________
    
    A: cross_section area
    E: elastic modulos
    L: length of beam
    I: second momento of area 
    G:shear modulus of material
    a:coefficient
    B:coefficient
    K: shear strain e x 
    y= parameter
    
    return: elemental sitffiness beam
    """
    


    def __init__(self,A, E , L, I, G, a, B, K, y):
        self.A = A
        self.E = E
        self.L = L
        self.I = I
        self.G - G
        self.a = a
        self.B = B
        self.K = K
        self.y = y
        

    def bar(self, E, A, L):
        self.result= np.array([K*G*J/L, 0],
                             [0, K*G*J/L])
        
        return self.result
        
        
    def beam(self, E, I, L, A, a):
        a=1/1+y
        B= y/4(1+y)
        y=12EI/KGAL
        
        self.result= np.array([[12*E*I*a/L**3, -6*E*I*a/L**2, -12*E*I*a/L**3, -6*E*L*a/L**2],
                               [-6*E*L*a/L**2, 4*E*I*(a+B)/L, 6*E*L*a/L**2, 2*E*I*(a-2B)/L],
                              [-12*E*I*a/L**3, -6*E*L*a/L**2, 12*E*I*a/L**3, 6*E*L*a/L**2],
                              [-6*E*I*a/L**2, 2*E*I*(a-2B)/L, 6*E*L*a/L**2, 4*E*I*(a+B)/L ]])
        
        return self.result
        


    def ElementMat(self, A, E, L, I):
        self.result=np.array([[K*G*J/L,0,0,0,0,0],
                             [0, 12*E*I*a/L**3, -6*E*I*a/L**2, -12*E*I*a/L**3, -6*E*L*a/L**2]])
        
        
        return self.result



    def GlobalMat(self):
        num_nodes=
        num_liber_per_node=
        num_liber= num_liber_per_node*num_node
        num_elements=
        element_conec=
        
        K_global= np.zeros((num_liber, num_liber))
        
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
