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
    
    
    return: elemental sitffiness beam
    """
    


    def __init__(self,A, E, L, I):
        self.A = A
        self.E = E
        self.L = L
        self.I = I
        

    def bar(self, E, A, L):
        self.result= np.array([[E*A/L, -E*A/L],
                               [-E*A/L, E*A/L]])
        
        return self.result
        
        
    def beam(self, E, I, L):
        self.result= np.array([[12*E*I/L**3, 6*E*I/L**2, -12*E*I/L**3, 6*E*L/L**2],
                               [6*E*L/L**2, 4*E*I/L, -6*E*L/L**2, 2*E*I/L],
                              [-12*E*I/L**3, -6*E*L/L**2, 12*E*I/L**3, -6*E*L/L**2],
                              [6*E*I/L**2, 2*E*I/L, -6*E*L/L**2, 4*E*I/L ]])
        
        return self.result
        


    def ElementMat(self, A, E, L, I):
        self.result=np.array([[A*E/L,0,0,-A*E/L,0,0],
                      [0,12*E*I/L**3,6*E*I/L**2,0,-12*E*I/L**3,6*E*I/L**2],
                      [0,6*E*I/L**2,4*E*I/L,0,-6*E*I/L**2,2*E*I/L],
                      [-A*E/L,0,0,A*E/L,0,0],
                      [0,-12*E*I/L**3,-6*E*I/L**2,0,12*E*I/L**3,-6*E*I/L**2],
                      [0,6*E*I/L**2,2*E*I/L,0,-6*E*I/L**2,4*E*I/L]])
        
        
        return self.result



    def GlobalMat(self):
        num_nodes=
        num_liber_per_node=
        num_liber= num_liber_per_node*num_node
        num_elements=
        element_conec=
        
        K_global= np.zeros((num_liber, num_liber))
