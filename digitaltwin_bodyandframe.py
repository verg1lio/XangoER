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
         self.A=A
         self.E=E
         self.L=L
         self.I=I
        
 def bar(self, E, A, L):
#matriz nodal de elemento de barra
#Ainda precisa integrar o angulo de torção na matriz
        self.result= np.array([[E*A/L, -E*A/L],
                               [-E*A/L, E*A/L]])
        
        return self.result 
        
        
    def beam(self, E, I, L):
#matriz nodal elemento de viga
        self.result= np.array([[12*E*I/L**3, 6*E*I/L**2, -12*E*I/L**3, 6*E*L/L**2],
                               [6*E*L/L**2, 4*E*I/L, -6*E*L/L**2, 2*E*I/L],
                              [-12*E*I/L**3, -6*E*L/L**2, 12*E*I/L**3, -6*E*L/L**2],
                              [6*E*I/L**2, 2*E*I/L, -6*E*L/L**2, 4*E*I/L ]])
        
        return self.result


