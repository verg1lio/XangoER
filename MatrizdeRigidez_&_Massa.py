import numpy as np
import matplotlib.pyplot as plt
import math

#Aço 1020
#densidade: 7,87 g/cm3, 7870 kg/m3 ou 0,284 lb/in3
#Modulo da elasticidade: 207 GPa = 207*10^9Pa, 30*10^6psi, razão de poisson 0,30

class Barra(): # para 2 graus de liberdade

    def __init__(self,area,elasticidade=1,comprimento=1,angulo=0,densidade=207,inercia=1):
        #Declaração dos atributos 

        self.area= area
        self.elasticidade = elasticidade
        self.comprimento = comprimento
        self.densidade = densidade
        self.angulo =  math.radians(angulo)
        self.inercia = inercia

        #Variveis para facilitar o calculo das matrizes
        self.cos_teta= math.cos(self.angulo)
        self.sen_teta= math.sin(self.angulo)
        self.c= (self.elasticidade*self.area)/(self.comprimento)
        self.m= (1/420)*self.densidade*self.area*self.comprimento

    def matriz_rigidez_local_barra(self): # Matriz de rigidez elementar da barra
        K=self.c*np.array(([1,0,-1,0],[0,0,0,0],[-1,0,1,0],[0,0,0,0]))
        return K  
        
    def matriz_rotacao(self): # Matriz de rotação 
        T = np.array([[self.cos_teta, self.sen_teta, 0, 0],
                    [-self.sen_teta, self.cos_teta, 0, 0],
                    [0, 0, self.cos_teta, self.sen_teta],
                    [0, 0, -self.sen_teta, self.cos_teta]])
        return T
    
    def matriz_rigidez_global_barra(self): # Matriz de rigidez global do elemento de barra 2D, Equação retirada do documento CILAMCE 2016
        K_e = self.c*np.array([[self.cos_teta**2, self.sen_teta*self.cos_teta, -self.cos_teta**2, -self.sen_teta * self.cos_teta],
                                   [0, self.cos_teta**2, -self.sen_teta*self.cos_teta, -self.sen_teta**2],
                                   [0, 0, self.cos_teta**2, self.sen_teta * self.cos_teta],
                                   [0, 0, 0, self.sen_teta**2]])
        K_e[np.tril_indices(4, k=-1)] = K_e.T[np.tril_indices(4, k=-1)]        # Preenchendo a parte inferior esquerda da matriz (simétrica)
        return K_e
    
    def matriz_rigidez_global_barra_livro(self): # Matriz de rigidez global do elemento de barra 2D, Equação retirada do livro que virgilio forneceu (The finite element method)
        K_e_livro = ((2*self.elasticidade*self.inercia)/self.comprimento**3)*np.array([[6, -3*self.comprimento, -6, -3*self.comprimento],
                                   [-3*self.comprimento, 2*self.comprimento**2, 3*self.comprimento, self.comprimento**2],
                                   [-6, 3*self.comprimento, 6, 3*self.comprimento],
                                   [-3*self.comprimento, self.comprimento**2, 3*self.comprimento, self.comprimento**2]])
        return K_e_livro

    def matriz_massa_elemento(self): # matriz de massa do elemento/ equação encontrada na pagina 311 do livro (The finite element method)
        M_e =self.m *np.array([[156, 22 * self.comprimento, 54, -13 * self.comprimento],
                                [22 * self.comprimento, 4 * self.comprimento**2, 13 * self.comprimento, -3 * self.comprimento**2],
                                [54, 13 * self.comprimento, 156, 22 * self.comprimento],
                                [-13 * self.comprimento, -3 * self.comprimento**2, 22 * self.comprimento, 4 * self.comprimento**2]])
        return M_e
    
#elemento de barra=Barra(Area,Elasticidade,Comprimento,angulo,densidade,inercia)

elemento1=Barra(0.008,207,2,60,78,1)

matriz_massa=elemento1.matriz_massa_elemento()
print(elemento1.matriz_rigidez_global_barra())
print(elemento1.matriz_massa_elemento())


# Plotando o gráfico da matriz de massa
plt.imshow(matriz_massa, cmap='viridis', interpolation='nearest')
plt.title('Matriz de Massa do Elemento')
plt.colorbar(label='Valores da Matriz')
plt.xlabel('Coluna')
plt.ylabel('Linha')
plt.show()