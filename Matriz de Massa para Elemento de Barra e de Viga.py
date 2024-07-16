#Inicializando as bibliotecas
import numpy as np
import sympy as sp

#Determinando o número de elementos
n_elem = int(input(print("Insira o número de elementos de barra: ")))

#Matriz de elemento de barra
def bar_mass_matrix(n_elem, rho, A, L):
    #Matriz de massa elementar
    m = (rho*A*L/6*n_elem) * np.array([[2,1], [1,2]])

    #Inicializando a matriz global de massa
    M = np.zeros((n_elem+1, n_elem+1))          #Matriz cheia de zeros

    #Juntando os elementos:
    for i in range(n_elem):
        M_montando = np.zeros((n_elem+1, n_elem+1))

        M_montando[i:i+2,i:i+2] = m

        M += M_montando

    return M

#Matriz de elemento de viga
#Determinando o número de elementos
n_elem = int(input(print("Insira o número de elementos de barra: ")))

def beam_mass_matrix(n_elem, rho, A, L):
    #Matriz de massa elementar
    m = (rho*A*L/420*n_elem) * np.array([[156, 22*L, 54, -13*L], [22*L, 4*L*L, 13*L, -3*L*L], [54, 13*L, 156, -22*L], [-13*L, -3*L*L, -22*L, 4*L*L]])

    #Inicializando a matriz global de massa
    M = np.zeros((2*n_elem+2, 2*n_elem+2))          #Matriz cheia de zeros

    #Juntando os elementos:
    for i in range(n_elem):
        M_montando = np.zeros((2*n_elem+2, 2*n_elem+2))

        M_montando[2*i:2*i+4,2*i:2*i+4] = m

        M += M_montando

    return M

#Parâmetros da barra
E = 210e9               #Módulo de Elasticidade lontigudinal, em Pa
A = 0.225               #Área da seção transversal, em m²
rho = 7850                #Massa específica do aço, em kg/m³
L = 1                   #Comprimento da barra, em m

#Criando a matriz global de massa para viga
M = beam_mass_matrix(n_elem, rho, A, L)

print(M)