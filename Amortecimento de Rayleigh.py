from numpy import pi
import sympy as sp

#Definindo as variáveis simbólicas
zeta1 = sp.Symbol('Z1')         #Coef de amortecimento do primeiro modo de vibração
zeta2 = sp.Symbol('Z2')         #Coef de amortecimento de segundo modo de vibração
freq_nat1 = sp.Symbol('W1')     #Primeiro modo de vibração
freq_nat2 = sp.Symbol('W2')     #Segundo modo de vibração
alfa = sp.Symbol('a')           #Coeficiente de proporcionalidade de Rayleigh alfa
beta = sp.Symbol('b')           #Coeficiente de proporcionalidade de Rayleigh beta

#Equações para determinação dos coeficientes alfa e beta
eq_zeta1 = sp.Eq(((1/2)*(alfa*freq_nat1 + (beta/freq_nat1))), zeta1)
eq_zeta2 = sp.Eq(((1/2)*(alfa*freq_nat2 + (beta/freq_nat2))), zeta2)

coeficientes = sp.solve((eq_zeta1, eq_zeta2), (alfa, beta))

#Matrizes de Rigidez e Massa
stiff_matrix = sp.MatrixSymbol('K', 2, 2)       #Inicializando uma matriz de símbolos

K = sp.Matrix(stiff_matrix)                     #Matriz simbólica de rigidez

K = K.subs([(stiff_matrix[0, 0], sp.Symbol('k11')),
            (stiff_matrix[0, 1], sp.Symbol('k12')),
            (stiff_matrix[1, 0], sp.Symbol('k21')),
            (stiff_matrix[1, 1], sp.Symbol('k22'))])

mass_matrix = sp.MatrixSymbol('M', 2, 2)       #Inicializando uma matriz de símbolos

M = sp.Matrix(mass_matrix)                     #Matriz simbólica de massa

M = M.subs([(mass_matrix[0, 0], sp.Symbol('m11')),
            (mass_matrix[0, 1], sp.Symbol('m12')),
            (mass_matrix[1, 0], sp.Symbol('m21')),
            (mass_matrix[1, 1], sp.Symbol('m22'))])

#Calculando uma matriz genérica de amortecimento de Rayleigh

C = coeficientes[alfa]*K + coeficientes[beta]*M

print("O coeficiente de proporcionalidade alfa é:")
print(coeficientes[alfa])
print("")
print("O coeficiente de proporcionalidade beta é:")
print(coeficientes[beta])
print("")
print("A matriz de rigidez 2x2 genérica é:")
print(K)
print("")
print("A matriz de massa 2x2 genérica é:")
print(M)
print("")
print("A matriz de amortecimendo 2x2 genérica é:")
print(C)

#Para o funcionamento apropriado do código acima, associar com as matrizes de massa e rigidez do chassi
#e realizar a análise modal para termos os coeficientes de amortecimento e as frequências naturais.