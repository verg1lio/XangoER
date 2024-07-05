import numpy as np
from scipy.linalg import eigh

#Função de análise modal a partir das matrizes de massa e rigidez globais
def modal_analysis(K, M, num_modes = 20):

    unsorted_eigenvalues, unsorted_eigenvectors = eigh(K, M)

    # Frequências naturais (raiz quadrada dos autovalores)
    unsorted_frequencies = np.sqrt(unsorted_eigenvalues) / (2*np.pi)            #Divisão por 2*pi para converter para hertz


    #Tratando os dados (tomando apenas as 20 primeiras frequências naturais)
    sorted_indices = np.argsort(unsorted_frequencies)                           #Ordena as frequências em ordem crescente
    top_indices = sorted_indices[:num_modes]                                    #Seleciona os índices dos primeiros n modos

    eigenvalues = np.array(unsorted_eigenvalues)[top_indices]                   #Filtra os primeiros n autovalores
    eigenvectors = np.array(unsorted_eigenvectors)[top_indices]                 #Filtra os primeiros n autovetores
    frequencies = np.array(unsorted_frequencies)[top_indices]                   #Filtra as primeiras n frequências

    return eigenvalues, eigenvectors, frequencies

#Inicializando matrizes aleatóras n*6 x n*6 para servirem como matrizes de massa e rigidez globais para teste
#Elas terão n_nós * 6 x n_nós * 6, pois cada nó terá 6 forças atuando nelas
n_nos = 40

valores_m = np.random.rand(n_nos*6)
M = np.diag(valores_m)                                          #Criando uma matriz diagonal de massa a partir de valores aleatórios de massa m

k_qualquer = np.random.randint(0, 1001, (n_nos*6, n_nos*6))     #Criando uma matriz aleatória de rigidez apenas com inteiros de 0 a 1000

K = (k_qualquer + k_qualquer.T) // 2                            #Forçando uma matriz simétrica de rigidez através do poder da álgebra

#Análise modal dos elementos da malha
autovalores, autovetores, frequências = modal_analysis(K, M)

#print("Autovalores (λ):")                                          #Autovalores não precisam ser exibidos
#print(autovalores)
#print("")

print("Frequências Naturais (ω):")
print(frequências)
print("")

print("Autovetores (Modos de Vibração):")
print(autovetores)

#COMENTÁRIOS FINAIS, PFV LEIA (PRINCIPALMENTE TU PATTY)
#No código real, é necessário que sejam colocadas as matrizes globais de massa e rigidez como entrada

#Cumprindo isso, dá certo e a gente só tava meio que se atrapalhando demais mesmo com teoria e o negócio era simples

#A versão complicada que faz barra a barra precisaria de mais coisa para funcionar, como por exemplo criar um código que que consulte um
#arquivo externo (como planilha do excel) para que ele possa pegar de elemento a elemento as características de inércia, comprimento,
#seção transversal, etc, para que a análise realmente seja feita sobre a estrutura que se quer representar.

#No final das contas, estávamos nos complicando demais sobre algo que realmente deveria ser simples e eu acho que boa parte disso foi
#eu induzindo o pessoal ao erro por causa de uma falha de entendimento teórico sobre coisa besta. Foi mal :(