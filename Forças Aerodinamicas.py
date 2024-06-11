import numpy as np
import matplotlib.pyplot as plt
import math

class forcas_aerodinamicas():
    def __init__(self,area_frontal,area_superior,densidade_ar,coeficient_drag,coeficient_lift,length,ni):
        self.af= area_frontal
        self.asup= area_superior
        self.rho= densidade_ar
        self.cd= coeficient_drag
        self.cl= coeficient_lift
        self.x= length
        self.ni=ni

    def aerodynamic_forces(self,velocidade):       # Função para o cálculo de forças aerodinamicas
        pressao_dinamica=0.5*self.rho*(velocidade**2)
        drag= self.cd * pressao_dinamica * self.af
        lift=self.cl * pressao_dinamica * self.asup
        return drag,lift
    
    def numero_Reynolds(self,velocidade):         # Função para encontrar o numero de reynolds
        return (self.x*velocidade)/self.ni


#Parâmetros físicos do ar:
p_atm= 101325                    #pressâo atmosférica (N/m^2)
rho= 1.184                       #densidade (kg/m^3)
mi= 1.849*10**-5                 #Viscosidade dinâmica (kg/m*s)
ni= (mi/rho)                     #Viscosidade Cinematica (m²/s)
#parâmetros do carro
length= 2                        #Comprimento do carro (m^2)
af= 1.5                          #Área Frontal do carro (m^2)
a_sup= 2                         #Área de Superfície do carro (m^2)
cd = 0.75                        #Coeficiente de arrasto por pressão do carro
cl= -0.3                         #Coeficiente de lift do carro

carro=forcas_aerodinamicas(af,a_sup,rho,cd,cl,length,ni) 

# Velocidades de 0 a 60 m/s
velocidades = np.linspace(0.1,60,50)
#calculo das forças
drags=[]
lifts=[]
for v in velocidades:
    drag,lift= carro.aerodynamic_forces(v)
    drags.append(drag)
    lifts.append(lift)

# Plotagem do gráfico 1
plt.figure(figsize=(10, 6))
plt.subplot(1,2,1)
plt.plot(velocidades, drags, label='Drag')
# Configurações do gráfico
plt.title('Gráfico de Arrasto em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)

# Plotagem do gráfico 2
plt.subplot(1,2,2)
plt.plot(velocidades, lifts, label='Downforce', linestyle='--')
# Configurações do gráfico2
plt.title('Gráfico de Downforce em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
plt.show()

