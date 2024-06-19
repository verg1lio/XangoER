import numpy as np
import matplotlib.pyplot as plt
import math

class forcas_aerodinamicas():
    def __init__(self,frontal_area,planform_area,air_density,coeficient_drag_p,coeficient_drag_f,coeficient_lift,length,ni,atm_pressure,attack_angle,span,chord):
        self.af= frontal_area
        self.asup= planform_area
        self.rho= air_density
        self.cd_p= coeficient_drag_p
        self.cd_f= coeficient_drag_f
        self.cl= coeficient_lift
        self.x= length
        self.ni= ni
        self.a1= a1
        self.a2= a2
        self.p_atm= atm_pressure
        self.alpha= attack_angle
        self.a_wing=span*chord
        self.AR= (span**2)/self.a_wing      
        

    def aerodynamic_forces(self,velocidade):                        # Função para o cálculo de forças aerodinâmicas
        pressao_dinamica=0.5*self.rho*(velocidade**2)
        drag= (self.cd_p * pressao_dinamica * self.af) + (self.cd_f* pressao_dinamica * self.asup) 
        lift=self.cl * pressao_dinamica * self.asup
        return drag,lift
    
#    def numero_Reynolds(self,velocidade):         # Função para encontrar o numero de reynolds
        return (self.x*velocidade)/self.ni
    
#    def bernoulli(self,velocidade):
        v2=velocidade*self.a1/self.a2                                   # Equação de conservação para encontrar a velocidade embaixo do carro
        p2=self.p_atm +0.5*rho*velocidade**2 - 0.5*v2**2                # Equação de bernoulli para encontrar a pressão embaixo do carro
        deltap=self.p_atm-p2                                            # Diferença de pressão entre a parte de baixo e a parte de cima do carro
        return deltap
    
    def aerodynamics_forces_wing(self,velocidade):                         #Função para calcular as forças aerodinamicas geradas pela asa
        pressao_dinamica=0.5*self.rho*(velocidade**2)
        cl_alpha=(2*math.pi)/(1+(2/self.AR))
        cl_wing=cl_alpha*self.alpha                              #Coeficiente de lift da asa
        cd_induced=(1/(math.pi*self.AR))*cl_wing**2    #Coeficiente de arrasto induzido
        cd_wing=cd_induced+0.05                             #Coeficiente de arrasto total da asa= Cd=Cd_induzido+Cd_forma
        drag_wing= (cd_wing * pressao_dinamica * self.a_wing)   #Calculo da força de arrasto
        lift_wing= (cl_wing * pressao_dinamica * self.a_wing)
        return drag_wing,lift_wing
    
#Parâmetros físicos do ar:
p_atm= 101325                    #pressâo atmosférica (N/m^2)
rho= 1.184                       #densidade (kg/m^3)
mi= 1.849*10**-5                 #Viscosidade dinâmica (kg/m*s)
ni= (mi/rho)                     #Viscosidade Cinematica (m²/s)
#parâmetros do carro
length= 2                        #Comprimento do carro (m^2)
af= 1.5                          #Área Frontal do carro (m^2)
a_sup= 2                         #Área de Superfície do carro (m^2)
cd_p = 0.45                      #Coeficiente de arrasto por pressão do carro
cd_f = 0.05                      #Coeficiente de arrasto por atrito do carro
cl= -0.3                         #Coeficiente de lift do carro
a1= 0.25                         #Área de entrada do ar embaixo do carro (m^2)
a2= 0.20                         #Área embaixo do carro (m^2)
#parâmetros da asa
chord=0.25                          #Comprimento da asa (m)
span=1                              #Largura da asa (m)
thickness=0.05                      #Expessura máxima da asa (m)
alpha=math.radians(3.75)            #Ângulo de incidencia do vento com a asa (radianos)

carro=forcas_aerodinamicas(af,a_sup,rho,cd_p,cd_f,cl,length,ni,p_atm,alpha,span,chord)    #Definição do objeto 

# Velocidades de 0 a 60 m/s
velocidades = np.linspace(0.1,60,50)
#calculo das forças
drags=[]
lifts=[]
drags_w=[]
lifts_w=[]
for v in velocidades:
    drag,lift= carro.aerodynamic_forces(v)
    drag_w,lift_w= carro.aerodynamics_forces_wing(v)
    drags.append(drag)
    lifts.append(lift)
    drags_w.append(drag_w)
    lifts_w.append(lift_w)

# Plotagem do gráfico 1
plt.figure(figsize=(12, 6))
plt.subplot(2,2,1)
plt.plot(velocidades, drags, label='Drag')
# Configurações do gráfico
plt.title('Gráfico de Arrasto do carro em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
# Plotagem do gráfico 2
plt.subplot(2,2,2)
plt.plot(velocidades, lifts, label='Downforce', linestyle='--')
# Configurações do gráfico2
plt.title('Gráfico de Downforce do carro em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
# Plotagem do gráfico 3
plt.subplot(2,2,3)
plt.plot(velocidades, drags_w, label='Drag')
# Configurações do gráfico
plt.title('Gráfico de Arrasto da Asa em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
# Plotagem do gráfico 4
plt.subplot(2,2,4)
plt.plot(velocidades, lifts_w, label='Lift', linestyle='--')
# Configurações do gráfico
plt.title('Gráfico de Lift da Asa em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()