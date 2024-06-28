import numpy as np
import matplotlib.pyplot as plt
import math

class Arrefecimento():
    def __init__(self,area_resfriamento,area_entrada,densidade_ar,calor_especifico_ar,temp_ar,temp_objeto,velocidade,viscosidade_cinematica_ar,condutividade_ar,viscosidade_dinamica_ar):
        self.a_res= area_resfriamento
        self.a_ent= area_entrada
        self.rho= densidade_ar
        self.c= calor_especifico_ar
        self.v=velocidade
        self.tf=temp_ar
        self.to=temp_objeto
        self.mi=viscosidade_dinamica_ar
        self.ni=viscosidade_cinematica_ar
        self.k_ar=condutividade_ar
        self.Vazao= velocidade*area_entrada*densidade_ar

    def arrefecimento_pwt(self,comprimento,massa,condutividade,temp_desejada):
        Pr=(self.mi*self.c)/self.k_ar                   # Calculo do numero de Prandlt
        Rey=comprimento*self.v/self.ni              # Calculo do numero de Reynolds
    
        if Rey < 200000:
            Nul= (0.664 * (Rey**0.5) * (Pr**(1/3)))
        else:
            X=200000*self.ni/self.v                 #Comprimento da placa onde o fluxo transiciona para turbulento
            Rey_X =X*self.v/self.ni                 # Calculo do numero de Reynolds para a distância X
            A= 0.037*(Rey_X**0.8) - 0.664*(Rey_X**0.5)
            Nul= (0.037 * ((Rey**0.8) - A) * (Pr**(1/3)))
        h_pwt= Nul * self.k_ar/comprimento
        qconv = h_pwt*self.a_res*(self.to - self.tf)
        quantidade_de_calor= massa*condutividade*(temp_desejada-self.to)
        return qconv,quantidade_de_calor
    
    def arrefecimento_freio(self,omega_inicial,diametro_freio,massa_carro,aceleração,rho_f,volume_f,c_f,t_parada,t_res,n_frenagem):
        rey_rotacional= (omega_inicial*(diametro_freio**2))/self.ni
        if rey_rotacional > 10**6:      
            h_freio=0.4*(self.k_ar/diametro_freio)*(rey_rotacional**0.8)     #Coeficiente de convecção do freio para fluxo turbulento
        else:
            h_freio=0.7*(self.k_ar/diametro_freio)*(rey_rotacional**0.55)     #Coeficiente de convecção do freio para laminar
        calor=h_freio*self.a_res*(self.to-self.tf)                            #Calor dissipado por convecção
        P_bav=massa_carro*aceleração*(self.v/2)                     #Potencia media após uma frenagem
        deltaT= (P_bav*t_parada)/(rho_f*c_f*volume_f)                         #Calculo da diferença de temperatura após a frenagem

        # Calcular os denominadores comuns
        denominator_common = (h_freio*self.a_res*t_res)/(rho_f*c_f*volume_f)

        # Calcular o numerador e o denominador separadamente
        numerator = (1 - math.exp(-n_frenagem * denominator_common)) * deltaT
        denominator = 1 - math.exp(-denominator_common)
        # Calcular a diferença de temperatura
        temperature_change = (numerator / denominator) + self.tf               #Calculo da temperatura do freio após n frenagens e tendo resfriamento a ar 

        return temperature_change

#Parâmetros físicos do ar:
rho= 1.184                       #densidade(kg/m^3)
c_ar= 1007                       #Calor especifico (J/(kg*K))
mi= 1.849*10**-5                 #Viscosidade dinâmica (kg/m*s)
k_ar=0.02551                     #Condutividade termica do ar (W/m*K)
ni= (mi/rho)                     #Viscosidade Cinematica (m²/s)
temp_ar=25                       #temperatura do ar (°C)

#Parametros do carro
a_ent=0.004                      #Área de admissão do vento para arrefecimento (m^2)
v=8.88                           #Velocidade do carro (m/s)
massa_c=230                      #massa do carro (kg)
a=1                              #Aceleração horizontal do carro (m/s^2)
t_parada=4                       #Tempo que leva até o carro parar após uma frenagem (s)
t_res=90                         #Tempo de resfriamento entre as frenagem (s)
#Parametros de PWT
    #Motor
massa_motor=70                                      #Massa do motor (kg)
temp_motor=100                                      #Temperatura do motor (°C) 
comprimento_motor=0.6                               #Comprimento do motor (m)
largura_motor=0.3                                   #Largura do motor (m)
a_res_motor= largura_motor*comprimento_motor        #Área de resfriamento do motor (m^2)
c_motor=420                                         #Calor especifico (J/(kg*K))
    #Bateria
massa_bateria=50                                    #Massa do Pack de baterias (kg)
temp_bateria=100                                    #Temperatura do pack de baterias(°C) 
comprimento_bateria=0.5                             #Comprimento do pack de baterias (m)
largura_bateria=0.5                                 #largura do pack de baterias (m)
a_res_bateria= largura_bateria*comprimento_bateria  #Área de resfriamento do motor (m^2)
c_bateria=420                                       #Calor especifico (J/(kg*K))   

#Parametros para arrefecimento dos freios, considerando que o material do freio é o aço 1020
omega_inicial=50                            #Velocidade angular inicial das rodas (rad/s)
diametro_freio=0.17                         #Diametro do disco de freio (m)                      
a_freio= math.pi*(diametro_freio/2)**2      #Área do disco de freio (m^2)
temp_freio=100                              #Temperatura do freio (°C)                   
volume_freio=a_freio*0.008                  #área do disco de freio * expessura dele
rho_freio=7900                              #densidade(kg/m^3)
c_freio=420                                 #Calor especifico (J/(kg*K))
k_freio=169200                              #Condutividade termica do ar (W/m*K)

motor=Arrefecimento(a_res_motor,a_ent,rho,c_ar,temp_ar,temp_motor,v,ni,k_ar,mi)             #Definição do objeto para motor
bateria=Arrefecimento(a_res_bateria,a_ent,rho,c_ar,temp_ar,temp_bateria,v,ni,k_ar,mi)       #Definição do objeto para bateria
freio_1=Arrefecimento(a_freio,a_ent,rho,c_ar,temp_ar,temp_freio,v,ni,k_ar,mi)               #Definição do objeto para freio

#Arrefecimento para PWT
calor_abs_motor,calor_necessario_motor=motor.arrefecimento_pwt(comprimento_motor,massa_motor,c_motor,80)
print(f"O arrefecimento do motor absorbe {calor_abs_motor:0.2f} J/s")
print(f"O motor precisa perder {calor_necessario_motor:0.2f} J para descer até a temperatura ideal")

calor_abs_bateria,calor_necessario_bateria=bateria.arrefecimento_pwt(comprimento_bateria,massa_bateria,c_bateria,80)
print(f"O arrefecimento da bateria absorbe {calor_abs_bateria:0.2f} J/s")
print(f"A bateria precisa perder {calor_necessario_bateria:0.2f} J para descer até a temperatura ideal")
# Calcular a temperatura para diferentes números de frenagens
n_frenagens = np.arange(1, 100)
temperaturas = [freio_1.arrefecimento_freio(omega_inicial, diametro_freio, massa_c, a,rho_freio, volume_freio, c_freio, t_parada, t_res, n) for n in n_frenagens]

# Plotar o gráfico
plt.plot(n_frenagens, temperaturas)
plt.xlabel('Número de frenagens')
plt.ylabel('Temperatura (°C)')
plt.title('Temperatura em função do número de frenagens')
plt.grid(True)
plt.show()