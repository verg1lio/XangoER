import numpy as np
import matplotlib.pyplot as plt
import math

class Arrefecimento():
    def __init__(self,area_resfriamento,area_entrada,densidade_ar,coeficiente_transferencia,calor_especifico_ar,temp_ar,temp_objeto,velocidade_ar,viscosidade_ar,condutividade_ar):
        self.a_res= area_resfriamento
        self.a_ent= area_entrada
        self.rho= densidade_ar
        self.h= coeficiente_transferencia
        self.c= calor_especifico_ar
        self.tf=temp_ar
        self.to=temp_objeto
        self.ni=viscosidade_ar
        self.k_ar=condutividade_ar
        self.Vazao= velocidade_ar*area_entrada*densidade_ar

    def arrefecimento_freio(self,omega_inicial,diametro_freio,massa_carro,aceleração,velocidade_carro,rho_f,volume_f,c_f,t_parada,t_res,n_frenagem):
        rey_rotacional= (omega_inicial*(diametro_freio**2))/self.ni
        if rey_rotacional > 10**6:      
            h_freio=0.4*(self.k_ar/diametro_freio)*(rey_rotacional**0.8)     #Coeficiente de convecção do freio para fluxo turbulento
        else:
            h_freio=0.7*(self.k_ar/diametro_freio)*(rey_rotacional**0.55)     #Coeficiente de convecção do freio para laminar
        calor=h_freio*self.a_res*(self.to-self.tf)                            #Calor dissipado por convecção
        P_bav=massa_carro*aceleração*(velocidade_carro/2)                     # Potencia media após uma frenagem
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
#parâmetros para a troca termica
temp_ar=30                       #temperatura do ar (°C)
a_ent=0.004                      #Área de admissão do vento para arrefecimento (m^2)
a_res = 1                        #Área de resfriamento da bateria (m^2)
h= 300                           #Coeficiente de transferência de calor (W/m2*K)
#Parametros do carro
v=8.88
massa_c=230
a=1
t_parada=4.5
t_res=90

#Parametros para arrefecimento dos freios, considerando que o material do freio é o aço 1020
omega_inicial=30
diametro_freio=0.15
a_freio= math.pi*(diametro_freio/2)**2
temp_freio=100
volume_freio=a_freio*0.004
rho_freio=7900                    #densidade(kg/m^3)
c_freio=420                       #Calor especifico (J/(kg*K))
k_freio=169200                    #Condutividade termica do ar (W/m*K)

freio_1=Arrefecimento(a_freio,a_ent,rho,h,c_ar,temp_ar,temp_freio,v,ni,k_ar)
# Calcular a temperatura para diferentes números de frenagens
n_frenagens = np.arange(1, 100)
temperaturas = [freio_1.arrefecimento_freio(omega_inicial, diametro_freio, massa_c, a, v, rho_freio, volume_freio, c_freio, t_parada, t_res, n) for n in n_frenagens]

# Plotar o gráfico
plt.plot(n_frenagens, temperaturas)
plt.xlabel('Número de frenagens')
plt.ylabel('Temperatura (°C)')
plt.title('Temperatura em função do número de frenagens')
plt.grid(True)
plt.show()