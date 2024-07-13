import numpy as np
import matplotlib.pyplot as plt
import math

class Arrefecimento():
    def __init__(self,area_resfriamento,densidade_ar,calor_especifico_ar,temp_ar,temp_objeto,velocidade,viscosidade_cinematica_ar,condutividade_ar,viscosidade_dinamica_ar,massa,temp_desejada,calor_especifico_objeto):
        self.a_res= area_resfriamento
        self.rho= densidade_ar
        self.c= calor_especifico_ar
        self.v=velocidade
        self.tf=temp_ar
        self.to=temp_objeto
        self.mi=viscosidade_dinamica_ar
        self.ni=viscosidade_cinematica_ar
        self.k_ar=condutividade_ar
        self.m=massa
        self.temp_desejada= temp_desejada
        self.c_objeto=calor_especifico_objeto
    def arrefecimento_pwt(self,comprimento):
        Pr=(self.mi*self.c)/self.k_ar                                               # Cálculo do número de Prandlt
        Rey=comprimento*self.v/self.ni                                              # Cálculo do número de Reynolds
    
        if Rey < 200000:                                                            #Fluxo laminar
            Nul= (0.664 * (Rey**0.5) * (Pr**(1/3)))
        else:                                                                       #Fluxo Turbulento
            X=200000*self.ni/self.v                                                 # Comprimento da placa onde o fluxo transiciona para turbulento
            Rey_X= X*self.v/self.ni                                                 # Cálculo do número de Reynolds para a distância X
            A= 0.037*(Rey_X**0.8) - 0.664*(Rey_X**0.5)
            Nul= (0.037 * ((Rey**0.8) - A) * (Pr**(1/3)))                           # Cálculo do número de Nusselt
        h_pwt= Nul * self.k_ar/comprimento                                          # Cálculo do coeficiente de convecção
        flag = 0
        flag_tempo = 0
        temp_atual = self.to
        quantidade_de_calor = self.m * self.c_objeto * (self.temp_desejada - temp_atual)  # Cálculo da quantidade de calor que o objeto precisa perder para alcançar a temperatura ideal
        temp_grafico = []
        tempo_grafico = []
        while flag == 0:
            qconv = h_pwt * self.a_res * (temp_atual - self.tf)  # Calor absorvido por convecção
            temp_final = (-qconv / (self.m * self.c_objeto)) + temp_atual
            temp_atual = temp_final
            temp_grafico.append(temp_atual)
            tempo_grafico.append(flag_tempo)
            flag_tempo += 1
            if temp_final <= self.temp_desejada:
                flag = 1
        return quantidade_de_calor, flag_tempo, temp_grafico, tempo_grafico
    
    def arrefecimento_freio(self,omega_inicial,diametro_freio,massa_carro,aceleração,rho_f,volume_f,t_parada,t_res,n_frenagem):
        rey_rotacional= (omega_inicial*(diametro_freio**2))/self.ni
        if rey_rotacional > 10**6:      
            h_freio=0.4*(self.k_ar/diametro_freio)*(rey_rotacional**0.8)                # Coeficiente de convecção do freio para fluxo turbulento
        else:
            h_freio=0.7*(self.k_ar/diametro_freio)*(rey_rotacional**0.55)               # Coeficiente de convecção do freio para laminar
        P_bav=massa_carro*aceleração*(self.v/2)                                         # Potência média após uma frenagem
        deltaT= (P_bav*t_parada)/(rho_f*self.c_objeto*volume_f)                         # Cálculo da diferença de temperatura após a frenagem
        # Calcular os denominadores comuns
        denominator_common = (h_freio*self.a_res*t_res)/(rho_f*self.c_objeto*volume_f)
        Temp_frenagem=[]
        for n in range (n_frenagem):
            # Calcular o numerador e o denominador separadamente
            numerator = (1 - math.exp(-n * denominator_common)) * deltaT
            denominator = 1 - math.exp(-denominator_common)
            # Calcular a diferença de temperatura
            temperature_change = (numerator / denominator) + self.tf               # Cálculo da temperatura do freio após n frenagens e tendo resfriamento a ar 
            Temp_frenagem.append(temperature_change)
        flag2 = 0
        flag_tempo2 = 0
        temp_atual_f = Temp_frenagem.pop()
        quantidade_de_calor_f = self.m * self.c_objeto * (self.temp_desejada - temp_atual_f)  # Cálculo da quantidade de calor que o objeto precisa perder para alcançar a temperatura ideal
        temp_grafico_f = []
        tempo_grafico_f = []
        while flag2 == 0:
            qconv_f = h_freio * self.a_res * (temp_atual_f - self.tf)  # Calor absorvido por convecção
            temp_final_f = (-qconv_f / (self.m * self.c_objeto)) + temp_atual_f
            temp_atual_f = temp_final_f
            temp_grafico_f.append(temp_atual_f)
            tempo_grafico_f.append(flag_tempo2)
            flag_tempo2 += 1
            if temp_final_f <= self.temp_desejada:
                flag2 = 1       
        return Temp_frenagem,quantidade_de_calor_f,flag_tempo2,temp_grafico_f,tempo_grafico_f
# Parâmetros físicos do ar:
rho = 1.184                        # Densidade(kg/m^3)
c_ar = 1007                        # Calor especifico (J/(kg*K))
mi = 1.849*10**-5                  # Viscosidade dinâmica (kg/m*s)
k_ar = 0.02551                     # Condutividade termica do ar (W/m*K)
ni = (mi/rho)                      # Viscosidade Cinematica (m²/s)
temp_ar = 25                       # Temperatura do ar (°C)

# Parametros do carro
a_ent = 0.004                      # Área de admissão do vento para arrefecimento (m^2)
v = 10                             # Velocidade do carro (m/s)
massa_c = 230                      # Massa do carro (kg)
a = 1                              # Aceleração horizontal do carro (m/s^2)
t_parada = 4                       # Tempo que leva até o carro parar após uma frenagem (s)
t_res = 90                         # Tempo de resfriamento entre as frenagem (s)

# Parametros de PWT
    # Motor
massa_motor = 70                                      # Massa do motor (kg)
temp_motor = 100                                      # Temperatura do motor (°C) 
comprimento_motor = 0.6                               # Comprimento do motor (m)
largura_motor = 0.3                                   # Largura do motor (m)
a_res_motor = largura_motor*comprimento_motor         # Área de resfriamento do motor (m^2)
c_motor = 420                                         # Calor especifico (J/(kg*K))
    # Bateria
massa_bateria = 50                                    # Massa do Pack de baterias (kg)
temp_bateria = 100                                    # Temperatura do pack de baterias(°C) 
comprimento_bateria = 0.5                             # Comprimento do pack de baterias (m)
largura_bateria = 0.5                                 # largura do pack de baterias (m)
a_res_bateria = largura_bateria*comprimento_bateria   # Área de resfriamento do motor (m^2)
c_bateria = 420                                       # Calor especifico (J/(kg*K))   
temp_ideal=40
# Parametros para arrefecimento dos freios, considerando que o material do freio é o aço 1020
omega_inicial = 50                            # Velocidade angular inicial das rodas (rad/s)
diametro_freio = 0.17                         # Diametro do disco de freio (m)                      
a_freio = math.pi*(diametro_freio/2)**2       # Área do disco de freio (m^2)                   
volume_freio = a_freio*0.008                  # Área do disco de freio * expessura dele
rho_freio = 7900                              # Densidade(kg/m^3)
c_freio = 420                                 # Calor especifico (J/(kg*K))
k_freio = 169200                              # Condutividade termica do freio (W/m*K)
massa_freio=volume_freio*rho_freio
temp_ideal_f=60
n_frenagens = 100
motor   = Arrefecimento(a_res_motor,rho,c_ar,temp_ar,temp_motor,v,ni,k_ar,mi,massa_motor,temp_ideal,c_motor)               # Definição do objeto para motor
bateria = Arrefecimento(a_res_bateria,rho,c_ar,temp_ar,temp_bateria,v,ni,k_ar,mi,massa_bateria,temp_ideal,c_bateria)       # Definição do objeto para bateria
freio_1 = Arrefecimento(a_freio,rho,c_ar,temp_ar,1,v,ni,k_ar,mi,massa_freio,temp_ideal_f,c_freio)                 # Definição do objeto para freio

# Arrefecimento para PWT
calor_necessario_motor, tempo_resfriamento_motor, temperatura_motor, tempo_motor = motor.arrefecimento_pwt(comprimento_motor)
print(f"O motor precisa perder {calor_necessario_motor:0.2f} J para esfriar até a temperatura ideal")
print(f"O motor leva {(tempo_resfriamento_motor)/60:0.2f} minutos esfriar até a temperatura ideal\n")

calor_necessario_bateria, tempo_resfriamento_bateria, temperatura_bateria, tempo_bateria=bateria.arrefecimento_pwt(comprimento_bateria)
print(f"A bateria precisa perder {calor_necessario_bateria:0.2f} J para esfriar até a temperatura ideal")
print(f"O motor leva {(tempo_resfriamento_bateria/60):0.2f} minutos esfriar até a temperatura ideal\n")
# Calcular a temperatura para diferentes números de frenagens
temperaturas,calor_necessario_freio,tempo_resfriamento_freio,temperatura_freio,tempo_freio= freio_1.arrefecimento_freio(omega_inicial, diametro_freio, massa_c, a,rho_freio, volume_freio, t_parada, t_res, n_frenagens)
print(f"O freio precisa perder {calor_necessario_freio:0.2f} J para esfriar até a temperatura ideal")
print(f"O freio leva {(tempo_resfriamento_freio)/60:0.2f} minutos esfriar até a temperatura ideal\n")

# Plotar o gráfico de freio
plt.figure(figsize=(12, 6))
plt.subplot(1,2,1)
plt.plot(temperaturas)
plt.xlabel('Número de frenagens')
plt.ylabel('Temperatura (°C)')
plt.title('Temperatura do freio em função do número de frenagens')
plt.grid(True)
# Plotar o gráfico do resfriamento do freio
plt.subplot(1,2,2)
plt.plot(tempo_freio, temperatura_freio)
plt.xlabel('Tempo(s)')
plt.ylabel('Temperatura (°C)')
plt.title('Temperatura do freio em função do Tempo')
plt.grid(True)
plt.show()

# Plotar o gráfico das Temperaturas de PWT
plt.figure(figsize=(12, 6))
plt.subplot(1,2,1)
plt.plot(tempo_motor, temperatura_motor)
plt.title('Gráfico de Temperatura do motor em função do tempo')
plt.xlabel('Tempo(s)')
plt.ylabel('Temperatura (°C)')
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(tempo_bateria, temperatura_bateria)
plt.title('Gráfico de Temperatura da bateria em função do tempo')
plt.xlabel('Tempo(s)')
plt.ylabel('Temperatura (°C)')
plt.grid(True)
plt.show()
