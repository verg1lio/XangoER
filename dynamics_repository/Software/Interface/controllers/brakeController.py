from Software.Systems.models.Dynamics import brake
from Software.Interface.controllers.envio_texto import envio
from Software.Interface.controllers.MainController import MainController
import numpy as np
from routes import route_map

class BrakeController(MainController):

    """
    Controlador responsável pela gestão dos modelos de freio (Brake)
    e pela navegação entre etapas da interface de execução do cálculo.

    Herdado de:
        MainController: Fornece funcionalidades comuns de gerenciamento de modelos.

    Atributos:
        envio (envio): Instância do controlador de envio de mensagens/textos.
    """

    def __init__(self):

        """
        Inicializa o controlador de freio.

        - Chama o construtor da classe base `MainController`.
        - Cria uma instância do controlador `envio` para mensagens de interface.
        """

        super().__init__()
        self.envio = envio()


    def instance_brake(self, dicionario, instanciado, falha):

        """
        Cria e armazena uma nova instância do modelo de freio (Brake).

        Fluxo:
        1. Recebe um dicionário com os parâmetros de entrada.
        2. Converte os valores para float (se possível).
        3. Popula os parâmetros do modelo com esses valores.
        4. Instancia a classe `brake` com os parâmetros.
        5. Salva o modelo no banco de dados (Excel) através de `self.model.createModel`.
        6. Chama a função `instanciado()` em caso de sucesso ou `falha()` em caso de erro.

        Args:
            dicionario (dict): Chaves correspondem aos nomes dos parâmetros, valores como strings numéricas.
            instanciado (callable): Função callback a ser executada se a criação for bem-sucedida.
            falha (callable): Função callback a ser executada se houver erro de conversão ou criação.

        Observações:
            - Se qualquer valor não puder ser convertido para float, a função `falha()` será chamada.
            - O número (`number`) do modelo é definido automaticamente com base na quantidade atual de modelos Brake.
        """

        lista = list(dicionario.values())

        params = [
            'RedP',  # Redução do pedal
            'a',  # Desaceleração [g]
            'psi',  # Distribuição de carga estática por eixo [adm]
            'μl',  # Coeficiente de atrito do contato pastilha/disco
            'HCG',  # Altura do centro de gravidade [m]
            'μ',  # Coeficiente de atrito do pneu/solo
            'FzF',  # Força de reação estática na dianteira [N]
            'FzR',  # Força de reação estática na traseira [N]
            'Rdp',  # Raio do pneu [m]
            'Dcm',  # Diâmetro do cilindro mestre em metros
            'Dwc',  # Diâmetro do cilindro da roda em metros
            'Npast',  # Número de pastilhas por disco
            'atrito_coeficiente',  # Coeficiente de atrito da pastilha
            'red',  # Raio efetivo do disco [m]
            'Mt',  # Massa total do veículo [kg]
            'm_wheel',  # Massa da roda [kg] <-- estimativa
            'm_tire',  # Massa do pneu [kg] <-- estimativa
            'L',  # Distância entre eixos [m]]
            'c_rr' # Coeficiente de resistência ao rolamento
        ]
            
        
        sucesso = 0

        try:
            

            for valor in range(len(lista)):
                lista[valor]= float(lista[valor])
                 

            sucesso = 1

                
            
        except ValueError:
            
            falha()
            
            
            
        if sucesso and len(lista):

            """for chave, valor in zip(params,lista):
                params[chave] = valor"""
                
            brake_model = dict(zip(params, lista))
            
            print(f"Chaves: {brake_model.keys()}\nValores: {brake_model.values()}")
            
            num = self.models.size("Brake") + 1
        
            brake_model = {"system": "Brake", "number": num, "model": brake_model}
            
            self.models.createModel(brake_model)
            
            instanciado()
   
    def result_data(self, dicionario, brake_model, func, page):

        """
        Prepara a navegação para a exibição dos resultados do cálculo de freio.

        - Lê os valores de força no pedal e velocidade do veículo do dicionário.
        - Registra no `route_map` a função que exibirá os resultados.
        - Navega para a rota `/run_brake/result_data`.

        Args:
            dicionario (dict): Contém as entradas 'Força do pedal' e 'Velocidade do veículo'.
            brake_model (object): Instância do modelo de freio a ser processada.
            func (callable): Função responsável por renderizar a página de resultados.
            page (ft.Page): Página Flet atual para manipulação da navegação.
        """

        brake_instance = brake(brake_model)
        
        pedal_forces = float(dicionario['Força do pedal'])

        vehicle_speed = float(dicionario['Velocidade do veículo'])

        tempo = float(dicionario['Tempo'])

        brake_instance.tempo = tempo
        brake_instance.initial_speed = vehicle_speed

        route_map.update({"/run_brake/result_data": lambda: func(brake_instance, pedal_forces, page)})

        page.go("/run_brake/result_data")
    
    def graph(self, brake_model, pedal_forces, func, page):
    
        """
        Prepara a navegação para a exibição do gráfico de resultados do freio.

        - Registra no `route_map` a função responsável por gerar o gráfico.
        - Navega para a rota `/run_brake/graph`.

        Args:
            brake_model (object): Instância do modelo de freio.
            pedal_forces (float): Força no pedal usada no cálculo.
            func (callable): Função responsável por renderizar o gráfico.
            page (ft.Page): Página Flet atual para manipulação da navegação.
        """

        route_map.update({"/run_brake/graph": lambda: func(brake_model, pedal_forces, page)})

        page.go("/run_brake/graph")

    def back_run(self, model, page, func):

        route_map.update({"/run_exec": lambda: func(model,page)})
        
        page.go("/{route}")