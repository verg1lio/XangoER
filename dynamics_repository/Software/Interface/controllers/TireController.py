from Software.Interface.controllers.MainController import MainController
from Software.Systems.models.Dynamics import tire
from Software.Interface.controllers.envio_texto import envio
import numpy as np
from routes import route_map
class TireController(MainController):
    """
    Controlador para gerenciar modelos do sistema de Pneus.

    Herda funcionalidades básicas de MainController e adiciona métodos específicos
    para criação e manipulação de modelos do tipo tire.
    """

    def __init__(self):
        """
        Inicializa a instância do controlador, incluindo o objeto de envio de dados.
        """
        super().__init__()
        self.envio = envio()


    def instance_tire(self, dicionario, instanciado, falha):

        """
        Cria uma instância do modelo Tire com base nos dados do dicionário fornecido.

        Converte os valores para float, cria o objeto de modelo,
        armazena-o na coleção e chama callbacks para sucesso ou falha.

        Args:
            dicionario (dict): Dados para criar o modelo, com chaves específicas.
            instanciado (callable): Função a ser chamada quando o modelo for criado com sucesso.
            falha (callable): Função a ser chamada em caso de erro de conversão ou falha.

        Nota:
            Espera que o dicionário tenha valores que correspondam às chaves:
            ["tire_Fz", "tire_Ca", "B0", "B1", "B2", "B3", "WB", "rear_axle_length",
             "track_y", "tire_k"] na mesma ordem.
        """
        
        lista = list(dicionario.values())
        
        sucesso = 0

        try:

            for valor in range(len(lista)):
                lista[valor]= float(lista[valor]) 

            sucesso = 1

                
            
        except ValueError:
            
            falha()
            
            
            
        if sucesso and len(lista):
            
            tire_keys = ["tire_Fz", "tire_Ca", "tire_friction_coef", "B0", "B1", "B2", "B3", "WB", "rear_axle_length",
                    "track_y", "tire_k"]
            
            tire_model = dict(zip(tire_keys,lista))

            num = self.models.size("Tire") + 1

            tire_model = {"system": "Tire", "number": num, "model": tire_model}
            
            self.models.createModel(tire_model)
            
            instanciado()
            
    def run_data(self, tire_model):

        tire_model = tire(**tire_model)

        #armazenando tire_Ca
        tire_Ca = tire_model.tire_Ca

        # Parâmetros fixos da simulação
        result = [(0.3336564873588197), (1.6271741344929977), (10), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]

         # Criação de arrays com valores para plotagem e análise
        slip_ratio = np.linspace(-1, 1, 1000)  # Razão de escorregamento variando de -1 a 1
        slip_angles = np.linspace(-9, 9, 1000)  # Ângulos de deslizamento variando de -9 a 9 graus

        # Dados experimentais
        ratio = np.linspace(-1, 1, 19)
        
        angles = np.array([-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
        
        tire_lateral_forces_1 = np.array([-2300, -2200, -2060, -1880, -1680, -1450, -1190, -850, -430, 60, 520, 890, 1170, 1390, 1580, 1730, 1890, 2000, 2090])
        
        tire_auto_align_moment_1 = np.array([-28.84, -28.68, -27.21, -26.41, -27.70, -24.21, -24.15, -15.88, -4.91, 14.72, 33.80, 43.79, 46.93, 49.09, 50.90, 50.10, 50.81, 48.12, 48.83])
        
        camber_experimental = [-2060.0, -1950.0, -1840.0, -1700.0, -1540.0, -1350.0, -1130.0, -860.0, -480.0, -30.0, 460.0, 880.0, 1230.0, 1490.0, 1720.0, 1910.0, 2090.0, 2230.0, 2310.0]
        
        camber_experimental_1 = [-1940.0, -1860.0, -1750.0, -1610.0, -1450.0, -1260.0, -1050.0, -760.0, -400.0, 60.0, 540.0, 970.0, 1290.0, 1550.0, 1780.0, 1980.0, 2150.0, 2280.0, 2370.0]
        
        camber_experimental_15 = [-1840.0, -1750.0, -1670.0, -1520.0, -1370.0, -1180.0, -960.0, -680.0, -310.0, 130.0, 610.0, 1020.0, 1360.0, 1630.0, 1850.0, 2040.0, 2220.0, 2360.0, 2430.0]


        tire_model.tire_Sa = slip_angles
        tire_model.tire_Ls = slip_ratio
        # Forças e gráficos...
        predicted_lat, predicted_align, predicted_long = tire_model.Tire_forces(result)

        # Testando cambers diferentes
        camber_values = [0.5, 1.0, 1.5]

        #Lista para armazenar valores
        forcas_camber = []

        #Gerando valores diferentes de forças de camber
        for camber in camber_values:
            tire_model.tire_Ca = camber
            forcas_camber.append(tire_model.Tire_forces(result)[0])

        #Plotagem
        camber_64 = tire_model.plot_camber_base64(
            forcas_camber[0], forcas_camber[1], forcas_camber[2],
            camber_experimental, camber_experimental_1, camber_experimental_15, angles
        )

        # Reset camber e plot principal
        tire_model.tire_Ca = tire_Ca

        #De boa
        graph_64 = tire_model.plot_graph_base_64(
            predicted_lat, predicted_align, predicted_long,
            tire_lateral_experimental=tire_lateral_forces_1,
            tire_auto_align_experimental=tire_auto_align_moment_1,
            angles=angles, ratio=ratio
        )


        #Tranquilo
        # Plot deformação e mecanismo
        deformacao_64 = tire_model.plotar_deformacao_64(track_y_values=np.linspace(0, 25, 1000))

        #Só depende de parâmetros próprios
        mechanism_64 = tire_model.plot_mechanism_base64()

        graphics = {"camber": camber_64, "graph": graph_64, "deformacao": deformacao_64, "mechanism": mechanism_64}

        return  graphics
    
    def run_deformacao(self, graph, func, page):

        route_map.update({"/run_tire/deformacao": lambda: func(graph, page)})

        page.go("/run_tire/deformacao")

    def run_graph(self, graph, func, page):

        route_map.update({"/run_tire/graph": lambda: func(graph, page)})

        page.go("/run_tire/graph")

    def run_mechanism(self, img_64, func, page):

        route_map.update({"/run_tire/run_mechanism": lambda: func(img_64, page)})

        page.go("/run_tire/run_mechanism")