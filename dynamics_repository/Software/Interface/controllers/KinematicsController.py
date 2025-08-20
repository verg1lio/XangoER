from Software.Interface.controllers.MainController import MainController
from Software.Systems.models.Dynamics import kmt
from Software.Interface.controllers.envio_texto import envio
import numpy as np
from routes import route_map


class KinematicsController(MainController):

    """
    Controlador para gerenciar modelos do sistema de Cinemática.

    Herda funcionalidades básicas de MainController e adiciona métodos específicos
    para criação e manipulação de modelos do tipo kmt.
    """

    def __init__(self):

        """
        Inicializa a instância do controlador, incluindo o objeto de envio de dados.
        """

        super().__init__()
        self.envio = envio()


    def instance_kmt(self, dicionario, instanciado, falha):

        """
        Cria uma instância do modelo Kinematics com base nos dados do dicionário fornecido.

        Converte os valores para os tipos esperados, cria o objeto de modelo,
        armazena-o na coleção e chama callbacks para sucesso ou falha.

        Args:
            dicionario (dict): Dados para criar o modelo, com chaves específicas.
            instanciado (callable): Função a ser chamada quando o modelo for criado com sucesso.
            falha (callable): Função a ser chamada em caso de erro de conversão ou falha.

        Nota:
            Espera as seguintes chaves no dicionario:
            "damper_type", "damper_F_static", "damper_K_friction", "L0", "L1", "L2", "L3", 
            "spring_type", "spring_k", "spring_x", "spring_non_lin_coef".
        """

        kmt_keys = ["damper_type", "damper_F_static", "damper_K_friction", "damper_F_viscous", "L0", "L1", "L2", "L3", 
                    "spring_type", "spring_k", "spring_x", "spring_non_lin_coef"]
        
        kmt_model = dict(zip(kmt_keys,dicionario.values()))
        
        sucesso = 0

        try:

            
            kmt_model["damper_type"] = str(kmt_model["damper_type"])
            kmt_model["damper_F_static"] = int(kmt_model["damper_F_static"])
            kmt_model["damper_K_friction"] = int(kmt_model["damper_K_friction"])
            kmt_model["spring_type"] = str(kmt_model["spring_type"])
            kmt_model["spring_k"] = int(kmt_model["spring_k"])
            kmt_model["spring_x"] = int(kmt_model["spring_x"])
            kmt_model["spring_non_lin_coef"] = int(kmt_model["spring_x"])
            kmt_model["L0"] = int(kmt_model["L0"])
            kmt_model["L1"] = int(kmt_model["L1"])
            kmt_model["L2"] = int(kmt_model["L2"])
            kmt_model["L3"] = int(kmt_model["L3"])

            sucesso = 1

                
            
        except ValueError:
            
            falha()
            
            
            
        if sucesso and len(kmt_model):
            
            """kmt_model = kmt(damper_F_static=kmt_model["damper_F_static"],
                            damper_K_friction=kmt_model["damper_K_friction"],
                            damper_type=kmt_model["damper_type"],
                            L0=kmt_model["L0"],
                            L1=kmt_model["L1"],
                            L2=kmt_model["L2"],
                            L3=kmt_model["L3"],
                            spring_k=kmt_model["spring_k"],
                            spring_type=kmt_model["spring_type"],
                            spring_non_lin_coef=kmt_model["spring_non_lin_coef"],
                            spring_x=kmt_model["spring_x"])"""
            
            
            num = self.models.size("Kinematics") + 1
        
            kmt_model = {"system": "Kinematics", "number": num, "model": kmt_model}
            
            self.models.createModel(kmt_model)
            
            instanciado()
            

    def run_kmt(self, model):

        model = kmt(**model)

        damper = model.plot_damper(damper_Vx_values=np.linspace(-10, 10, 100), damper_Vy_values=np.linspace(-10, 10, 100))

        spring = model.plot_spring(spring_x_values=np.linspace(-10, 10, 100) , spring_y_values=np.linspace(-10, 10, 100))

        cinematica, text = model.plotar_cinematica_base64()

        graph = [damper, spring, cinematica]

        return graph, text
    
    def routes(self):

        print(route_map)
    
    def run_graphs(self, graphics, run_plotagem, page):


        route_map.update({"/run_kinematics/graph": lambda: run_plotagem(graphics,page)})
        
        page.go("/run_kinematics/graph")    
    
    def dados_cinematica(self, graphics, text, func, page):

        route_map.update({"/run_kinematics/dados_cinematica": lambda: func(graphics, text, page)})
        
        page.go("/run_kinematics/dados_cinematica")    