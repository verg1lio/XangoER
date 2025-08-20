from Software.Systems.models.Dynamics import dt, tire
from Software.Interface.controllers.envio_texto import envio
from routes import route_map
import numpy as np
from Software.Interface.controllers.MainController import MainController

class DrivetrainController(MainController):

    """
    Controlador responsável pela gestão dos modelos de transmissão (Drivetrain)
    e pela navegação entre etapas da interface de execução do cálculo.

    Herdado de:
        MainController: Fornece funcionalidades comuns de gerenciamento de modelos.

    Atributos:
        envio (envio): Instância do controlador de envio de mensagens/textos.
    """

    def __init__(self):

        """
        Inicializa o controlador, herdando funcionalidades do MainController
        e instanciando o módulo de envio de texto.
        """

        super().__init__()
        self.envio = envio()


    def instance_dt(self, dicionario, instanciado, falha):

        """
        Cria e registra um novo modelo de trem de força a partir dos parâmetros fornecidos.

        Args:
            dicionario (dict): Dicionário contendo os parâmetros do modelo, 
                onde cada valor deve ser conversível para float.
            instanciado (callable): Função callback chamada quando o modelo é criado com sucesso.
            falha (callable): Função callback chamada quando ocorre erro de conversão de dados.
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

            dt_keys = [
                "gravitational_center_x", "gravitational_center_y", "mass", "axis_distance",
                "friction_coef", "tire_radius", "ideal_acceleration", "primary_reduction",
                "reducao_unica", "rpm", "torque", "cp"]
            
            dt_model = dict(zip(dt_keys,lista))

            num = self.models.size("Drivetrain") + 1
        
            dt_model = {"system": "Drivetrain", "number": num, "model": dt_model}
            
            self.models.createModel(dt_model)
            
            instanciado()
            
    def dados(self, dados):

        """
        Organiza os dados de simulação da transmissão em strings formatadas.

        Args:
            dados (list[dict]): Lista de dicionários contendo dados de simulação com chaves:
                'forca_trativa', 'va', 'velocidade_linear', 'fa', 'rr', 'forca_final'.

        Returns:
            tuple[str]: Tupla contendo strings concatenadas (uma por linha) para:
                força trativa, velocidade angular, velocidade linear,
                força de arrasto, resistência ao rolamento e força final.
        """

        forca_trativa = ""
        velocidade_angular = ""
        velocidade_linear = ""
        forca_arrasto = ""
        resistencia_rolamento = ""
        forca_final = ""

        for dado in dados:
            forca_trativa += f"{dado["forca_trativa"]}\n"
            velocidade_angular += f"{dado["va"]}\n" 
            velocidade_linear += f"{dado["velocidade_linear"]}\n"
            forca_arrasto += f"{dado["fa"]}\n"
            resistencia_rolamento += f"{dado["rr"]}\n"
            forca_final += f"{dado["forca_final"]}\n"

        return forca_trativa, velocidade_angular, velocidade_linear, forca_arrasto, resistencia_rolamento, forca_final
    

    def info_func(self, dt_model, dicionario, func, page):

        """
        Atualiza o modelo com o RPM final e redireciona para a rota de exibição dos resultados.

        Args:
            dt_model: Instância do modelo de trem de força.
            dicionario (dict): Contém o valor de 'RPM Final'.
            func (callable): Função responsável por renderizar a página.
            page: Objeto de página do Flet para navegação.
        """
        
        dt_model = dt(**dt_model)

        dt_model.new_rpm = int(dicionario["RPM Final"])

        route_map.update({"/run_dt/result_data": lambda: func(dt_model, page)})

        page.go("/run_dt/result_data")


    def performance_run(self, dt_model, func, page):

        """
        Redireciona para a rota de execução de performance do modelo.

        Args:
            dt_model: Instância do modelo de trem de força.
            func (callable): Função responsável por executar e renderizar a performance.
            page: Objeto de página do Flet.
        """

        route_map.update({"/run_dt/performance_run": lambda: func(dt_model, page)})

        page.go("/run_dt/performance_run")

    def slip_ratio(self, dt_model, func, velocidade_angular, rpm_faixa, page):

        """
        Redireciona para a rota de cálculo do slip ratio (escorregamento do pneu).

        Args:
            dt_model: Instância do modelo de trem de força.
            func (callable): Função que processa o slip ratio.
            velocidade_angular: Valores de velocidade angular.
            rpm_faixa: Faixa de valores de RPM.
            page: Objeto de página do Flet.
        """

        route_map.update({"/run_dt/slip_ratio": lambda: func(dt_model, velocidade_angular, rpm_faixa, page)})

        page.go("/run_dt/slip_ratio")

    def calcular_slip(self, dados,velocidade_linear,raio_pneu, rpm_faixa):

        """
        Calcula e plota o slip ratio a partir de dados experimentais e parâmetros do pneu.

        Args:
            dados (list[dict]): Lista de dicionários contendo valores de 'va' (velocidade angular).
            velocidade_linear (float): Velocidade linear do veículo.
            raio_pneu (float): Raio do pneu.
            rpm_faixa (list[float]): Lista de valores de RPM para plotagem.

        Returns:
            list[str]: Lista de imagens em Base64 com os gráficos gerados.
        """

        velocidade_angular = np.array([dado["va"] for dado in dados])

        slip_ratio = tire.slip_ratio_1(velocidade_angular=velocidade_angular,velocidade_linear=velocidade_linear,raio_pneu=raio_pneu)

        #Transformando num array  
        rpm_faixa = np.array(rpm_faixa)

        #Plotagem de gráfico do slip ratio e saídas de seus dados no terminal
        tire.show_slip_ratio(rpm_faixa, slip_ratio, velocidade_angular)

        #Salvando os dados como array para cálculo de força longitudinal
        slip_ratios = np.array(slip_ratio)

        #Criando instância da classe Tire
        Slip_model = tire(tire_Fz=1500, tire_Sa=0, tire_Ls=slip_ratios, tire_friction_coef=1.45, tire_Ca=0)

        #Dados experimentais para instância em Tire
        result = [(0.3336564873588197), (1.6271741344929977), (1), (4.3961693695846655), (931.4055775279057), (366.4936818126405)]

        #Recebendo valores de força lateral, torque auto allinhante e força longitudinal
        tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces = Slip_model.Tire_forces(result)

        #Plotagem de gráficos
        graph = Slip_model.plot_graph_base_64( tire_lateral_forces, tire_auto_align_moment, tire_longitudinal_forces)

        return graph

    def halfshafts(self, dt_model, func, page):

        """
        Redireciona para a rota de simulação dos semi-eixos (halfshafts).

        Args:
            dt_model: Instância do modelo de trem de força.
            func (callable): Função que executa e renderiza a simulação.
            page: Objeto de página do Flet.
        """
        
        route_map.update({"/run_dt/halfshafts": lambda: func(dt_model, page)})

        page.go("/run_dt/halfshafts")