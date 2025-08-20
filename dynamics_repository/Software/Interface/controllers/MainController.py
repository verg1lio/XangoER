from Software.Systems.management.Models import Models
from routes import route_map
import pandas as pd
    

class MainController:
    """
    Controlador principal para gerenciar a comunicação entre a interface e os modelos
    do sistema, incluindo navegação de rotas, manipulação de DataFrames e execução
    de modelos.
    """

    def __init__(self):
        """Inicializa o controlador principal e a instância de Models."""
        self.models = Models()

    def models_saved(self, system):
        """
        Obtém o número de modelos salvos para um sistema específico e uma lista de identificadores.

        Args:
            system (str): Nome do sistema (ex.: "Freio", "Transmissão").

        Returns:
            tuple: (quantidade de modelos salvos, lista com números sequenciais dos modelos)
        """
        self.models.atualizar()
        size = self.models.size(system)
        lista = [n + 1 for n in range(size)]
        return size, lista
    
    def run(self, system, num, route, func, page):
        """
        Executa um modelo salvo e redireciona para a rota correspondente.

        Args:
            system (str): Nome do sistema.
            num (int): Número do modelo salvo.
            route (str): Nome da rota de destino.
            func (callable): Função que será executada com o modelo.
            page (ft.Page): Página do Flet para navegação.
        """
        
        system, model = self.models.getModel(system=system, number=num)
        
        self.change_route(model=model, route=route, func=func, page=page)

    def change_route(self, model, route, func, page):
        """
        Atualiza o mapeamento de rotas e redireciona para a rota especificada.

        Args:
            model (object): Modelo que será passado para a função.
            route (str): Nome da rota.
            func (callable): Função que receberá o modelo e a página.
            page (ft.Page): Página do Flet para navegação.
        """
        route_map.update({f"/{route}": lambda: func(model, page)})
        page.go(f"/{route}")

    def get_data_set(self):
        """
        Obtém o conjunto de dados atualizado dos modelos.

        Returns:
            tuple: (DataFrame de transmissão, freio, pneu e cinemática)
        """
        dt, brake, tire, kmt = self.models.atualizar_lista()
        return dt, brake, tire, kmt 
    
    def get_all_dataframe_names(self):
        """
        Retorna a lista de nomes de DataFrames disponíveis.

        Returns:
            list: Lista com nomes como ["Freio", "Pneu", "Cinemática", "Transmissão"].
        """
        return ["Freio", "Pneu", "Cinemática", "Transmissão"]

    def get_dataframe(self, nome_df: str) -> pd.DataFrame:
        """
        Obtém um DataFrame específico pelo seu nome.

        Args:
            nome_df (str): Nome do DataFrame desejado.

        Returns:
            pd.DataFrame: DataFrame correspondente ou vazio se não encontrado.
        """
        dt, brake, tire, kmt = self.models.atualizar_lista()

        map_nome_para_df = {
            "Freio": brake,
            "Pneu": tire,
            "Cinemática": kmt,
            "Transmissão": dt,
        }

        df = map_nome_para_df.get(nome_df, pd.DataFrame())

        """if nome_df == "Freio" and not df.empty and df.shape[1] > 0:
            df = df.drop(df.columns[0], axis=1)"""

        return df

    def atualizar_celula(self, nome_df: str, linha: int, coluna: str, valor):
        """
        Atualiza o valor de uma célula específica em um DataFrame.

        Args:
            nome_df (str): Nome do DataFrame.
            linha (int): Índice da linha.
            coluna (str): Nome da coluna.
            valor: Novo valor a ser definido.
        """
        system = self.models.convertGUI(nome_df)
        self.models.atualizar_celula(system, linha, coluna, valor)

    def excluir_linha(self, nome_df: str, linha: int):
        """
        Remove uma linha de um DataFrame.

        Args:
            nome_df (str): Nome do DataFrame.
            linha (int): Índice da linha a ser removida.
        """
        system = self.models.convertGUI(nome_df)
        self.models.excluir_linha(system, linha)
