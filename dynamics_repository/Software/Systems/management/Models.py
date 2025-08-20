import pandas as pd

class Models:
    """
    Classe Models

    Responsável por carregar, armazenar, criar, listar, buscar, atualizar e excluir modelos
    dos sistemas de um carro (Drivetrain, Brake, Tire, Kinematics), utilizando
    dados salvos em um arquivo Excel.

    Estrutura geral:
    - Cada modelo é armazenado como um dicionário no formato:
        {
            "system": <str>,   # Nome do sistema
            "number": <int>,   # Número sequencial do modelo
            "model": <obj>     # Instância da classe correspondente
        }
    - Atributo principal:
        models (list): Lista contendo todos os modelos carregados.

    Métodos principais:
        __init__():
            Inicializa a classe e carrega os modelos do arquivo Excel.

        createModel(model):
            Adiciona um novo modelo à lista e salva no Excel.

        size(system):
            Retorna a quantidade de modelos salvos para o sistema especificado.

        convert(opt):
            Converte entrada numérica (string) do terminal para nome do sistema.

        convertGUI(opt):
            Converte rótulos usados na interface gráfica para nomes internos do sistema.

        getModel(system, number):
            Busca um modelo específico pelo nome do sistema e número do modelo.

        getDataFrames():
            Lista todos os modelos carregados e retorna DataFrames separados
            para cada sistema.

        saveModels():
            Salva todos os modelos no arquivo Excel, em abas separadas.

        getDicts(system, model):
            Converte dicionários carregados do Excel em instâncias de classes
            específicas de cada sistema.

        bancoAcesso():
            Carrega todos os modelos a partir do arquivo Excel, convertendo-os
            para instâncias das classes correspondentes.

        atualizar_lista():
            Recarrega todos os modelos do Excel e retorna os DataFrames.

        atualizar():
            Recarrega todos os modelos do Excel sem retornar DataFrames.

        atualizar_celula(system, index, coluna, valor):
            Atualiza o valor de um atributo ou chave específica em um modelo,
            identificando-o pelo sistema e índice, e salva no Excel.

        excluir_linha(system, index):
            Remove um modelo específico de um sistema e salva no Excel.
    """

    

    def __init__(self) -> None:
        """
         Inicializa a classe Models.

        - Cria a lista interna `self.models` vazia.
        - Carrega todos os modelos salvos no arquivo Excel, convertendo-os em
          instâncias das classes correspondentes e armazenando-os em `self.models`.
        """
        self.models = []      

        self.bancoAcesso()

    def createModel(self, model):
        """
        Adiciona um novo modelo à lista e salva no Excel.

        Args:
            model (dict): Dicionário contendo informações do modelo no formato:
                {
                    "system": str,    # Nome do sistema
                    "number": int,    # Número identificador
                    "model": object   # Instância da classe do modelo
                }
        """
        self.models.append(model)

        self.saveModels()
 
    def size(self, system):
        """
        Retorna a quantidade de modelos salvos para um sistema específico.

        Args:
            system (str): Nome interno do sistema (ex: "Brake", "Drivetrain").

        Returns:
            int: Quantidade de modelos encontrados.
        """
        saved= []


        for model in self.models:

            if model["system"] == system:

                saved.append(model)

        return len(saved)

    def getModel(self, system=None, number=None):
        """
        Busca um modelo específico pelo nome do sistema e número.

        Args:
            system (str): Nome interno do sistema.
            number (int): Número identificador do modelo.

        Returns:
            tuple: (system, objeto_modelo) ou None se não encontrado.
        """


        for model in self.models:

            if model["system"] == system and model["number"] == number:

                return system, model["model"]

    def getDataFrames(self):
        """
        Retorna todos os modelos separados por sistema em formato DataFrame.

        Returns:
            tuple: DataFrames na ordem (Drivetrain, Brake, Tire, Kinematics).
        """
        
        dt = []
        brake = []
        tire = []
        kmt = []

        for m in self.models:

            if m["system"] == "Drivetrain":

                dt.append(m["model"])

            elif m["system"] == "Brake":

                brake.append(m["model"])

            elif m["system"] == "Tire":
                
                tire.append(m['model'])

            elif m["system"] == "Kinematics":

                kmt.append(m["model"])


        dt = pd.DataFrame(dt)

        brake = pd.DataFrame(brake)

        tire = pd.DataFrame(tire)

        kmt = pd.DataFrame(kmt)


        return dt, brake, tire, kmt

    def saveModels(self):
        """
        Salva todos os modelos atuais no arquivo Excel.

        - Cada sistema é salvo em uma aba diferente.
        """
        dt, brake, tire, kinematics = self.getDataFrames()

        with pd.ExcelWriter("Software/Systems/management/saved_models.xlsx", engine='openpyxl') as writer:

            dt.to_excel(writer, sheet_name="Drivetrain", index=False)

            brake.to_excel(writer, sheet_name="Brake", index=False)

            tire.to_excel(writer, sheet_name="Tire", index=False)  

            kinematics.to_excel(writer, sheet_name="Kinematics", index=False)

    def convertGUI(self, opt):
            """
            Converte rótulos usados na GUI para nomes internos de sistema.

            Args:
                opt (str): Nome do sistema como exibido na interface gráfica.

            Returns:
                str: Nome interno padrão do sistema.
            """
            if opt == "Cinemática":
                return "Kinematics"
            if opt == "Pneu":
                return "Tire"
            if opt == "Transmissão":
                return "Drivetrain"
            if opt == "Freio":
                return "Brake"
            
    def getDicts(self, system: str, model: dict):
        """
         Converte dicionários de dados em dicionários das classes correspondentes.

        Args:
            system (str): Nome interno do sistema.
            model (list[dict]): Lista de dicionários com atributos do modelo.

        Returns:
            list: Lista de dicionários no formato:
                {
                    "system": str,
                    "number": int,
                    "model": dict
                }
        """
        if system == "Drivetrain":

            lists = []

            for num in range(len(model)):

                list_model = {"system":system, "number": num+1, "model":model[num]}


                lists.append(list_model)


            return lists

        if system == "Brake":

            lists = []

            for num in range(len(model)):

                list_model = {"system":system, "number": num+1, "model":model[num]}


                lists.append(list_model)

        
            return lists

        if system == "Tire":

            lists = []
            
            for num in range(len(model)):
                
                list_model = {"system":system, "number": num+1, "model":model[num]}


                lists.append(list_model)


            return lists

        if system == "Kinematics":

            lists = []

            for num in range(len(model)):

                list_model = {"system":system, "number": num+1, "model":model[num]}


                lists.append(list_model)


            return lists
        
    def bancoAcesso(self):

        """
        Carrega os modelos salvos no Excel para a lista `self.models`.

        - Lê os dados de cada aba do arquivo.
        - Converte cada linha em instância de classe.
        - Adiciona todos os modelos na lista interna.
        """

        dt_df = pd.read_excel("Software/Systems/management/saved_models.xlsx", sheet_name="Drivetrain")

        brake_df = pd.read_excel("Software/Systems/management/saved_models.xlsx", sheet_name="Brake")

        tire_df = pd.read_excel("Software/Systems/management/saved_models.xlsx", sheet_name="Tire")

        kinematics_df = pd.read_excel("Software/Systems/management/saved_models.xlsx", sheet_name="Kinematics")
        
        dt_df = dt_df.to_dict(orient="records")

        brake_df = brake_df.to_dict(orient="records")

        tire_df = tire_df.to_dict(orient="records")

        kinematics_df = kinematics_df.to_dict(orient="records")

        dt_model = self.getDicts(system = "Drivetrain", model = dt_df)

        brake_model = self.getDicts(system='Brake', model = brake_df)

        tire_model = self.getDicts(system="Tire", model = tire_df)

        kinematics_model = self.getDicts(system="Kinematics", model = kinematics_df)


        for model in dt_model:
            self.models.append(model)

        for model in brake_model:
            self.models.append(model)

        for model in tire_model:
            self.models.append(model)

        for model in kinematics_model:
            self.models.append(model)          

    def atualizar_lista(self):

        """
        Recarrega os modelos a partir do Excel e retorna os DataFrames.

        Returns:
            tuple: DataFrames na ordem (Drivetrain, Brake, Tire, Kinematics).
        """

        self.models = []
        self.bancoAcesso()
        dt, brake, tire, kmt = self.getDataFrames()
        return dt, brake, tire, kmt
    
    def atualizar(self):
        """
        Recarrega os modelos a partir do Excel sem retornar DataFrames.
        """
        
        self.models = []
        self.bancoAcesso()

    
    def atualizar_celula(self, system: str, index: int, coluna: str, valor):
        """
        Atualiza o valor de uma célula no modelo identificado pelo sistema e índice.

        Args:
            system (str): Nome do sistema ("Brake", "Drivetrain", etc).
            index (int): Índice (linha) do modelo na lista.
            coluna (str): Nome da coluna a ser alterada.
            valor: Novo valor a ser atribuído.
        """
        # Encontra o modelo correto
        modelos = [m for m in self.models if m["system"] == system]
        if index < 0 or index >= len(modelos):
            return  # índice inválido

        modelo = modelos[index]
        obj = modelo["model"]

        # Atualiza o atributo do objeto, se existir
        if hasattr(obj, coluna):
            try:
                setattr(obj, coluna, valor)
            except Exception as e:
                print(f"Erro ao atualizar {coluna}: {e}")
        else:
            # Se não for atributo, tenta atualizar dict interno (caso exista)
            if isinstance(obj, dict) and coluna in obj:
                obj[coluna] = valor

        # Atualizar o dicionário do modelo para refletir a mudança (opcional)
        # Pode ser necessário se for usado para salvar/exibir
        if hasattr(obj, "__dict__"):
            modelo["model"] = obj  # já está atualizado

        # Após alteração, salvar no Excel
        self.saveModels()

    def excluir_linha(self, system: str, index: int):
        """
        Remove o modelo (linha) do sistema dado pelo índice.

        Args:
            system (str): Nome do sistema.
            index (int): Índice do modelo na lista.
        """
        modelos = [m for m in self.models if m["system"] == system]
        if index < 0 or index >= len(modelos):
            return

        modelo_remover = modelos[index]
        try:
            self.models.remove(modelo_remover)
        except ValueError:
            pass

        # Após remoção, salvar no Excel
        self.saveModels()