class envio:

    """
    Método para receber as entradas das caixas de texto e armazenar os valores após o envio
    Estruture a ordem de chamada dos parâmetros conforme a ordem dos parâmetros de seu código
    """

    def __init__(self):
        self.value = None
        self.text = None
        self.dicionario = {}



    def atualizar_dict(self):
        self.dicionario[self.text] = self.value

    def dicionarios(self):
        return self.dicionario
