class Peso:
    """
    Classe que calcula o peso total de um sistema de powertrain, considerando os pesos 
    individuais da bateria, do inversor, do motor e do chicote elétrico.
    """
    
    def __init__(self, peso_bateria, peso_inversor, peso_motor, peso_chicote):
        self.peso_bateria = peso_bateria
        self.peso_inversor = peso_inversor
        self.peso_motor = peso_motor
        self.peso_chicote = peso_chicote

    def peso_total(self):
        """
        Calcula o peso total dos componentes do sistema de powertrain.

        Soma o peso da bateria, do inversor, do motor e do chicote elétrico.

        Returns
        -------
        float
            Somatório dos pesos dos componentes do sistema.

        Examples
        --------
        >>> peso_pwt = Peso(10, 10, 65, 5)
        >>> peso_pwt.peso_total()
        90
        """
        return self.peso_bateria + self.peso_inversor + self.peso_motor + self.peso_chicote

    @staticmethod
    def example():
        """
        Método de exemplo para demonstrar o uso da classe Peso.

        Exibe o peso total calculado para um conjunto de valores fornecidos.
        """
        peso_pwt = Peso(10, 10, 65, 5)
        total = peso_pwt.peso_total()
        print(f"O peso total é {total} kg")

# Exemplo de uso da classe Peso
Peso.example()
