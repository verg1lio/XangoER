import numpy as np
class Pedal:
    """A Pedal object representing an accelerator input for modulation control.

    This class models a simple accelerator pedal, producing a normalized output
    between 0 and 1 that can be used as a modulation reference for an inverter.

    Parameters
    ----------
    ganho : float, optional
        Gain factor to scale the pedal reference (default is 1.0).

    Attributes
    ----------
    posicao : float
        Current pedal position, normalized between 0.0 (released) and 1.0 (fully pressed).
    ganho : float
        Gain applied to the pedal position when generating the reference.

    Methods
    -------
    set_posicao(valor):
        Sets the pedal position, ensuring saturation between 0.0 and 1.0.
    get_referencia():
        Returns the modulation reference scaled by the gain.
    """

    def __init__(self, ganho=1.0):
        self.posicao = 0.0
        self.ganho = ganho

    def set_posicao(self, valor):
        """Sets the pedal position, with saturation between 0.0 and 1.0."""
        self.posicao = np.clip(valor, 0.0, 1.0)

    def get_referencia(self):
        """Returns the pedal reference scaled by the gain."""
        return self.posicao * self.ganho