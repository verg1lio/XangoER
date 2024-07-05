class LagCompensator:
  def __init__(self, Kc: float, beta: float) -> None:
    self.Kc = Kc
    self.beta = beta
    self.gain = Kc*beta
    print("we on babyy!")

  def Gc_of_s(self, s: float):
    """
    Runs a signl through the transfer function
    
    Parameters
    ----------
    s : float
      Signal input
    """
    K = self.gain
    T = self.T
    beta = self.beta

    return K*(T*s + 1)/(beta*T*s + 1)

  def G_of_s(self, s: float):
    """
    Runs a signal through the open-loop transfer function

    Parameters
    ----------
    s : float
      Signal input
    """

    return 1 / (s*(s + 1)*(0.5*s + 1))
  
  def G1_of_s(self, s: float):
    """
    Runs a signl through the G1 function
    
    Parameters
    ----------
    s : float
      Signal input
    """

    return self.gain * self.G_of_s(s)