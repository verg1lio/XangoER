import numpy as np
import control.matlab as ml
import matplotlib.pyplot as plt
from typing import Callable

class LagCompensator:
  def __init__(self, Kc: float, beta: float, timeConstant: float) -> None:
    self.Kc = Kc # no ideia wtf this is
    self.ß = beta
    self.K = Kc*beta # gain
    self.T = timeConstant # T

  def Gc_of_s(self, s: complex):
    """
    Runs a signl through the transfer function
    
    Parameters
    ----------
    s : complex
      Signal input
    """

    return self.K*(self.T*s + 1)/(self.ß*self.T*s + 1)

  def G_of_s(self, s: complex):
    """
    Runs a signal through the open-loop transfer function

    Parameters
    ----------
    s : complex
      Signal input
    """

    return 1 / (s*(s + 1)*(0.5*s + 1))
  
  def G1_of_s(self, s: complex):
    """
    Runs a signl through the G1 function
    
    Parameters
    ----------
    s : complex
      Signal input
    """

    return self.K * self.G_of_s(s)


def ogataExample(signals: np.ndarray):
  G = LagCompensator(1, 10, 0.2).Gc_of_s # transfer function

  magnitude = 20*np.log10(np.abs(G(signals)))
  phase = np.arctan2(np.imag(G(signals)), np.real(G(signals))) * 180/np.pi

  return magnitude, phase

def plotExample(example: Callable):
  frequencyInHz = np.logspace(-4, 4, 1000)
  frequencyInRadPerSec = 2*np.pi*frequencyInHz
  signals = 1j*frequencyInRadPerSec

  magninute, phase = example(signals)

  plt.figure()
  plt.subplot(211)
  plt.semilogx(frequencyInRadPerSec, magninute)
  plt.ylabel('dB')
  plt.subplot(212)
  plt.semilogx(frequencyInHz, phase)
  plt.ylabel('phase')
  plt.xlabel('rad/s')
  plt.xlim([np.power(10.0, -4), np.power(10, 3)])
  plt.savefig('bode-diagram.png')

plotExample(ogataExample)