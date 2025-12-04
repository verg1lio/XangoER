import sys
import os

# Adiciona o diretório pai ao path do Python
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from Constants.constants import PI23, SQRT3_2, TWO_THIRDS


class Inversor:
    """A three-phase sinusoidal inverter model for DC-AC conversion.

    This class models a simplified three-phase inverter that converts a DC input
    voltage (from a battery pack) into sinusoidal AC phase voltages. It supports
    synchronous reference frame (Park) transformation for motor control.

    Parameters
    ----------
    eficiencia : float, optional
        Conversion efficiency of the inverter (default is 0.95).
    freq_chaveamento : float, optional
        Switching frequency [Hz] (default is 10000).

    Attributes
    ----------
    eficiencia : float
        Efficiency of the inverter.
    freq_chaveamento : float
        Switching frequency [Hz].
    Vdc : float
        DC input voltage [V].
    tete : float
        Electrical angle of the stator [rad].
    ws : float
        Synchronous electrical angular velocity [rad/s].
    h : float
        Switching period [s].
    pi23 : float
        Constant = 2π/3.
    sqrt3_2 : float
        Constant = √3/2.
    two_thirds : float
        Constant = 2/3.

    Methods
    -------
    set_Vdc(Vdc):
        Sets the DC input voltage.
    set_ws(ws):
        Sets the synchronous electrical angular velocity.
    source_voltage(modulacao):
        Computes the three-phase sinusoidal voltages based on DC input and modulation index.
    park_transform(va, vb, vc, theta_e):
        Applies Clarke and Park transforms to compute (vq, vd).
    """

    def __init__(self, eficiencia=0.95, freq_chaveamento=1000):
        self.eficiencia = eficiencia
        self.freq_chaveamento = freq_chaveamento
        self.Vdc = 0.0
        self.tete = 0.0
        self.ws = 0.0
        self.h = 1.0 / freq_chaveamento

        self.pi23 = PI23
        self.sqrt3_2 = SQRT3_2
        self.two_thirds = TWO_THIRDS

    def set_Vdc(self, Vdc):
        """Sets the DC input voltage [V]."""
        self.Vdc = Vdc

    def set_ws(self, ws):
        """Sets the synchronous electrical angular velocity [rad/s]."""
        self.ws = ws

    def source_voltage(self, modulacao):
        """Computes three-phase sinusoidal voltages given DC input and modulation index.

        Parameters
        ----------
        modulacao : float
            Modulation index (0.0–1.0) typically obtained from a Pedal object.

        Returns
        -------
        tuple of float
            Phase voltages (vs1, vs2, vs3).
        """
        self.tete += self.h * self.ws
        if self.tete >= 2 * np.pi:
            self.tete -= 2 * np.pi

        cos_tete = np.cos(self.tete)
        cos_tete_pi23 = np.cos(self.tete - self.pi23)
        cos_tete_pi23_2 = np.cos(self.tete + self.pi23)

        V_amp = self.Vdc * modulacao * self.eficiencia

        vs1 = V_amp * cos_tete
        vs2 = V_amp * cos_tete_pi23
        vs3 = V_amp * cos_tete_pi23_2

        return vs1, vs2, vs3

    def park_transform(self, va, vb, vc, theta_e):
        """Applies Clarke and Park transforms to compute dq-axis voltages.

        Parameters
        ----------
        va : float
            Phase-a voltage [V].
        vb : float
            Phase-b voltage [V].
        vc : float
            Phase-c voltage [V].
        theta_e : float
            Electrical angle [rad].

        Returns
        -------
        tuple of float
            dq-axis voltages (vq, vd).
        """
        cos_theta = np.cos(theta_e)
        sin_theta = np.sin(theta_e)

        valpha = self.two_thirds * (va - 0.5 * vb - 0.5 * vc)
        vbeta = self.two_thirds * (self.sqrt3_2 * vb - self.sqrt3_2 * vc)

        vd = valpha * cos_theta + vbeta * sin_theta
        vq = -valpha * sin_theta + vbeta * cos_theta

        return vq, vd
