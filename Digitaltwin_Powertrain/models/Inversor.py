import numpy as np

# ============================
# Constantes globais
# ============================
PI = np.pi
PI23 = 2 * np.pi / 3       # 120°
SQRT3_2 = np.sqrt(3) / 2   # √3/2 (usado na Clarke)
TWO_THIRDS = 2 / 3         # Fator da Clarke


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

    def __init__(self, eficiencia=0.95, freq_chaveamento=10000):
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
        if self.tete >= 2 * PI:
            self.tete -= 2 * PI

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
