import numpy as np

class Tire:
    """A Tire object based on a simplified Pacejka Magic Formula model.

    This class models the longitudinal tire force as a function of slip ratio,
    vertical load, and tire-road friction coefficient. It uses a simplified
    Pacejka formulation with configurable parameters.

    Parameters
    ----------
    pacejka_params : tuple
        Pacejka parameters (E, Cy, Cx, Cz, c1, c2) for the tire model.
    tire_friction_coef : float
        Tire-road friction coefficient (μ).

    Attributes
    ----------
    pacejka_params : tuple
        Tuple containing Pacejka model coefficients (E, Cy, Cx, Cz, c1, c2).
    tire_friction_coef : float
        Friction coefficient (μ) used in the longitudinal force computation.

    Methods
    -------
    Tire_forces(Fz, Ls):
        Computes the longitudinal tire force given vertical load (Fz) and slip ratio (Ls).
    SlipRatio(velocidade_angular, raio_pneu, velocidade_linear, eps=0.1):
        Computes the slip ratio based on wheel angular speed, tire radius, and vehicle speed.
    """

    def __init__(self, pacejka_params, tire_friction_coef):
        self.pacejka_params = pacejka_params
        self.tire_friction_coef = tire_friction_coef

    def Tire_forces(self, Fz, s):
        """
        Improved longitudinal Pacejka MF model using structured parameters.
        """

        E, Cy, Cx, Cz, c1, c2 = self.pacejka_params

        # B vem diretamente do parâmetro c1
        B = c1

        # Cx já é o shape factor C
        C = Cx

        # D = mu * Fz
        mu = self.tire_friction_coef
        D = mu * Fz

        # Magic Formula
        arg = B * s
        Fx = D * np.sin(
            C * np.arctan(arg - E * (arg - np.arctan(arg)))
        )

        return Fx

    @staticmethod
    def SlipRatio(velocidade_angular, raio_pneu, velocidade_linear, eps=0.1):
        """Computes slip ratio given wheel angular velocity, tire radius, and vehicle speed.

        Parameters
        ----------
        velocidade_angular : float or array_like
            Angular velocity of the wheel [rad/s].
        raio_pneu : float
            Tire effective rolling radius [m].
        velocidade_linear : float or array_like
            Linear velocity of the vehicle [m/s].
        eps : float, optional
            Small threshold to avoid division by zero. Default is 0.1.

        Returns
        -------
        float or ndarray
            Slip ratio (dimensionless). Positive values indicate wheel spin,
            negative values indicate wheel lock.
        """
        v = np.array(velocidade_linear, dtype=float)
        omega = np.array(velocidade_angular, dtype=float)
        denom = np.where(np.abs(v) < eps, eps, v)  # prevent division by near-zero
        return (omega * raio_pneu - v) / denom
