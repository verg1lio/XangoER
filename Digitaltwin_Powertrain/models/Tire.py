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

    def Tire_forces(self, Fz, Ls):
        """Computes the longitudinal tire force using a simplified Pacejka model.

        Parameters
        ----------
        Fz : float
            Vertical load on the tire [N].
        Ls : float
            Slip ratio (dimensionless).

        Returns
        -------
        float
            Longitudinal tire force [N].
        """
        E, Cy, Cx, Cz, c1, c2 = self.pacejka_params

        # Slip stiffness coefficient
        Cs = c1 * np.sin(2 * np.arctan(Fz / c2))

        # Peak factor (friction coefficient × vertical load)
        D = self.tire_friction_coef * Fz

        # Avoid division by zero
        if (Cx * D) == 0:
            return 0.0

        # Stiffness factor
        Bx = Cs / (Cx * D)

        # Pacejka argument
        arg_arctan = Bx * Ls

        # Longitudinal force (simplified Pacejka "Magic Formula")
        tire_longitudinal_force = D * np.sin(
            Cx * np.arctan(arg_arctan - E * (arg_arctan - np.arctan(arg_arctan)))
        )
        return tire_longitudinal_force

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
