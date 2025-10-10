import sys
import os

# Adiciona o diretório pai ao path do Python
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from models import PIDController as PID
from Constants.constants import PI23, RQ23, SQRT3_2
class Motor:
    """Permanent Magnet Synchronous Motor (PMSM) with Field-Oriented Control (FOC).

    This class models a PMSM including its electrical and mechanical
    dynamics, along with Field-Oriented Control (FOC) for current and speed
    regulation. It also supports the application of external torque profiles.

    Parameters
    ----------
    rs : float
        Stator resistance [Ω].
    ld : float
        Direct-axis inductance [H].
    lq : float
        Quadrature-axis inductance [H].
    jm : float
        Rotor inertia [kg·m²].
    kf : float
        Friction coefficient [N·m·s].
    lambda_m : float
        Permanent magnet flux linkage [Wb].
    p : int
        Number of pole pairs.
    valor_mu : float
        Modulation index for inverter voltage scaling.
    TL : bool, optional
        Enable external load torque. Default is False.
    torque : float or callable, optional
        Constant external torque value [N·m] or a function torque(t) that
        defines time-dependent torque profile. Default is 0.0.

    Attributes
    ----------
    pi23 : float
        Constant π/3 used for phase transformations.
    rq23 : float
        Constant 2/3 used for Clarke transformation.
    one_point_five_p : float
        Precomputed constant 1.5 * p.
    torque_constant : float
        Torque constant Kt = 1.5 * p * λm [N·m/A].
    inv_ld : float
        Inverse of d-axis inductance [1/H].
    inv_lq : float
        Inverse of q-axis inductance [1/H].
    inv_jm : float
        Inverse of rotor inertia [1/(kg·m²)].
    m : float
        Thermal mass parameter.
    C : float
        Thermal capacity parameter.
    inv_mC : float
        Inverse of m*C for thermal calculations.
    id_controller : PIDController
        PI controller for d-axis current regulation.
    iq_controller : PIDController
        PI controller for q-axis current regulation.
    speed_controller : PIDController
        PI controller for rotor speed regulation.
    TL : bool
        Whether external torque is enabled.
    _external_torque : callable
        Function defining external torque as a function of time.
    max_current : float
        Maximum allowed phase current [A].
    Vdc : float
        DC-link voltage [V].
    Vs : float
        Initial stator voltage [V].
    Vlimit : float
        Maximum voltage magnitude [V].
    id_ref : float
        Reference value for d-axis current [A].
    iq_ref : float
        Reference value for q-axis current [A].
    speed_ref : float
        Reference rotor speed [rad/s].

    Methods
    -------
    set_external_torque(torque):
        Defines an external torque as a constant or function of time.
    enable_external_torque(enable):
        Enables or disables application of external torque.
    set_load(t):
        Returns the load torque at time t [N·m].
    field_oriented_control(isd, isq, wm, dt):
        Executes the FOC control loop, computing voltage commands vd, vq.
    inverse_park_transform(vd, vq, theta_e):
        Transforms dq voltages into abc phase voltages.
    abc_currents_from_dq(isd, isq, theta_e, flux_d):
        Computes abc phase currents from dq currents.
    """

    def __init__(self, rs, ld, lq, jm, kf, lambda_m, p, valor_mu, TL=False, torque=0.0, speed_ref= 0.0):
        self.pi23 = PI23
        self.rq23 = RQ23
        self.rs = rs
        self.ld = ld
        self.lq = lq
        self.jm = jm
        self.kf = kf
        self.lambda_m = lambda_m
        self.p = p
        self.valor_mu = valor_mu

        # Precomputed constants
        self.one_point_five_p = 1.5 * p
        self.torque_constant = self.one_point_five_p * lambda_m
        self.inv_ld = 1.0 / ld
        self.inv_lq = 1.0 / lq
        self.inv_jm = 1.0 / jm

        # Thermal model parameters
        self.m = 22.0
        self.C = 0.385
        self.inv_mC = 1.0 / (self.m * self.C)

        # Controllers
        self.id_controller = PID.Controller(kp=0.5, ki=100.0, kd = 0, limit=600.0)
        self.iq_controller = PID.Controller(kp=0.5, ki=100.0, kd = 0 ,limit=600.0)
        self.speed_controller = PID.Controller(kp=10.0, ki=5,kd = 0, limit=600.0)

        # External torque setup
        self.TL = bool(TL)
        if callable(torque):
            self._external_torque = torque
        else:
            v = float(torque)
            self._external_torque = (lambda t, v=v: v)

        # Limits
        self.max_current = 220.0 * np.sqrt(2)
        self.Vdc = 600
        self.Vs = (self.valor_mu * self.Vdc) / np.sqrt(3)
        self.Vlimit = 600

        # References
        self.id_ref = 0.0
        self.iq_ref = 0.0
        self.speed_ref = speed_ref

    def set_external_torque(self, torque):
        """Defines the external torque source.

        Parameters
        ----------
        torque : float or callable
            Constant torque [N·m] or function torque(t).
        """
        if callable(torque):
            self._external_torque = torque
        else:
            v = float(torque)
            self._external_torque = (lambda t, v=v: v)

    def enable_external_torque(self, enable: bool):
        """Enable or disable the application of external load torque.

        Parameters
        ----------
        enable : bool
            True to enable, False to disable.
        """
        self.TL = bool(enable)

    def set_load(self, t):
        """Returns load torque at time t.

        Parameters
        ----------
        t : float
            Time [s].

        Returns
        -------
        float
            Load torque [N·m].
        """
        if self.TL and (self._external_torque is not None):
            try:
                return float(self._external_torque(t))
            except Exception:
                pass
        if t < 0.5:
            return 0.0
        elif t < 1.0:
            return 100.0
        else:
            return -200.0

    def field_oriented_control(self, isd, isq, wm, dt):
        """Executes the Field-Oriented Control (FOC) algorithm.

        Parameters
        ----------
        isd : float
            Direct-axis stator current [A].
        isq : float
            Quadrature-axis stator current [A].
        wm : float
            Mechanical rotor speed [rad/s].
        dt : float
            Time step [s].

        Returns
        -------
        vd : float
            Direct-axis voltage command [V].
        vq : float
            Quadrature-axis voltage command [V].
        id_ref : float
            Reference direct-axis current [A].
        iq_ref : float
            Reference quadrature-axis current [A].
        speed_error : float
            Speed regulation error [rad/s].
        """
        speed_error = self.speed_ref - wm
        torque_ref = self.speed_controller.update(speed_error, dt)

        iq_ref = torque_ref / self.torque_constant if self.torque_constant != 0 else 0.0
        id_ref = 0.0

        if abs(iq_ref) > self.max_current:
            iq_ref = np.sign(iq_ref) * self.max_current

        error_d = id_ref - isd
        error_q = iq_ref - isq

        we = self.p * wm
        decoupling_d = -we * self.lq * isq
        decoupling_q = we * (self.ld * isd + self.lambda_m)

        vd = self.id_controller.update(error_d, dt) + decoupling_d
        vq = self.iq_controller.update(error_q, dt) + decoupling_q

        vmag = np.sqrt(vd**2 + vq**2)
        if vmag > self.Vlimit:
            v_scale = self.Vlimit / vmag
            vd *= v_scale
            vq *= v_scale

        return vd, vq, id_ref, iq_ref, speed_error

    def inverse_park_transform(self, vd, vq, theta_e):
        """Inverse Park transform to obtain abc voltages.

        Parameters
        ----------
        vd : float
            Direct-axis voltage [V].
        vq : float
            Quadrature-axis voltage [V].
        theta_e : float
            Electrical rotor angle [rad].

        Returns
        -------
        vs1 : float
            Phase A voltage [V].
        vs2 : float
            Phase B voltage [V].
        vs3 : float
            Phase C voltage [V].
        0.0 : float
            Dummy return for compatibility.
        """
        cos_theta = np.cos(theta_e)
        sin_theta = np.sin(theta_e)
        valpha = vd * cos_theta - vq * sin_theta
        vbeta = vd * sin_theta + vq * cos_theta
        vs1 = valpha
        vs2 = -0.5 * valpha + SQRT3_2 * vbeta
        vs3 = -0.5 * valpha - SQRT3_2 * vbeta
        return vs1, vs2, vs3, 0.0

    def abc_currents_from_dq(self, isd, isq, theta_e, flux_d):
        """Compute abc phase currents from dq currents.

        Parameters
        ----------
        isd : float
            Direct-axis current [A].
        isq : float
            Quadrature-axis current [A].
        theta_e : float
            Electrical rotor angle [rad].
        flux_d : float
            Flux linkage along d-axis [Wb].

        Returns
        -------
        is1 : float
            Phase A current [A].
        is2 : float
            Phase B current [A].
        is3 : float
            Phase C current [A].
        fs_val : float
            Flux contribution from d-axis [Wb].
        """
        cos_theta = np.cos(theta_e)
        sin_theta = np.sin(theta_e)
        is1 = self.rq23 * (isd * cos_theta - isq * sin_theta)
        is2 = self.rq23 * (isd * np.cos(theta_e - self.pi23) - isq * np.sin(theta_e - self.pi23))
        is3 = self.rq23 * (isd * np.cos(theta_e + self.pi23) - isq * np.sin(theta_e + self.pi23))
        fs_val = self.rq23 * flux_d
        return is1, is2, is3, fs_val, fs_val, fs_val
