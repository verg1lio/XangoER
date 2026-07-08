import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from Constants.constants import PI23, SQRT3_2


class Motor:
    """Permanent Magnet Synchronous Motor (PMSM) — parameter & coordinate-transform model.

    This class is a stateless container of PMSM electrical/mechanical/thermal
    parameters plus the Park/Clarke transformations used by the simulation
    logger. The actual control loops (FOC, current PIs, speed PI) and the
    plant ODE live in :class:`Simulation.Simulation` — Motor does not own any
    controller state or dynamics.

    Parameters
    ----------
    rs : float
        Stator resistance [Ω].
    ld, lq : float
        Direct- and quadrature-axis inductances [H].
    jm : float
        Rotor inertia [kg·m²].
    kf : float
        Viscous friction coefficient [N·m·s].
    lambda_m : float
        Permanent magnet flux linkage [Wb].
    p : int
        Number of pole pairs.
    valor_mu : float
        Modulation index for inverter voltage scaling (0–1).
    speed_ref : float, optional
        Reference rotor speed [rad/s] for the over-speed soft limiter.
        Default 0.0.
    k_h_iron : float, optional
        Hysteresis iron-loss coefficient [W/Hz]. Loss ≈ k_h·|f_elec|.
        Default 0.9 — calibrated so iron loss ≈ 600 W at 670 Hz (EMRAX 228
        @ 4000 RPM, 60/40 hysteresis/eddy split per datasheet).
    k_e_iron : float, optional
        Eddy-current iron-loss coefficient [W/Hz²]. Loss ≈ k_e·f_elec².
        Default 9e-4 — pairs with k_h above to give ~1 kW total iron loss
        at 4000 RPM. Tune lower for low-iron motors (PCB stators, etc.).
    alpha_cu : float, optional
        Copper temperature coefficient [1/K] for Rs(T) = Rs0·[1+α·(T−T_ref)].
        Default 0.00393 (annealed copper).
    T_ref : float, optional
        Reference temperature [°C] at which Rs0 was measured. Default 25.
    T_alarm, T_max : float, optional
        Thermal derating thresholds [°C]. Below T_alarm: no derating.
        Above T_max: zero current. Linear ramp in between. Defaults
        130 / 160 °C (typical winding insulation class H limits).

    Attributes
    ----------
    pi23 : float
        Constant 2π/3 used for phase transformations.
    m, C : float
        Thermal mass [kg] and copper specific heat [J/(kg·K)] used by the
        winding thermal ODE (defaults EMRAX 228 LC: 13.5 kg, 385 J/(kg·K)).
    max_current : float
        Maximum phase current amplitude [A] (peak, = √2·I_rms) used by
        Simulation as the upper bound on iq_ref. Configurable via the
        ``max_current`` kwarg; default 300·√2 A.
    Vdc : float
        DC-link voltage [V] used as the default DC voltage when no battery
        model is provided. Default 600 V.
    """

    def __init__(self, rs, ld, lq, jm, kf, lambda_m, p, valor_mu, speed_ref=0.0,
                 k_h_iron=0.9, k_e_iron=9e-4,
                 alpha_cu=0.00393, T_ref=25.0,
                 T_alarm=130.0, T_max=160.0,
                 h_cool=10.0, A_cool=0.15, T_ambient=25.0,
                 max_current=None):
        self.pi23 = PI23

        # Electrical / mechanical parameters
        # rs0 = nominal Rs at T_ref (immutable reference for Rs(T) law).
        # rs    = effective Rs at current temperature (updated by Simulation
        # each step). Tools that need the cold-resistance use rs0.
        self.rs0 = rs
        self.rs = rs
        self.ld = ld
        self.lq = lq
        self.jm = jm
        self.kf = kf
        self.lambda_m = lambda_m
        self.p = p
        self.valor_mu = valor_mu

        # Thermal model — EMRAX 228 LC defaults
        self.m = 13.5      # thermal mass [kg]
        self.C = 385.0     # specific heat (copper) [J/(kg·K)]

        # Operational limits
        # max_current é a AMPLITUDE de pico da corrente de fase [A]
        # (frame dq amplitude-invariante: iq = √2·I_rms).  O teto de torque
        # resultante é T = 1.5·p·λm·max_current — escolha o valor que
        # reproduz o torque de pico do datasheet (EMRAX 228: 230 N·m →
        # max_current ≈ 323 A com λm=0.04748, p=10).
        if max_current is not None:
            self.max_current = float(max_current)
        else:
            self.max_current = 300.0 * np.sqrt(2)
        self.Vdc = 600.0

        # Reference speed (used by the simulation's soft over-speed limiter)
        self.speed_ref = speed_ref

        # Iron losses (Steinmetz form, separated terms)
        self.k_h_iron = float(k_h_iron)
        self.k_e_iron = float(k_e_iron)

        # Temperature compensation
        self.alpha_cu = float(alpha_cu)
        self.T_ref    = float(T_ref)

        # Thermal derating thresholds
        self.T_alarm = float(T_alarm)
        self.T_max   = float(T_max)

        # Cooling model — convective: P_cool = h_cool · A_cool · (T − T_ambient)
        # EMRAX 228 LC defaults: natural convection (h≈10), surface ≈0.15 m²
        # For forced liquid: h ≈ 500–2000 W/(m²·K); for forced air: h ≈ 50–200
        self.h_cool    = float(h_cool)
        self.A_cool    = float(A_cool)
        self.T_ambient = float(T_ambient)

    def inverse_park_transform(self, vd, vq, theta_e):
        """Inverse Park transform: dq voltages → abc phase voltages.

        Parameters
        ----------
        vd, vq : float
            Direct- and quadrature-axis voltages [V].
        theta_e : float
            Electrical rotor angle [rad].

        Returns
        -------
        (vs1, vs2, vs3, 0.0) : tuple of float
            Three-phase voltages [V]. The trailing 0.0 is kept for
            backwards-compatible unpacking in :mod:`Simulation`.
        """
        cos_theta = np.cos(theta_e)
        sin_theta = np.sin(theta_e)
        valpha = vd * cos_theta - vq * sin_theta
        vbeta  = vd * sin_theta + vq * cos_theta
        vs1 = valpha
        vs2 = -0.5 * valpha + SQRT3_2 * vbeta
        vs3 = -0.5 * valpha - SQRT3_2 * vbeta
        return vs1, vs2, vs3, 0.0

    def abc_currents_from_dq(self, isd, isq, theta_e, flux_d):
        """Compute abc phase currents from dq currents (amplitude-invariant).

        The torque model uses ``Kt = 1.5·p·λm``, which is consistent with the
        amplitude-invariant inverse Park used here (no √(2/3) scaling).

        Parameters
        ----------
        isd, isq : float
            Direct- and quadrature-axis currents [A].
        theta_e : float
            Electrical rotor angle [rad].
        flux_d : float
            Flux linkage along d-axis [Wb] (passed through unchanged for logging).

        Returns
        -------
        (is1, is2, is3, flux_d, flux_d, flux_d) : tuple of float
        """
        is1 = isd * np.cos(theta_e)             - isq * np.sin(theta_e)
        is2 = isd * np.cos(theta_e - self.pi23) - isq * np.sin(theta_e - self.pi23)
        is3 = isd * np.cos(theta_e + self.pi23) - isq * np.sin(theta_e + self.pi23)
        return is1, is2, is3, flux_d, flux_d, flux_d