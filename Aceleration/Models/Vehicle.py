"""
Vehicle.py — Modelo longitudinal do veículo Fórmula SAE

Histórico de correções
-----------------------
v2 (este):
  Bug latente corrigido em calculate_reflected_inertia():
    ANTES: j_translation = mass × r²  depois divide tudo por N²
           → resulta em mass × r²/N² = mass × (r/N)²  ← correto por acidente
           → mas a função retorna J_wrong se chamada isoladamente
           → Simulation.py compensava com uma subtração separada
    DEPOIS: formula diretamente correta, sem compensação externa necessária:
           j_translation = mass × (r/N)²  (reflexão direta para o motor)
           j_rot_roda+tx refletidas = J_slow / N²
           j_total = j_motor_side (não incluído aqui) + j_rot_refl + j_transl

  Adicionado: J_sprocket (inércia do sprocket/coroa no lado do motor).
    Componentes no lado de ALTA VELOCIDADE (motor) somam-se DIRETAMENTE
    ao rotor — sem reflexão por N².
    Exemplo: sprocket de alumínio FSAE (m≈0.3 kg, r_ext≈55 mm)
             J ≈ 0.5 × 0.3 × 0.055² ≈ 0.0005 kg·m²
    Se desconhecido, use 0.0 e o modelo ainda é conservador.

Equação de inércia refletida (dois corpos separados — ver Simulation.py):
  ─────────────────────────────────────────────────────────────────────
  Para a equação do MOTOR (J_eff):
    J_eff = J_motor + J_sprocket + J_rot_refletida
  Para a equação do VEÍCULO (Newton translacional):
    m × dv/dt = F_tração − F_resistência
  ─────────────────────────────────────────────────────────────────────
  Isso evita dupla contagem da massa translacional m×(r/N)².
  calculate_reflected_inertia() retorna J_rot_refletida + J_translacional
  para uso em modelos de CORPO ÚNICO (acoplamento rígido). Para o modelo
  de dois corpos de Simulation.py, apenas J_rot_refletida é usada —
  o campo j_rot_reflected é disponibilizado separadamente.
"""

import numpy as np


class Vehicle:
    """Vehicle model for longitudinal dynamic simulation.

    Models aerodynamic drag, rolling resistance, road grade, load
    transfer, and drivetrain inertia reflected to the motor shaft.

    Parameters
    ----------
    mass : float
        Total vehicle mass (includes driver) [kg].
    wheel_radius : float
        Effective rolling radius of the tyre [m].
    wheel_mass : float
        Mass of one wheel+tyre assembly [kg].
    drag_coeff : float
        Aerodynamic drag coefficient (Cd).
    frontal_area : float
        Vehicle frontal area [m²].
    rolling_resistance : float
        Rolling resistance coefficient (Cr).
    road_grade : float, optional
        Road slope angle [rad]. Default 0.
    environment_density : float, optional
        Air density [kg/m³]. Default 1.225.
    L : float, optional
        Wheelbase [m].
    h : float, optional
        Height of centre of gravity [m].
    dist_cg : float, optional
        Longitudinal distance from the FRONT (non-driven) axle to the CG [m].
        Used as ``a`` in Fz_rear = m·g·a/L (moment about front axle → rear load).
        For typical FSAE RWD with dist_cg=0.6 m, L=1.5 m: 40 % weight on rear.
    n_driven_wheels : int, optional
        Number of driven (powered) wheels. Default 2 (rear-wheel drive).
    lift_coeff : float, optional
        Aerodynamic vertical-force coefficient (Cl). Positive Cl = downforce
        (force pointing down, increasing tyre load). Default 2.0 for a
        full-aero FSAE package; use 0.0 for cars without wings.
    lift_area : float, optional
        Reference area for Cl [m²]. Default 0.7 m² (typical FSAE wing area).
    downforce_balance_rear : float, optional
        Fraction of total downforce applied to the rear axle (0–1).
        Default 0.6 (slight rear bias typical of FSAE acceleration setups).
    """

    def __init__(self,
                 mass: float,
                 wheel_radius: float,
                 wheel_mass: float,
                 drag_coeff: float,
                 frontal_area: float,
                 rolling_resistance: float,
                 road_grade: float = 0.0,
                 environment_density: float = 1.225,
                 L: float = None,
                 h: float = None,
                 dist_cg: float = None,
                 n_driven_wheels: int = 2,
                 rolling_resistance_v_ref: float = 150.0,
                 has_wing: bool = True,
                 lift_coeff_front: float = 1.0,
                 area_front: float = 0.35,
                 lift_coeff_rear: float = 1.0,
                 area_rear: float = 0.35,
                 **kwargs):

        self.mass               = float(mass)
        self.wheel_radius       = float(wheel_radius)
        self.wheel_mass         = float(wheel_mass)
        self.drag_coeff         = float(drag_coeff)
        self.frontal_area       = float(frontal_area)
        self.rolling_resistance = float(rolling_resistance)
        self.road_grade         = float(road_grade)
        self.environment_density = float(environment_density)
        self.L                  = L
        self.h                  = h
        self.dist_cg            = dist_cg
        self.n_driven_wheels    = int(n_driven_wheels)

        # Cr(v) = Cr0 · (1 + v / v_ref) — velocidade de referência para
        # variação da resistência ao rolamento. 150 m/s ≈ +13% a 72 km/h.
        self.rolling_resistance_v_ref = float(rolling_resistance_v_ref)

        # Aerodinâmica vertical — asas dianteira e traseira separadas
        self.has_wing         = bool(has_wing)
        self.lift_coeff_front = float(lift_coeff_front)
        self.area_front       = float(area_front)
        self.lift_coeff_rear  = float(lift_coeff_rear)
        self.area_rear        = float(area_rear)

        # Atributos agregados para retrocompatibilidade (calculados dos componentes)
        self.lift_coeff   = self.lift_coeff_front + self.lift_coeff_rear
        self.lift_area    = self.area_front + self.area_rear
        _df_front = self.lift_coeff_front * self.area_front
        _df_rear  = self.lift_coeff_rear  * self.area_rear
        _df_total = _df_front + _df_rear
        self.downforce_balance_rear = _df_rear / _df_total if _df_total > 0 else 0.5

        self.g        = 9.81
        self.half_rho = 0.5 * self.environment_density

    # ─────────────────────────────────────────────────────────────────────────
    # Forças de resistência
    # ─────────────────────────────────────────────────────────────────────────

    def calculate_resistance_forces(self, velocity: float) -> float:
        """Total longitudinal resistance force [N].

        Includes aerodynamic drag, velocity-dependent rolling resistance,
        and road grade.

        Rolling resistance: Cr(v) = Cr0 · (1 + v / v_ref)
        At v_ref=150 m/s: +6.7% @ 10 m/s, +13% @ 20 m/s, +20% @ 30 m/s.
        """
        v      = abs(float(velocity))
        Cr_v   = self.rolling_resistance * (1.0 + v / self.rolling_resistance_v_ref)
        Faero  = self.half_rho * self.drag_coeff * self.frontal_area * v ** 2
        Froll  = Cr_v * self.mass * self.g * np.cos(self.road_grade)
        Fgrade = self.mass * self.g * np.sin(self.road_grade)
        return Faero + Froll + Fgrade

    def calculate_load_torque(self, velocity: float, transmission) -> float:
        """Motor-equivalent torque to overcome vehicle resistance [N·m]."""
        F_resist     = self.calculate_resistance_forces(velocity)
        wheel_torque = F_resist * self.wheel_radius
        return transmission.wheel_to_motor_torque(wheel_torque)

    # ─────────────────────────────────────────────────────────────────────────
    # Transferência de carga dinâmica
    # ─────────────────────────────────────────────────────────────────────────

    def calculate_downforce(self, velocity: float) -> float:
        """Aerodynamic downforce [N] — TOTAL (front + rear), positive = down.

        F_down = ½ρ·(Cl_front·A_front + Cl_rear·A_rear)·v²

        Returns 0 when has_wing=False (car without wings).
        """
        if not self.has_wing:
            return 0.0
        v = abs(float(velocity))
        F_front = self.half_rho * self.lift_coeff_front * self.area_front * v ** 2
        F_rear  = self.half_rho * self.lift_coeff_rear  * self.area_rear  * v ** 2
        return F_front + F_rear

    def calculate_load_transfer(self, a: float, velocity: float = 0.0) -> float:
        """Dynamic normal force on the driven (rear) axle, PER WHEEL [N].

        Aggregates three contributions:
          • Static weight on rear axle (from CG location)
          • Longitudinal load transfer (m·a·h/L) during acceleration
          • Aerodynamic downforce on rear (½ρCl·A·v² · balance_rear)

        Parameters
        ----------
        a : float
            Longitudinal vehicle acceleration [m/s²].
        velocity : float, optional
            Vehicle longitudinal speed [m/s]. Default 0 (static).

        Returns
        -------
        float
            Normal force per rear wheel [N] (axle total / 2).
        """
        if self.L is None or self.h is None or self.dist_cg is None:
            base = self.mass * self.g * 0.5
            F_down_rear = (self.half_rho * self.lift_coeff_rear * self.area_rear
                           * abs(velocity) ** 2) if self.has_wing else 0.0
            return base + F_down_rear * 0.5

        Fz_static     = (self.mass * self.g * self.dist_cg) / self.L
        load_transfer = (self.mass * a * self.h) / self.L
        F_down_rear   = ((self.half_rho * self.lift_coeff_rear * self.area_rear
                          * abs(velocity) ** 2) if self.has_wing else 0.0)
        return (Fz_static + load_transfer + F_down_rear) / 2.0

    # ─────────────────────────────────────────────────────────────────────────
    # Inércia refletida — FÓRMULA CORRETA
    # ─────────────────────────────────────────────────────────────────────────

    def calculate_reflected_inertia(self, transmission) -> float:
        """Total drivetrain inertia reflected to the motor shaft [kg·m²].

        Computes the equivalent rotational inertia seen by the motor,
        following the kinetic-energy equivalence method:

            J_reflected = J_motor_side
                        + J_slow_side / N²
                        + m_vehicle × (r/N)²

        where:
          J_motor_side   = sprocket_inertia  (rotates at motor speed)
          J_slow_side    = J_wheels + J_axles + J_differential
                           (rotate at wheel speed ω_wheel = ω_motor / N)
          m × (r/N)²     = translational mass reflected to motor shaft
                           (r/N is wheel radius as seen from the motor)

        Note on the two-body model (used in Simulation.py)
        ---------------------------------------------------
        When motor and vehicle are modelled as separate dynamic bodies
        (Newton eq. for each), the translational term m×(r/N)² belongs
        to the vehicle equation, NOT to the motor equation.  In that
        case use only:
            J_eff_motor = J_motor (from Motor.jm) + J_motor_side (sprocket)
                        + J_slow_side / N²

        This value is available as ``j_rot_reflected`` attribute after
        calling this method.

        Parameters
        ----------
        transmission : Transmission
            Must have attributes: final_drive_ratio, efficiency,
            axle_inertia, diff_inertia.

        Returns
        -------
        float
            Total inertia reflected to motor shaft [kg·m²], including
            the translational mass term (for rigid single-body models).
        """
        N = transmission.final_drive_ratio
        r = self.wheel_radius

        # ── Lado lento (velocidade da roda) ───────────────────────────────
        # Fator de forma k para conjunto pneu+aro:
        #   k = 0.80 para rodas com raio interno moderado (fórmula student típica)
        #   k → 1.0  para anel fino; k → 0.5 para disco sólido
        k = 0.80
        j_wheel_single     = k * self.wheel_mass * r**2
        j_all_wheels       = 4 * j_wheel_single          # 4 rodas (todas giram)

        # Componentes da transmissão no lado lento (2 semi-eixos + diferencial)
        j_axles            = self.n_driven_wheels * transmission.axle_inertia
        j_differential     = transmission.diff_inertia

        # Coroa (lado do diferencial — baixa velocidade) refletida por N²
        j_coroa            = getattr(transmission, 'coroa_inertia', 0.0)
        j_slow_total       = j_all_wheels + j_axles + j_differential + j_coroa

        # ── Reflexão: lado lento → lado do motor ──────────────────────────
        # Pelo princípio de conservação de energia cinética:
        #   ½ J_slow × ω_roda²  =  ½ (J_slow/N²) × ω_motor²
        j_slow_reflected   = j_slow_total / (N ** 2)

        # ── Lado rápido (velocidade do motor) ────────────────────────────
        # Sprocket/coroa — parâmetro da transmissão, soma DIRETAMENTE (sem reflexão)
        j_motor_side       = getattr(transmission, 'sprocket_inertia', 0.0)

        # ── Massa translacional do veículo refletida ──────────────────────
        # A roda "vê" ω_roda = v / r, e o motor vê ω_motor = v / (r/N)
        # Energia cinética translacional: ½ m v² = ½ m (r/N)² ω_motor²
        # Portanto a inércia equivalente no motor é:
        #   J_transl = m × (r/N)²     ← CORRETO
        # (erro clássico: usar m × r², que é N² vezes maior)
        j_translational    = self.mass * (r / N) ** 2

        # ── Total refletido (modelo de corpo único) ───────────────────────
        j_total_reflected  = j_motor_side + j_slow_reflected + j_translational

        # Disponibiliza componentes separados para Simulation.py
        # (dois corpos: motor + veículo)
        self.j_rot_reflected  = j_motor_side + j_slow_reflected  # sem massa transl.
        self.j_translational  = j_translational

        return j_total_reflected


  