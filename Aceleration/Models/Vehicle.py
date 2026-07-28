# -*- coding: utf-8 -*-
"""
Vehicle.py — Modelo longitudinal do veículo Fórmula SAE
========================================================

Modela o carro como um corpo em translação pura, respondendo a três
perguntas para a Simulation a cada passo:

  1. Quais forças RESISTEM ao movimento?  (calculate_resistance_forces:
     arrasto aerodinâmico + rolamento dependente de velocidade + rampa)
  2. Quanta carga vertical Fz há na roda motriz?  (calculate_load_transfer:
     peso estático + transferência longitudinal m·a·h/L + downforce traseiro)
     — é o Fz que alimenta a Magic Formula do pneu.
  3. Qual a inércia do trem de força vista pelo eixo do motor?
     (calculate_reflected_inertia — usada para montar J_eff da equação
     rotacional em Simulation.__init__)

Modelo de DOIS CORPOS (contexto essencial para a inércia refletida)
-------------------------------------------------------------------
Motor e veículo são corpos dinâmicos SEPARADOS na Simulation, acoplados
apenas pelo escorregamento do pneu:

    J_eff × dω/dt = T_em − T_carga − kf·ω − T_ferro     (rotacional, motor)
    m     × dv/dt = F_tração − F_resistência             (translacional, carro)

Para a equação do MOTOR, J_eff deve conter SOMENTE inércias rotacionais:
    J_eff = J_rotor + J_pinhão + (J_rodas + J_semi-eixos + J_diff + J_coroa)/N²

A massa translacional refletida m·(r/N)² pertence à equação do VEÍCULO e
NÃO pode entrar em J_eff (seria dupla contagem). Porém,
calculate_reflected_inertia() RETORNA o total incluindo m·(r/N)² — para
compatibilidade com modelos de corpo único (acoplamento rígido). A
Simulation compensa subtraindo o termo translacional, e este módulo também
disponibiliza as parcelas separadas nos atributos ``j_rot_reflected`` e
``j_translational`` após cada chamada.

Histórico de correções
-----------------------
v2 (este):
  Bug latente corrigido em calculate_reflected_inertia():
    ANTES: j_translation = mass × r² e depois dividia tudo por N²
           → resultava em mass × (r/N)²  ← correto por acidente
           → mas a função retornava valor errado se chamada isoladamente
    DEPOIS: fórmula diretamente correta j_translation = mass × (r/N)².
  Adicionado: sprocket_inertia (pinhão no lado do motor) somado DIRETO,
  sem reflexão por N², pois gira em velocidade de motor.
"""

import numpy as np


class Vehicle:
    """Modelo do veículo para simulação dinâmica longitudinal.

    Modela arrasto aerodinâmico, resistência ao rolamento dependente de
    velocidade, inclinação de pista, transferência de carga dinâmica,
    downforce com asas dianteira/traseira separadas e a inércia do trem
    de força refletida ao eixo do motor.

    Parameters
    ----------
    mass : float
        Massa total do veículo (piloto incluso) [kg].
    wheel_radius : float
        Raio efetivo de rolamento do pneu [m].
    wheel_mass : float
        Massa de UM conjunto roda+pneu [kg] (usada na inércia das rodas).
    drag_coeff : float
        Coeficiente de arrasto aerodinâmico (Cd).
    frontal_area : float
        Área frontal do veículo [m²].
    rolling_resistance : float
        Coeficiente base de resistência ao rolamento Cr0 — o efetivo
        cresce com a velocidade: Cr(v) = Cr0·(1 + v/v_ref).
    road_grade : float, optional
        Ângulo de inclinação da pista [rad]. Default 0 (pista plana).
    environment_density : float, optional
        Densidade do ar [kg/m³]. Default 1.225 (nível do mar, 15 °C).
    L : float, optional
        Entre-eixos [m]. Se None (junto com h e dist_cg), a transferência
        de carga usa um fallback de 50 % do peso por eixo.
    h : float, optional
        Altura do centro de gravidade [m] — controla a transferência
        longitudinal m·a·h/L.
    dist_cg : float, optional
        Distância longitudinal do eixo DIANTEIRO (não-motriz) ao CG [m].
        Usada como ``a`` em Fz_traseiro = m·g·a/L (momento em torno do
        eixo dianteiro → carga no traseiro). Para um FSAE RWD típico com
        dist_cg=0.6 m e L=1.5 m: 40 % do peso no eixo traseiro.
    n_driven_wheels : int, optional
        Número de rodas motrizes. Default 2 (tração traseira).
    rolling_resistance_v_ref : float, optional
        Velocidade de referência [m/s] da lei Cr(v). Default 150
        (→ +13 % de rolamento a 72 km/h).
    has_wing : bool, optional
        Se False, todo o downforce é zerado (carro sem asas).
    lift_coeff_front, area_front : float, optional
        Cl e área de referência [m²] da asa DIANTEIRA. O downforce
        dianteiro carrega o eixo dianteiro (não-motriz) — aumenta arrasto
        induzido implícito no Cd, mas não a tração.
    lift_coeff_rear, area_rear : float, optional
        Cl e área de referência [m²] da asa TRASEIRA. O downforce traseiro
        soma-se ao Fz da roda motriz → mais tração em alta velocidade.

    Attributes
    ----------
    lift_coeff, lift_area, downforce_balance_rear : float
        Agregados de retrocompatibilidade calculados a partir das asas
        individuais (Cl total, área total e fração do downforce no
        traseiro, respectivamente).
    j_rot_reflected, j_translational : float
        Preenchidos por calculate_reflected_inertia(): parcela rotacional
        (para o modelo de dois corpos) e parcela translacional m·(r/N)².
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
        # variação da resistência ao rolamento. Com v_ref=150 m/s o efeito
        # é suave: +13 % de rolamento a 72 km/h.
        self.rolling_resistance_v_ref = float(rolling_resistance_v_ref)

        # Aerodinâmica vertical — asas dianteira e traseira separadas.
        # Só o downforce TRASEIRO entra no Fz da roda motriz (RWD).
        self.has_wing         = bool(has_wing)
        self.lift_coeff_front = float(lift_coeff_front)
        self.area_front       = float(area_front)
        self.lift_coeff_rear  = float(lift_coeff_rear)
        self.area_rear        = float(area_rear)

        # Atributos agregados para retrocompatibilidade com código que
        # usava asa única (calculados a partir dos componentes).
        self.lift_coeff   = self.lift_coeff_front + self.lift_coeff_rear
        self.lift_area    = self.area_front + self.area_rear
        _df_front = self.lift_coeff_front * self.area_front
        _df_rear  = self.lift_coeff_rear  * self.area_rear
        _df_total = _df_front + _df_rear
        self.downforce_balance_rear = _df_rear / _df_total if _df_total > 0 else 0.5

        self.g        = 9.81
        # ½ρ pré-computado — aparece em todas as fórmulas aerodinâmicas
        self.half_rho = 0.5 * self.environment_density

    # ─────────────────────────────────────────────────────────────────────────
    # Forças de resistência
    # ─────────────────────────────────────────────────────────────────────────

    def calculate_resistance_forces(self, velocity: float) -> float:
        """Força de resistência longitudinal total [N].

        Soma três parcelas:
          • Arrasto aerodinâmico:  F_aero  = ½ρ·Cd·A·v²
          • Rolamento:             F_roll  = Cr(v)·m·g·cos(θ), com
                                   Cr(v) = Cr0·(1 + v/v_ref)
                                   (+6.7 % @ 10 m/s, +13 % @ 20 m/s,
                                    +20 % @ 30 m/s para v_ref=150)
          • Rampa:                 F_grade = m·g·sin(θ)

        O valor retornado é subtraído da força de tração na equação
        translacional do veículo (m·dv/dt = F_tração − F_resist).
        """
        v      = abs(float(velocity))
        Cr_v   = self.rolling_resistance * (1.0 + v / self.rolling_resistance_v_ref)
        Faero  = self.half_rho * self.drag_coeff * self.frontal_area * v ** 2
        Froll  = Cr_v * self.mass * self.g * np.cos(self.road_grade)
        Fgrade = self.mass * self.g * np.sin(self.road_grade)
        return Faero + Froll + Fgrade

    def calculate_load_torque(self, velocity: float, transmission) -> float:
        """Torque equivalente no MOTOR para vencer a resistência do veículo [N·m].

        Converte a força resistiva em torque de roda (F·r) e reflete ao
        eixo do motor dividindo por N·η. Função de conveniência para
        análises isoladas — a Simulation calcula o torque de carga pelo
        caminho completo do pneu (com saturação de tração).
        """
        F_resist     = self.calculate_resistance_forces(velocity)
        wheel_torque = F_resist * self.wheel_radius
        return transmission.wheel_to_motor_torque(wheel_torque)

    # ─────────────────────────────────────────────────────────────────────────
    # Transferência de carga dinâmica
    # ─────────────────────────────────────────────────────────────────────────

    def calculate_downforce(self, velocity: float) -> float:
        """Downforce aerodinâmico TOTAL [N] (dianteira + traseira), positivo = para baixo.

        F_down = ½ρ·(Cl_front·A_front + Cl_rear·A_rear)·v²

        Retorna 0 quando has_wing=False (carro sem asas). Usada para
        diagnóstico/exibição; a Simulation usa apenas a parcela traseira,
        via calculate_load_transfer.
        """
        if not self.has_wing:
            return 0.0
        v = abs(float(velocity))
        F_front = self.half_rho * self.lift_coeff_front * self.area_front * v ** 2
        F_rear  = self.half_rho * self.lift_coeff_rear  * self.area_rear  * v ** 2
        return F_front + F_rear

    def calculate_load_transfer(self, a: float, velocity: float = 0.0) -> float:
        """Força normal dinâmica no eixo motriz (traseiro), POR RODA [N].

        Soma três contribuições e divide o total do eixo por 2:

          • Peso estático no traseiro:  m·g·dist_cg/L
            (momento em torno do eixo dianteiro)
          • Transferência longitudinal: m·a·h/L
            (acelerar "senta" o carro na traseira → mais Fz → mais tração;
            é por isso que o Fz cresce durante a largada)
          • Downforce traseiro:         ½ρ·Cl_rear·A_rear·v²
            (cresce com v² → recupera tração em alta velocidade)

        O resultado alimenta a Magic Formula do pneu (que trabalha POR
        pneu) — a Simulation multiplica a força do pneu por
        n_driven_wheels para obter o teto de tração do eixo.

        Parameters
        ----------
        a : float
            Aceleração longitudinal do veículo [m/s²]. Como Fz depende de
            ``a`` e ``a`` depende de Fz, a Simulation resolve esse loop
            algébrico com 2 iterações de ponto fixo.
        velocity : float, optional
            Velocidade longitudinal [m/s]. Default 0 (estático).

        Returns
        -------
        float
            Força normal por roda traseira [N] (total do eixo / 2).
        """
        # Fallback sem geometria: 50 % do peso + downforce traseiro,
        # dividido por roda (sem transferência longitudinal).
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
    # Inércia refletida ao eixo do motor
    # ─────────────────────────────────────────────────────────────────────────

    def calculate_reflected_inertia(self, transmission) -> float:
        """Inércia total do trem de força refletida ao eixo do motor [kg·m²].

        Calcula a inércia rotacional equivalente vista pelo motor pelo
        método de equivalência de energia cinética:

            J_reflected = J_motor_side
                        + J_slow_side / N²
                        + m_vehicle × (r/N)²

        onde:
          J_motor_side   = sprocket_inertia (pinhão — gira em vel. de motor,
                           soma direto, SEM reflexão)
          J_slow_side    = J_rodas + J_semi-eixos + J_diferencial + J_coroa
                           (giram em vel. de roda ω_roda = ω_motor/N →
                           refletidas por 1/N²)
          m × (r/N)²     = massa translacional refletida ao eixo do motor
                           (r/N é o raio da roda "visto" pelo motor)

        Nota sobre o modelo de DOIS CORPOS (usado em Simulation.py)
        -----------------------------------------------------------
        Quando motor e veículo são corpos dinâmicos separados (uma equação
        de Newton para cada), o termo translacional m×(r/N)² pertence à
        equação do VEÍCULO, NÃO à do motor. Nesse caso use apenas:

            J_eff_motor = J_rotor (Motor.jm) + J_motor_side (pinhão)
                        + J_slow_side / N²

        Essa parcela fica disponível no atributo ``j_rot_reflected`` após
        chamar este método (e ``j_translational`` guarda m×(r/N)²).
        A Simulation usa o retorno completo e subtrai o termo
        translacional — as duas rotas dão o mesmo J_eff.

        Parameters
        ----------
        transmission : Transmission
            Deve ter os atributos: final_drive_ratio, axle_inertia,
            diff_inertia (sprocket_inertia e coroa_inertia são lidos com
            getattr e default 0).

        Returns
        -------
        float
            Inércia total refletida ao eixo do motor [kg·m²], INCLUINDO o
            termo translacional (para modelos rígidos de corpo único).
        """
        N = transmission.final_drive_ratio
        r = self.wheel_radius

        # ── Lado lento (velocidade da roda) ───────────────────────────────
        # Fator de forma k para o conjunto pneu+aro (J = k·m·r²):
        #   k = 0.80 para rodas com raio interno moderado (Formula Student típica)
        #   k → 1.0  para anel fino; k → 0.5 para disco sólido
        k = 0.80
        j_wheel_single     = k * self.wheel_mass * r**2
        j_all_wheels       = 4 * j_wheel_single          # 4 rodas (todas giram)

        # Componentes da transmissão no lado lento (2 semi-eixos + diferencial)
        j_axles            = self.n_driven_wheels * transmission.axle_inertia
        j_differential     = transmission.diff_inertia

        # Coroa (lado do diferencial — baixa velocidade) — também refletida por N²
        j_coroa            = getattr(transmission, 'coroa_inertia', 0.0)
        j_slow_total       = j_all_wheels + j_axles + j_differential + j_coroa

        # ── Reflexão: lado lento → lado do motor ──────────────────────────
        # Pelo princípio de conservação de energia cinética:
        #   ½ J_slow × ω_roda²  =  ½ (J_slow/N²) × ω_motor²
        j_slow_reflected   = j_slow_total / (N ** 2)

        # ── Lado rápido (velocidade do motor) ────────────────────────────
        # Pinhão/sprocket — parâmetro da transmissão, soma DIRETAMENTE
        # (já gira na velocidade do motor, não há reflexão)
        j_motor_side       = getattr(transmission, 'sprocket_inertia', 0.0)

        # ── Massa translacional do veículo refletida ──────────────────────
        # A roda "vê" ω_roda = v/r, e o motor vê ω_motor = v/(r/N).
        # Energia cinética translacional: ½ m v² = ½ m (r/N)² ω_motor²
        # Portanto a inércia equivalente no eixo do motor é:
        #   J_transl = m × (r/N)²     ← CORRETO
        # (erro clássico: usar m × r², que é N² vezes maior)
        j_translational    = self.mass * (r / N) ** 2

        # ── Total refletido (modelo de corpo único) ───────────────────────
        j_total_reflected  = j_motor_side + j_slow_reflected + j_translational

        # Disponibiliza as parcelas separadas para Simulation.py
        # (modelo de dois corpos: motor + veículo)
        self.j_rot_reflected  = j_motor_side + j_slow_reflected  # sem massa transl.
        self.j_translational  = j_translational

        return j_total_reflected
