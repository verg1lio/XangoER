# -*- coding: utf-8 -*-
"""
Simulation.py — Motor de simulação do DigitalTwin Aceleration (Xangô E-Racing)
===============================================================================

Este módulo é o CORAÇÃO do digital twin: ele conecta todos os modelos
(Motor, Vehicle, Transmission, BatteryPack, Tire, Pedal, PIDController),
possui TODO o estado dinâmico e o avança no tempo com um integrador RK4 de
passo fixo a 10 kHz (hp = 1e-4 s).

Visão geral da arquitetura
--------------------------
Os modelos em ``Models/`` são contêineres de parâmetros sem estado próprio.
A Simulation lê os parâmetros deles no __init__ e monta:

  • o VETOR DE ESTADO de 11 variáveis:
        [isd, isq, iso,   → correntes dq0 do motor [A]
         wm, theta_m,     → velocidade [rad/s] e ângulo mecânico [rad]
         temp,            → temperatura do enrolamento [°C]
         vv, vp,          → velocidade [m/s] e posição [m] do veículo
         soc, Iast, tacc] → estados da bateria
  • a ODE física pura (_physics_ode) — sem controladores, sem side-effects;
  • o loop de controle FOC executado UMA vez por passo, ANTES do RK4, com
    as saídas (vd, vq, I_bat) congeladas durante as 4 avaliações do RK4.

Modelo de dois corpos (motor ↔ veículo)
---------------------------------------
Motor e veículo são corpos dinâmicos separados, acoplados só pelo pneu:

    J_eff × dω/dt = T_em − T_carga − kf·ω − T_ferro     (rotacional)
    m     × dv/dt = F_tração − F_resistência             (translacional)

onde T_carga = F_tração_limitada · r / (N·η) é a reação da força de tração
no eixo do motor, e F_tração é SATURADA pelo máximo que o pneu transmite
(Magic Formula) — é daí que surge o comportamento de wheelspin.

J_eff contém APENAS inércias rotacionais (rotor + pinhão + lado lento/N²);
a massa translacional fica na equação do veículo (sem dupla contagem).

Histórico de correções
-----------------------
v1  (original)  — PIDs dentro da ODE + cl errado + double-append + post-proc PIDs
v2  (anterior)  — PIDs separados da ODE via RK4 fixo, cl corrigido
v3  (este)      — Corrente da bateria calculada corretamente (BUG RAIZ)

BUG RAIZ — corrente incorreta para o modelo de bateria (NÃO REGREDIR)
----------------------------------------------------------------------
Versões anteriores passavam `abs(isq)` como corrente para calcular_tensao() e
calcular_derivadas(). isq é a corrente no eixo q do referencial dq rotativo —
ela NÃO é a corrente DC do banco de baterias.

Consequência com os parâmetros reais do projeto:
  R_interno_total = n_serie × R_célula = 264 × 0.02 = 5.28 Ω
  isq ≈ 700 A (saturado pelo controlador de velocidade)
  Queda = 5.28 × 700 = 3,696 V  >>  V_nominal = 976 V  →  Vdc < 0

Com Vdc < 0 → Vclamp = max(Vdc, 1.0) = 1 V. A tensão disponível (1 V) não
cancela o back-EMF de desacoplamento no eixo d (≈ ωe × Lq × isq ≈ 473 V):
    d_isd/dt = (vd - rs×isd + ωe×Lq×isq) / Ld
           = (±1 - 0 + 473) / 9.65e-5  →  isd → -∞

Isso provoca em cascata:
  • isd → -1000 A → torque eletromagnético oscila ±600 Nm
  • Perdas cobre = 1.5 × rs × (isd² + isq²) ≈ 15,800 W → temperatura → 4000 °C
  • Tensão da bateria oscila entre -4000 V e +100 V

CORREÇÃO (mantida nesta versão):
  A corrente DC extraída da bateria relaciona-se com a potência elétrica consumida:
      P_ac = potência mecânica / η_drive          [W]
      I_bat = P_ac / V_nominal                    [A] — corrente DC no barramento

  I_bat é calculada UMA VEZ por passo a partir dos estados e comandos congelados,
  e passada como parâmetro fixo tanto para calcular_tensao() quanto para a ODE
  (calcular_derivadas). Isso mantém a consistência entre o cálculo da tensão
  terminal e a evolução do SoC.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from Models.BatteryPack import BatteryPack
from Models import PIDController as PID
from Models.Tire import Tire
from Models.Transmission import Transmission
from Models.Vehicle import Vehicle
from Models.Motor import Motor
from Models.Pedal import Pedal
from typing import Optional, Dict, Any


class Simulation:
    """Digital-twin do powertrain — PMSM + veículo + bateria + pneu.

    Arquitetura de controle (sequência por passo)
    ---------------------------------------------
    Loop RK4 de passo fixo hp. A cada passo k:
      1. Desempacota o estado x[k].
      2. Calcula a corrente DC da bateria por balanço de potência.
      3. Obtém a tensão terminal Vdc da bateria com a corrente correta.
      4. Monta iq_ref (pedal → limitadores) e id_ref (field weakening),
         executa os PIs de corrente (1 vez) → vd, vq congelados.
      5. Integra _physics_ode com vd, vq, I_bat fixos (RK4, 4 avaliações puras).
      6. Loga x[k] e os sinais de controle.

    _physics_ode é função pura: sem estado interno, sem PIDs — isso é o que
    torna o RK4 matematicamente correto (as 4 avaliações veem a mesma
    entrada de controle, como um ZOH real de 10 kHz).

    Controle de torque em malha aberta (decisão de projeto)
    -------------------------------------------------------
    Não há PI de velocidade: iq_ref vem direto do pedal
    (iq_ref = max_current × posição), passando por três protetores:
      • limitador SUAVE de sobre-rotação (smoothstep 1→0 entre 90 % e
        100 % de speed_ref);
      • teto de potência DC (FSAE EV.4.1: 80 kW);
      • derate térmico (rampa linear entre T_alarm e T_max).
    Um PI de velocidade saturaria o integrador durante toda a aceleração
    (erro grande) e oscilaria ao alcançar a referência — o pedal em malha
    aberta reproduz melhor o carro real e evita esse windup.

    Acoplamento motor–veículo (dois corpos dinâmicos separados)
    -----------------------------------------------------------
        J_eff  × dω/dt = Tₑ  − Tcarga − kf·ω
        m      × dv/dt = Ftração − Fresistência

    onde Tcarga = Ftração_limitada × r / (N·η)

    J_eff = J_rotor + J_rodas_refletido + J_transmissão_refletido
    (sem massa translacional do veículo — evita dupla contagem).
    """

    # Passo de integração E de controle [s] — 1e-4 s = 10 kHz. Fixo: passo
    # variável quebraria a premissa de controladores discretos com Ts fixo.
    HP: float = 1e-4

    def __init__(self,
                 motor: Optional[Motor] = None,
                 vehicle: Optional[Vehicle] = None,
                 transmission: Optional[Transmission] = None,
                 battery: Optional[BatteryPack] = None,
                 tire: Optional[Tire] = None,
                 pedal: Optional[Pedal] = None,
                 tmax: float = 60.0,
                 dmax: Optional[float] = 75.0,
                 steps: int = 12000,
                 p_max_dc: Optional[float] = 80e3):

        # Referências aos modelos (todos opcionais — subsistemas ausentes
        # são simplesmente pulados na ODE e no vetor de estado)
        self.motor        = motor
        self.vehicle      = vehicle
        self.transmission = transmission
        self.battery      = battery
        self.tire         = tire
        self.pedal        = pedal      # se None → throttle = 100% durante toda simulação
        self.tmax         = float(tmax)   # limite de segurança [s]
        self.dmax         = float(dmax) if dmax is not None else None  # critério de parada [m]
        self.hp           = self.HP

        # Limite de potência no barramento DC (FSAE EV.4.1: 80 kW).
        # None → sem limite.  Aplicado como teto de iq_ref no loop de controle:
        #     iq_max = P_max_dc·η / (Kt·ω)
        # → torque constante até a velocidade-base, potência constante acima.
        self.p_max_dc = float(p_max_dc) if p_max_dc is not None else None
        # Eficiência estimada motor+inversor — usada no balanço de potência
        # DC↔mecânica (corrente de bateria e teto de potência).
        self.eta_drive = 0.90

        # ── Parâmetros do motor ──────────────────────────────────────────────
        # Copiados do objeto Motor via getattr (com defaults seguros) para
        # acesso rápido como atributos locais no caminho quente do loop.
        self.p        = getattr(motor, 'p',        1)
        self.ld       = getattr(motor, 'ld',       1.0)
        self.lq       = getattr(motor, 'lq',       1.0)
        # rs0 = referência fria; self.rs é re-derivado a cada passo
        # (Rs(T) = rs0 · [1 + α·(T − T_ref)])  para refletir aquecimento
        # do cobre durante a prova — afeta perdas, queda de tensão e
        # convergência do PI de corrente.
        self.rs0      = getattr(motor, 'rs0', getattr(motor, 'rs', 0.0))
        self.rs       = self.rs0
        self.lambda_m = getattr(motor, 'lambda_m', 0.0)
        self.jm       = getattr(motor, 'jm',       1.0)
        self.kf       = getattr(motor, 'kf',       0.0)

        # Perdas no ferro (Steinmetz concentrado: P_iron = k_h·|f_e| + k_e·f_e²)
        self.k_h_iron = getattr(motor, 'k_h_iron', 0.0)
        self.k_e_iron = getattr(motor, 'k_e_iron', 0.0)
        # Coeficientes para Rs(T)
        self.alpha_cu = getattr(motor, 'alpha_cu', 0.0)
        self.T_ref    = getattr(motor, 'T_ref',   25.0)
        # Derate térmico do iq_pedal (rampa linear entre T_alarm e T_max)
        self.T_alarm  = getattr(motor, 'T_alarm', 1e6)   # default = sem derate
        self.T_max    = getattr(motor, 'T_max',   1e6)
        # Piso de velocidade para converter P_iron → torque sem singularidade
        self.wm_floor = 1.0   # [rad/s]

        # Constante de torque (convenção amplitude-invariante) e limites
        self.Kt          = 1.5 * self.p * self.lambda_m       # [N·m/A]
        self.max_current = getattr(motor, 'max_current', np.inf)
        self.Vdc_default = getattr(motor, 'Vdc',         600.0)
        self.mod_index   = getattr(motor, 'valor_mu',    1.0)
        self.speed_ref   = getattr(motor, 'speed_ref',   700.23)  # [rad/s]

        # ── Inércia efetiva no eixo do motor ─────────────────────────────────
        # J_total de calculate_reflected_inertia inclui J_translação = m(r/N)².
        # Esse termo pertence à equação do veículo (corpo separado) → removemos
        # para evitar dupla contagem.  O max(..., jm) garante J_eff ≥ J_rotor
        # mesmo com parâmetros degenerados.
        if vehicle is not None and transmission is not None:
            j_all         = vehicle.calculate_reflected_inertia(transmission)
            r, N          = vehicle.wheel_radius, transmission.final_drive_ratio
            j_translation = vehicle.mass * (r / N) ** 2
            self.j_eff    = max(self.jm + j_all - j_translation, self.jm)
        else:
            self.j_eff = self.jm

        # Inversos pré-computados — troca divisão por multiplicação no loop
        # de 10 kHz (a ODE é chamada 4× por passo).
        self.inv_ld  = 1.0 / self.ld  if self.ld  > 0 else 0.0
        self.inv_lq  = 1.0 / self.lq  if self.lq  > 0 else 0.0
        self.inv_jm  = 1.0 / self.j_eff if self.j_eff > 0 else 0.0
        # Indutância de sequência zero (não fornecida pelo datasheet):
        # aproximada como 10 % da média de Ld/Lq — iso decai e não afeta torque.
        self.L0      = max(0.1 * (self.ld + self.lq) / 2.0, 1e-6)

        # Tensão nominal da bateria (denominador do balanço de potência)
        if battery is not None:
            self.Vdc_nominal = battery.calcular_tensao_nominal()
        else:
            self.Vdc_nominal = self.Vdc_default

        # Piso de tensão: protege contra Vdc negativo por parâmetros extremos
        # (o colapso de barramento descrito no BUG RAIZ do cabeçalho)
        self.Vdc_floor = 0.3 * self.Vdc_nominal

        # Filtro passa-baixa em wm para o limitador suave de sobre-rotação.
        # τ = 5 ms → corner ≈ 31.8 Hz, bem acima do BW do loop externo de
        # velocidade (~5 Hz) e bem abaixo da frequência de ripple numérico.
        # Atenua a realimentação ruído-do-torque → wm → iq_ref que dispara
        # ciclo-limite na transição do limitador.
        self.tau_wm = 5e-3

        # ── Térmico ──────────────────────────────────────────────────────────
        # Modelo de massa concentrada: dT/dt = (P_ger − P_dissip)/(m·C)
        m_th = getattr(motor, 'm', 22.0)
        C_th = getattr(motor, 'C', 385.0)   # J/(kg·K) — calor específico do cobre
        self.inv_mC = 1.0 / (m_th * C_th) if m_th * C_th > 0 else 0.0

        # Resfriamento convectivo: P_cool = h_cool · A_cool · (T − T_ambient)
        # h=0 → adiabático (comportamento anterior); h>0 → dissipa calor
        self.h_cool    = getattr(motor, 'h_cool',    0.0)
        self.A_cool    = getattr(motor, 'A_cool',    0.0)
        self.T_ambient = getattr(motor, 'T_ambient', 25.0)

        # ── Controladores PI de corrente ─────────────────────────────────────
        #
        # Sintonia por cancelamento de polo (bandwidth alvo 500 Hz):
        #   a planta elétrica de cada eixo é 1ª ordem com τe = L/Rs;
        #   escolhendo kp = ωc·L e ki = ωc·Rs, o zero do PI (ki/kp = Rs/L)
        #   cancela o polo da planta e a malha fechada vira 1ª ordem com
        #   banda ωc.  Os ganhos são calculados dos parâmetros REAIS do
        #   motor, garantindo tuning consistente quando o usuário troca o
        #   motor no dashboard.
        wc_curr = 2.0 * np.pi * 500.0     # [rad/s]
        kp_curr = wc_curr * self.ld
        ki_curr = wc_curr * self.rs

        self.id_ctrl = PID.Controller(
            kp=kp_curr, ki=ki_curr, kd=0.0,
            limit=self.Vdc_nominal, Ts=self.hp)
        self.iq_ctrl = PID.Controller(
            kp=kp_curr, ki=ki_curr, kd=0.0,
            limit=self.Vdc_nominal, Ts=self.hp)

        # Não há PI de velocidade: o controle de rotação é open-loop via
        # pedal + limitador suave smoothstep em wm (ver loop principal).
        # Mantemos a sintonia que SERIA usada por um PI de velocidade só para
        # diagnóstico no print abaixo (ωcv = 2π·5 rad/s, ζ ≈ 0.9).
        wcv      = 2.0 * np.pi * 5.0
        kp_spd_d = wcv * self.j_eff / self.Kt if self.Kt > 0 else 0.0
        ki_spd_d = 0.30 * kp_spd_d
        zeta_d   = (kp_spd_d / (2.0 * np.sqrt(max(ki_spd_d * self.j_eff, 1e-9))))

        print(f"[Simulation] Kt={self.Kt:.4f} N·m/A  "
              f"J_eff={self.j_eff:.4f} kg·m²  "
              f"Vdc_nom={self.Vdc_nominal:.1f} V  "
              f"PI corrente kp={kp_curr:.3f} ki={ki_curr:.3f}  "
              f"(ref. PI velocidade kp={kp_spd_d:.3f} ki={ki_spd_d:.3f} ζ≈{zeta_d:.2f})")

        self._reset_ic()
        self._init_storage()

    # ─────────────────────────────────────────────────────────────────────────
    # Condições iniciais e armazenamento
    # ─────────────────────────────────────────────────────────────────────────

    def _reset_ic(self):
        """Define as condições iniciais de todos os estados.

        Carro parado, motor parado, correntes nulas, enrolamento a 25 °C.
        Os estados da bateria são herdados do objeto BatteryPack (permite
        iniciar com SoC < 1 para estudos de fim de prova).
        """
        self.isd0 = self.isq0 = self.iso0 = 0.0
        self.wm0  = self.theta_m0 = 0.0
        self.temp0 = 25.0
        self.vv0 = self.vp0 = 0.0
        if self.battery is not None:
            self.soc0  = getattr(self.battery, 'soc',             1.0)
            self.Iast0 = getattr(self.battery, 'Iast',            0.0)
            self.tacc0 = getattr(self.battery, 'tempo_acumulado', 0.0)
        else:
            self.soc0 = 1.0; self.Iast0 = self.tacc0 = 0.0
        # Estado auxiliar do filtro de wm (não faz parte do vetor de estado x)
        self.wm_filt = 0.0

    # Lista canônica das variáveis temporais armazenadas pelo loop.
    # Centralizar evita drift entre _init_storage, o slice final e o dict
    # de retorno em simulate().  ATENÇÃO: os nomes internos self.* diferem
    # das chaves do dict de retorno (ex.: 'corrented' → chave 'isd') — o
    # dashboard usa os self.*; consumidores externos devem usar o dict.
    _LOG_FIELDS = (
        'tempo', 'corrented', 'correnteq',
        'corrente1', 'corrente2', 'corrente3',
        'tensaosd', 'tensaosq', 'tensao1', 'tensao2', 'tensao3',
        'fluxosd', 'fluxosq', 'conjugado', 'velocidade',
        'conjcarga', 'torque_mecanico', 'temperatura',
        'vd_control', 'vq_control', 'vd_real', 'vq_real',
        'speed_error', 'iq_ref_trace', 'id_ref_trace',
        'vehicle_velocity', 'vehicle_position', 'vehicle_acceleration',
        'wheel_torque', 'tractive_force_hist', 'resistive_force_hist',
        'slip_ratio_hist', 'longitudinal_force_hist', 'fz_hist',
        'soc_hist', 'Iast_hist', 'tempo_acumulado_hist',
        'battery_voltage_hist', 'battery_current_hist',
        'pedal_position_hist', 'effective_speed_ref_hist',
    )

    def _init_storage(self, N: int = 0):
        """Pré-aloca arrays numpy de tamanho ``N`` para o log do loop.

        Pré-alocação evita o overhead de ``list.append`` em ~40 variáveis ×
        ~60k passos. ``N=0`` cria arrays vazios (uso pelo ``__init__``
        antes de ``simulate()`` conhecer a duração efetiva).
        """
        for name in self._LOG_FIELDS:
            setattr(self, name, np.zeros(N))
        self.k_log = 0     # próximo índice livre para escrita

    # ─────────────────────────────────────────────────────────────────────────
    # Vetor de estado
    # ─────────────────────────────────────────────────────────────────────────

    def _build_x0(self) -> np.ndarray:
        """Monta o vetor de estado inicial x0.

        O vetor é montado por blocos condicionais — subsistemas ausentes
        (sem veículo ou sem bateria) simplesmente não entram no vetor:
            [isd, isq, iso, wm, theta_m, temp]           sempre (6)
            + [vv, vp]            se veículo+transmissão  (8)
            + [soc, Iast, tacc]   se bateria             (11)
        """
        x = [self.isd0, self.isq0, self.iso0,
             self.wm0, self.theta_m0, self.temp0]
        if self.vehicle and self.transmission:
            x += [self.vv0, self.vp0]
        if self.battery:
            x += [self.soc0, self.Iast0, self.tacc0]
        return np.array(x, dtype=float)

    def _unpack(self, x: np.ndarray):
        """Desempacota o vetor de estado na mesma ordem de _build_x0.

        Devolve SEMPRE a tupla completa de 11 valores — estados de
        subsistemas ausentes recebem defaults neutros (vv=vp=0, soc=1).
        Isso mantém a ODE e o logger independentes da composição do vetor.
        """
        isd, isq, iso, wm, theta_m, temp = x[0], x[1], x[2], x[3], x[4], x[5]
        i = 6
        if self.vehicle and self.transmission:
            vv, vp = x[i], x[i+1]; i += 2
        else:
            vv = vp = 0.0
        if self.battery:
            soc, Iast, tacc = x[i], x[i+1], x[i+2]
        else:
            soc = 1.0; Iast = tacc = 0.0
        return isd, isq, iso, wm, theta_m, temp, vv, vp, soc, Iast, tacc

    # ─────────────────────────────────────────────────────────────────────────
    # ODE física pura — sem PIDs, sem side-effects
    # vd, vq, I_bat são entradas CONGELADAS para o passo atual
    # ─────────────────────────────────────────────────────────────────────────

    def _physics_ode(self, _t: float, x: np.ndarray,
                     vd: float, vq: float, I_bat: float,
                     telemetry: Optional[dict] = None) -> np.ndarray:
        """Derivadas físicas dx/dt com entradas de controle congeladas.

        Blocos avaliados em sequência:
          1. Rs(T) — resistência do cobre na temperatura DESTE estado;
          2. Dinâmica elétrica dq (equações do PMSM em referencial rotórico);
          3. Acoplamento pneu–veículo (slip → Magic Formula → tração
             saturada → torque de carga + aceleração do carro);
          4. Perdas no ferro (Steinmetz) → torque de freio + calor;
          5. Dinâmica mecânica do eixo do motor;
          6. Dinâmica térmica do enrolamento;
          7. Derivadas da bateria (SoC/Iast) com I_bat congelado.

        I_bat é a corrente DC estimada do banco de baterias, calculada por
        balanço de potência antes desta chamada.  Ela é fixa durante as 4
        avaliações do RK4 para manter consistência física.

        Quando ``telemetry`` é um dict, intermediários (Ft, Fz, slip, Fresist,
        Tcarga, d_vv) são gravados ali para que o logger reaproveite o cálculo
        de pneu — evita executar a Magic Formula 2× por passo. Passar
        ``telemetry`` apenas na 1ª avaliação do RK4 (k1) mantém o log
        consistente com x[k] (estado antes do passo).
        """
        isd, isq, iso, wm, theta_m, temp, vv, vp, soc, Iast, tacc = (
            self._unpack(x))

        # Velocidade elétrica = pares de polos × velocidade mecânica
        we = self.p * wm

        # ── Rs(T): cobre aquece, resistência sobe (~0.4 %/K) ────────────────
        # Avaliado dentro da ODE para que cada estágio do RK4 use o Rs
        # consistente com o ``temp`` do estado correspondente.  Para α=0
        # (compatibilidade), reduz-se a Rs constante = rs0.
        rs_eff = self.rs0 * (1.0 + self.alpha_cu * (temp - self.T_ref))

        # ── Dinâmica elétrica dq ─────────────────────────────────────────────
        # Equações clássicas do PMSM em referencial síncrono:
        #   Ld·d(isd)/dt = vd − Rs·isd + ωe·Lq·isq          (acoplamento +q→d)
        #   Lq·d(isq)/dt = vq − Rs·isq − ωe·(Ld·isd + λm)   (acopl. −d→q + bEMF)
        #   L0·d(iso)/dt = −Rs·iso                          (seq. zero decai)
        d_isd = (vd - rs_eff*isd + we*self.lq*isq) * self.inv_ld
        d_isq = (vq - rs_eff*isq - we*(self.ld*isd + self.lambda_m)) * self.inv_lq
        d_iso = (-rs_eff * iso) / self.L0

        # ── Acoplamento pneu–veículo ─────────────────────────────────────────
        Tcarga = 0.0
        d_vv   = 0.0
        Ft = Fz = Fx_max = Fresist = slip = 0.0

        if self.vehicle and self.transmission and self.tire:
            try:
                N   = self.transmission.final_drive_ratio
                eta = self.transmission.efficiency
                r   = self.vehicle.wheel_radius
                m   = self.vehicle.mass

                # Força que o MOTOR pede no contato pneu-solo
                # (torque em roda / raio) — ainda sem limite de aderência
                Tmotor_roda   = self.transmission.motor_to_wheel_torque(self.Kt * isq)
                Ftração_ideal = Tmotor_roda / r

                # Slip ratio entre a roda (ω_motor/N) e o carro (vv)
                omega_roda = self.transmission.motor_to_wheel_speed(wm)
                slip       = Tire.SlipRatio(omega_roda, r, vv)
                Fresist    = self.vehicle.calculate_resistance_forces(vv)

                # Loop algébrico de 2 iterações (ponto fixo): Fz depende da
                # aceleração (transferência de carga) e a aceleração depende
                # de Fz (teto de tração).  vv é passado a load_transfer para
                # que o termo aerodinâmico ½ρCl·A·v² seja incluído na carga
                # vertical da roda traseira (aumenta Fz_max → mais tração
                # disponível em alta velocidade).
                #
                # Fz é POR RODA (load_transfer divide o eixo por 2) e a Magic
                # Formula retorna Fx de UM pneu — o teto de tração do eixo é
                # n_driven × Fx(Fz_por_roda).  Em linha reta com tração
                # simétrica as duas rodas motrizes contribuem igualmente.
                n_drv = getattr(self.vehicle, 'n_driven_wheels', 2)

                # Iteração 0: Fz com a=0 → tração saturada → aceleração a0
                Fz0 = self.vehicle.calculate_load_transfer(0.0, vv)
                Fx0 = self.tire.Tire_forces(Fz0, slip) * n_drv
                Ft0 = np.sign(Ftração_ideal) * min(abs(Ftração_ideal), abs(Fx0))
                a0  = (Ft0 - Fresist) / m

                # Iteração 1: Fz refinado com a0 → tração final Ft
                Fz     = self.vehicle.calculate_load_transfer(a0, vv)
                Fx_max = self.tire.Tire_forces(Fz, slip) * n_drv
                Ft     = np.sign(Ftração_ideal) * min(abs(Ftração_ideal), abs(Fx_max))

                # Torque de carga no motor = reação da força de tração
                # refletida pela transmissão (r/N com perda η)
                Tcarga = Ft * r / (N * eta)

                # Equação translacional do veículo; a condição vv<0.01
                # impede o carro de "andar para trás" parado no grid
                # quando a resistência excede a tração inicial.
                Fr   = Ft - Fresist
                d_vv = 0.0 if (vv < 0.01 and Fr < 0) else Fr / m

            except Exception as e:
                print(f"[ODE] pneu: {e}")

        # ── Perdas no ferro (Steinmetz) — torque de freio + calor ────────────
        # f_elec = we / (2π).  P_iron dissipa energia do estator/airgap; o
        # equivalente mecânico é um torque opondo movimento:
        #     T_iron = P_iron / max(|wm|, wm_floor)   (sinal = sign(wm))
        # A energia perdida vai integralmente para o aquecimento (somada ao
        # cobre na ODE térmica).  wm_floor evita singularidade em wm → 0.
        f_elec = abs(we) * (1.0 / (2.0 * np.pi))
        P_iron = self.k_h_iron * f_elec + self.k_e_iron * f_elec * f_elec
        wm_safe = max(abs(wm), self.wm_floor)
        T_iron_brake = (P_iron / wm_safe) * np.sign(wm) if wm != 0.0 else 0.0

        # ── Dinâmica mecânica do eixo do motor ───────────────────────────────
        # J_eff·dω/dt = T_em − T_carga − atrito viscoso − freio de ferro
        Ce      = self.Kt * isq
        d_wm    = (Ce - Tcarga - self.kf * wm - T_iron_brake) * self.inv_jm
        d_theta = wm

        # ── Dinâmica térmica ─────────────────────────────────────────────────
        # Balanço: dT/dt = (P_gerado − P_dissipado) / (m·C)
        # P_cu = 1.5·Rs·(isd²+isq²) — perdas Joule trifásicas na convenção
        #        amplitude-invariante (o 1.5 é o mesmo do Kt)
        # P_cool = h_cool · A_cool · (T − T_ambient)  — convecção (natural ou forçada)
        # h_cool = 0 → adiabático (padrão sem resfriamento configurado)
        P_cu   = 1.5 * rs_eff * (isd**2 + isq**2)
        P_cool = self.h_cool * self.A_cool * max(temp - self.T_ambient, 0.0)
        d_temp = (P_cu + P_iron - P_cool) * self.inv_mC

        # ── Bateria — usa I_bat congelado (correto: corrente DC, não isq) ────
        if self.battery:
            try:
                dsoc, dIast, dtacc = self.battery.calcular_derivadas(I_bat)
            except Exception:
                dsoc = dIast = dtacc = 0.0
        else:
            dsoc = dIast = dtacc = 0.0

        # Telemetria opcional para o logger reaproveitar o cálculo de pneu.
        if telemetry is not None:
            telemetry['Ft']      = float(Ft)
            telemetry['Fz']      = float(Fz)
            telemetry['Fx_max']  = float(Fx_max)
            telemetry['slip']    = float(slip)
            telemetry['Fresist'] = float(Fresist)
            telemetry['Tcarga']  = float(Tcarga)
            telemetry['d_vv']    = float(d_vv)

        # ── Monta vetor de derivadas (mesma ordem de _build_x0) ─────────────
        # Nota: a derivada da POSIÇÃO é a velocidade vv (estado), por isso
        # o par [d_vv, vv].
        dx = [d_isd, d_isq, d_iso, d_wm, d_theta, d_temp]
        if self.vehicle and self.transmission:
            dx += [d_vv, vv]
        if self.battery:
            dx += [dsoc, dIast, dtacc]
        return np.array(dx, dtype=float)

    # ─────────────────────────────────────────────────────────────────────────
    # RK4 de passo fixo com entradas congeladas
    # ─────────────────────────────────────────────────────────────────────────

    def _rk4(self, t: float, x: np.ndarray,
              vd: float, vq: float, I_bat: float, dt: float,
              telemetry: Optional[dict] = None) -> np.ndarray:
        """Um passo de Runge-Kutta clássico de 4ª ordem.

        x[k+1] = x[k] + (dt/6)·(k1 + 2·k2 + 2·k3 + k4)

        vd, vq, I_bat ficam FIXOS nas 4 avaliações — fisicamente isso
        representa o zero-order-hold do controlador digital de 10 kHz, e
        matematicamente garante que a ODE seja pura (mesma entrada em
        todos os estágios).

        ``telemetry``, se fornecido, é repassado APENAS ao k1 — captura os
        intermediários físicos avaliados em x[k] (estado pré-passo) para o
        logger evitar recomputar pneu/transmissão.
        """
        k1 = self._physics_ode(t,        x,              vd, vq, I_bat, telemetry)
        k2 = self._physics_ode(t + dt/2, x + (dt/2)*k1, vd, vq, I_bat)
        k3 = self._physics_ode(t + dt/2, x + (dt/2)*k2, vd, vq, I_bat)
        k4 = self._physics_ode(t + dt,   x + dt*k3,     vd, vq, I_bat)
        return x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

    # ─────────────────────────────────────────────────────────────────────────
    # Loop principal de simulação
    # ─────────────────────────────────────────────────────────────────────────

    def simulate(self, t0: float = 0.0,
                 tf: Optional[float] = None) -> Dict[str, Any]:
        """Executa a simulação completa com RK4 de passo fixo.

        Critério de parada (o que ocorrer primeiro):
          • posição do veículo >= dmax  [m]  — critério principal
            (75 m = distância da prova de Acceleration da FSAE)
          • tempo >= tmax               [s]  — limite de segurança

        Sequência por passo (detalhada nos comentários numerados abaixo):
          1. Desempacota o estado atual x[k].
          2. Calcula I_bat via balanço de potência (correção do bug raiz).
          3. Obtém Vdc terminal da bateria com I_bat correto.
          4. Monta as referências de corrente (pedal → iq_ref com
             limitadores; field weakening → id_ref).
          5. Executa os PIs de corrente com desacoplamento (1 chamada cada).
          6. Integra a física com RK4 usando vd, vq, I_bat congelados.
          7. Loga o estado ANTERIOR ao passo + sinais de controle.
          8. Verifica o critério de parada por distância.

        Returns
        -------
        dict
            28 séries temporais com chaves públicas ('t', 'isd', 'soc',
            'vehicle_velocity', ...). Consumidores externos devem usar
            este dict; o dashboard usa os atributos self.* diretamente.
        """
        if tf is None:
            tf = self.tmax

        dt = self.hp
        N  = int((tf - t0) / dt)

        criterio = (f"dmax={self.dmax:.1f} m" if self.dmax is not None
                    else f"tmax={tf:.1f} s")
        print(f"🚀 Simulação: [{t0:.2f}, {tf:.2f}]s  "
              f"dt={dt:.1e}s  N={N}  critério={criterio}  "
              f"J_eff={self.j_eff:.4f} kg·m²  "
              f"Vdc_nom={self.Vdc_nominal:.1f} V")

        # Prepara armazenamento e zera controladores — execuções repetidas
        # de simulate() partem sempre do mesmo estado.
        self._init_storage(N)
        self.id_ctrl.reset()
        self.iq_ctrl.reset()
        self.wm_filt = 0.0   # filtro LP do limitador de velocidade

        x = self._build_x0()
        t = float(t0)

        for _ in range(N):

            # ── 1. Estado atual ─────────────────────────────────────────────
            (isd, isq, iso, wm, theta_m, temp,
             vv, vp, soc, Iast, tacc) = self._unpack(x)

            theta_e = self.p * theta_m     # ângulo elétrico
            we      = self.p * wm          # velocidade elétrica

            # ── 2. Corrente DC por balanço de potência (CORREÇÃO PRINCIPAL) ──
            #
            # P_mec = Ce × wm            [W - potência mecânica no eixo]
            # P_ac  = P_mec / η_drive    [W - potência elétrica consumida]
            # I_bat = P_ac / V_nominal   [A - corrente DC no banco]
            #
            # Usa V_nominal (circuito aberto) no denominador para evitar
            # referência circular (a tensão terminal depende de I_bat).
            # Apenas potência positiva é extraída da bateria (sem regeneração
            # neste modelo simplificado).
            #
            # Na iteração inicial (isd=isq=0), I_bat=0 e Vdc = Voc nominal.
            Ce_now = self.Kt * isq
            P_mec  = Ce_now * wm
            P_ac   = max(P_mec / self.eta_drive, 0.0)
            I_bat  = P_ac / max(self.Vdc_nominal, 1.0)

            # ── 3. Tensão terminal da bateria ────────────────────────────────
            # Vdc define o limite de tensão ±Vlim dos comandos do inversor
            # neste passo — sob carga pesada a bateria "afunda" e o FOC
            # perde margem de tensão (motivando o field weakening abaixo).
            if self.battery is not None and 0.0 < soc <= 1.0:
                try:
                    Vdc = self.battery.calcular_tensao(I_bat, soc, Iast, tacc)
                    # Protege contra colapso de tensão por parâmetros extremos
                    Vdc = max(Vdc, self.Vdc_floor)
                except Exception:
                    Vdc = self.Vdc_nominal
            else:
                Vdc = self.Vdc_nominal
            Vlim = Vdc   # limite de tensão para os comandos do inversor

            # ── 4. Controle de torque direto com limitador suave de sobre-rotação ─
            #
            # Arquitetura:
            #   • Pedal → iq_ref = max_current × pedal_pos  (controle de torque)
            #   • Limitador proporcional suave quando wm > speed_ref (sem integral)
            #
            # Não usa PI de velocidade para gerar iq_ref: evita o problema de
            # windup zero — o anti-windup congela o integrador durante a aceleração
            # inteira (erro grande → saturação), de modo que quando wm alcança
            # speed_ref a saída do PI cai abruptamente para zero, causando oscilação.
            if self.pedal is not None:
                pedal_pos = self.pedal.update(t)
            else:
                pedal_pos = 1.0
            effective_speed_ref = self.speed_ref

            # Derate térmico: reduz iq_max linearmente entre T_alarm e T_max.
            # Acima de T_max → corrente zero (proteção de isolamento).  Para
            # T_max = T_alarm = 1e6 (defaults) o fator é sempre 1 e nada muda.
            if temp >= self.T_max:
                thermal_derate = 0.0
            elif temp > self.T_alarm:
                thermal_derate = (self.T_max - temp) / max(self.T_max - self.T_alarm, 1e-6)
            else:
                thermal_derate = 1.0

            # Comando de corrente do pedal — saturado pelo derate térmico
            iq_pedal = self.max_current * pedal_pos * thermal_derate

            # Filtro LP de 1ª ordem em wm (τ=tau_wm) — usado APENAS pelo
            # limitador suave e pelo teto de potência. Atenua o ripple
            # numérico de torque que realimentava em iq_ref e disparava
            # ciclo-limite no joelho da rampa. Não filtra wm na dinâmica
            # do motor (que continua exata).
            _alpha_wm = dt / (self.tau_wm + dt)
            self.wm_filt += _alpha_wm * (wm - self.wm_filt)

            # Limitador suave de sobre-rotação: SMOOTHSTEP 3s²−2s³ na faixa
            # 90 %→100 % de speed_ref.  Derivada nula nos dois extremos (vs.
            # rampa linear, que tem "joelho" de ganho descontínuo em 0.9·sref
            # e amplifica ruído numérico via o loop wm → iq_ref → Tem → wm).
            _lo = 0.90 * self.speed_ref
            _hi = self.speed_ref
            _prog = float(np.clip(
                (self.wm_filt - _lo) / max(_hi - _lo, 1e-12), 0.0, 1.0))
            _s     = 1.0 - _prog                       # 1 em _lo, 0 em _hi
            _fator = _s * _s * (3.0 - 2.0 * _s)        # smoothstep
            iq_from_speed = self.max_current * _fator

            # Limitador de potência DC (FSAE EV.4.1: P_dc ≤ p_max_dc).
            #   P_dc ≈ P_mec/η = Kt·iq·wm/η  →  iq_max = P_max_dc·η/(Kt·wm)
            # Produz o perfil clássico: torque constante até a velocidade-
            # base, potência constante (torque ∝ 1/ω) acima dela.
            # Usa wm_filt pela mesma razão do limitador de rotação: evita
            # realimentar ripple numérico de torque no comando de corrente.
            if self.p_max_dc is not None and self.Kt > 0.0:
                wm_pwr        = max(self.wm_filt, self.wm_floor)
                iq_from_power = (self.p_max_dc * self.eta_drive
                                 / (self.Kt * wm_pwr))
            else:
                iq_from_power = self.max_current

            # iq_ref final = o MAIS RESTRITIVO dos três limites
            iq_ref = float(np.clip(min(iq_pedal, iq_from_speed, iq_from_power),
                                   0.0, self.max_current))
            sp_err = effective_speed_ref - wm   # logado apenas para diagnóstico

            # ── Field weakening (enfraquecimento de campo) ────────────────────
            # Velocidade base: ωe_base = Vlim / λm  (ponto onde o back-EMF
            # no eixo q (ωe·λm) atinge o limite de tensão com id=0).
            # Acima da velocidade base, injetamos id_ref < 0 para reduzir
            # o fluxo efetivo λd = λm + Ld·id, liberando margem de tensão
            # para controlar iq e manter torque.
            #
            # Fórmula: id_fw = (fw_limit / ωe − λm) / Ld   (negativo acima da base)
            # fw_limit = 0.90·Vlim — operamos dentro de 90% do limite para
            # deixar margem ao PI de corrente.
            # id_fw é limitado a −0.5·max_current (proteção contra
            # desmagnetização e sobrecorrente de eixo d).
            we_abs = abs(we)
            id_ref = 0.0
            if (we_abs > self.wm_floor * self.p
                    and self.lambda_m > 0.0
                    and self.ld > 0.0):
                fw_limit = 0.90 * Vlim
                if we_abs * self.lambda_m > fw_limit:
                    id_fw  = (fw_limit / we_abs - self.lambda_m) / self.ld
                    id_ref = float(np.clip(id_fw, -self.max_current * 0.5, 0.0))

            # ── 5. PIs de corrente → comandos de tensão ──────────────────────
            #
            # vd = PI_d(id_ref − isd) − ωe·Lq·isq            (desacoplamento)
            # vq = PI_q(iq_ref − isq) + ωe·(Ld·isd + λm)     (desac. + bEMF)
            #
            # Os termos de desacoplamento (feedforward) cancelam o
            # acoplamento cruzado das equações dq, deixando para cada PI
            # uma planta de 1ª ordem independente.
            #
            # Anti-windup em DOIS estágios:
            #  (a) Interno ao PID: saturação ao próprio limit (Vdc_nominal) com
            #      back-calculation automática em update().
            #  (b) Externo: a soma PI+desacoplamento é cortada a ±Vlim, que pode
            #      ser MENOR que Vdc_nominal sob carga (queda interna na bateria).
            #      Se a soma exceder ±Vlim, devolvemos o excesso ao integrador
            #      via back_calculate() para evitar windup silencioso.
            dec_d = -we * self.lq * isq
            dec_q =  we * (self.ld * isd + self.lambda_m)

            vd_pid = self.id_ctrl.update(id_ref - isd, dt=dt)
            vq_pid = self.iq_ctrl.update(iq_ref - isq, dt=dt)

            vd_total = vd_pid + dec_d
            vq_total = vq_pid + dec_q

            vd = float(np.clip(vd_total, -Vlim, Vlim))
            vq = float(np.clip(vq_total, -Vlim, Vlim))

            if vd_total != vd:
                self.id_ctrl.back_calculate(vd_total - vd)
            if vq_total != vq:
                self.iq_ctrl.back_calculate(vq_total - vq)

            # ── 6. RK4 com vd, vq, I_bat CONGELADOS ─────────────────────────
            # telem captura intermediários de pneu/transmissão em x[k] (k1
            # do RK4) — reaproveitados pelo logger abaixo sem recomputar.
            telem = {}
            x = self._rk4(t, x, vd, vq, I_bat, dt, telemetry=telem)

            # ── 7. Logging — estado ANTES do passo ──────────────────────────
            Ft_log      = telem.get('Ft',      0.0)
            Fz_log      = telem.get('Fz',      0.0)
            Fx_log      = telem.get('Fx_max',  0.0)
            slip_log    = telem.get('slip',    0.0)
            Fresist_log = telem.get('Fresist', 0.0)
            Tcarga_log  = telem.get('Tcarga',  0.0)
            vv_acc_log  = telem.get('d_vv',    0.0)
            wt_log      = (self.transmission.motor_to_wheel_torque(Ce_now)
                           if self.transmission is not None else 0.0)

            # Correntes de fase abc para o dashboard (best-effort — se o
            # motor não fornecer a transformação, loga zeros)
            try:
                is1, is2, is3, *_ = self.motor.abc_currents_from_dq(
                    isd, isq, theta_e, self.ld*isd + self.lambda_m)
            except Exception:
                is1 = is2 = is3 = 0.0

            # Tensões de fase reais: transformada inversa de Park dos comandos vd, vq
            # (sincronizada com theta_e do rotor).
            try:
                va, vb, vc, _ = self.motor.inverse_park_transform(vd, vq, theta_e)
                vd_r, vq_r = vd, vq
            except Exception:
                va = vb = vc = vd_r = vq_r = 0.0

            # Escrita indexada em arrays pré-alocados (sem .append).
            k = self.k_log
            self.tempo[k]                    = t
            self.corrented[k]                = isd
            self.correnteq[k]                = isq
            self.corrente1[k]                = is1
            self.corrente2[k]                = is2
            self.corrente3[k]                = is3
            self.tensaosd[k]                 = vd
            self.tensaosq[k]                 = vq
            self.tensao1[k]                  = va
            self.tensao2[k]                  = vb
            self.tensao3[k]                  = vc
            self.fluxosd[k]                  = self.ld*isd + self.lambda_m
            self.fluxosq[k]                  = self.lq*isq
            self.conjugado[k]                = Ce_now
            self.velocidade[k]               = wm * 60.0 / (2*np.pi)
            self.torque_mecanico[k]          = Ce_now - Tcarga_log
            self.conjcarga[k]                = Tcarga_log
            self.temperatura[k]              = temp
            self.vd_control[k]               = vd
            self.vq_control[k]               = vq
            self.vd_real[k]                  = vd_r
            self.vq_real[k]                  = vq_r
            self.speed_error[k]              = sp_err
            self.iq_ref_trace[k]             = iq_ref
            self.id_ref_trace[k]             = id_ref
            self.vehicle_velocity[k]         = vv
            self.vehicle_position[k]         = vp
            self.vehicle_acceleration[k]     = vv_acc_log
            self.wheel_torque[k]             = wt_log
            self.tractive_force_hist[k]      = Ft_log
            self.resistive_force_hist[k]     = Fresist_log
            self.slip_ratio_hist[k]          = slip_log
            self.longitudinal_force_hist[k]  = Fx_log
            self.fz_hist[k]                  = Fz_log
            self.soc_hist[k]                 = soc
            self.Iast_hist[k]                = Iast
            self.tempo_acumulado_hist[k]     = tacc
            self.battery_voltage_hist[k]     = Vdc
            self.battery_current_hist[k]     = I_bat   # ← corrente DC real
            self.pedal_position_hist[k]      = pedal_pos
            self.effective_speed_ref_hist[k] = effective_speed_ref
            self.k_log = k + 1

            t += dt

            # ── 8. Critério de parada por distância ──────────────────────────
            if self.dmax is not None and vp >= self.dmax:
                break

        # Trunca arrays para o tamanho efetivamente preenchido (a parada
        # por distância normalmente encerra antes de N passos).
        n_used = self.k_log
        for _name in self._LOG_FIELDS:
            setattr(self, _name, getattr(self, _name)[:n_used])

        print(f"✅ Concluído — {n_used} pontos  "
              f"t_final={self.tempo[-1]:.3f} s  "
              f"dist_final={self.vehicle_position[-1]:.2f} m  "
              f"SoC final={self.soc_hist[-1]:.3f}  "
              f"Vdc final={self.battery_voltage_hist[-1]:.1f} V")

        # Dict de retorno com CHAVES PÚBLICAS (nomes ≠ atributos internos).
        # Consumidores externos (otimizador, scripts) devem usar estas chaves.
        return {
            't':                      np.array(self.tempo),
            'isd':                    np.array(self.corrented),
            'isq':                    np.array(self.correnteq),
            'is1':                    np.array(self.corrente1),
            'is2':                    np.array(self.corrente2),
            'is3':                    np.array(self.corrente3),
            'va':                     np.array(self.tensao1),
            'vb':                     np.array(self.tensao2),
            'vc':                     np.array(self.tensao3),
            'soc':                    np.array(self.soc_hist),
            'battery_voltage':        np.array(self.battery_voltage_hist),
            'battery_current':        np.array(self.battery_current_hist),
            'vehicle_velocity':       np.array(self.vehicle_velocity),
            'vehicle_acceleration':   np.array(self.vehicle_acceleration),
            'vehicle_position':       np.array(self.vehicle_position),
            'traction_force':         np.array(self.tractive_force_hist),
            'resistance_force':       np.array(self.resistive_force_hist),
            'slip_ratio':             np.array(self.slip_ratio_hist),
            'Fz':                     np.array(self.fz_hist),
            'rpm':                    np.array(self.velocidade),
            'torque_electromagnetic': np.array(self.conjugado),
            'torque_mechanical':      np.array(self.torque_mecanico),
            'vd_control':             np.array(self.vd_control),
            'vq_control':             np.array(self.vq_control),
            'iq_ref':                 np.array(self.iq_ref_trace),
            'id_ref':                 np.array(self.id_ref_trace),
            'pedal_position':         np.array(self.pedal_position_hist),
            'effective_speed_ref':    np.array(self.effective_speed_ref_hist),
        }
