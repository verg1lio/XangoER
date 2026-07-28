# -*- coding: utf-8 -*-
"""
Motor.py — Motor síncrono de ímãs permanentes (PMSM) — EMRAX 228
=================================================================

Contêiner de parâmetros elétricos, mecânicos e térmicos do motor, mais as
transformações de coordenadas (Park inversa e dq→abc) usadas para logging.

DIVISÃO DE RESPONSABILIDADES (importante):
  • Motor.py NÃO possui dinâmica nem estado de controlador. Ele apenas
    ARMAZENA parâmetros (Rs, Ld, Lq, λm, p, limites térmicos, ...) e
    fornece funções puras de transformação de coordenadas.
  • Toda a dinâmica (equações dq, temperatura, RK4) e todo o controle
    (PIs de corrente, field weakening, limitadores) vivem em
    Simulation/Simulation.py, que LÊ os parâmetros deste objeto no
    __init__ via getattr.

Convenção de referencial (crítica para consistência):
  O projeto usa a transformada de Park INVARIANTE EM AMPLITUDE — a
  amplitude de pico da corrente de fase é igual a |i_dq|. Consequências:
    • Torque: T = 1.5 · p · λm · iq   (o fator 1.5 vem dessa convenção);
    • max_current é AMPLITUDE DE PICO de fase (= √2 · I_rms);
    • As transformações abaixo não usam o fator √(2/3) da forma
      invariante em potência.
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from Constants.constants import PI23, SQRT3_2


class Motor:
    """PMSM — contêiner de parâmetros + transformações de coordenadas.

    Esta classe é um contêiner SEM ESTADO de parâmetros elétricos,
    mecânicos e térmicos do PMSM, mais as transformações Park/Clarke
    usadas pelo logger da simulação. As malhas de controle (FOC, PIs de
    corrente) e a ODE da planta vivem em :class:`Simulation.Simulation` —
    Motor não possui nenhum estado de controlador ou dinâmica.

    Parameters
    ----------
    rs : float
        Resistência de estator por fase [Ω], medida a T_ref (a simulação
        aplica a lei Rs(T) = rs·[1 + α·(T − T_ref)] a partir deste valor).
    ld, lq : float
        Indutâncias de eixo direto e de quadratura [H]. Para o EMRAX 228
        (rotor de ímãs superficiais) Ld ≈ Lq — sem torque de relutância.
    jm : float
        Inércia do rotor [kg·m²].
    kf : float
        Coeficiente de atrito viscoso [N·m·s] (rolamentos + ventilação).
        Atenção: valores altos dupla-contariam as perdas no ferro, que já
        são modeladas separadamente (calibrado em 0.005 no projeto).
    lambda_m : float
        Fluxo concatenado dos ímãs permanentes λm [Wb]. Define a constante
        de torque Kt = 1.5·p·λm e o back-EMF ωe·λm.
    p : int
        Número de PARES de polos (EMRAX 228: 10).
    valor_mu : float
        Índice de modulação do inversor (0–1) para escala de tensão.
    speed_ref : float, optional
        Velocidade de referência do rotor [rad/s] usada pelo limitador
        suave de sobre-rotação da simulação (fade de iq entre 90 % e
        100 % deste valor). Default 0.0.
    k_h_iron : float, optional
        Coeficiente de perdas por HISTERESE [W/Hz]: P_h ≈ k_h·|f_elec|.
        Default 0.9 — calibrado para perda ≈ 600 W a 670 Hz (EMRAX 228
        @ 4000 RPM, divisão 60/40 histerese/Foucault do datasheet).
    k_e_iron : float, optional
        Coeficiente de perdas por CORRENTES DE FOUCAULT [W/Hz²]:
        P_e ≈ k_e·f_elec². Default 9e-4 — junto com k_h dá ~1 kW total
        de perdas no ferro a 4000 RPM. Reduzir para motores de baixo
        ferro (estator PCB etc.).
    alpha_cu : float, optional
        Coeficiente térmico do cobre [1/K] para Rs(T) = Rs0·[1+α·(T−T_ref)].
        Default 0.00393 (cobre recozido — resistência sobe ~0.4 %/K).
    T_ref : float, optional
        Temperatura de referência [°C] na qual rs foi medida. Default 25.
    T_alarm, T_max : float, optional
        Limiares do derate térmico [°C]. Abaixo de T_alarm: sem redução.
        Acima de T_max: corrente zero. Rampa linear no intervalo.
        Defaults 130 / 160 °C (limites típicos de isolamento classe H).
    h_cool : float, optional
        Coeficiente de convecção [W/(m²·K)] do modelo de resfriamento
        P_cool = h·A·(T − T_amb). Default 10 (convecção natural).
        Referências: ar forçado ≈ 50–200; líquido forçado ≈ 500–2000.
    A_cool : float, optional
        Área efetiva de troca térmica [m²]. Default 0.15 (carcaça EMRAX).
    T_ambient : float, optional
        Temperatura ambiente [°C] do resfriamento. Default 25.
    max_current : float, optional
        Amplitude de PICO máxima da corrente de fase [A]. Se None, usa
        300·√2 ≈ 424 A. O projeto passa 323 A → torque de pico
        1.5·10·0.04748·323 ≈ 230 N·m (datasheet EMRAX 228).

    Attributes
    ----------
    pi23 : float
        Constante 2π/3 (defasagem de 120° entre fases).
    rs0 : float
        Rs nominal a T_ref — referência IMUTÁVEL da lei Rs(T).
    rs : float
        Rs efetiva na temperatura corrente (atualizada pela simulação).
    m, C : float
        Massa térmica [kg] e calor específico do cobre [J/(kg·K)] usados
        pela ODE térmica do enrolamento (defaults EMRAX 228 LC: 13.5 kg,
        385 J/(kg·K)).
    max_current : float
        Amplitude máxima de pico da corrente de fase [A] — teto de iq_ref
        na simulação.
    Vdc : float
        Tensão DC de fallback [V] usada quando nenhum modelo de bateria é
        fornecido à simulação. Default 600 V.
    """

    def __init__(self, rs, ld, lq, jm, kf, lambda_m, p, valor_mu, speed_ref=0.0,
                 k_h_iron=0.9, k_e_iron=9e-4,
                 alpha_cu=0.00393, T_ref=25.0,
                 T_alarm=130.0, T_max=160.0,
                 h_cool=10.0, A_cool=0.15, T_ambient=25.0,
                 max_current=None):
        self.pi23 = PI23

        # ── Parâmetros elétricos / mecânicos ─────────────────────────────
        # rs0 = Rs nominal a T_ref (referência imutável da lei Rs(T)).
        # rs  = Rs efetiva na temperatura corrente (atualizada pela
        # Simulation a cada passo). Ferramentas que precisam da
        # resistência "fria" devem usar rs0.
        self.rs0 = rs
        self.rs = rs
        self.ld = ld
        self.lq = lq
        self.jm = jm
        self.kf = kf
        self.lambda_m = lambda_m
        self.p = p
        self.valor_mu = valor_mu

        # ── Modelo térmico — defaults EMRAX 228 LC ───────────────────────
        # Massa térmica concentrada (lumped): todo o calor gerado aquece
        # uma única massa equivalente m com calor específico C.
        self.m = 13.5      # massa térmica [kg]
        self.C = 385.0     # calor específico (cobre) [J/(kg·K)]

        # ── Limites operacionais ─────────────────────────────────────────
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

        # Velocidade de referência (usada pelo limitador suave de
        # sobre-rotação da simulação — NÃO é setpoint de malha fechada)
        self.speed_ref = speed_ref

        # ── Perdas no ferro (forma de Steinmetz, termos separados) ───────
        # P_iron = k_h·|f_e| (histerese) + k_e·f_e² (Foucault),
        # com f_e = frequência elétrica em Hz. Avaliadas na ODE da simulação.
        self.k_h_iron = float(k_h_iron)
        self.k_e_iron = float(k_e_iron)

        # ── Compensação de temperatura Rs(T) ─────────────────────────────
        self.alpha_cu = float(alpha_cu)
        self.T_ref    = float(T_ref)

        # ── Limiares de derate térmico ───────────────────────────────────
        self.T_alarm = float(T_alarm)
        self.T_max   = float(T_max)

        # ── Modelo de resfriamento convectivo ────────────────────────────
        # P_cool = h_cool · A_cool · (T − T_ambient)
        # Defaults EMRAX 228 LC: convecção natural (h≈10), superfície ≈0.15 m²
        # Para líquido forçado: h ≈ 500–2000 W/(m²·K); ar forçado: h ≈ 50–200
        self.h_cool    = float(h_cool)
        self.A_cool    = float(A_cool)
        self.T_ambient = float(T_ambient)

    def inverse_park_transform(self, vd, vq, theta_e):
        """Transformada inversa de Park: tensões dq → tensões de fase abc.

        Duas etapas encadeadas (ambas invariantes em amplitude):
          1. Rotação dq → αβ pelo ângulo elétrico θe (Park inversa);
          2. Projeção αβ → abc (Clarke inversa): fase a alinhada com α,
             fases b/c defasadas ±120° via os fatores −1/2 e ±√3/2.

        Usada apenas para LOGGING das tensões de fase que o inversor
        aplicaria — não realimenta a dinâmica (a ODE trabalha direto em dq).

        Parameters
        ----------
        vd, vq : float
            Tensões de eixo direto e de quadratura [V] (comandos do FOC).
        theta_e : float
            Ângulo elétrico do rotor [rad] (= p · θ_mecânico).

        Returns
        -------
        (vs1, vs2, vs3, 0.0) : tuple of float
            Tensões trifásicas [V]. O 0.0 final é mantido para
            desempacotamento retrocompatível em :mod:`Simulation`.
        """
        cos_theta = np.cos(theta_e)
        sin_theta = np.sin(theta_e)
        # Park inversa: rotaciona o vetor (vd, vq) do referencial girante
        # para o referencial estacionário αβ
        valpha = vd * cos_theta - vq * sin_theta
        vbeta  = vd * sin_theta + vq * cos_theta
        # Clarke inversa (amplitude-invariante): αβ → abc
        vs1 = valpha
        vs2 = -0.5 * valpha + SQRT3_2 * vbeta
        vs3 = -0.5 * valpha - SQRT3_2 * vbeta
        return vs1, vs2, vs3, 0.0

    def abc_currents_from_dq(self, isd, isq, theta_e, flux_d):
        """Correntes de fase abc a partir das correntes dq (amplitude-invariante).

        Projeção direta de cada fase: i_k = isd·cos(θe − φk) − isq·sin(θe − φk),
        com φk ∈ {0, +2π/3, −2π/3}. O modelo de torque usa Kt = 1.5·p·λm,
        consistente com esta forma amplitude-invariante (sem fator √(2/3)).

        Também é função exclusiva de LOGGING (aba "Correntes" do dashboard) —
        a dinâmica elétrica evolui inteiramente no referencial dq.

        Parameters
        ----------
        isd, isq : float
            Correntes de eixo direto e de quadratura [A].
        theta_e : float
            Ângulo elétrico do rotor [rad].
        flux_d : float
            Fluxo concatenado no eixo d [Wb] — repassado sem alteração no
            retorno, apenas por conveniência do logger.

        Returns
        -------
        (is1, is2, is3, flux_d, flux_d, flux_d) : tuple of float
        """
        is1 = isd * np.cos(theta_e)             - isq * np.sin(theta_e)
        is2 = isd * np.cos(theta_e - self.pi23) - isq * np.sin(theta_e - self.pi23)
        is3 = isd * np.cos(theta_e + self.pi23) - isq * np.sin(theta_e + self.pi23)
        return is1, is2, is3, flux_d, flux_d, flux_d
