# -*- coding: utf-8 -*-
"""
Tire.py — Modelo de pneu Pacejka PAC2002 (Magic Formula 5.2) longitudinal
==========================================================================

Responde à pergunta central do acoplamento motor↔veículo: "quanta força
longitudinal Fx o pneu consegue transmitir ao solo, dada a carga vertical
Fz e o escorregamento κ (slip ratio)?"

O fluxo dentro da Simulation é:
  1. ``SlipRatio()`` calcula κ a partir de ω_roda, raio e velocidade do carro;
  2. ``Tire_forces(Fz, κ)`` avalia a Magic Formula e devolve Fx de UM pneu;
  3. A Simulation multiplica por n_driven_wheels e usa esse valor como TETO
     da força de tração (a força pedida pelo motor é saturada nesse limite).

É esse teto que produz o comportamento realista de wheelspin: pedir mais
torque do que o pneu transmite não acelera mais o carro — apenas aumenta o
slip e (na curva calibrada, pós-pico) REDUZ Fx.
"""

import numpy as np


class Tire:
    """Pacejka MF 5.2 / PAC2002 — força longitudinal Fx de escorregamento puro.

    Implementa a Magic Formula longitudinal canônica do PAC2002 com
    sensibilidade a carga (dfz) e camber (γ).  Forma geral (Pacejka,
    "Tire and Vehicle Dynamics", cap. 4):

        Fx0 = Dx · sin( Cx · arctan( Bx·κ − Ex·(Bx·κ − arctan(Bx·κ)) ) )

    onde, com  dfz = (Fz − Fz0) / Fz0  (desvio relativo da carga nominal):

        Cx  = PCX1                                  (fator de forma)
        μx  = (PDX1 + PDX2·dfz) · (1 − PDX3·γ²)     (atrito de pico)
        Dx  = |μx| · λμx · Fz                       (valor de pico de Fx)
        Kxκ = Fz · (PKX1 + PKX2·dfz) · exp(PKX3·dfz) (rigidez longitudinal)
        Bx  = Kxκ / (Cx · Dx)                       (fator de rigidez)
        Ex  = (PEX1 + PEX2·dfz + PEX3·dfz²) · (1 − PEX4·sgn(κ))  (curvatura)

    Papel intuitivo de cada fator: Kxκ é a inclinação da curva na origem
    (região linear), Dx é a altura do pico, Cx controla quanto a curva cai
    após o pico (assíntota Dx·sin(Cx·π/2)) e Ex ajusta a posição/formato
    do pico.

    Parameters
    ----------
    pacejka_params : dict
        Coeficientes PAC2002 longitudinais.  Chaves esperadas:
        PCX1, PDX1, PDX2, PDX3, PEX1, PEX2, PEX3, PEX4, PKX1, PKX2.
        PKX3 é opcional (default 0 → exp(0)=1).
    tire_friction_coef : float
        Fator de correção de atrito λμx (scrub/Calspan factor).
        Para datasets ajustados em Calspan/TTC o pico |μx| sai ~50%
        inflado pela esteira de papel-lixa; valores realistas em pista
        ficam em **λμx ≈ 0.55–0.70**.  NÃO usar > 1.0 com o dataset
        `HOOSIER_FSAE_LONG` (já tem |PDX1| ≈ 3.5).
    Fz0 : float, optional
        Carga nominal de referência [N] do ajuste PAC2002.  Default
        654 N (~150 lbf), padrão TTC para slick FSAE 13".  Convém
        aproximar da carga estática real (mass·g / n_rodas) para
        reduzir extrapolação no termo dfz, sobretudo porque PDX2 > 0
        no dataset Hoosier amplifica μx em carga abaixo de Fz0.
    camber : float, optional
        Camber γ em radianos.  Default 0 (FSAE com 0° estático).

    Caveats
    -------
    * **Convenção ISO de Fz.** O ajuste PAC2002 da Hoosier FSAE usa
      Fz < 0 (apontando para baixo, ISO 8855).  Com PDX1 ≈ −3.51 e Fz
      negativo, D = μx·Fz sai positivo naturalmente.  Aqui usamos Fz
      positivo (magnitude) e tomamos `|μx|`, o que é numericamente
      idêntico — não é workaround, é apenas a outra metade da convenção.
    * **Pico Calspan inflado.** |PDX1| ≈ 3.5 vem de bancada TTC e é
      tipicamente ~50% maior que o atrito real em asfalto.  Compensar
      via λμx ≈ 0.6, não via reajuste de PDX1.
    * **Baixa carga FSAE.** O dataset Hoosier responde mal em Fz baixo
      (reportado pela comunidade Racer/overtake.gg).  PDX2 = +0.633 faz
      μx subir quando dfz < 0, contra a tendência real.  Mitigar
      escolhendo Fz0 próximo da carga estática operacional do carro.
    * Ex é clampado em ≤ 1 (Pacejka 2002, eq. 4.E12 — requisito do MF).
    * O slip aplicado à MF não é limitado artificialmente — o limite físico
      é κ ∈ [-1, +1] aplicado em SlipRatio(). Isso permite operação na cauda
      descendente da curva, que é fisicamente real em patinagem total.
    """

    # Parâmetros brutos do ajuste TTC/Calspan — APENAS REFERÊNCIA.
    # NÃO usar diretamente: geram curva monotonicamente crescente (sem pico),
    # porque PEX1=1.346 é clampado para Ex=1.0 pela MF, eliminando o pico,
    # e PKX1=58.52 coloca o pico teórico em κ≈0.13 — irrealisticamente baixo.
    # Mantidos aqui como registro da origem dos dados (bancada TTC).
    HOOSIER_FSAE_LONG_CALSPAN = {
        'PCX1':  1.27872550e+00,
        'PDX1': -3.51415830e+00,
        'PDX2':  6.33383390e-01,
        'PDX3':  7.11629130e+00,
        'PEX1':  1.34635150e+00,
        'PEX2':  3.84490190e-17,
        'PEX3': -3.92323870e-17,
        'PEX4':  7.42710500e-03,
        'PKX1':  5.85243740e+01,
        'PKX2':  5.47564740e+00,
        'PKX3':  0.0,
    }

    # Parâmetros OPERACIONAIS — calibrados para comportamento físico realista.
    # É este o dicionário usado pelo dashboard e pelo otimizador.
    #
    #  Problema do conjunto Calspan bruto (acima):
    #    PEX1 = 1.346 → clampado para Ex = 1.0 (máximo MF) → curva sem pico,
    #    monotonicamente crescente até κ=1 (derrapagem total = maior tração!).
    #    PCX1 = 1.279 → assintota Dx·sin(Cx·π/2) = 0.906·Dx → queda máx. 9.4%
    #    PKX1 = 58.52 → Bx=21.7 → κ_peak(E=0) ≈ 0.13.
    #
    #  Correções aplicadas (sem alterar o NÍVEL de atrito PDX1/PDX2/PDX3):
    #    PCX1 1.279 → 1.50  : assintota cai para 0.707·Dx → queda máx. 29.3%
    #                          cria pico real e degradação pós-pico visível
    #    PKX1 58.52 → 23.0  : Bx = 8.53 → κ_peak ≈ 0.20 (FSAE slick típico)
    #    PKX2  5.48 →  2.15 : escalonado proporcionalmente (mesma sensibilidade à carga)
    #    PEX1  1.35 → -0.50 : curva de curvatura — de platô plano para queda progressiva;
    #                          PEX1<0 é fisicamente correto para pneus racing
    #                          (ver Pacejka 2002, §4.3.2)
    #
    #  Resultado verificado numericamente:
    #    κ_peak = 0.204  |  queda a κ=0.5: -11.2%  |  queda a κ=1.0: -19.7%
    #
    #  É essa queda pós-pico que permite ao otimizador dispensar penalidade
    #  explícita de slip: patinar já custa tempo naturalmente.
    HOOSIER_FSAE_LONG = {
        'PCX1':  1.50,
        'PDX1': -3.51415830e+00,
        'PDX2':  6.33383390e-01,
        'PDX3':  7.11629130e+00,
        'PEX1': -0.50,
        'PEX2':  3.84490190e-17,
        'PEX3': -3.92323870e-17,
        'PEX4':  7.42710500e-03,
        'PKX1':  23.0,
        'PKX2':  2.15,
        'PKX3':  0.0,
    }

    def __init__(self, pacejka_params, tire_friction_coef,
                 Fz0=654.0, camber=0.0):
        # Cópia defensiva: o dashboard passa o dict do store e poderia
        # mutá-lo depois; PKX3 ganha default 0 se ausente.
        self.pacejka_params = dict(pacejka_params)
        self.pacejka_params.setdefault('PKX3', 0.0)
        self.tire_friction_coef = float(tire_friction_coef)
        self.Fz0 = float(Fz0)
        self.camber = float(camber)

    def Tire_forces(self, Fz, s):
        """Força longitudinal Fx(Fz, κ) [N] pela Magic Formula PAC2002.

        Parameters
        ----------
        Fz : float
            Carga vertical no pneu [N] (magnitude; ver caveat de convenção
            ISO na docstring da classe). Valores ≤ 0 são protegidos por um
            piso de 1e-6 N.
        s : float
            Slip ratio κ (adimensional, tipicamente vindo de ``SlipRatio``).

        Returns
        -------
        float
            Força longitudinal de UM pneu [N]. Positiva em tração
            (κ > 0), negativa em frenagem (κ < 0).
        """
        p = self.pacejka_params
        # Proteções numéricas: Fz=0 (roda no ar) e Dx=0 anulariam divisões.
        Fz_safe = max(abs(float(Fz)), 1e-6)
        dfz = (Fz_safe - self.Fz0) / self.Fz0

        # Fator de forma e atrito de pico (com sensibilidade a carga e camber)
        Cx = p['PCX1']
        mux = (p['PDX1'] + p['PDX2'] * dfz) * (1.0 - p['PDX3'] * self.camber ** 2)
        Dx = abs(mux) * self.tire_friction_coef * Fz_safe
        D_safe = max(Dx, 1e-6)

        # Rigidez longitudinal e fator de rigidez Bx = Kxκ/(Cx·Dx)
        Kxk = Fz_safe * (p['PKX1'] + p['PKX2'] * dfz) * np.exp(p['PKX3'] * dfz)
        Bx = Kxk / (Cx * D_safe)
        B_safe = max(abs(Bx), 1e-6)

        # Fator de curvatura Ex — assimétrico em tração/frenagem via PEX4;
        # clampado em ≤ 1 (requisito da MF, Pacejka 2002 eq. 4.E12).
        sgn_s = float(np.sign(s)) if np.sign(s) != 0 else 1.0
        Ex = (p['PEX1'] + p['PEX2'] * dfz + p['PEX3'] * dfz ** 2) \
             * (1.0 - p['PEX4'] * sgn_s)
        Ex = min(Ex, 1.0)

        # Fórmula Mágica propriamente dita
        arg = Bx * float(s)
        Fx = Dx * np.sin(
            Cx * np.arctan(arg - Ex * (arg - np.arctan(arg)))
        )
        return Fx

    @staticmethod
    def SlipRatio(velocidade_angular, raio_pneu, velocidade_linear,
                  eps=0.1, slip_max=1.0):
        """Slip ratio saturado em [-slip_max, +slip_max] (convenção SAE J670).

        Forma robusta para baixa velocidade:
            s = (ω·r − v) / max(|ω·r|, |v|, eps)

        Esta normalização (em vez da divisão clássica por v) é naturalmente
        limitada em [-1, +1]: quando a roda gira no vazio (v=0, ω>0) o
        denominador vira |ω·r| e s → +1 (não infinito); quando v cresce
        acompanhando a roda, s → 0.  Isso evita o slip explosivo na partida
        (v≈0), que colocava a Magic Formula na cauda descendente e fazia o
        carro "perder" tração fictícia no primeiro instante da largada.

        Parameters
        ----------
        velocidade_angular : float or array_like
            Velocidade angular da RODA [rad/s] (já convertida do motor
            pela relação N).
        raio_pneu : float
            Raio efetivo de rolamento do pneu [m].
        velocidade_linear : float or array_like
            Velocidade linear do veículo [m/s].
        eps : float, optional
            Piso do denominador para evitar divisão por zero quando ambas
            as velocidades são ≈ 0 (carro parado). Default 0.1.
        slip_max : float, optional
            Saturação dura do slip retornado.  Default 1.0 — slip > 1
            sinaliza patinagem total e não tem significado físico adicional
            no modelo simplificado de Pacejka.

        Returns
        -------
        float or ndarray
            Slip ratio adimensional em [-slip_max, +slip_max].  Positivo =
            patinagem (aceleração), negativo = travamento (frenagem).
        """
        v = np.asarray(velocidade_linear, dtype=float)
        omega = np.asarray(velocidade_angular, dtype=float)
        v_wheel = omega * raio_pneu           # velocidade periférica da roda
        # Denominador robusto: maior entre |ω·r|, |v| e o piso eps
        denom = np.maximum.reduce([np.abs(v_wheel), np.abs(v), np.full_like(v, eps)])
        slip = (v_wheel - v) / denom
        return np.clip(slip, -slip_max, slip_max)
