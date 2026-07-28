# -*- coding: utf-8 -*-
"""
Transmission.py — Modelo do trem de força (pinhão → corrente → coroa →
diferencial LSD → semi-eixos)
======================================================================

Representa a transmissão de relação única do carro. É um contêiner de
parâmetros SEM estado dinâmico próprio: fornece a relação N, as eficiências
e as inércias que a Vehicle e a Simulation usam para acoplar motor e rodas.

Papel de cada grandeza no restante do sistema:
  • ``final_drive_ratio`` (N)  — converte velocidade e torque entre motor e
    roda; é a principal variável de decisão do otimizador (Otimizador/).
  • ``efficiency`` (η)         — multiplica o torque no sentido motor→roda e
    divide no sentido roda→motor (perdas por atrito de corrente e engrenagens).
  • Inércias — separadas por LADO da transmissão, porque a reflexão para o
    eixo do motor depende da velocidade a que cada peça gira:
      - ``sprocket_inertia`` (pinhão): gira NA VELOCIDADE DO MOTOR → soma
        direto em J_eff, sem reflexão;
      - ``coroa_inertia``, ``diff_inertia``, ``axle_inertia``: giram na
        velocidade da RODA → são refletidas por 1/N² em
        Vehicle.calculate_reflected_inertia().
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class Transmission:
    """Modelo de transmissão pinhão/corrente/coroa/diferencial LSD/semi-eixos.

    Cadeia cinemática completa:
        Motor → Pinhão (Z₁) → Corrente → Coroa (Z₂) → Diferencial LSD
              → Semi-eixos → Rodas

    A relação final pode ser definida de duas formas:
      1. Por número de dentes: N = Z_coroa / Z_pinhao (preferida — reflete
         o hardware real e é o formato que o otimizador entrega);
      2. Por valor manual ``final_drive_ratio`` (usada quando Z₁/Z₂ ainda
         não foram escolhidos).

    Parameters
    ----------
    final_drive_ratio : float, optional
        Relação de transmissão final (N = n_motor / n_roda). IGNORADO se
        Z_pinhao e Z_coroa forem fornecidos. Default 5.0.
    Z_pinhao : int, optional
        Número de dentes do pinhão (lado do motor). Se fornecido junto com
        Z_coroa, a relação é calculada como N = Z_coroa / Z_pinhao.
    Z_coroa : int, optional
        Número de dentes da coroa (lado do diferencial).
    chain_efficiency : float, optional
        Eficiência do conjunto pinhão/corrente/coroa (0–1). Default 0.98
        (corrente lubrificada típica FSAE).
    diff_efficiency : float, optional
        Eficiência do diferencial LSD (0–1). Default 0.97.
    axle_inertia : float, optional
        Inércia de CADA semi-eixo [kg·m²] (multiplicada por
        n_driven_wheels na Vehicle). Default 0.0015.
    diff_inertia : float, optional
        Inércia do corpo do diferencial [kg·m²]. Default 0.012.
    sprocket_inertia : float, optional
        Inércia do pinhão no eixo do motor [kg·m²] — gira em velocidade de
        motor, portanto soma DIRETO em J_eff, sem reflexão por N².
        Default 0.0002 (pinhão de aço pequeno).
    coroa_inertia : float, optional
        Inércia da coroa no lado do diferencial [kg·m²] — gira em velocidade
        de roda, refletida por 1/N² em Vehicle.calculate_reflected_inertia().
        Default 0.0015.
    lsd_bias_ratio : float, optional
        Razão de torque bias do diferencial LSD (roda externa / interna no
        limite de bloqueio). 1.0 = diferencial aberto, ∞ = bloqueado.
        Torsen típico: 3.0–5.0; clutch-pack FSAE: 2.0–3.0. Default 3.0.
        Nota: em trajetória reta com tração simétrica (o caso da prova de
        aceleração simulada aqui) o LSD se comporta como diferencial aberto —
        este parâmetro fica registrado para futuras análises de curva/
        manobrabilidade, mas NÃO afeta a simulação longitudinal atual.
    efficiency : float, optional
        Eficiência combinada total (override). Se fornecido, substitui
        chain_efficiency × diff_efficiency. Mantido para retrocompatibilidade
        com código que passava um único valor de eficiência.

    Attributes
    ----------
    final_drive_ratio : float
        Relação N efetiva (calculada dos dentes ou fornecida manualmente).
    efficiency : float
        Eficiência total η = chain_efficiency × diff_efficiency (ou override).
    Z_pinhao, Z_coroa : int or None
        Dentes do pinhão e da coroa. None se a relação foi manual.
    """

    def __init__(self,
                 final_drive_ratio=5.0,
                 axle_inertia=0.0015,
                 diff_inertia=0.012,
                 sprocket_inertia=0.0002,
                 efficiency=None,
                 Z_pinhao=None,
                 Z_coroa=None,
                 chain_efficiency=0.98,
                 diff_efficiency=0.97,
                 coroa_inertia=0.0015,
                 lsd_bias_ratio=3.0,
                 **kwargs):

        # ── Relação de transmissão ────────────────────────────────────────────
        # Prioridade: par de dentes (Z₂/Z₁) > valor manual > default 5.0.
        # O teste int(Z_pinhao) > 0 protege contra divisão por zero se o
        # usuário digitar 0 no dashboard.
        if Z_pinhao is not None and Z_coroa is not None and int(Z_pinhao) > 0:
            self.Z_pinhao = int(Z_pinhao)
            self.Z_coroa  = int(Z_coroa)
            self.final_drive_ratio = float(Z_coroa) / float(Z_pinhao)
        else:
            self.Z_pinhao = None
            self.Z_coroa  = None
            self.final_drive_ratio = float(final_drive_ratio) if final_drive_ratio else 5.0

        # ── Eficiências ───────────────────────────────────────────────────────
        # η_total = η_corrente × η_diferencial, salvo override explícito.
        self.chain_efficiency = float(chain_efficiency)
        self.diff_efficiency  = float(diff_efficiency)
        # Override explícito (retrocompatibilidade com código que passa efficiency=x)
        if efficiency is not None:
            self.efficiency = float(efficiency)
        else:
            self.efficiency = self.chain_efficiency * self.diff_efficiency

        # ── Inércias ──────────────────────────────────────────────────────────
        # Separadas por lado: o pinhão gira em velocidade de MOTOR (alta) e
        # os demais componentes em velocidade de RODA (baixa) — a distinção
        # importa na reflexão para J_eff (ver Vehicle.calculate_reflected_inertia).
        self.axle_inertia     = float(axle_inertia)
        self.diff_inertia     = float(diff_inertia)
        self.sprocket_inertia = float(sprocket_inertia)  # pinhão — lado motor (alta vel.)
        self.coroa_inertia    = float(coroa_inertia)     # coroa  — lado diff  (baixa vel.)

        # ── LSD ───────────────────────────────────────────────────────────────
        self.lsd_bias_ratio = float(lsd_bias_ratio)

    # ── Conversões torque ─────────────────────────────────────────────────────
    # Convenção: a eficiência SEMPRE reduz o torque útil, independente do
    # sentido — por isso multiplica em motor→roda e divide em roda→motor.

    def motor_to_wheel_torque(self, motor_torque):
        """Torque do motor → torque na roda [N·m]: T_roda = T_motor · N · η."""
        return motor_torque * self.final_drive_ratio * self.efficiency

    def wheel_to_motor_torque(self, wheel_torque):
        """Torque na roda → torque equivalente no motor [N·m]: T_motor = T_roda / (N·η)."""
        return wheel_torque / (self.final_drive_ratio * self.efficiency)

    # ── Conversões velocidade ─────────────────────────────────────────────────
    # Velocidade não sofre perda: apenas a razão cinemática N se aplica.

    def motor_to_wheel_speed(self, motor_speed):
        """Velocidade angular do motor → roda [rad/s]: ω_roda = ω_motor / N."""
        return motor_speed / self.final_drive_ratio

    def wheel_to_motor_speed(self, wheel_speed):
        """Velocidade angular da roda → motor [rad/s]: ω_motor = ω_roda · N."""
        return wheel_speed * self.final_drive_ratio
