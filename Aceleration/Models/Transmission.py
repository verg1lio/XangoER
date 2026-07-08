import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class Transmission:
    """Modelo de transmissão pinhão/corrente/coroa/diferencial LSD/semi-eixos.

    Representa o trem de força completo entre o eixo do motor e as rodas.
    Cadeia cinemática:
        Motor → Pinhão → Corrente → Coroa → Diferencial LSD → Semi-eixos → Rodas

    Parameters
    ----------
    final_drive_ratio : float, optional
        Relação de transmissão final (N = n_motor / n_roda). Ignorado se
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
        Inércia de cada semi-eixo [kg·m²]. Default 0.0015.
    diff_inertia : float, optional
        Inércia do corpo do diferencial [kg·m²]. Default 0.012.
    sprocket_inertia : float, optional
        Inércia do pinhão no eixo do motor [kg·m²] — gira em velocidade de
        motor, sem reflexão por N². Default 0.0002.
    coroa_inertia : float, optional
        Inércia da coroa no lado do diferencial [kg·m²] — gira em velocidade
        de roda, refletida por N² em Vehicle.calculate_reflected_inertia().
        Default 0.0015.
    lsd_bias_ratio : float, optional
        Razão de torque bias do diferencial LSD (roda externa / interna no
        limite de bloqueio). 1.0 = diferencial aberto, ∞ = bloqueado.
        Torsen típico: 3.0–5.0; clutch-pack FSAE: 2.0–3.0. Default 3.0.
        Nota: em trajetória reta com tração simétrica o LSD age como diff
        aberto — este parâmetro é relevante para análise de manobrabilidade.
    efficiency : float, optional
        Eficiência combinada total (override). Se fornecido, substitui
        chain_efficiency × diff_efficiency. Mantido para retrocompatibilidade.

    Attributes
    ----------
    final_drive_ratio : float
        Relação N efetiva (calculada ou fornecida).
    efficiency : float
        Eficiência total η = chain_efficiency × diff_efficiency (ou override).
    Z_pinhao, Z_coroa : int or None
        Dentes do pinhão e coroa. None se não especificados.
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
        if Z_pinhao is not None and Z_coroa is not None and int(Z_pinhao) > 0:
            self.Z_pinhao = int(Z_pinhao)
            self.Z_coroa  = int(Z_coroa)
            self.final_drive_ratio = float(Z_coroa) / float(Z_pinhao)
        else:
            self.Z_pinhao = None
            self.Z_coroa  = None
            self.final_drive_ratio = float(final_drive_ratio) if final_drive_ratio else 5.0

        # ── Eficiências ───────────────────────────────────────────────────────
        self.chain_efficiency = float(chain_efficiency)
        self.diff_efficiency  = float(diff_efficiency)
        # Override explícito (retrocompatibilidade com código que passa efficiency=x)
        if efficiency is not None:
            self.efficiency = float(efficiency)
        else:
            self.efficiency = self.chain_efficiency * self.diff_efficiency

        # ── Inércias ──────────────────────────────────────────────────────────
        self.axle_inertia     = float(axle_inertia)
        self.diff_inertia     = float(diff_inertia)
        self.sprocket_inertia = float(sprocket_inertia)  # pinhão — lado motor (alta vel.)
        self.coroa_inertia    = float(coroa_inertia)     # coroa  — lado diff  (baixa vel.)

        # ── LSD ───────────────────────────────────────────────────────────────
        self.lsd_bias_ratio = float(lsd_bias_ratio)

    # ── Conversões torque ─────────────────────────────────────────────────────

    def motor_to_wheel_torque(self, motor_torque):
        """Torque do motor → torque na roda [N·m]."""
        return motor_torque * self.final_drive_ratio * self.efficiency

    def wheel_to_motor_torque(self, wheel_torque):
        """Torque na roda → torque equivalente no motor [N·m]."""
        return wheel_torque / (self.final_drive_ratio * self.efficiency)

    # ── Conversões velocidade ─────────────────────────────────────────────────

    def motor_to_wheel_speed(self, motor_speed):
        """Velocidade angular do motor → velocidade angular da roda [rad/s]."""
        return motor_speed / self.final_drive_ratio

    def wheel_to_motor_speed(self, wheel_speed):
        """Velocidade angular da roda → velocidade angular do motor [rad/s]."""
        return wheel_speed * self.final_drive_ratio
