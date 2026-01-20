import sys
import os

# Adiciona o diretório pai ao path do Python
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
class Transmission:
    """Modelo de transmissão do veículo.

    Esta classe representa a transmissão responsável por transferir
    torque e velocidade entre o motor e as rodas, considerando a
    relação de transmissão final (final drive ratio) e a eficiência mecânica.

    Parameters
    ----------
    final_drive_ratio : float
        Relação de transmissão final (razão entre motor e roda).
    efficiency : float, optional
        Eficiência da transmissão (0 < efficiency ≤ 1). Default é 0.95.

    Attributes
    ----------
    final_drive_ratio : float
        Relação de transmissão final (motor → roda).
    efficiency : float
        Eficiência da transmissão.
    
    Methods
    -------
    motor_to_wheel_torque(motor_torque):
        Converte torque do motor em torque na roda.
    wheel_to_motor_torque(wheel_torque):
        Converte torque da roda em torque no motor.
    motor_to_wheel_speed(motor_speed):
        Converte velocidade do motor em velocidade da roda.
    wheel_to_motor_speed(wheel_speed):
        Converte velocidade da roda em velocidade do motor.
    """

    def __init__(self, final_drive_ratio, efficiency=0.95):
        self.final_drive_ratio = final_drive_ratio
        self.efficiency = efficiency

    def motor_to_wheel_torque(self, motor_torque):
        """Converte torque do motor para torque na roda.

        Parameters
        ----------
        motor_torque : float
            Torque gerado pelo motor [Nm].

        Returns
        -------
        float
            Torque resultante na roda [Nm].
        """
        return motor_torque * self.final_drive_ratio * self.efficiency

    def wheel_to_motor_torque(self, wheel_torque):
        """Converte torque na roda para torque no motor.

        Parameters
        ----------
        wheel_torque : float
            Torque aplicado na roda [Nm].

        Returns
        -------
        float
            Torque equivalente no motor [Nm].
        """
        return wheel_torque / (self.final_drive_ratio * self.efficiency)

    def motor_to_wheel_speed(self, motor_speed):
        """Converte velocidade angular do motor para velocidade angular da roda.

        Parameters
        ----------
        motor_speed : float
            Velocidade angular do motor [rad/s].

        Returns
        -------
        float
            Velocidade angular da roda [rad/s].
        """
        return motor_speed / self.final_drive_ratio

    def wheel_to_motor_speed(self, wheel_speed):
        """Converte velocidade angular da roda para velocidade angular do motor.

        Parameters
        ----------
        wheel_speed : float
            Velocidade angular da roda [rad/s].

        Returns
        -------
        float
            Velocidade angular do motor [rad/s].
        """
        return wheel_speed * self.final_drive_ratio
