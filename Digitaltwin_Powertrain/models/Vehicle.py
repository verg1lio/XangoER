import numpy as np
class Vehicle:
    """Vehicle model for dynamic simulation.

    This class represents a vehicle's longitudinal dynamics,
    including aerodynamic drag, rolling resistance, and road grade effects.

    Parameters
    ----------
    mass : float
        Vehicle mass [kg].
    wheel_radius : float
        Radius of the wheels [m].
    drag_coeff : float
        Aerodynamic drag coefficient.
    frontal_area : float
        Frontal area of the vehicle [m²].
    rolling_resistance : float
        Rolling resistance coefficient.
    road_grade : float, optional
        Road slope angle [radians]. Default is 0.
    environment_density : float, optional
        Air density [kg/m³]. Default is 1.225.

    Attributes
    ----------
    mass : float
        Vehicle mass [kg].
    wheel_radius : float
        Wheel radius [m].
    drag_coeff : float
        Aerodynamic drag coefficient.
    frontal_area : float
        Frontal area [m²].
    rolling_resistance : float
        Rolling resistance coefficient.
    road_grade : float
        Road slope angle [radians].
    environment_density : float
        Air density [kg/m³].
    g : float
        Gravitational acceleration [m/s²].
    half_rho : float
        Precomputed 0.5 * environment_density for aerodynamic calculations.

    Methods
    -------
    calculate_resistance_forces(velocity):
        Computes the total resistance force acting on the vehicle.
    calculate_load_torque(velocity, transmission):
        Computes the torque reflected to the motor due to vehicle resistance.
    """

    def __init__(self, mass, wheel_radius, drag_coeff, frontal_area, rolling_resistance,
                 road_grade=0, environment_density=1.225):
        self.mass = mass  # kg
        self.wheel_radius = wheel_radius  # m
        self.drag_coeff = drag_coeff
        self.frontal_area = frontal_area  # m²
        self.rolling_resistance = rolling_resistance
        self.road_grade = road_grade  # radians
        self.environment_density = environment_density  # kg/m³
        self.g = 9.81  # gravitational acceleration [m/s²]
        self.half_rho = 0.5 * environment_density  # precompute for aerodynamic force

    def calculate_resistance_forces(self, velocity):
        """Computes the total resistance forces acting on the vehicle.

        Parameters
        ----------
        velocity : float
            Vehicle longitudinal speed [m/s].

        Returns
        -------
        float
            Total resistance force [N], including aerodynamic drag,
            rolling resistance, and road grade effect.
        """
        v2 = velocity**2  # precompute velocity squared
        aerodynamic_force = self.half_rho * self.drag_coeff * self.frontal_area * v2
        rolling_force = self.rolling_resistance * self.mass * self.g * np.cos(self.road_grade)
        grade_force = self.mass * self.g * np.sin(self.road_grade)

        return aerodynamic_force + rolling_force + grade_force

    def calculate_load_torque(self, velocity, transmission):
        """Computes the torque reflected to the motor from vehicle resistance.

        Parameters
        ----------
        velocity : float
            Vehicle longitudinal speed [m/s].
        transmission : Transmission
            Transmission object for converting wheel torque to motor torque.

        Returns
        -------
        float
            Equivalent motor torque required to overcome resistance [Nm].
        """
        resistance_force = self.calculate_resistance_forces(velocity)
        wheel_torque = resistance_force * self.wheel_radius
        motor_torque = transmission.wheel_to_motor_torque(wheel_torque)
        return motor_torque
