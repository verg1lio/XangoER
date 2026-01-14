import numpy as np

class Vehicle:
    """Vehicle model for dynamic simulation.

    This class represents a vehicle's longitudinal dynamics,
    including aerodynamic drag, rolling resistance, and road grade effects.
    
    Updated to handle geometric parameters (L, h, dist_cg).

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
    L : float, optional
        Wheelbase [m].
    h : float, optional
        Height of Center of Gravity [m].
    dist_cg : float, optional
        Distance from front/rear axle to CG [m].
    **kwargs : dict
        catches any extra arguments to prevent TypeErrors.

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
    L : float
        Wheelbase [m].
    h : float
        Height of Center of Gravity [m].
    dist_cg : float
        Distance from axle to CG [m].
    """

    def __init__(self, mass, wheel_radius,wheel_mass, drag_coeff, frontal_area, rolling_resistance,
                 road_grade=0, environment_density=1.225, 
                 L=None, h=None, dist_cg=None):
        
        # Standard dynamics parameters
        self.mass = mass  # kg
        self.wheel_radius = wheel_radius  # m
        self.wheel_mass = wheel_mass
        self.drag_coeff = drag_coeff
        self.frontal_area = frontal_area  # m²
        self.rolling_resistance = rolling_resistance
        self.road_grade = road_grade  # radians
        self.environment_density = environment_density  # kg/m³
        
        # Geometric parameters 
        self.L = L
        self.h = h
        self.dist_cg = dist_cg
        
        # Constants and Precomputations
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
    
    def calculate_load_transfer(self,a):
        """
        Calcula a força vertical dinâmica (Fz) considerando a transferência de carga longitudinal.

        A equação baseia-se em:
        Fz = (m * g * dist_cg / L) +/- (m * a * h / L)

        Parameters
        ----------
        m : float
            Massa total do veículo [kg].
        g : float
            Aceleração da gravidade [m/s²].
        L : float
            Distância entre-eixos (Wheelbase) [m].
        h : float
            Altura do Centro de Gravidade (CG) [m].
        a : float
            Aceleração longitudinal do veículo [m/s²].
        dist_cg : float
            Distância do CG até o eixo OPOSTO ao que está sendo calculado.
            - Para Fz Traseira: use 'b' (ou a distância que define a alavanca de peso).
            - Para Fz Dianteira: use 'a'.

        Returns
        -------
        float
            Força normal vertical (Fz) no eixo [N].
        """
        # Componente estática (distribuição de peso parada)
        # Nota: Se dist_cg for 'b' e estamos calculando traseira, certifique-se
        # que sua geometria veicular define 'b' corretamente como a alavanca para o eixo traseiro.
        Fz_static = (self.mass * self.g * self.dist_cg) / self.L

        # Componente dinâmica (Transferência de carga devido à aceleração)
        load_transfer = (self.mass * a * self.h) / self.L

        return (Fz_static + load_transfer)/2
        
    def calculate_reflected_inertia(self, transmission):
        """
        Calcula a inércia total da carga (veículo + rodas + transmissão) refletida no eixo do motor.
        
        Aplica o conceito de conservação de energia cinética, onde a massa translacional
        do veículo é convertida em uma inércia rotacional equivalente.

        Fórmula:
        J_refletido = (J_translacao + J_rotacional_rodas + J_transmissao_lenta) / N^2

        Onde:
        - J_translacao = Massa_Total * Raio_Pneu^2
        - N = Relação de transmissão total (Gear Ratio)

        Parameters
        ----------
        transmission : Transmission
            Objeto contendo 'final_drive_ratio', 'diff_inertia' e 'axle_inertia'.
        
        Returns
        -------
        float
            Inércia equivalente total vista pelo eixo do motor [kg.m²].
        """

        N = transmission.final_drive_ratio
        
        # Inércia Translacional Equivalente 
        # Representa a massa do carro como se fosse um grande volante de inércia acoplado às rodas.
      
        j_translation = self.mass * (self.wheel_radius ** 2)

        # Inércia Rotacional das Rodas
        # k = fator de forma (0.8 é uma boa aproximação para conjunto pneu+roda)

        k = 0.8
        j_wheel_single = k * self.wheel_mass * (self.wheel_radius ** 2)
        
        # Consideramos as 4 rodas do veículo para a inércia rotacional
        j_wheels_rotational = 4 * j_wheel_single

        # Inércia dos Componentes da Transmissão 

        j_drivetrain_low_speed = transmission.diff_inertia + (2 * transmission.axle_inertia)

        # Soma Total na Roda 
        # Aqui somamos todas as inércias que giram na velocidade da roda
        j_total_at_wheel = j_translation + j_wheels_rotational + j_drivetrain_low_speed

        # Reflexão para o Motor 
       
        j_reflected = j_total_at_wheel / (N ** 2)

        return j_reflected