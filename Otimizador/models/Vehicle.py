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
    
    def LoadTransfer(self,a):
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
    
    def calculate_single_wheel_inertia(self, k_factor=0.8):
        """
        Método auxiliar para estimar a inércia de uma única roda (J_roda) 
        usando uma massa estimada e um fator de forma (k).
        
        Parameters
        ----------
        estimated_wheel_mass : float
            Massa estimada de uma única roda e pneu [kg].
        k_factor : float
            Fator de forma (k) para a inércia. 0.5 (sólido) a 1.0 (anel). 0.8 é um bom meio termo.

        Returns
        -------
        float
            Inércia rotacional de uma roda [kg.m²].
        """
        # I = k * m * r^2
        return k_factor * self.wheel_mass * (self.wheel_radius ** 2)
    
    def calculate_reflected_inertia(self, transmission):
        """
        Calcula a inércia equivalente do veículo refletida no eixo do motor.

        J_refletido = (J_translação + J_rodas) / (i^2)

        Parameters
        ----------
        transmission : Transmission
            Objeto da transmissão contendo 'final_drive_ratio'.
        wheel_inertia : float, optional
            Momento de inércia de UMA roda/pneu [kg.m²]. 
            Se não fornecido, considera apenas a massa do veículo.

        Returns
        -------
        float
            Inércia da carga vista pelo motor [kg.m²].
        """
        # Converte Massa Translacional em Inércia Rotacional (J = m * r²)
        j_translation = self.mass * (self.wheel_radius ** 2)

        # Calcula J_roda usando a massa estimada
        j_single_wheel = self.calculate_single_wheel_inertia() 
        j_wheels = j_single_wheel * 4.0 # Inércia de 4 rodas
        
        # Inércia total no eixo da roda (Low speed side)
        j_load_total = j_translation + j_wheels

        # Reflete para o eixo do motor (High speed side)
        # DIVIDIMOS pelo quadrado da relação, pois o motor gira mais rápido
        ratio = transmission.final_drive_ratio
        
        if ratio == 0:
            return 0.0
            
        return j_load_total / (ratio ** 2)
