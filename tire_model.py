import numpy as np
class Dynamics:
    
    def __init__(self, spring_type, spring_k, spring_F, spring_non_lin_coef, tire_Fz, tire_Sa, tire_Ls, tire_type):
        # Modelo de mola
        self.spring_type = spring_type # Hooke, Softening
        self.spring_k = spring_k # rigidez da mola [N/m]
        self.spring_F = spring_F # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef # coeficiente de ganho não-linear
        # Modelo de pneu
        self.tire_Fz = tire_Fz # carga vertical no pneu [N]
        self.tire_Sa = tire_Sa # slip angle do pneu [deg]
        self.tire_Ls = tire_Ls # longitudinal slip do pneu [Admensional]
        self.tire_type = tire_type # Default, Admensional
    
    def Tire(self):
        
        self.tire_Sa = np.radians(self.tire_Sa)
        
        # Pacejka parâmetros

        Cy = 1.30 # C para força lateral
        Cx = 1.65 # C para força longitudinal
        Cz = 2.40 # C para momento de torque auto-alinhante
        if self.tire_type == 'Default':
            D = 0.85*self.tire_Fz
            tire_lateral_force = 
