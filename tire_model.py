import numpy as np


class Dynamics:
    def __init__(self, spring_type, spring_k, spring_F, spring_non_lin_coef, tire_Fz, tire_Sa, tire_Ls, tire_type):
        # Modelo de mola
        self.spring_type = spring_type  # Hooke, Softening
        self.spring_k = spring_k  # rigidez da mola [N/m]
        self.spring_F = spring_F  # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef  # coeficiente de ganho não-linear
        # Modelo de pneu
        self.tire_Fz = tire_Fz  # carga vertical no pneu [N]
        self.tire_Sa = np.radians(tire_Sa)  # slip angle do pneu [rad]
        self.tire_Ls = tire_Ls  # longitudinal slip do pneu [Admensional]
        self.tire_type = tire_type  # Default, Admensional

    def Tire(self):
        # Pacejka parâmetros
        B = 1.00
        E = 1.00
        Cy = 1.30  # C para força lateral
        Cx = 1.65  # C para força longitudinal
        Cz = 2.40  # C para momento de torque auto-alinhante

        if self.tire_type == 'Default':
            D = 0.85 * self.tire_Fz
            tire_lateral_force = D * np.sin(Cy * np.arctan(B * self.tire_Sa - E * (B * self.tire_Sa - np.arctan(B * self.tire_Sa))))
            tire_longitudinal_force = D * np.sin(Cx * np.arctan(B * self.tire_Sa - E * (B * self.tire_Sa - np.arctan(B * self.tire_Sa))))
            tire_auto_align_moment = D * np.sin(Cz * np.arctan(B * self.tire_Sa - E * (B * self.tire_Sa - np.arctan(B * self.tire_Sa))))

        print("Tire Lateral Force:", tire_lateral_force)
        print("Tire Longitudinal Force:", tire_longitudinal_force)
        print("Tire Auto-Aligning Moment:", tire_auto_align_moment)


# Exemplo de uso:
dynamics = Dynamics(spring_type="Hooke", spring_k=1000, spring_F=500, spring_non_lin_coef=0.1,
                    tire_Fz=5000, tire_Sa=5, tire_Ls=0.1, tire_type="Default")
dynamics.Tire()





