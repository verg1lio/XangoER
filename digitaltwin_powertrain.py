import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

class PWT:
    """A PWT (Power Train) object.

    This class will create a Power Train digital twin with the battery, power conversor, motor and controller elements provided.

    Parameters AINDA A ALTERAR, ISSO É UM EXEMPLO
    ----------
    shaft_elements : list
        List with the shaft elements.
    disk_elements : list
        List with the disk elements.
    bearing_elements : list
        List with the bearing elements.
    automeshing : boolean
        Set it True to use the automeshing method. Default is False.
        If automeshing is True, the previous shaft_elements parameter is now defined by:
            shaft_elements : list
                List with the length, inner and outter diameter and material of each element, as follows:
                    [[length, inner, outter, material],[length, inner, outter, material],...]
        For the other parameters please check the respective bearing and disk classes for more information.
    **kwargs : dict, optional
        If automeshing is True, these parameters needs to be informed.
        The possible arguments are:
            alpha : float
                Proportional damping coefficient, associated to the element Mass matrix
            beta : float
                Proportional damping coefficient, associated to the element Stiffness matrix

    Returns
    -------
    A rotor object.

    Attributes
    ----------
    MM : array
        Global mass matrix.
    KK : array
        Global stiffness matrix.
    CCgyros: array
        Global gyroscopic matrix.
    CCtotal: array
        Global damping matrix

    Examples
    --------
    >>> import lmest_rotor as lm
    >>> rotor = lm.rotor_example()
    >>> rotor.MM
    array(30x30)
    """

    def __init__(self, time, bat_capacity, bat_voltage, bat_current, bat_p_stacks, bat_s_stacks, cnv_input_voltage, 
                 cnv_output_voltage, cnv_efficiency, mot_type, mot_K, mot_T):
        self.bat_capacity = bat_capacity # in Ah
        self.bat_voltage = bat_voltage # in V
        self.bat_current = bat_current # in A
        self.bat_p_stacks = bat_p_stacks # battery parallel stacks
        self.bat_s_stacks = bat_s_stacks # battery series stacks
        self.bat_energy = self.bat_capacity * self.bat_voltage # in Wh
        self.cnv_input_voltage = cnv_input_voltage # in V
        self.cnv_output_voltage = cnv_output_voltage # in V
        self.cnv_efficiency = cnv_efficiency # in %
        self.mot_type = mot_type # Types = DC, Triphase, XXX, XXX
        self.mot_K = mot_K  # Motor gain
        self.mot_T = mot_T  # Time constant
        self.mot_numerator = [K]
        self.mot_denominator = [T, 1]
        self.mot_tf = signal.TransferFunction(self.mot_numerator, self.mot_denominator)
        


    def Battery(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        if current <= 0:
            return 0
        else:
            discharge_time = self.bat_capacity / self.bat_current # in hours
            if discharge_time > time:
                discharge_time = time
            discharge_energy = self.bat_current * self.bat_voltage * discharge_time # in Wh
            self.bat_energy -= discharge_energy
            self.bat_capacity = self.bat_energy / self.bat_voltage # in Ah
        return discharge_energy



    def Conversor(self, discharge_energy):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        input_power = discharge_energy / self.cnv_efficiency
        cnv_input_current = input_power / self.cnv_input_voltage
        cnv_output_current = discharge_energy / self.cnv_output_voltage
        return cnv_input_current, cnv_output_current



    def Motor(self, input_signal, time):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        time, response = signal.step(self.mot_tf, T=time, X0=0.0, input=input_signal)
        return time, response



    def Controller(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return Cont


def example():
    # Motor parameters
    K = 1.0       # Motor gain
    T = 0.1       # Time constant

    # Create the DC motor model
    mot_tf = signal.TransferFunction([K], [T, 1])

    # Time vector for simulation
    t = np.linspace(0, 5, 1000)

    # Input voltage (step input)
    u = np.ones_like(t)

    # Simulate the motor response to the input
    time, response = mot_tf.simulate(u, t)

    # Plot the step response of the motor
    plt.figure()
    plt.plot(time, response)
    plt.xlabel('Time (s)')
    plt.ylabel('Angular Velocity')
    plt.title('Step Response of the DC Motor')
    plt.grid(True)
    plt.show()


