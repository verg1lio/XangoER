import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

class PWT:
    """A PWT (Power Train) object.
222
    This class will create a Power Train digital twin with the battery, power conversor, motor and controller elements provided.

    Parameters AINDA A ALTERAR, ISSO É UM EXEMPLO
    ----------
    bat_type : string
    String descevendo modelo usado para a bateria. No momento temos implementado:
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

    def __init__(self, bat_type, bat_Idis, bat_Ich, bat_t, bat_Es, bat_K, bat_Q, bat_A, bat_B, bat_Iast, bat_R,
    cnv_input_voltage, cnv_output_voltage, cnv_efficiency, 
    mot_type, mot_K, mot_T, 
    time=None):
    # time):
        model = input('Select model: S for Shepherd or R for Rint: ').strip().upper()
        if bat_type == 'S':
            self.bat_type = 'Modified Shepherd'
        elif bat_type == 'R':
            self.bat_type = 'Rint model'
        else:
            self.bat_type = bat_type 
        self.bat_Idis = bat_Idis # Idis is the current discharge
        self.bat_Ich = bat_Ich # Ich is the current charge
        self.bat_t = bat_t # Depicts the sampling time
        self.bat_Es = bat_Es # Es represents the open-circuit voltage
        self.bat_K = bat_K # K indicates the change in polarization resistance factor in (mΩ A/h)
        self.bat_Q = bat_Q # Q: is the battery nominal capacity (Ah)
        self.bat_A = bat_A # A: voltage factor
        self.bat_B = bat_B # B: capacity factor.
        self.bat_Iast = bat_Iast # depicts the filtered current
        self.bat_R = bat_R # R is the internal resistance 
        self.cnv_input_voltage = cnv_input_voltage # in V
        self.cnv_output_voltage = cnv_output_voltage # in V
        self.cnv_efficiency = cnv_efficiency # in %
        self.mot_type = mot_type # Types = DC, Triphase, XXX, XXX
        self.mot_K = mot_K # Motor gain
        self.mot_T = mot_T # Time constant
        self.mot_numerator = [mot_K]
        self.mot_denominator = [mot_T, 1] 
        self.mot_tf = signal.TransferFunction(self.mot_numerator, self.mot_denominator) # Usado na função de exemplo

    #Podemos entrar com a variável time, ou utilizar valores da capacidade/corrente da bateria.
        if time == None:
            self.time = self.bat_capacity / self.bat_current #A*h/A = h
        else:
            self.time = time #Talvez esse seja a única versão "correta"

    def PWT(self, ):
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
        if self.bat_type == 'Modified Shepherd':
            V_dis = self.bat_Es - self.bat_R * self.bat_Idis - self.bat_K * (self.bat_Q / (self.bat_Q - (self.bat_Idis * self.bat_t))) * ((self.bat_Idis * self.bat_t) + self.bat_Iast) + self.bat_A  * np.exp(- self.bat_B * self.bat_Idis * self.bat_t) 
            V_ch = self.bat_Es - self.bat_R * self.bat_Ich - self.bat_K * (self.bat_Q / ( (self.bat_Ich * self.bat_t) - 0.1*self.bat_Q )) * self.bat_Iast - self.bat_K*(((self.bat_Q)/(self.bat_Q-(self.bat_Ich*self.bat_t))) * (self.bat_Ich* self.bat_t)) + self.bat_A  * np.exp(- self.bat_B * self.bat_Ich * self.bat_t)
        elif self.bat_type == 'Rint model':
            V_dis = self.bat_Es - self.bat_R * self.bat_Idis 
            V_ch =  self.bat_Es - self.bat_R * self.bat_Ich 
        return V_dis, V_ch



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



    # def Motor(self, input_signal): Original
    def Motor(self):
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


        time, response = signal.step(self.mot_tf, T = np.linspace(0,self.time,10000), X0=0.0)
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

        #return Cont
    
def example():
    # Motor parameters
    K = 1.0 # Motor gain
    T = 0.1 # Time constant

    # Create the DC motor model
    mot_tf = signal.TransferFunction([K], [T, 100])

    # Time vector for simulation
    t = np.linspace(0, 5, 1000)

    # Input voltage (step input)
    u = np.ones_like(t)

    # Simulate the motor response to the input
    time, response = signal.step(mot_tf, T=t)

    # Plot the step response of the motor
    plt.figure()
    plt.plot(time, response)
    plt.xlabel('Time (s)')
    plt.ylabel('Angular Velocity')
    plt.title('Step Response of the DC Motor')
    plt.grid(True)
    plt.show()
    powertrain = PWT(bat_capacity = 10, 
    bat_voltage = 2, 
    bat_current = 1, 
    bat_p_stacks = 1, 
    bat_s_stacks = 1, 
    cnv_input_voltage = 1, 
    cnv_output_voltage = 1, 
    cnv_efficiency = 1, 
    mot_type = "ACREDITO SER UMA STRING", #Tipo do motor, DC/Triphase/Etc (Corrente Direta/Trifásico)
    mot_K = 20, # Velocidade Angular do Motor 
    mot_T = 10) # Constante de tempo

    tempo_motor = powertrain.Motor()[0]
    resposta_motor = powertrain.Motor()[1]

    plt.figure()
    plt.plot(tempo_motor, resposta_motor)
    plt.xlabel('Time (s)')
    plt.ylabel('Angular Velocity')
    plt.title('Step Response of the DC Motor')
    plt.grid(True)
    plt.show()

example()
