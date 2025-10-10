import numpy as np
class Controller:
    """A PID controller with anti-windup and output saturation.

    This class implements a Proportional-Integral-Derivative (PID) controller 
    with optional anti-windup compensation and output limitation.

    Parameters
    ----------
    kp : float
        Proportional gain.
    ki : float
        Integral gain.
    kd : float
        Derivative gain.
    limit : float
        Output saturation limit (absolute value).
    Ts : float, optional
        Default sampling time [s]. If no `dt` is passed to `update()`, this 
        value will be used. Default is 0.001.

    Attributes
    ----------
    kp : float
        Proportional gain.
    ki : float
        Integral gain.
    kd : float
        Derivative gain.
    limit : float
        Output saturation limit.
    Ts : float
        Default sampling period [s].
    integral : float
        Current integral term accumulator.
    prev_error : float
        Previous control error, used for derivative calculation.
    prev_output : float
        Previous output value.
    anti_windup_enabled : bool
        Flag to enable or disable anti-windup correction.

    Methods
    -------
    update(error, dt=None):
        Computes the control signal given the current error.
    reset():
        Resets controller states (integral, error, and output).
    set_parameters(kp, ki, kd):
        Updates PID gains.
    """

    def __init__(self, kp, ki, kd, limit, Ts=0.1):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.limit = limit
        self.Ts = Ts
        
        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_output = 0.0
        self.anti_windup_enabled = True
        
    def update(self, error, dt=None):
        if dt is None:
            dt = self.Ts
            
        # Proportional term
        proportional = self.kp * error
        
        # Integral term
        self.integral += error * dt
        
        # Derivative term
        derivative = 0.0
        if dt > 0:
            derivative = self.kd * (error - self.prev_error) / dt
        
        # Control signal
        output = proportional + self.ki * self.integral + derivative
        
        # Saturation and anti-windup
        if output > self.limit:
            output = self.limit
            if self.anti_windup_enabled:
                self.integral -= error * dt
        elif output < -self.limit:
            output = -self.limit
            if self.anti_windup_enabled:
                self.integral -= error * dt
        
        self.prev_error = error
        self.prev_output = output
        
        return output
    
    def reset(self):
        """Resets controller internal states."""
        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_output = 0.0
        
    def set_parameters(self, kp, ki, kd):
        """Updates PID controller parameters."""
        self.kp = kp
        self.ki = ki
        self.kd = kd
