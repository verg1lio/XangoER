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

    def __init__(self, kp, ki, kd, limit, Ts=0.01):
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

        # Integrate provisionally; back-calculation reverses any excess that
        # gets clipped by saturation below.
        self.integral += error * dt

        derivative = 0.0
        if dt > 0 and self.kd != 0.0:
            derivative = self.kd * (error - self.prev_error) / dt

        output_raw = self.kp * error + self.ki * self.integral + derivative

        if output_raw > self.limit:
            output = self.limit
        elif output_raw < -self.limit:
            output = -self.limit
        else:
            output = output_raw

        # Back-calculation anti-windup: when the output is clipped, subtract
        # the exact excess from the integrator. The next step then starts
        # from a state consistent with the actually-applied output, which
        # avoids the asymmetric limit-cycle behaviour of plain clamping
        # anti-windup and converges cleanly when saturation releases.
        if self.anti_windup_enabled and self.ki > 0.0 and output_raw != output:
            self.integral -= (output_raw - output) / self.ki

        self.prev_error = error
        self.prev_output = output

        return output

    def back_calculate(self, excess):
        """External anti-windup hook: corrects the integrator by ``excess/ki``.

        Use when the controller output is further saturated downstream — e.g.
        after adding a feedforward decoupling term and clipping the sum to a
        dynamic limit that the PID itself does not know about. ``excess`` is
        the post-clip *shortfall*: ``output_raw - output_clipped``. Positive
        excess shrinks the integrator (pulls toward less positive output);
        negative excess grows it.
        """
        if self.anti_windup_enabled and self.ki > 0.0 and excess != 0.0:
            self.integral -= excess / self.ki

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
