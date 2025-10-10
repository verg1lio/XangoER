"""
Digitaltwin_Powertrain - Simulação de sistema de propulsão para veículo elétrico
"""

__version__ = "0.1.0"
__author__ = "Leonardo"
__description__ = "Digital Twin para sistema de propulsão"

# Importações principais para facilitar o uso
from models.BatteryPack import BatteryPack
from models.Inversor import Inversor
from models.PIDController import Controller
from models.Tire import Tire
from models.Transmission import Transmission
from models.Vehicle import Vehicle
from models.Motor import Motor
from models.Pedal import Pedal
from simulation import Simulation

__all__ = ['BatteryPack', 'Simulation', 'Inversor', 'Controller', 'Pedal', 'Tire', 'Transmission', 'Vehicle', 'Motor']