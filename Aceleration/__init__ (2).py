"""
Digitaltwin_Powertrain - Simulação de sistema de propulsão para veículo elétrico
"""

__version__ = "0.1.5"
__author__ = "Leonardo, Marco Affonso"
__description__ = "Digital Twin para sistema de propulsão"

# Importações principais para facilitar o uso
from Models.BatteryPack import BatteryPack
from Models.Inversor import Inversor
from Models.PIDController import Controller
from Models.Tire import Tire
from Models.Transmission import Transmission
from Models.Vehicle import Vehicle
from Models.Motor import Motor
from Models.Pedal import Pedal
from Simulation.Simulation import Simulation

__all__ = ['BatteryPack', 'Simulation', 'Inversor', 'Controller', 'Pedal', 'Tire', 'Transmission', 'Vehicle', 'Motor']