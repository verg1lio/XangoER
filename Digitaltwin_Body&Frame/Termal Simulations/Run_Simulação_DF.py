from Simulation import Motor, Simulation
from Simulação_bateria_DF import SimulacaoBateria
from models.Motor import Motor
from models.Vehicle import Vehicle
from models.BatteryPack import BatteryPack
from models.Transmission import Transmission
from models.Tire import Tire
from models import Inversor

'''
# =========================
# CRIAÇÃO DOS OBJETOS
# =========================

# ---- MOTOR ----
motor = Motor(
    rs=0.03,
    ld=1.2e-3,
    lq=1.3e-3,
    lambda_m=0.07,
    p=4,
    jm=0.01,          # inércia do rotor
    kf=1e-3,          # atrito viscoso
    valor_mu=0.01 )    # atrito seco

# ---- VEÍCULO ----
vehicle = Vehicle(
    mass=300.0,                 # kg
    wheel_radius=0.25,           # m
    drag_coeff=0.9,              # coef. de arrasto
    frontal_area=1.2,            # m²
    rolling_resistance=0.015,    # rolamento
    road_grade=0.0,              # plano
    environment_density=1.225    # ar ao nível do mar
)

# ---- BATERIA ----
battery = BatteryPack(
    tipo_celula='Li-ion',
    n_serie=96,
    n_paralelo=4,
    soc_inicial=1.0
)

# ---- TRANSMISSÃO ----
transmission = Transmission(
    final_drive_ratio=8.0,
    efficiency=0.95
)


# ---- PNEU ----
pacejka_params = (
    0.97,   # E   → fator de curvatura
    1.3,    # Cy  → rigidez lateral
    1.2,    # Cx  → rigidez longitudinal
    1.0,    # Cz  → carga normal de referência
    10.0,   # c1  → ganho slip
    1.9     # c2  → forma da curva
)

tire = Tire(
    pacejka_params=pacejka_params,
    tire_friction_coef=1.0
)

# ---- INVERSOR ----
inversor = Inversor(
) 
'''
# cria e roda a simulação dinâmica
sim = Simulation(
    motor=motor,
    vehicle=vehicle,
    transmission=transmission,
    battery=battery,
    tire=tire,
    tmax=600,
    steps=600
)

results = sim.simulate()

velocidade = results['vehicle_velocity']
correntes = results['battery_current']

sim_bat = SimulacaoBateria(velocidade, correntes, 96,9)
sim_bat.simular()
sim_bat.plot1()
sim_bat.plot3D()
sim_bat.plot_planos_centrais()
