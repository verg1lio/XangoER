import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
import numpy as np
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import plotly.express as px
import dash_bootstrap_components as dbc

# ---------------------
# Constantes pré-calculadas
# ---------------------
SQRT3 = np.sqrt(3)
SQRT3_2 = SQRT3 / 2
PI = np.pi
PI23 = 2 * PI / 3
RQ23 = np.sqrt(2 / 3)
TWO_THIRDS = 2 / 3

# ---------------------
# Classe BatteryPack (modelo Modified Shepherd simplificado)
# ---------------------
class BatteryPack:
    def __init__(self, tipo_celula, n_serie, n_paralelo, soc_inicial=1.0):
        self.parametros = self.definir_celula(tipo_celula)
        if self.parametros is None:
            raise ValueError(f"Tipo de célula desconhecido: {tipo_celula}")

        # Extrair parâmetros
        params = self.parametros
        self.E0 = params['E0']
        self.K = params['K']
        self.Q = params['Q']
        self.A = params['A']
        self.B = params['B']
        self.R = params['R']
        self.peso = params['peso']
        self.volume = params['volume']

        self.n_serie = n_serie
        self.n_paralelo = n_paralelo

        self.soc = soc_inicial
        self.carga_total = self.Q * n_paralelo  # Ah
        self.inv_carga_total = 1.0 / (self.carga_total * 3600.0)  # Pré-calcular inverso

        # Variáveis de estado para integração
        self.Iast = 0.0  # Corrente acumulada (Coulombs)
        self.tempo_acumulado = 0.0  # Tempo total (s)

        # Histórico
        self.tensao_hist = []
        self.corrente_hist = []
        self.soc_hist = []
        self.tempo = []

    def definir_celula(self, tipo):
        catalogo = {
            'Li-ion': {
                'E0': 3.7,
                'K': 0.005,
                'Q': 3.0,
                'A': 0.1,
                'B': 0.01,
                'R': 0.02,
                'peso': 0.045,
                'volume': 0.025
            },
            'LiFePO4': {
                'E0': 3.2,
                'K': 0.003,
                'Q': 2.8,
                'A': 0.05,
                'B': 0.015,
                'R': 0.008,
                'peso': 0.07,
                'volume': 0.03
            }
        }
        return catalogo.get(tipo)

    def calcular_derivadas(self, corrente):
        """Calcula as derivadas para integração numérica"""
        dsoc_dt = -corrente * self.inv_carga_total  # d(SOC)/dt
        dIast_dt = corrente  # d(Iast)/dt
        dtime_dt = 1.0  # d(tempo_acumulado)/dt

        return dsoc_dt, dIast_dt, dtime_dt

    def calcular_tensao(self, corrente, soc, Iast, tempo_acumulado):
        """Calcula a tensão da bateria baseada no estado atual"""
        carga_utilizada = (1 - soc) * self.carga_total
        denom1 = max(self.carga_total - carga_utilizada, 1e-6)
        inv_denom1 = 1.0 / denom1  # Pré-calcular inverso

        # Pré-calcular termos comuns
        term_exp = self.A * np.exp(-self.B * carga_utilizada)
        term_K = self.K * self.carga_total * inv_denom1
        term_carga = term_K * carga_utilizada
        term_Iast = term_K * (Iast / 3600)

        if corrente >= 0:  # Descarga
            if corrente < 100 and soc >= 0.0:
                Vt = self.E0 - self.R * corrente - term_Iast - term_carga + term_exp
            else:  # Limitar corrente máxima de descarga
                Vt = self.E0 - self.R * corrente * 0.2 - term_Iast - term_carga + term_exp
        else:  # Carga
            corrente_abs = abs(corrente)
            denom2 = max(corrente_abs * tempo_acumulado - 0.1 * self.carga_total, 1e-6)
            term_K_charge = self.K * self.carga_total / denom2 * (Iast / 3600)
            Vt = self.E0 - self.R * corrente_abs - term_K_charge - term_carga + term_exp

        return Vt * self.n_serie

    def calcular_peso_total(self):
        return self.n_serie * self.n_paralelo * self.peso

    def calcular_volume_total(self):
        return self.n_serie * self.n_paralelo * self.volume

    def calcular_tensao_nominal(self):
        return self.E0 * self.n_serie

# ---------------------
# Modelo de Inversor Modificado
# ---------------------
class Inversor:
    """Modelo de inversor trifásico para conversão DC-AC usando abordagem senoidal"""
    def __init__(self, eficiencia=0.95, freq_chaveamento=10000):
        self.eficiencia = eficiencia
        self.freq_chaveamento = freq_chaveamento
        self.Vdc = 0  # Tensão de entrada DC (será definida pela bateria)
        self.Vs = 0  # Tensão de fase (será calculada a partir de Vdc)
        self.tete = 0  # Ângulo elétrico do estator
        self.ws = 0  # Velocidade síncrona (rad/s)
        self.h = 1.0 / freq_chaveamento  # Período de chaveamento
        self.pi23 = PI23
        self.sqrt3_2 = SQRT3_2
        self.two_thirds = TWO_THIRDS

    def set_Vdc(self, Vdc):
        """Define a tensão DC de entrada e calcula a tensão de fase"""
        self.Vdc = Vdc

    def set_ws(self, ws):
        """Define a velocidade síncrona"""
        self.ws = ws

    def source_voltage(self):
        """Tensão da fonte"""
        # Atualiza o ângulo elétrico do estator
        self.tete += self.h * self.ws
        if self.tete >= 2 * PI:
            self.tete -= 2 * PI

        # Pré-calcular funções trigonométricas
        cos_tete = np.cos(self.tete)
        cos_tete_pi23 = np.cos(self.tete - self.pi23)
        cos_tete_pi23_2 = np.cos(self.tete + self.pi23)

        # Calcula as tensões de cada fase
        vs1 = self.Vdc * cos_tete
        vs2 = self.Vdc * cos_tete_pi23
        vs3 = self.Vdc * cos_tete_pi23_2

        return vs1, vs2, vs3

    def park_transform(self, va, vb, vc, theta_e):
        """Calcula Vq e Vd a partir das tensões de fase"""
        # Pré-calcular funções trigonométricas
        cos_theta = np.cos(theta_e)
        sin_theta = np.sin(theta_e)

        # Transformada de Clarke
        valpha = self.two_thirds * (va - 0.5 * vb - 0.5 * vc)
        vbeta = self.two_thirds * (self.sqrt3_2 * vb - self.sqrt3_2 * vc)

        # Transformada de Park
        vd = valpha * cos_theta + vbeta * sin_theta
        vq = -valpha * sin_theta + vbeta * cos_theta

        return vq, vd

# ---------------------
# Controladores, Transmissão, Veículo, Motor
# ---------------------
class PIController:
    """Controlador PI para FOC"""
    def __init__(self, kp, ki, limit):
        self.kp = kp
        self.ki = ki
        self.limit = limit
        self.integral = 0.0
        self.prev_error = 0.0

    def update(self, error, dt):
        self.integral += error * dt
        self.prev_error = error

        output = self.kp * error + self.ki * self.integral

        if output > self.limit:
            output = self.limit
            self.integral -= error * dt
        elif output < -self.limit:
            output = -self.limit
            self.integral -= error * dt

        return output

    def reset(self):
        self.integral = 0.0
        self.prev_error = 0.0

class Tire:
    def __init__(self, pacejka_params, tire_friction_coef):
        self.pacejka_params = pacejka_params
        self.tire_friction_coef = tire_friction_coef

    # CORRIJA O MÉTODO para receber Fz e Ls como argumentos
    def Tire_forces(self, Fz, Ls):
        E, Cy, Cx, Cz, c1, c2 = self.pacejka_params

        Cs = c1 * np.sin(2 * np.arctan(Fz / c2))
        D = self.tire_friction_coef * Fz

        # Adicionar verificação para evitar divisão por zero
        if (Cx * D) == 0:
            return 0.0

        Bx = Cs / (Cx * D)

        arg_arctan = Bx * Ls
        tire_longitudinal_force = D * np.sin(Cx * np.arctan(arg_arctan - E * (arg_arctan - np.arctan(arg_arctan))))
        return tire_longitudinal_force 

    @staticmethod
    def SlipRatio(velocidade_angular, raio_pneu, velocidade_linear, eps=0.1):
        # velocidade_linear em m/s, raio_pneu em m
        v = np.array(velocidade_linear, dtype=float)
        omega = np.array(velocidade_angular, dtype=float)
        denom = np.where(np.abs(v) < eps, eps, v)
        return (omega * raio_pneu - v) / denom

class Transmission:
    """Modelo de transmissão do veículo"""
    def __init__(self, final_drive_ratio, efficiency=0.95):
        self.final_drive_ratio = final_drive_ratio
        self.efficiency = efficiency

    def motor_to_wheel_torque(self, motor_torque):
        """Converte torque do motor para torque na roda"""
        return motor_torque * self.final_drive_ratio * self.efficiency

    def wheel_to_motor_torque(self, wheel_torque):
        """Converte torque na roda para torque no motor"""
        return wheel_torque / (self.final_drive_ratio * self.efficiency)

    def motor_to_wheel_speed(self, motor_speed):
        """Converte velocidade do motor para velocidade da roda"""
        return motor_speed / self.final_drive_ratio

    def wheel_to_motor_speed(self, wheel_speed):
        """Converte velocidade da roda para velocidade do motor"""
        return wheel_speed * self.final_drive_ratio

class Vehicle:
    """Modelo de veículo para simulação de dinâmica"""
    def __init__(self, mass, wheel_radius, drag_coeff, frontal_area, rolling_resistance,
                 road_grade=0, environment_density=1.225):
        self.mass = mass  # kg
        self.wheel_radius = wheel_radius  # m
        self.drag_coeff = drag_coeff
        self.frontal_area = frontal_area  # m²
        self.rolling_resistance = rolling_resistance
        self.road_grade = road_grade  # radianos
        self.environment_density = environment_density  # kg/m³ (densidade do ar)
        self.g = 9.81  # Aceleração da gravidade
        self.half_rho = 0.5 * environment_density  # Pré-calcular

    def calculate_resistance_forces(self, velocity):
        """Calcula as forças de resistência ao movimento"""
        v2 = velocity**2  # Pré-calcular
        aerodynamic_force = self.half_rho * self.drag_coeff * self.frontal_area * v2
        rolling_force = self.rolling_resistance * self.mass * self.g * np.cos(self.road_grade)
        grade_force = self.mass * self.g * np.sin(self.road_grade)

        return aerodynamic_force + rolling_force + grade_force

    def calculate_load_torque(self, velocity, transmission):
        """Calcula o torque de carga refletido no motor"""
        resistance_force = self.calculate_resistance_forces(velocity)
        wheel_torque = self.calculate_resistance_forces(velocity) * self.wheel_radius
        motor_torque = transmission.wheel_to_motor_torque(wheel_torque)
        return motor_torque

class Motor:
    """Motor PMSM com FOC usando solve_ivp. Agora com integração à bateria através do inversor."""
    def __init__(self, rs, ld, lq, jm, kf, lambda_m, p, valor_mu, TL=False, torque=0.0,
                 vehicle=None, transmission=None, battery: BatteryPack = None, tire: Tire = None):
        self.pi23 = PI23
        self.rq23 = RQ23
        self.rs = rs
        self.ld = ld
        self.lq = lq
        self.jm = jm
        self.kf = kf
        self.lambda_m = lambda_m
        self.p = p
        self.valor_mu = valor_mu

        # Pré-calcular constantes
        self.one_point_five_p = 1.5 * p
        self.torque_constant = self.one_point_five_p * lambda_m
        self.inv_ld = 1.0 / ld
        self.inv_lq = 1.0 / lq
        self.inv_jm = 1.0 / jm

        self.m = 22.0
        self.C = 0.385
        self.inv_mC = 1.0 / (self.m * self.C)

        self.id_controller = PIController(kp=0.5, ki=100.0, limit=600.0)
        self.iq_controller = PIController(kp=0.5, ki=100.0, limit=600.0)
        self.speed_controller = PIController(kp=10.0, ki=0, limit=600.0)

        self.TL = bool(TL)
        if callable(torque):
            self._external_torque = torque
        else:
            v = float(torque)
            self._external_torque = (lambda t, v=v: v)

        self.max_current = 220.0 * np.sqrt(2)
        self.Vdc = 600  # Valor inicial, será atualizado pela bateria
        self.Vs = (self.valor_mu * self.Vdc) / np.sqrt(3)
        self.Vlimit = 600

        self.id_ref = 0.0
        self.iq_ref = 0.0
        self.speed_ref = 471.23

        self.tmax = 11
        self.hp = self.tmax / 2000.0

        # Modelos de transmissão, veículo, bateria e inversor
        self.transmission = transmission
        self.vehicle = vehicle
        self.battery = battery
        self.tire = tire
        self.inversor = Inversor()  # Novo modelo de inversor

        self.reset_initial_conditions()
        self.initialize_storage()

    def reset_initial_conditions(self):
        self.cl = 0.0
        self.wm0 = 0.0
        self.theta_m0 = 0.0
        self.isd0 = 0.0
        self.isq0 = 0.0
        self.iso0 = 0.0
        self.temp0 = 25.0
        self.int_id0 = 0.0
        self.int_iq0 = 0.0
        self.int_speed0 = 0.0
        self.vehicle_velocity0 = 0.0
        self.vehicle_position0 = 0.0

        # Estados iniciais da bateria
        if self.battery:
            self.soc0 = self.battery.soc
            self.Iast0 = self.battery.Iast
            self.tempo_acumulado0 = self.battery.tempo_acumulado
        else:
            self.soc0 = 1.0
            self.Iast0 = 0.0
            self.tempo_acumulado0 = 0.0

    def initialize_storage(self):
        self.tempo = []
        self.corrented = []
        self.correnteq = []
        self.corrente1 = []
        self.corrente2 = []
        self.corrente3 = []
        self.tensaosd = []
        self.tensaosq = []
        self.tensao1 = []
        self.tensao2 = []
        self.tensao3 = []
        self.fluxosd = []
        self.fluxosq = []
        self.conjugado = []
        self.velocidade = []
        self.conjcarga = []
        self.torque_mecanico = []
        self.temperatura = []
        self.vd_control = []
        self.vq_control = []
        self.vd_real = []
        self.vq_real = []
        self.speed_error = []
        self.iq_ref_trace = []
        self.id_ref_trace = []

        # Armazenamento para dados do veículo
        self.vehicle_velocity = []  # m/s
        self.vehicle_position = []  # m
        self.vehicle_acceleration = []  # m/s²
        self.wheel_torque = []  # Nm

        #Armazenamento para dados de dinâmica
        self.tractive_force_hist = []   # Força de tração real (N)
        self.resistive_force_hist = []  # Soma das forças de resistência (N)
        self.slip_ratio_hist = []       # Slip ratio (adimensional)
        self.longitudinal_force_hist = [] # Força longitudinal máxima (Fx_max) (N)
        
        # Armazenamento para dados da bateria
        self.soc_hist = []
        self.Iast_hist = []
        self.tempo_acumulado_hist = []
        self.battery_voltage_hist = []
        self.battery_current_hist = []

    def set_external_torque(self, torque):
        if callable(torque):
            self._external_torque = torque
        else:
            v = float(torque)
            self._external_torque = (lambda t, v=v: v)

    def enable_external_torque(self, enable: bool):
        self.TL = bool(enable)

    def set_load(self, t):
        if self.TL and (self._external_torque is not None):
            try:
                return float(self._external_torque(t))
            except Exception:
                pass
        if t < 0.5:
            return 0.0
        elif t < 1.0:
            return 100.0
        else:
            return -200.0

    def field_oriented_control(self, isd, isq, wm, dt):
        speed_error = self.speed_ref - wm
        torque_ref = self.speed_controller.update(speed_error, dt)

        iq_ref = torque_ref / self.torque_constant if self.torque_constant != 0 else 0.0
        id_ref = 0.0

        if abs(iq_ref) > self.max_current:
            iq_ref = np.sign(iq_ref) * self.max_current

        error_d = id_ref - isd
        error_q = iq_ref - isq

        we = self.p * wm
        decoupling_d = -we * self.lq * isq
        decoupling_q = we * (self.ld * isd + self.lambda_m)

        vd = self.id_controller.update(error_d, dt) + decoupling_d
        vq = self.iq_controller.update(error_q, dt) + decoupling_q

        self.vd_control.append(vd)
        self.vq_control.append(vq)

        vmag = np.sqrt(vd**2 + vq**2)
        if vmag > self.Vlimit:
            v_scale = self.Vlimit / vmag
            vd *= v_scale
            vq *= v_scale

        return vd, vq, id_ref, iq_ref, speed_error

    def inverse_park_transform(self, vd, vq, theta_e):
        cos_theta = np.cos(theta_e)
        sin_theta = np.sin(theta_e)
        valpha = vd * cos_theta - vq * sin_theta
        vbeta = vd * sin_theta + vq * cos_theta
        vs1 = valpha
        vs2 = -0.5 * valpha + SQRT3_2 * vbeta
        vs3 = -0.5 * valpha - SQRT3_2 * vbeta
        return vs1, vs2, vs3, 0.0

    def abc_currents_from_dq(self, isd, isq, theta_e, flux_d):
        cos_theta = np.cos(theta_e)
        sin_theta = np.sin(theta_e)
        is1 = self.rq23 * (isd * cos_theta - isq * sin_theta)
        is2 = self.rq23 * (isd * np.cos(theta_e - self.pi23) - isq * np.sin(theta_e - self.pi23))
        is3 = self.rq23 * (isd * np.cos(theta_e + self.pi23) - isq * np.sin(theta_e + self.pi23))
        fs_val = self.rq23 * flux_d
        return is1, is2, is3, fs_val, fs_val, fs_val

    def combined_edos(self, t, x):
        """Equações diferenciais combinadas do motor, veículo e bateria"""
        # Estados do motor
        isd, isq, iso, wm, theta_m, temp, int_id, int_iq, int_speed = x[:9]

        # Estados do veículo (se existirem)
        if self.vehicle and self.transmission:
            vehicle_velocity, vehicle_position = x[9:11]
            # Estados da bateria (se existir)
            if self.battery:
                soc, Iast, tempo_acumulado = x[11:14]
            else:
                soc, Iast, tempo_acumulado = 1.0, 0.0, 0.0
        else:
            vehicle_velocity, vehicle_position = 0.0, 0.0
            # Estados da bateria (se existir)
            if self.battery:
                soc, Iast, tempo_acumulado = x[9:12]
            else:
                soc, Iast, tempo_acumulado = 1.0, 0.0, 0.0

        theta_e = self.p * theta_m
        we = self.p * wm

        # Cálculo do torque de carga baseado na dinâmica do veículo
        if self.vehicle and self.transmission:
            cl = self.vehicle.calculate_load_torque(vehicle_velocity, self.transmission)

            # 1. Calcule a velocidade angular da roda
            velocidade_angular_roda = self.transmission.motor_to_wheel_speed(wm)

            # 2. Calcule o Slip Ratio dinamicamente
            # Note que agora chamamos Tire.SlipRatio, pois é um método estático
            slip = Tire.SlipRatio(velocidade_angular_roda, self.vehicle.wheel_radius, vehicle_velocity)

            # 3. Calcule a força vertical (Fz)
            Fz = self.vehicle.mass * self.vehicle.g 

            # 4. Calcule a força longitudinal máxima (Fx_max) que o pneu pode fornecer
            Fx_max = self.tire.Tire_forces(Fz, slip)

            # 5. Calcule a força de tração ideal (sem limite de aderência)
            wheel_torque_ideal = self.transmission.motor_to_wheel_torque(self.torque_constant * isq)
            traction_force_ideal = wheel_torque_ideal / self.vehicle.wheel_radius

            # 6. Limite a força de tração pela aderência máxima do pneu
            traction_force_limitada = np.sign(traction_force_ideal) * min(abs(traction_force_ideal), abs(Fx_max))

            # 7. Calcule a aceleração real com a condição física correta
            resistance_force = self.vehicle.calculate_resistance_forces(vehicle_velocity)
            forca_resultante = traction_force_limitada - resistance_force

            # CONDIÇÃO FÍSICA: Se o carro está praticamente parado e a força resultante é negativa,
            # a aceleração real é zero (ele não anda para trás).
            # Usamos uma pequena tolerância (ex: 0.01 m/s) para segurança com números de ponto flutuante.
            if vehicle_velocity < 0.01 and forca_resultante < 0:
                vehicle_acceleration = 0.0
            else:
                vehicle_acceleration = forca_resultante / self.vehicle.mass


        speed_error = self.speed_ref - wm
        torque_ref_unsat = self.speed_controller.kp * speed_error + self.speed_controller.ki * int_speed
        int_speed_dot = speed_error
        torque_ref = torque_ref_unsat

        iq_ref = torque_ref / self.torque_constant if self.torque_constant != 0 else 0.0
        id_ref = 0.0
        if abs(iq_ref) > self.max_current:
            iq_ref = np.sign(iq_ref) * self.max_current

        error_d = id_ref - isd
        error_q = iq_ref - isq

        vd_unclamped = self.id_controller.kp * error_d + self.id_controller.ki * int_id + (-we * self.lq * isq)
        vq_unclamped = self.iq_controller.kp * error_q + self.iq_controller.ki * int_iq + (we * (self.ld * isd + self.lambda_m))

        # Atualiza a tensão DC do inversor com a tensão da bateria
        if self.battery:
            self.Vdc1 = 1500#self.battery.calcular_tensao(-isq, soc, Iast, tempo_acumulado)
            self.inversor.set_Vdc(self.Vdc1)
        else:
            self.Vdc = 1600

        vd = np.clip(vd_unclamped, -self.Vdc1, self.Vdc1)
        vq = np.clip(vq_unclamped, -self.Vdc1, self.Vdc1)

        int_id_dot = 0.0 if abs(vd_unclamped) > self.Vlimit else error_d
        int_iq_dot = 0.0 if abs(vq_unclamped) > self.Vlimit else error_q

        d_isd = (vd - self.rs * isd + we * self.lq * isq) * self.inv_ld
        d_isq = (vq - self.rs * isq - we * (self.ld * isd + self.lambda_m)) * self.inv_lq
        L0 = 0.1 * (self.ld + self.lq) / 2.0
        d_iso = (0.0 - self.rs * iso) / L0

        ce = self.torque_constant * isq

        d_wm = (ce - cl - self.kf * wm) * self.inv_jm
        d_theta_m = wm

        i_rms_sq = isd**2 + isq**2
        copper_losses = 3.0 * self.rs * (i_rms_sq / 2)
        d_temp = copper_losses * self.inv_mC

        if abs(iq_ref) >= self.max_current and np.sign(int_speed) == np.sign(speed_error):
            int_speed_dot = 0.0

        # Derivadas da bateria
        if self.battery:
            corrente_battery = abs(isq)
            dsoc_dt, dIast_dt, dtime_dt = self.battery.calcular_derivadas(corrente_battery)
        else:
            dsoc_dt, dIast_dt, dtime_dt = 0.0, 0.0, 0.0

        # Derivadas dos estados do veículo
        if self.vehicle and self.transmission:
            d_vehicle_velocity = vehicle_acceleration
            d_vehicle_position = vehicle_velocity

            # Montar vetor de derivadas completo (motor + veículo + bateria)
            if self.battery:
                dxdt = np.array([
                    d_isd, d_isq, d_iso,
                    d_wm, d_theta_m, d_temp,
                    int_id_dot, int_iq_dot, int_speed_dot,
                    d_vehicle_velocity, d_vehicle_position,
                    dsoc_dt, dIast_dt, dtime_dt
                ], dtype=float)
            else:
                dxdt = np.array([
                    d_isd, d_isq, d_iso,
                    d_wm, d_theta_m, d_temp,
                    int_id_dot, int_iq_dot, int_speed_dot,
                    d_vehicle_velocity, d_vehicle_position
                ], dtype=float)
        else:
            # Montar vetor de derivadas (motor + bateria)
            if self.battery:
                dxdt = np.array([
                    d_isd, d_isq, d_iso,
                    d_wm, d_theta_m, d_temp,
                    int_id_dot, int_iq_dot, int_speed_dot,
                    dsoc_dt, dIast_dt, dtime_dt
                ], dtype=float)
            else:
                dxdt = np.array([
                    d_isd, d_isq, d_iso,
                    d_wm, d_theta_m, d_temp,
                    int_id_dot, int_iq_dot, int_speed_dot
                ], dtype=float)

        return dxdt

    def simulate(self, t0=0.0, tf=None):
        if tf is None:
            tf = self.tmax

        # Garantir que t_eval esteja dentro do t_span
        num_points = int((tf - t0) / self.hp) + 1
        t_eval = np.linspace(t0, tf, num_points)

        # Resto do código permanece igual
        # Configuração inicial dos estados
        if self.vehicle and self.transmission:
            if self.battery:
                x0 = np.array([
                    self.isd0, self.isq0, self.iso0,
                    self.wm0, self.theta_m0, self.temp0,
                    self.int_id0, self.int_iq0, self.int_speed0,
                    self.vehicle_velocity0, self.vehicle_position0,
                    self.soc0, self.Iast0, self.tempo_acumulado0
                ], dtype=float)
            else:
                x0 = np.array([
                    self.isd0, self.isq0, self.iso0,
                    self.wm0, self.theta_m0, self.temp0,
                    self.int_id0, self.int_iq0, self.int_speed0,
                    self.vehicle_velocity0, self.vehicle_position0
                ], dtype=float)
        else:
            if self.battery:
                x0 = np.array([
                    self.isd0, self.isq0, self.iso0,
                    self.wm0, self.theta_m0, self.temp0,
                    self.int_id0, self.int_iq0, self.int_speed0,
                    self.soc0, self.Iast0, self.tempo_acumulado0
                ], dtype=float)
            else:
                x0 = np.array([
                    self.isd0, self.isq0, self.iso0,
                    self.wm0, self.theta_m0, self.temp0,
                    self.int_id0, self.int_iq0, self.int_speed0
                ], dtype=float)

        self.initialize_storage()

        self.id_controller.reset()
        self.iq_controller.reset()
        self.speed_controller.reset()

        try:
            # Usa RK45 para resolver as equações diferenciais combinadas
            sol = solve_ivp(
                fun=self.combined_edos,
                t_span=(t0, tf),
                y0=x0,
                method='RK45',
                t_eval=t_eval,
                atol=1e-6,
                rtol=1e-6
            )
        except Exception as e:
            print(f"Erro na simulação: {e}")
            sol = type('obj', (object,), {
                't': np.array([t0, tf]),
                'y': np.zeros((len(x0), 2))
            })

        dt = self.hp

        for idx, t in enumerate(sol.t):
            # Estados do motor
            isd = sol.y[0, idx]
            isq = sol.y[1, idx]
            iso = sol.y[2, idx]
            wm = sol.y[3, idx]
            theta_m = sol.y[4, idx]
            temp = sol.y[5, idx]
            int_id = sol.y[6, idx]
            int_iq = sol.y[7, idx]
            int_speed = sol.y[8, idx]

            # Estados do veículo (se existirem)
            if self.vehicle and self.transmission:
                vehicle_velocity = sol.y[9, idx]
                vehicle_position = sol.y[10, idx]
                # Estados da bateria (se existir)
                if self.battery:
                    soc = max(sol.y[11, idx], 0.000)
                    Iast = sol.y[12, idx]
                    tempo_acumulado = sol.y[13, idx]
                else:
                    soc, Iast, tempo_acumulado = 1.0, 0.0, 0.0
            else:
                vehicle_velocity, vehicle_position = 0.0, 0.0
                # Estados da bateria (se existir)
                if self.battery:
                    soc = max(sol.y[9, idx], 0.0000)
                    Iast = sol.y[10, idx]
                    tempo_acumulado = sol.y[11, idx]
                else:
                    soc, Iast, tempo_acumulado = 1.0, 0.0, 0.0

            theta_e = self.p * theta_m
            we = self.p * wm

            vd_control, vq_control, id_ref, iq_ref, speed_error = self.field_oriented_control(isd, isq, wm, dt)

            flux_d = self.ld * isd + self.lambda_m
            flux_q = self.lq * isq
            flux_o = 0.1 * (self.ld + self.lq) / 2.0 * iso

            ce = self.torque_constant * isq

            # Cálculo do torque de carga para armazenamento
            if self.vehicle and self.transmission:
                cl = self.vehicle.calculate_load_torque(vehicle_velocity, self.transmission)
            else:
                cl = self.set_load(t)

            cm = ce - cl

            # Integração com a bateria e inversor
            if self.battery:
                # Calcular corrente total estimada
                corrente_battery = abs(isq)

                # Calcular tensão da bateria
                battery_voltage = self.battery.calcular_tensao(corrente_battery, soc, Iast, tempo_acumulado) if soc >= 0.0 else 0.0

                # Configurar tensão DC no inversor
                self.inversor.set_Vdc(battery_voltage)
                self.Vdc = battery_voltage

                # Configurar velocidade síncrona no inversor
                self.inversor.set_ws(we)

                # Gerar tensões de fase através do inversor
                va, vb, vc = self.inversor.source_voltage()

                vq_r, vd_r = self.inversor.park_transform(va, vb, vc, theta_e)

            else:
                # Sem bateria, usar valor padrão
                self.inversor.set_Vdc(self.Vdc)
                self.inversor.set_ws(we)
                battery_voltage = self.Vdc
                corrente_battery = 0.0
                soc, Iast, tempo_acumulado = 1.0, 0.0, 0.0
                va, vb, vc = 0.0, 0.0, 0.0
                vq_r, vd_r = 0.0, 0.0

            # Calcular correntes de fase
            is1, is2, is3, fs1, fs2, fs3 = self.abc_currents_from_dq(isd, isq, theta_e, flux_d)

            # Armazenar dados
            self.tempo.append(t)
            self.corrented.append(isd)
            self.correnteq.append(isq)
            self.corrente1.append(is1)
            self.corrente2.append(is2)
            self.corrente3.append(is3)
            self.tensaosd.append(vd_control)
            self.tensaosq.append(vq_control)
            self.tensao1.append(va)
            self.tensao2.append(vb)
            self.tensao3.append(vc)
            self.fluxosd.append(flux_d)
            self.fluxosq.append(flux_q)
            self.conjugado.append(ce)
            self.velocidade.append(wm * 60 / (2 * PI))  # Convert rad/s to RPM
            self.torque_mecanico.append(cm)
            self.conjcarga.append(cl)
            self.temperatura.append(temp)
            self.vd_control.append(vd_control)
            self.vq_real.append(vq_r)
            self.vd_real.append(vd_r)
            self.vq_control.append(vq_control)
            self.speed_error.append(speed_error)
            self.iq_ref_trace.append(iq_ref)
            self.id_ref_trace.append(id_ref)

            # Armazena dados do veículo
            self.vehicle_velocity.append(vehicle_velocity)
            self.vehicle_position.append(vehicle_position)

            # Armazena dados da bateria
            if self.battery:
                self.soc_hist.append(soc)
                self.Iast_hist.append(Iast)
                self.tempo_acumulado_hist.append(tempo_acumulado)
                self.battery_voltage_hist.append(battery_voltage)
                self.battery_current_hist.append(corrente_battery)

            # Calcula e armazena aceleração do veículo
            if self.vehicle and self.transmission and self.tire:
                wheel_torque = self.transmission.motor_to_wheel_torque(ce)
                
                # Recalcula as condições para obter a aceleração correta para salvar
                velocidade_angular_roda = self.transmission.motor_to_wheel_speed(wm)
                slip = Tire.SlipRatio(velocidade_angular_roda, self.vehicle.wheel_radius, vehicle_velocity)
                Fz = self.vehicle.mass * self.vehicle.g
                Fx_max = self.tire.Tire_forces(Fz, slip)

                traction_force_ideal = wheel_torque / self.vehicle.wheel_radius
                traction_force_limitada = np.sign(traction_force_ideal) * min(abs(traction_force_ideal), abs(Fx_max))
                resistance_force = self.vehicle.calculate_resistance_forces(vehicle_velocity)
                forca_resultante = traction_force_limitada - resistance_force
                
                if vehicle_velocity < 0.01 and forca_resultante < 0:
                    vehicle_acceleration = 0.0
                else:
                    vehicle_acceleration = forca_resultante / self.vehicle.mass

                self.vehicle_acceleration.append(vehicle_acceleration)
                self.wheel_torque.append(wheel_torque) # Armazena o torque ideal na roda
                self.tractive_force_hist.append(traction_force_limitada)
                self.resistive_force_hist.append(resistance_force)
                self.slip_ratio_hist.append(slip)
                self.longitudinal_force_hist.append(Fx_max)

            else:
                self.vehicle_acceleration.append(0)
                self.wheel_torque.append(0)
                self.tractive_force_hist.append(0)
                self.resistive_force_hist.append(0)
                self.slip_ratio_hist.append(0)
                self.longitudinal_force_hist.append(0)

# ---------------------
# Instâncias iniciais (padrão)
# ---------------------
transmission_default = Transmission(
    final_drive_ratio=4.0,
    efficiency=0.95
)

vehicle_default = Vehicle(
    mass=230,  # kg (exemplo)
    wheel_radius=0.16,  # m
    drag_coeff=0.3,
    frontal_area=1,  # m²
    rolling_resistance=0.015,
    road_grade=0  # radianos (0 = plano)
)

battery_default = BatteryPack(tipo_celula='Li-ion', n_serie=162, n_paralelo=1, soc_inicial=1.0)

# ---------------------
# Dash app com layout modificado
# ---------------------
external_stylesheets = [
    'https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css',
    'https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap',
    'https://codepen.io/chriddyp/pen/bWLwgP.css',
    dbc.themes.BOOTSTRAP
]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

# Layout principal com sidebar
app.layout = html.Div([
    dcc.Store(id='simulation-data'),
    dcc.Store(id='parameters-store', data={
        'vehicle': {'mass': 300, 'wheel_radius': 0.16, 'drag_coeff': 0.3,
                    'frontal_area': 1, 'rolling_resistance': 0.015, 'road_grade': 0},
        'transmission': {'final_drive_ratio': 4.0, 'efficiency': 0.95},
        'battery': {'tipo_celula': 'Li-ion', 'n_serie': 162, 'n_paralelo': 1, 'soc_inicial': 1.0},
        'motor': {'rs': 0.04585, 'ld': 0.00067, 'lq': 0.00067, 'jm': 0.05769,
                  'kf': 0.1, 'lambda_m': 0.13849, 'p': 10, 'valor_mu': 0.95},
        'inversor': {'eficiencia': 0.95, 'freq_chaveamento': 10000},
        'simulacao': {'tmax': 11, 'velocidade_ref': 471.23},
        'tire': {
            'pacejka_params': [0.333, 1.627, 1, 4.396, 931.4, 366.4],
            'tire_friction_coef': 1.45
        }
    }),
    dcc.Store(id='sidebar-state', data={'visible': True}),
    dcc.Store(id='stored-figures', data={}),
    dcc.Store(id='simulation-run-flag', data={'run': False}),

    html.Div([
        # Botão vertical fino
        html.Button(
            id='sidebar-toggle',
            children=html.I(className='fas fa-chevron-right'),
            style={
                'position': 'fixed',
                'top': '10px',
                'left': '0px',
                'zIndex': 1000,
                'width': '40px',
                'height': '40px',
                'borderRadius': '0 4px 4px 0',
                'backgroundColor': "#830B0B",
                'color': 'white',
                'border': 'none',
                'cursor': 'pointer',
                'padding': '0',
                'display': 'flex',
                'alignItems': 'center',
                'justifyContent': 'center',
                'fontSize': '12px',
                'boxShadow': '1px 0 2px rgba(0,0,0,0.2)',
                'transition': 'all 0.2s ease'
            }
        ),

        # Sidebar com controles
        html.Div([
            html.H3("Parâmetros de Simulação", style={'marginLeft': '30px', 'marginBottom': '20px'}),

            # Controles do Veículo
            dbc.Card([
                dbc.CardHeader(html.H5("Veículo")),
                dbc.CardBody([
                    dbc.Label("Massa (kg)"),
                    dbc.Input(id='vehicle-mass', type='number', value=230, step=1),
                    dbc.Label("Raio da Roda (m)"),
                    dbc.Input(id='wheel-radius', type='number', value=0.16, step=0.01),
                    dbc.Label("Coef. de Arrasto"),
                    dbc.Input(id='drag-coeff', type='number', value=0.78, step=0.01),
                    dbc.Label("Área Frontal (m²)"),
                    dbc.Input(id='frontal-area', type='number', value=0.68, step=0.1),
                    dbc.Label("Resistência Rolamento"),
                    dbc.Input(id='rolling-resistance', type='number', value=0.015, step=0.001),
                ])
            ], style={'marginBottom': '15px'}),

            # Controles da Transmissão
            dbc.Card([
                dbc.CardHeader(html.H5("Transmissão")),
                dbc.CardBody([
                    dbc.Label("Relação Diferencial"),
                    dbc.Input(id='final-drive-ratio', type='number', value=4.0, step=0.1),
                    dbc.Label("Eficiência"),
                    dbc.Input(id='transmission-efficiency', type='number', value=0.95, step=0.01, min=0, max=1),
                ])
            ], style={'marginBottom': '15px'}),

            # Controles da Bateria
            dbc.Card([
                dbc.CardHeader(html.H5("Bateria")),
                dbc.CardBody([
                    dbc.Label("Tipo de Célula"),
                    dbc.Select(
                        id='battery-type',
                        options=[{'label': 'Li-ion', 'value': 'Li-ion'},
                                 {'label': 'LiFePO4', 'value': 'LiFePO4'}],
                        value='Li-ion'
                    ),
                    dbc.Label("Células em Série"),
                    dbc.Input(id='n-serie', type='number', value=162, step=1),
                    dbc.Label("Células em Paralelo"),
                    dbc.Input(id='n-paralelo', type='number', value=1, step=1),
                    dbc.Label("SoC Inicial"),
                    dbc.Input(id='soc-inicial', type='number', value=1.0, step=0.01, min=0, max=1),
                ])
            ], style={'marginBottom': '15px'}),

            # Controles do Motor
            dbc.Card([
                dbc.CardHeader(html.H5("Motor")),
                dbc.CardBody([
                    dbc.Label("Resistência (Ω)"),
                    dbc.Input(id='motor-rs', type='number', value=0.04585, step=0.0001),
                    dbc.Label("Indutância d (H)"),
                    dbc.Input(id='motor-ld', type='number', value=0.00067, step=0.00001),
                    dbc.Label("Indutância q (H)"),
                    dbc.Input(id='motor-lq', type='number', value=0.00067, step=0.00001),
                    dbc.Label("Inércia (kg.m²)"),
                    dbc.Input(id='motor-jm', type='number', value=0.05769, step=0.0001),
                    dbc.Label("Atrito (N.m.s)"),
                    dbc.Input(id='motor-kf', type='number', value=0.1, step=0.01),
                    dbc.Label("Fluxo (Wb)"),
                    dbc.Input(id='motor-lambda', type='number', value=0.13849, step=0.0001),
                    dbc.Label("Polos"),
                    dbc.Input(id='motor-poles', type='number', value=10, step=1),
                    dbc.Label("Modulação"),
                    dbc.Input(id='motor-modulation', type='number', value=0.95, step=0.01, min=0, max=1),
                ])
            ], style={'marginBottom': '15px'}),

            # Controles de Simulação
            dbc.Card([
                dbc.CardHeader(html.H5("Simulação")),
                dbc.CardBody([
                    dbc.Label("Tempo Máximo (s)"),
                    dbc.Input(id='sim-time', type='number', value=11, step=1),
                    dbc.Label("Velocidade Ref (rad/s)"),
                    dbc.Input(id='speed-ref', type='number', value=471.23, step=1),
                ])
            ], style={'marginBottom': '15px'}),

            html.Button("Executar Simulação", id="run-simulation", n_clicks=0,
                        className="btn btn-primary", style={'width': '100%', 'marginBottom': '10px'}),
            html.Button("Aplicar Parâmetros", id="apply-params", n_clicks=0,
                        className="btn btn-secondary", style={'width': '100%'}),

        ], id='sidebar', style={'width': '300px', 'padding': '15px', 'backgroundColor': '#f8f9fa',
                               'height': '100vh', 'overflowY': 'auto', 'position': 'fixed', 'transition': 'left 0.3s ease'}),

        # Área principal com botões e gráfico
        html.Div([
            html.H1("Digital-Twin de Powertrain Xangô-E.Racing",
                    style={'textAlign': 'center', 'color': '#2c3e50', 'marginBottom': '20px'}),

            # Botões de visualização
            html.Div([
                dbc.Button([html.I(className="fas fa-tachometer-alt"), " Velocidade e Torque"],
                           id='btn-velocity-torque', color="#830B0B", className="me-1 mb-1"),
                dbc.Button([html.I(className="fas fa-bolt"), " Correntes"],
                           id='btn-currents', color="#830B0B", className="me-1 mb-1"),
                dbc.Button([html.I(className="fas fa-plug"), " Tensões"],
                           id='btn-voltages', color="#830B0B", className="me-1 mb-1"),
                dbc.Button([html.I(className="fas fa-magnet"), " Fluxos e Temperatura"],
                           id='btn-flux-temp', color="#830B0B", className="me-1 mb-1"),
                dbc.Button([html.I(className="fas fa-sliders-h"), " Sinais de Controle"],
                           id='btn-control', color="#830B0B", className="me-1 mb-1"),
                dbc.Button([html.I(className="fas fa-chart-bar"), " Visão Completa"],
                           id='btn-complete', color="#830B0B", className="me-1 mb-1"),
                dbc.Button([html.I(className="fas fa-car"), " Desempenho do Veículo"],
                           id='btn-vehicle', color="#830B0B", className="me-1 mb-1"),
                dbc.Button([html.I(className="fas fa-dot-circle"), " Desempenho do Pneu"],
                           id='btn-tire', color="#830B0B", className="me-1 mb-1"),
                dbc.Button([html.I(className="fas fa-battery-three-quarters"), " Bateria"],
                           id='btn-battery', color="#830B0B", className="me-1 mb-1"),
            ], style={'textAlign': 'center', 'marginBottom': '15px'}),

            # Gráfico
            dcc.Graph(id='motor-graph', style={'height': '120vh', 'borderRadius': '8px',
                                              'boxShadow': '0 4px 6px rgba(0,0,0,0.1)'}),

            # Loading component
            dcc.Loading(
                id="loading-simulation",
                type="dot",
                color="#9B0909",
                children=html.Div(id="loading-output"),
                style={'position': 'fixed', 'top': '50%', 'left': '50%', 'transform': 'translate(-50%, -50%)',
                       'zIndex': 2000, 'display': 'none'}  # Inicialmente oculto
            ),

        ], id='main-content', style={'transition': 'margin-left 0.3s ease'}),
    ])
])

# ---------------------
# Funções de plotagem (definidas globalmente)
# ---------------------
def create_velocity_torque_plot(motor):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Velocidade Mecânica", "Torque do Motor"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.velocidade, name='Velocidade'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.conjugado, name='Torque Elétrico'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.conjcarga, name='Torque de Carga', line=dict(dash='dash')), row=2, col=1)
    fig.update_layout(height=900, title_text="Velocidade e Torque", template="plotly_white")
    fig.update_xaxes(title_text="Tempo (s)", row=2, col=1)
    fig.update_yaxes(title_text="Velocidade (RPM)", row=1, col=1)
    fig.update_yaxes(title_text="Torque (Nm)", row=2, col=1)
    return fig

def create_currents_plot(motor):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Correntes dq", "Correntes de Fase"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrented, name='Id'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.correnteq, name='Iq'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente1, name='Fase 1'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente2, name='Fase 2'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente3, name='Fase 3'), row=2, col=1)
    fig.update_layout(height=900, title_text="Correntes", template="plotly_white")
    return fig

def create_voltages_plot(motor):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Tensões de Controle dq", "Tensões de Fase"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensaosd, name='Vd'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vd_real, name='Vd real'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensaosq, name='Vq'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vq_real, name='Vq real'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao1, name='Fase 1'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao2, name='Fase 2'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao3, name='Fase 3'), row=2, col=1)
    fig.update_layout(height=900, title_text="Tensões", template="plotly_white")
    return fig

def create_flux_temp_plot(motor):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Fluxos Magnéticos", "Temperatura do Motor"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.fluxosd, name='Fluxo d'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.fluxosq, name='Fluxo q'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.temperatura, name='Temperatura'), row=2, col=1)
    fig.update_layout(height=900, title_text="Fluxos e Temperatura", template="plotly_white")
    return fig

def create_control_plot(motor):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Sinais de Controle FOC", "Erro de Velocidade"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vd_control, name='Vd control'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vq_control, name='Vq control'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.speed_error, name='Erro de Velocidade'), row=2, col=1)
    fig.update_layout(height=900, title_text="Sinais de Controle", template="plotly_white")
    return fig

def create_complete_plot(motor):
    fig = make_subplots(rows=3, cols=2, subplot_titles=(
        "Velocidade Mecânica", "Torque do Motor", "Correntes", "Correntes de Fase", "Tensões de Controle", "Tensões de Fase"
    ), vertical_spacing=0.08, horizontal_spacing=0.1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.velocidade, name='Velocidade'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.conjugado, name='Torque Elétrico'), row=1, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.conjcarga, name='Torque de Carga', line=dict(dash='dash')), row=1, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrented, name='Id'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.correnteq, name='Iq'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente1, name='Fase 1'), row=2, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente2, name='Fase 2'), row=2, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente3, name='Fase 3'), row=2, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensaosd, name='Vd'), row=3, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vd_real, name='Vd real'), row=3, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensaosq, name='Vq'), row=3, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vq_real, name='Vq real'), row=3, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao1, name='Fase 1'), row=3, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao2, name='Fase 2'), row=3, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao3, name='Fase 3'), row=3, col=2)
    fig.update_layout(height=900, title_text="Visão Completa da Simulação", template="plotly_white", hovermode="x unified", showlegend=False)
    return fig

def create_vehicle_plot(motor):
    fig = make_subplots(
        rows=2, cols=2, 
        subplot_titles=("Velocidade (km/h)", "Aceleração (m/s²)", "Forças (N)", "Torque na Roda (Nm)"), 
        vertical_spacing=0.15, horizontal_spacing=0.1
    )
    velocity_kmh = [v * 3.6 for v in motor.vehicle_velocity]
    fig.add_trace(go.Scatter(x=motor.tempo, y=velocity_kmh, name='Velocidade'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vehicle_acceleration, name='Aceleração'), row=1, col=2)
    
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tractive_force_hist, name='Força Trativa'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.resistive_force_hist, name='Forças Resistivas'), row=2, col=1)
    
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.wheel_torque, name='Torque na Roda'), row=2, col=2)
    
    fig.update_layout(height=900, title_text="Desempenho do Veículo", template="plotly_white", legend_tracegroupgap=180)
    fig.update_yaxes(title_text="Força (N)", row=2, col=1)
    return fig

def create_tire_plot(motor):
    fig = make_subplots(
        rows=1, cols=2, 
        subplot_titles=("Slip Ratio x Tempo", "Força Longitudinal x Slip Ratio")
    )

    fig.add_trace(
        go.Scatter(x=motor.tempo, y=motor.slip_ratio_hist, name='Slip Ratio', line=dict(color='blue')),
        row=1, col=1
    )
    fig.update_xaxes(title_text="Tempo [s]", row=1, col=1)
    fig.update_yaxes(title_text="Slip Ratio", row=1, col=1)

    if len(motor.slip_ratio_hist) > 0:
        df_tire = pd.DataFrame({
            'slip': motor.slip_ratio_hist,
            'fx': motor.longitudinal_force_hist
        })
        df_tire_sorted = df_tire.sort_values(by='slip').drop_duplicates(subset='slip')

        fig.add_trace(
            go.Scatter(
                x=df_tire_sorted['slip'], 
                y=df_tire_sorted['fx'], 
                name='Força Longitudinal', 
                mode='lines', 
                line=dict(color='green')
            ),
            row=1, col=2
        )
    
    fig.update_xaxes(title_text="Slip Ratio", row=1, col=2)
    fig.update_yaxes(title_text="Força Longitudinal (N)", row=1, col=2)

    fig.update_layout(height=500, title_text="Análise de Desempenho do Pneu", template="plotly_white")
    
    return fig

def create_battery_plot(motor):
    fig = make_subplots(rows=3, cols=1, subplot_titles=("Tensão do Banco (V)", "Corrente (A)", "SoC"), vertical_spacing=0.12)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.battery_voltage_hist, name='Tensão Banco'), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.battery_current_hist, name='Corrente'), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.soc_hist, name='SoC'), row=3, col=1)
    fig.update_layout(height=900, title_text="Estado da Bateria", template="plotly_white")
    fig.update_yaxes(title_text="Tensão (V)", row=1, col=1)
    fig.update_yaxes(title_text="Corrente (A)", row=2, col=1)
    fig.update_yaxes(title_text="SoC (-)", row=3, col=1)
    return fig

# ---------------------
# Callback para mostrar/esconder la sidebar
# ---------------------
@app.callback(
    [Output('sidebar', 'style'),
     Output('main-content', 'style'),
     Output('sidebar-toggle', 'children'),
     Output('sidebar-state', 'data')],
    [Input('sidebar-toggle', 'n_clicks')],
    [State('sidebar-state', 'data')]
)
def toggle_sidebar(n_clicks, sidebar_state):
    if n_clicks is None:
        sidebar_style = {'width': '300px', 'padding': '15px', 'backgroundColor': '#f8f9fa',
                         'height': '100vh', 'overflowY': 'auto', 'position': 'fixed',
                         'transition': 'left 0.3s ease', 'left': '0'}
        content_style = {'transition': 'margin-left 0.3s ease', 'marginLeft': '300px'}
        icon = html.I(className='fas fa-chevron-right')
        return sidebar_style, content_style, icon, {'visible': True}

    visible = not sidebar_state.get('visible', True)

    if visible:
        sidebar_style = {'width': '300px', 'padding': '15px', 'backgroundColor': '#f8f9fa',
                         'height': '100vh', 'overflowY': 'auto', 'position': 'fixed',
                         'transition': 'left 0.3s ease', 'left': '0'}
        content_style = {'transition': 'margin-left 0.3s ease', 'marginLeft': '300px'}
        icon = html.I(className='fas fa-chevron-right')
    else:
        sidebar_style = {'width': '300px', 'padding': '15px', 'backgroundColor': '#f8f9fa',
                         'height': '100vh', 'overflowY': 'auto', 'position': 'fixed',
                         'transition': 'left 0.3s ease', 'left': '-300px'}
        content_style = {'transition': 'margin-left 0.3s ease', 'marginLeft': '0'}
        icon = html.I(className='fas fa-chevron-right', style={'transform': 'rotate(180deg)'})

    return sidebar_style, content_style, icon, {'visible': visible}

# ---------------------
# Callbacks para atualizar parâmetros
# ---------------------
@app.callback(
    Output('parameters-store', 'data'),
    [Input('apply-params', 'n_clicks')],
    [State('vehicle-mass', 'value'),
     State('wheel-radius', 'value'),
     State('drag-coeff', 'value'),
     State('frontal-area', 'value'),
     State('rolling-resistance', 'value'),
     State('final-drive-ratio', 'value'),
     State('transmission-efficiency', 'value'),
     State('battery-type', 'value'),
     State('n-serie', 'value'),
     State('n-paralelo', 'value'),
     State('soc-inicial', 'value'),
     State('motor-rs', 'value'),
     State('motor-ld', 'value'),
     State('motor-lq', 'value'),
     State('motor-jm', 'value'),
     State('motor-kf', 'value'),
     State('motor-lambda', 'value'),
     State('motor-poles', 'value'),
     State('motor-modulation', 'value'),
     State('sim-time', 'value'),
     State('speed-ref', 'value'),
     State('parameters-store', 'data')]
)
def update_parameters(n_clicks, mass, wheel_radius, drag_coeff, frontal_area, rolling_resistance,
                      final_drive_ratio, transmission_efficiency,
                      battery_type, n_serie, n_paralelo, soc_inicial,
                      motor_rs, motor_ld, motor_lq, motor_jm, motor_kf, motor_lambda, motor_poles, motor_modulation,
                      sim_time, speed_ref, current_params):
    if n_clicks is None or n_clicks == 0:
        return current_params

    current_params['vehicle'] = {
        'mass': mass,
        'wheel_radius': wheel_radius,
        'drag_coeff': drag_coeff,
        'frontal_area': frontal_area,
        'rolling_resistance': rolling_resistance,
        'road_grade': 0
    }

    current_params['transmission'] = {
        'final_drive_ratio': final_drive_ratio,
        'efficiency': transmission_efficiency
    }

    current_params['battery'] = {
        'tipo_celula': battery_type,
        'n_serie': n_serie,
        'n_paralelo': n_paralelo,
        'soc_inicial': soc_inicial
    }

    current_params['motor'] = {
        'rs': motor_rs,
        'ld': motor_ld,
        'lq': motor_lq,
        'jm': motor_jm,
        'kf': motor_kf,
        'lambda_m': motor_lambda,
        'p': motor_poles,
        'valor_mu': motor_modulation
    }

    current_params['simulacao'] = {
        'tmax': sim_time,
        'velocidade_ref': speed_ref
    }

    return current_params

# ---------------------
# Callback para executar simulação e armazenar figuras
# ---------------------
@app.callback(
    [Output('stored-figures', 'data'),
     Output('simulation-run-flag', 'data'),
     Output('loading-output', 'children')],
    [Input('run-simulation', 'n_clicks')],
    [State('parameters-store', 'data')]
)
def run_simulation_once(n_clicks, parameters):
    if n_clicks is None or n_clicks == 0:
        return dash.no_update, {'run': False}, ""
    
    # Extrair parâmetros e executar simulação
    vehicle_params = parameters['vehicle']
    transmission_params = parameters['transmission']
    battery_params = parameters['battery']
    motor_params = parameters['motor']
    simulacao_params = parameters['simulacao']
    tire_params = parameters['tire']

    transmission = Transmission(**transmission_params)
    vehicle = Vehicle(**vehicle_params)
    battery = BatteryPack(**battery_params)
    tire = Tire(**tire_params)

    motor = Motor(
        rs=motor_params['rs'],
        ld=motor_params['ld'],
        lq=motor_params['lq'],
        jm=motor_params['jm'],
        kf=motor_params['kf'],
        lambda_m=motor_params['lambda_m'],
        p=motor_params['p'],
        valor_mu=motor_params['valor_mu'],
        TL=False,
        torque=0.0,
        vehicle=vehicle,
        transmission=transmission,
        battery=battery,
        tire=tire
    )
    
    motor.tmax = simulacao_params['tmax']
    motor.speed_ref = simulacao_params['velocidade_ref']
    
    motor.simulate()
    
    # Gerar todas as figuras
    figures = {
        'velocity_torque': create_velocity_torque_plot(motor),
        'currents': create_currents_plot(motor),
        'voltages': create_voltages_plot(motor),
        'flux_temp': create_flux_temp_plot(motor),
        'control': create_control_plot(motor),
        'complete': create_complete_plot(motor),
        'vehicle': create_vehicle_plot(motor),
        'tire': create_tire_plot(motor),
        'battery': create_battery_plot(motor)
    }
    
    stored_figures = {}
    for key, fig in figures.items():
        stored_figures[key] = fig.to_dict()
    
    return stored_figures, {'run': True}, ""

# ---------------------
# Callback para atualizar gráfico
# ---------------------
@app.callback(
    Output('motor-graph', 'figure'),
    [Input('btn-velocity-torque', 'n_clicks'),
     Input('btn-currents', 'n_clicks'),
     Input('btn-voltages', 'n_clicks'),
     Input('btn-flux-temp', 'n_clicks'),
     Input('btn-control', 'n_clicks'),
     Input('btn-complete', 'n_clicks'),
     Input('btn-vehicle', 'n_clicks'),
     Input('btn-tire', 'n_clicks'),
     Input('btn-battery', 'n_clicks')],
    [State('stored-figures', 'data')]
)
def update_graph_from_stored(btn_velocity, btn_currents, btn_voltages, btn_flux, btn_control,
                            btn_complete, btn_vehicle, btn_tire, btn_battery, stored_figures):
    ctx = dash.callback_context
    if not ctx.triggered or not stored_figures:
        return {}
    
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    figure_map = {
        'btn-velocity-torque': 'velocity_torque',
        'btn-currents': 'currents',
        'btn-voltages': 'voltages',
        'btn-flux-temp': 'flux_temp',
        'btn-control': 'control',
        'btn-complete': 'complete',
        'btn-vehicle': 'vehicle',
        'btn-tire': 'tire',
        'btn-battery': 'battery'
    }
    
    figure_key = figure_map.get(button_id)
    if figure_key and figure_key in stored_figures:
        return stored_figures[figure_key]
    
    return {}

if __name__ == '__main__':
    app.run(debug=True)