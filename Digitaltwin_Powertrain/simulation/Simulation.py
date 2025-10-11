import sys
import os

# Adiciona o diretório pai ao path do Python
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.integrate import solve_ivp
from models.BatteryPack import BatteryPack
from models import Inversor
from models import PIDController as PID
from models.Tire import Tire
from models.Transmission import Transmission
from models.Vehicle import Vehicle
from models.Motor import Motor
from typing import Optional, Dict, Any


class Simulation:
    """Vehicle + motor + battery + inverter simulation orchestrator using PID controllers."""

    def __init__(self,
                 motor: Optional[Motor] = None,
                 vehicle: Optional[Vehicle] = None,
                 transmission: Optional[Transmission] = None,
                 battery: Optional[BatteryPack] = None,
                 tire: Optional[Tire] = None,
                 inversor: Optional[Inversor] = None,
                 tmax: float = 10.0,
                 steps: int = 2000):
        
        # store models
        self.motor = motor
        self.vehicle = vehicle
        self.transmission = transmission
        self.battery = battery
        self.tire = tire
        self.inversor = inversor if inversor is not None else Inversor()

        # simulation time configuration - REDUZIDO PARA TESTE
        self.tmax = float(tmax)
        self.steps = int(steps)
        # conservative controller update timestep used inside combined_edos
        self.hp = self.tmax / float(self.steps) if self.steps > 0 else 1e-3

        # extract motor parameters safely with defaults
        self.p = getattr(self.motor, 'p', 1)
        self.ld = getattr(self.motor, 'ld', 1.0)
        self.lq = getattr(self.motor, 'lq', 1.0)
        self.rs = getattr(self.motor, 'rs', 0.0)
        self.lambda_m = getattr(self.motor, 'lambda_m', 0.0)
        self.jm = getattr(self.motor, 'jm', 1.0)
        self.kf = getattr(self.motor, 'kf', 0.0)
        
        # CORREÇÃO: torque_constant calculado corretamente
        self.torque_constant = (3.0/2.0) * self.p * self.lambda_m
        
        self.max_current = getattr(self.motor, 'max_current', np.inf)
        self.Vdc_default = getattr(self.motor, 'Vdc', 600.0)
        self.modulation_index = getattr(self.motor, 'valor_mu', 1.0)
        self.speed_ref = getattr(self.motor, 'speed_ref', 471.23)  # rad/s

        # precompute inverses used in ODEs
        self.inv_ld = 1.0 / self.ld if self.ld != 0 else 0.0
        self.inv_lq = 1.0 / self.lq if self.lq != 0 else 0.0
        self.inv_jm = 1.0 / self.jm if self.jm != 0 else 0.0

        # thermal params
        self.m_thermal = getattr(self.motor, 'm', 22.0)
        self.C_thermal = getattr(self.motor, 'C', 0.385)
        self.inv_mC = 1.0 / (self.m_thermal * self.C_thermal) if (self.m_thermal * self.C_thermal) != 0 else 0.0

        # PIDController gains: try to read from motor's controller prototypes
        id_proto = getattr(self.motor, 'id_controller', None)
        iq_proto = getattr(self.motor, 'iq_controller', None)
        sp_proto = getattr(self.motor, 'speed_controller', None)

        # CORREÇÃO: Ganhos mais conservadores para evitar instabilidade
        id_kp = getattr(id_proto, 'kp', 0.1) if id_proto is not None else 0.1
        id_ki = getattr(id_proto, 'ki', 10.0) if id_proto is not None else 10.0
        id_kd = getattr(id_proto, 'kd', 0.0) if id_proto is not None else 0.0
        id_limit = getattr(id_proto, 'limit', 600.0) if id_proto is not None else 600.0

        iq_kp = getattr(iq_proto, 'kp', 0.1) if iq_proto is not None else 0.1
        iq_ki = getattr(iq_proto, 'ki', 10.0) if iq_proto is not None else 10.0
        iq_kd = getattr(iq_proto, 'kd', 0.0) if iq_proto is not None else 0.0
        iq_limit = getattr(iq_proto, 'limit', 600.0) if iq_proto is not None else 600.0

        sp_kp = getattr(sp_proto, 'kp', 10.0) if sp_proto is not None else 5.0
        sp_ki = getattr(sp_proto, 'ki', 2.0) if sp_proto is not None else 2.0
        sp_kd = getattr(sp_proto, 'kd', 0.1) if sp_proto is not None else 0.1
        sp_limit = getattr(sp_proto, 'limit', 600.0) if sp_proto is not None else 600.0

        # instantiate internal PID controllers used by the integrator
        self.id_controller = PID.Controller(kp=id_kp, ki=id_ki, kd=id_kd, limit=id_limit, Ts=self.hp)
        self.iq_controller = PID.Controller(kp=iq_kp, ki=iq_ki, kd=iq_kd, limit=iq_limit, Ts=self.hp)
        self.speed_controller = PID.Controller(kp=sp_kp, ki=sp_ki, kd=sp_kd, limit=sp_limit, Ts=self.hp)

        # initial conditions and storage
        self.reset_initial_conditions()
        self.initialize_storage()

    def reset_initial_conditions(self) -> None:
        """Set sensible default initial conditions used to build x0 for integration."""
        # motor electrical/mechanical states
        self.isd0 = 0.0
        self.isq0 = 0.0
        self.iso0 = 0.0
        self.wm0 = 0.0
        self.theta_m0 = 0.0
        self.temp0 = 25.0

        # vehicle initial states
        self.vehicle_velocity0 = 0.0
        self.vehicle_position0 = 0.0

        # battery initial states
        if self.battery is not None:
            self.soc0 = getattr(self.battery, 'soc', 1.0)
            self.Iast0 = getattr(self.battery, 'Iast', 0.0)
            self.tempo_acumulado0 = getattr(self.battery, 'tempo_acumulado', 0.0)
        else:
            self.soc0 = 1.0
            self.Iast0 = 0.0
            self.tempo_acumulado0 = 0.0

    def initialize_storage(self) -> None:
        """Initialize lists to store simulation outputs."""
        # motor signals
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

        # vehicle
        self.vehicle_velocity = []
        self.vehicle_position = []
        self.vehicle_acceleration = []
        self.wheel_torque = []

        # tyre/adhesion
        self.tractive_force_hist = []
        self.resistive_force_hist = []
        self.slip_ratio_hist = []
        self.longitudinal_force_hist = []

        # battery
        self.soc_hist = []
        self.Iast_hist = []
        self.tempo_acumulado_hist = []
        self.battery_voltage_hist = []
        self.battery_current_hist = []

    def combined_edos(self, t: float, x: np.ndarray) -> np.ndarray:
        """Assemble combined ODEs (motor electrical/mechanical + vehicle + battery)."""
        try:
            # Unpack mandatory motor states
            isd, isq, iso, wm, theta_m, temp = x[:6]
            idx = 6

            # vehicle states
            if (self.vehicle is not None) and (self.transmission is not None):
                vehicle_velocity = x[idx]; vehicle_position = x[idx + 1]
                idx += 2
            else:
                vehicle_velocity = 0.0
                vehicle_position = 0.0

            # battery states
            if self.battery is not None:
                soc = x[idx]; Iast = x[idx + 1]; tempo_acumulado = x[idx + 2]
            else:
                soc = 1.0; Iast = 0.0; tempo_acumulado = 0.0

            # electrical transforms
            theta_e = self.p * theta_m
            we = self.p * wm

            # --- Vehicle/tyre related computations ---
            if (self.vehicle is not None) and (self.transmission is not None) and (self.tire is not None):
                try:
                    cl = self.vehicle.calculate_load_torque(vehicle_velocity, self.transmission)
                    vel_ang_roda = self.transmission.motor_to_wheel_speed(wm)
                    slip = self.tire.SlipRatio(vel_ang_roda, self.vehicle.wheel_radius, vehicle_velocity)
                    Fz = self.vehicle.mass * self.vehicle.g
                    Fx_max = self.tire.Tire_forces(Fz, slip)
                    wheel_torque_ideal = self.transmission.motor_to_wheel_torque(self.torque_constant * isq)
                    traction_force_ideal = wheel_torque_ideal / self.vehicle.wheel_radius
                    traction_force_limited = np.sign(traction_force_ideal) * min(abs(traction_force_ideal), abs(Fx_max))
                    resistance_force = self.vehicle.calculate_resistance_forces(vehicle_velocity)
                    resultant_force = traction_force_limited - resistance_force
                    if vehicle_velocity < 0.01 and resultant_force < 0:
                        vehicle_acceleration = 0.0
                    else:
                        vehicle_acceleration = resultant_force / self.vehicle.mass
                except Exception as e:
                    print(f"Erro em cálculos veiculares: {e}")
                    cl = 0.0
                    vehicle_acceleration = 0.0
                    traction_force_limited = 0.0
                    resistance_force = 0.0
                    Fx_max = 0.0
                    slip = 0.0
            else:
                cl = 0.0
                vehicle_acceleration = 0.0
                traction_force_limited = 0.0
                resistance_force = 0.0
                Fx_max = 0.0
                slip = 0.0

            # --- Speed controller (PID) ---
            speed_error = (self.speed_ref - wm) if (self.motor is not None) else 0.0
            torque_ref = self.speed_controller.update(speed_error, dt=self.hp)

            iq_ref = torque_ref / self.torque_constant if self.torque_constant != 0 else 0.0
            id_ref = 0.0
            
            # Limitar corrente de referência
            if abs(iq_ref) > self.max_current:
                iq_ref = np.sign(iq_ref) * self.max_current

            # current errors
            error_d = id_ref - isd
            error_q = iq_ref - isq

            # decoupling
            decoupling_d = -we * self.lq * isq
            decoupling_q = we * (self.ld * isd + self.lambda_m)

            # PID controllers for d/q voltages
            vd_pid = self.id_controller.update(error_d, dt=self.hp)
            vq_pid = self.iq_controller.update(error_q, dt=self.hp)

            vd_unclamped = vd_pid + decoupling_d
            vq_unclamped = vq_pid + decoupling_q

            # DC link from battery
            if self.battery is not None:
                try:
                    batt_current_est = abs(isq)
                    Vdc_now = self.battery.calcular_tensao(batt_current_est, soc, Iast, tempo_acumulado) if soc >= 0.0 else self.Vdc_default
                    if hasattr(self.inversor, 'set_Vdc'):
                        self.inversor.set_Vdc(Vdc_now)
                        # va,vb,vc = self.inversor.source_voltage(self.modulation_index)
                        # vq_r, vd_r = self.inversor.park_transform(va, vb, vc, theta_e)
                        
                except Exception as e:
                    print(f"Erro em bateria: {e}")
                    Vdc_now = self.Vdc_default
            else:
                Vdc_now = self.Vdc_default

            Vclamp = Vdc_now if Vdc_now > 0 else self.Vdc_default
            
            # CORREÇÃO: Limitação mais conservadora
            vd = np.clip(vd_unclamped, -Vclamp, Vclamp)
            vq = np.clip(vq_unclamped, -Vclamp, Vclamp)

            # electrical dynamics (linearized RL)
            d_isd = (vd - self.rs * isd + we * self.lq * isq) * self.inv_ld
            d_isq = (vq - self.rs * isq - we * (self.ld * isd + self.lambda_m)) * self.inv_lq
            
            # CORREÇÃO: Evitar divisão por zero em iso
            L0 = max(0.1 * (self.ld + self.lq) / 2.0, 1e-6)
            d_iso = (0.0 - self.rs * iso) / L0

            # mechanical dynamics
            ce = self.torque_constant * isq
            d_wm = (ce - cl - self.kf * wm) * self.inv_jm
            d_theta_m = wm

            # thermal dynamics (simple copper loss heating)
            i_rms_sq = isd**2 + isq**2
            copper_losses = 3.0 * self.rs * (i_rms_sq / 2.0)
            d_temp = copper_losses * self.inv_mC 

            # battery derivatives
            if self.battery is not None:
                try:
                    corrente_batt = abs(isq)
                    dsoc_dt, dIast_dt, dtime_dt = self.battery.calcular_derivadas(corrente_batt)
                except Exception as e:
                    print(f"Erro em derivadas da bateria: {e}")
                    dsoc_dt, dIast_dt, dtime_dt = 0.0, 0.0, 0.0
            else:
                dsoc_dt, dIast_dt, dtime_dt = 0.0, 0.0, 0.0

            # assemble derivative vector consistent with simulate() x0 layout
            dx = [d_isd, d_isq, d_iso, d_wm, d_theta_m, d_temp]
            if (self.vehicle is not None) and (self.transmission is not None):
                dx.extend([vehicle_acceleration, vehicle_velocity])
                if self.battery is not None:
                    dx.extend([dsoc_dt, dIast_dt, dtime_dt])
            else:
                if self.battery is not None:
                    dx.extend([dsoc_dt, dIast_dt, dtime_dt])
            #print(f"t={t:.3f}, wm={wm:.2f}")
            return np.array(dx, dtype=float)
            
        except Exception as e:
            print(f"Erro crítico em combined_edos: {e}")
            # Retornar derivadas zeradas em caso de erro
            n_states = len(x)
            return np.zeros(n_states, dtype=float)

    def simulate(self, t0: float = 0.0, tf: Optional[float] = None) -> Dict[str, Any]:
        """Run the simulation and return results."""
        if tf is None:
            tf = self.tmax

        print(f"🚀 Iniciando simulação: t0={t0}, tf={tf}, steps={self.steps}")

        # CORREÇÃO: Reduzir número de pontos para evitar travamento
        num_points = min(self.steps, 2000)  # Máximo 500 pontos
        t_eval = np.linspace(t0, tf, num_points)

        # build initial state vector
        if (self.vehicle is not None) and (self.transmission is not None):
            if self.battery is not None:
                x0 = np.array([
                    self.isd0, self.isq0, self.iso0,
                    self.wm0, self.theta_m0, self.temp0,
                    self.vehicle_velocity0, self.vehicle_position0,
                    self.soc0, self.Iast0, self.tempo_acumulado0
                ], dtype=float)
            else:
                x0 = np.array([
                    self.isd0, self.isq0, self.iso0,
                    self.wm0, self.theta_m0, self.temp0,
                    self.vehicle_velocity0, self.vehicle_position0
                ], dtype=float)
        else:
            if self.battery is not None:
                x0 = np.array([
                    self.isd0, self.isq0, self.iso0,
                    self.wm0, self.theta_m0, self.temp0,
                    self.soc0, self.Iast0, self.tempo_acumulado0
                ], dtype=float)
            else:
                x0 = np.array([
                    self.isd0, self.isq0, self.iso0,
                    self.wm0, self.theta_m0, self.temp0
                ], dtype=float)

        # reset storage and controllers
        self.initialize_storage()
        
        # reset internal PIDs states
        try:
            self.id_controller.reset()
            self.iq_controller.reset()
            self.speed_controller.reset()
        except Exception as e:
            print(f"Aviso ao resetar controladores: {e}")

        # integrate
        try:
            print("⏳ Executando integração numérica...")
            
            # CORREÇÃO: Método mais robusto e tolerâncias relaxadas
            sol = solve_ivp(
                fun=self.combined_edos, 
                t_span=(t0, tf), 
                y0=x0,
                method='RK23',  # Métoddo mais rapido de integração RK23
                t_eval=t_eval, 
                atol=1e-4,  # Relaxado
                rtol=1e-3,  # Relaxado  
                max_step=0.01  # Limitar passo máximo
            )
            
            print(f"✅ Integração concluída! {len(sol.t)} pontos calculados")
            
        except Exception as e:
            print(f"❌ Erro na integração: {e}")
            # Criar solução fallback
            t_fallback = np.linspace(t0, tf, 100)
            y_fallback = np.zeros((len(x0), len(t_fallback)))
            for i in range(len(x0)):
                y_fallback[i, :] = x0[i]  # Estado constante
            sol = type('obj', (object,), {'t': t_fallback, 'y': y_fallback})()

        # post-process: reconstruir sinais para logging
        print("📊 Processando resultados...")
        # post-process: reconstruct controller outputs deterministically for logging
        pid_id_log = PID.Controller(kp=self.id_controller.kp, ki=self.id_controller.ki,
                                   kd=getattr(self.id_controller, 'kd', 0.0), limit=self.id_controller.limit, Ts=self.hp)
        pid_iq_log = PID.Controller(kp=self.iq_controller.kp, ki=self.iq_controller.ki,
                                   kd=getattr(self.iq_controller, 'kd', 0.0), limit=self.iq_controller.limit, Ts=self.hp)
        pid_sp_log = PID.Controller(kp=self.speed_controller.kp, ki=self.speed_controller.ki,
                                   kd=getattr(self.speed_controller, 'kd', 0.0), limit=self.speed_controller.limit, Ts=self.hp)
        
        for idx, t in enumerate(sol.t):

            y = sol.y[:, idx]
            # motor states
            isd = y[0]; isq = y[1]; iso = y[2]
            wm = y[3]; theta_m = y[4]; temp = y[5]
            read_idx = 6

            if (self.vehicle is not None) and (self.transmission is not None):
                vehicle_velocity = y[read_idx]; vehicle_position = y[read_idx + 1]
                read_idx += 2
            else:
                vehicle_velocity = 0.0; vehicle_position = 0.0

            if self.battery is not None:
                soc = y[read_idx]; Iast = y[read_idx + 1]; tempo_acumulado = y[read_idx + 2]
            else:
                soc = 1.0; Iast = 0.0; tempo_acumulado = 0.0

            # time step for PID reconstruction
            if idx == 0:
                dt = self.hp
            else:
                dt = sol.t[idx] - sol.t[idx - 1]
                if dt <= 0:
                    dt = self.hp

            # speed PID -> torque reference
            speed_error = (self.speed_ref - wm) if (self.motor is not None) else 0.0
            torque_ref = pid_sp_log.update(speed_error, dt=dt)

            iq_ref = torque_ref / self.torque_constant if self.torque_constant != 0 else 0.0
            id_ref = 0.0

            # current controllers
            error_d = id_ref - isd
            error_q = iq_ref - isq

            # electrical freq and decoupling
            theta_e = self.p * theta_m
            we = self.p * wm
            decoupling_d = -we * self.lq * isq
            decoupling_q = we * (self.ld * isd + self.lambda_m)

            # DC link from battery
            if self.battery is not None:
                try:
                    batt_current_est = abs(isq)
                    Vdc_now = self.battery.calcular_tensao(batt_current_est, soc, Iast, tempo_acumulado) if soc >= 0.0 else self.Vdc_default
                    if hasattr(self.inversor, 'set_Vdc'):
                        self.inversor.set_Vdc(Vdc_now)
                except Exception as e:
                    print(f"Erro em bateria: {e}")
                    Vdc_now = self.Vdc_default
            else:
                Vdc_now = self.Vdc_default

            Vclamp = Vdc_now if Vdc_now > 0 else self.Vdc_default

            vd_pid = pid_id_log.update(error_d, dt=dt)
            vq_pid = pid_iq_log.update(error_q, dt=dt)

            vd_control = np.clip(vd_pid + decoupling_d, -Vclamp, Vclamp)
            vq_control = np.clip(vq_pid + decoupling_q, -Vclamp, Vclamp)

            # battery/inverter interaction for logging
            if self.battery is not None:
                battery_current_est = abs(isq)
                battery_voltage = self.battery.calcular_tensao(battery_current_est, soc, Iast, tempo_acumulado) if soc >= 0.0 else 0.0
            else:
                battery_voltage = self.Vdc_default

            # inverter voltages (best-effort)
            try:
                if hasattr(self.inversor, 'set_Vdc'):
                    self.inversor.set_Vdc(battery_voltage)
                if hasattr(self.inversor, 'set_ws'):
                    self.inversor.set_ws(self.p * wm)
                modulation = getattr(self.motor, 'valor_mu', self.modulation_index)
                try:
                    va, vb, vc = self.inversor.source_voltage(modulation)
                except TypeError:
                    va, vb, vc = self.inversor.source_voltage()
                if hasattr(self.inversor, 'park_transform'):
                    vq_r, vd_r = self.inversor.park_transform(va, vb, vc, theta_e)
                else:
                    vq_r, vd_r = 0.0, 0.0
            except Exception:
                va, vb, vc, vq_r, vd_r = 0.0, 0.0, 0.0, 0.0, 0.0

            # abc currents
            try:
                is1, is2, is3, fs1, fs2, fs3 = self.motor.abc_currents_from_dq(isd, isq, theta_e, self.ld * isd + self.lambda_m)
            except Exception:
                is1, is2, is3, fs1, fs2, fs3 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

            cl = self.vehicle.calculate_load_torque(vehicle_velocity, self.transmission) if (self.vehicle is not None and self.transmission is not None) else 0.0
            ce = self.torque_constant * isq

            # append logs
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
            self.fluxosd.append(self.ld * isd + self.lambda_m)
            self.fluxosq.append(self.lq * isq)
            self.conjugado.append(ce)
            self.velocidade.append(wm * 60.0 / (2.0 * np.pi))
            self.torque_mecanico.append(ce - cl)
            self.conjcarga.append(cl)
            self.temperatura.append(temp)
            self.vd_real.append(vd_r)
            self.vq_real.append(vq_r)
            self.vd_control.append(vd_control)
            self.vq_control.append(vq_control)
            self.speed_error.append(speed_error)
            self.iq_ref_trace.append(iq_ref)
            self.id_ref_trace.append(id_ref)

            self.vehicle_velocity.append(vehicle_velocity)
            self.vehicle_position.append(vehicle_position)

            self.soc_hist.append(soc)
            self.Iast_hist.append(Iast)
            self.tempo_acumulado_hist.append(tempo_acumulado)
            self.battery_voltage_hist.append(battery_voltage)
            self.battery_current_hist.append(abs(isq))

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
                self.wheel_torque.append(wheel_torque) 
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

        # assemble results
        results = {
            't': np.array(self.tempo),
            'isd': np.array(self.corrented),
            'isq': np.array(self.correnteq),
            'is1': np.array(self.corrente1),
            'is2': np.array(self.corrente2),
            'is3': np.array(self.corrente3),
            'va': np.array(self.tensao1),
            'vb': np.array(self.tensao2),
            'vc': np.array(self.tensao3),
            'soc': np.array(self.soc_hist),
            'battery_voltage': np.array(self.battery_voltage_hist),
            'vehicle_velocity': np.array(self.vehicle_velocity),
            'vehicle_acceleration': np.array(self.vehicle_acceleration),
            'traction_force': np.array(self.tractive_force_hist),
            'resistance_force': np.array(self.resistive_force_hist),
            'slip_ratio': np.array(self.slip_ratio_hist),
            'rpm': np.array(self.velocidade),
            'torque_electromagnetic': np.array(self.conjugado),
            'torque_mechanical': np.array(self.torque_mecanico),
            'vd_control': np.array(self.vd_control),
            'vq_control': np.array(self.vq_control),
            'iq_ref': np.array(self.iq_ref_trace),
            'id_ref': np.array(self.id_ref_trace)
        }

        print(f"✅ Simulação finalizada! {len(results['t'])} pontos coletados")
        return results