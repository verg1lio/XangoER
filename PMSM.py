import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
import numpy as np
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import plotly.express as px

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

class Motor:
    """Motor PMSM com FOC usando solve_ivp"""
    def __init__(self, rs, ld, lq, jm, kf, lambda_m, p, valor_mu, TL=False, torque=0.0):
        self.pi23 = 2 * np.pi / 3
        self.rq23 = np.sqrt(2 / 3)
        self.rq3 = np.sqrt(3)
        self.rs = rs
        self.ld = ld
        self.lq = lq
        self.jm = jm
        self.kf = kf
        self.lambda_m = lambda_m
        self.p = p
        self.valor_mu = valor_mu

        self.m = 22.0
        self.C = 0.385

        self.id_controller = PIController(kp=0.5, ki=100.0, limit=1000.0)
        self.iq_controller = PIController(kp=0.5, ki=100.0, limit=1000.0)
        self.speed_controller = PIController(kp=1.0, ki=5, limit=500.0)

        self.TL = bool(TL)
        if callable(torque):
            self._external_torque = torque
        else:
            self._external_torque = (lambda t, v=float(torque): v)

        self.max_current = 220.0 * np.sqrt(2)
        self.Vdc = 600.0
        self.Vs = (self.valor_mu * self.Vdc) / np.sqrt(3)
        self.Vlimit = self.Vs * 3

        self.id_ref = 0.0
        self.iq_ref = 0.0
        self.speed_ref = 471.23

        self.tmax = 2.0
        self.hp = self.tmax / 2000.0

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
        self.speed_error = []
        self.iq_ref_trace = []
        self.id_ref_trace = []

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

        iq_ref = torque_ref / (1.5 * self.p * self.lambda_m) if (1.5 * self.p * self.lambda_m) != 0 else 0.0
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
        vmax = self.Vlimit
        if vmag > vmax:
            vd *= vmax / vmag
            vq *= vmax / vmag

        return vd, vq, id_ref, iq_ref, speed_error

    def inverse_park_transform(self, vd, vq, theta_e):
        valpha = vd * np.cos(theta_e) - vq * np.sin(theta_e)
        vbeta = vd * np.sin(theta_e) + vq * np.cos(theta_e)
        v0 = 0.0
        vs1 = valpha
        vs2 = -0.5 * valpha + (np.sqrt(3) / 2.0) * vbeta
        vs3 = -0.5 * valpha - (np.sqrt(3) / 2.0) * vbeta
        return vs1, vs2, vs3, v0

    def abc_currents_from_dq(self, isd, isq, theta_e, flux_d):
        is1 = self.rq23 * (isd * np.cos(theta_e) - isq * np.sin(theta_e))
        is2 = self.rq23 * (isd * np.cos(theta_e - self.pi23) - isq * np.sin(theta_e - self.pi23))
        is3 = self.rq23 * (isd * np.cos(theta_e + self.pi23) - isq * np.sin(theta_e + self.pi23))
        fs1 = self.rq23 * flux_d
        fs2 = self.rq23 * flux_d
        fs3 = self.rq23 * flux_d
        return is1, is2, is3, fs1, fs2, fs3

    def edos(self, t, x):
        isd, isq, iso, wm, theta_m, temp, int_id, int_iq, int_speed = x
        theta_e = self.p * theta_m
        we = self.p * wm

        cl = self.set_load(t)

        speed_error = self.speed_ref - wm
        torque_ref_unsat = self.speed_controller.kp * speed_error + self.speed_controller.ki * int_speed
        int_speed_dot = speed_error
        torque_ref = torque_ref_unsat

        denom = (1.5 * self.p * self.lambda_m)
        iq_ref = torque_ref / denom if denom != 0 else 0.0
        id_ref = 0.0
        if abs(iq_ref) > self.max_current:
            iq_ref = np.sign(iq_ref) * self.max_current

        error_d = id_ref - isd
        error_q = iq_ref - isq

        vd_unclamped = self.id_controller.kp * error_d + self.id_controller.ki * int_id + (- we * self.lq * isq)
        vq_unclamped = self.iq_controller.kp * error_q + self.iq_controller.ki * int_iq + ( we * (self.ld * isd + self.lambda_m) )

        vd = np.clip(vd_unclamped, -self.Vlimit, self.Vlimit)
        vq = np.clip(vq_unclamped, -self.Vlimit, self.Vlimit)

        int_id_dot = 0.0 if abs(vd_unclamped) > self.Vlimit else error_d
        int_iq_dot = 0.0 if abs(vq_unclamped) > self.Vlimit else error_q

        d_isd = (vd - self.rs * isd + we * self.lq * isq) / self.ld
        d_isq = (vq - self.rs * isq - we * (self.ld * isd + self.lambda_m)) / self.lq
        L0 = 0.1 * (self.ld + self.lq) / 2.0
        d_iso = (0.0 - self.rs * iso) / L0

        ce = 1.5 * self.p * (self.lambda_m * isq)
        d_wm = (ce - cl - self.kf * wm) / self.jm
        d_theta_m = wm

        i_rms = np.sqrt(isd**2 + isq**2) / np.sqrt(2)
        copper_losses = 3.0 * self.rs * (i_rms**2)
        d_temp = copper_losses / (self.m * self.C)

        if abs(iq_ref) >= self.max_current and np.sign(int_speed) == np.sign(speed_error):
            int_speed_dot = 0.0

        dxdt = np.array([
            d_isd, d_isq, d_iso,
            d_wm, d_theta_m, d_temp,
            int_id_dot, int_iq_dot, int_speed_dot
        ], dtype=float)

        return dxdt

    def simulate(self, t0=0.0, tf=None):
        if tf is None:
            tf = self.tmax
        t_eval = np.arange(t0, tf + self.hp, self.hp)

        x0 = np.array([
            self.isd0, self.isq0, self.iso0,
            self.wm0, self.theta_m0, self.temp0,
            self.int_id0, self.int_iq0, self.int_speed0
        ], dtype=float)

        sol = solve_ivp(fun=self.edos, t_span=(t0, tf), y0=x0, method='RK45', t_eval=t_eval, vectorized=False, atol=1e-6, rtol=1e-6)

        self.initialize_storage()

        self.id_controller.reset()
        self.iq_controller.reset()
        self.speed_controller.reset()

        dt = self.hp
        for idx, t in enumerate(sol.t):
            isd = sol.y[0, idx]
            isq = sol.y[1, idx]
            iso = sol.y[2, idx]
            wm = sol.y[3, idx]
            theta_m = sol.y[4, idx]
            temp = sol.y[5, idx]
            int_id = sol.y[6, idx]
            int_iq = sol.y[7, idx]
            int_speed = sol.y[8, idx]

            theta_e = self.p * theta_m
            we = self.p * wm

            vd_control, vq_control, id_ref, iq_ref, speed_error = self.field_oriented_control(isd, isq, wm, dt)

            flux_d = self.ld * isd + self.lambda_m
            flux_q = self.lq * isq
            flux_o = 0.1 * (self.ld + self.lq) / 2.0 * iso

            ce = 1.5 * self.p * (self.lambda_m * isq)

            cl = self.set_load(t)
            cm = ce - cl

            vs1, vs2, vs3, v0 = self.inverse_park_transform(vd_control, vq_control, theta_e)

            is1, is2, is3, fs1, fs2, fs3 = self.abc_currents_from_dq(isd, isq, theta_e, flux_d)

            self.tempo.append(t)
            self.corrented.append(isd)
            self.correnteq.append(isq)
            self.corrente1.append(is1)
            self.corrente2.append(is2)
            self.corrente3.append(is3)
            self.tensaosd.append(vd_control)
            self.tensaosq.append(vq_control)
            self.tensao1.append(is1 * self.rs)
            self.tensao2.append(is2 * self.rs)
            self.tensao3.append(is3 * self.rs)
            self.fluxosd.append(flux_d)
            self.fluxosq.append(flux_q)
            self.conjugado.append(ce)
            #self.velocidade.append(wm)
            self.velocidade.append(wm * 60 / (2 * np.pi))  # Convert rad/s to RPM
            self.torque_mecanico.append(cm)
            self.conjcarga.append(cl)
            self.temperatura.append(temp)
            self.vd_control.append(vd_control)
            self.vq_control.append(vq_control)
            self.speed_error.append(speed_error)
            self.iq_ref_trace.append(iq_ref)
            self.id_ref_trace.append(id_ref)

# Criar instância do motor e simular
motor = Motor(
    rs=0.04585,
    ld=0.00067,
    lq=0.00067,
    jm=0.05769,
    kf=0.1,
    lambda_m=0.13849,
    p=10,
    valor_mu=0.95,
    TL=False,
    torque=0.0
)

motor.speed_ref = 471.23
motor.simulate()

# Criar aplicação Dash
external_stylesheets = [
    'https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css',
    'https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap',
    'https://codepen.io/chriddyp/pen/bWLwgP.css'
]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

# Layout da aplicação
app.layout = html.Div([
    html.Div([
        html.H1("Sistema de Controle FOC para Motor PMSM", 
                style={'textAlign': 'center', 'color': '#2c3e50', 'marginBottom': '30px'}),
        
        html.Div([
            html.Button([
                html.I(className="fas fa-tachometer-alt", style={'marginRight': '8px'}),
                "Velocidade e Torque"
            ], id='btn-velocity-torque', className='btn btn-primary'),
            
            html.Button([
                html.I(className="fas fa-bolt", style={'marginRight': '8px'}),
                "Correntes"
            ], id='btn-currents', className='btn btn-success'),
            
            html.Button([
                html.I(className="fas fa-plug", style={'marginRight': '8px'}),
                "Tensões"
            ], id='btn-voltages', className='btn btn-warning'),
            
            html.Button([
                html.I(className="fas fa-magnet", style={'marginRight': '8px'}),
                "Fluxos e Temperatura"
            ], id='btn-flux-temp', className='btn btn-info'),
            
            html.Button([
                html.I(className="fas fa-sliders-h", style={'marginRight': '8px'}),
                "Sinais de Controle"
            ], id='btn-control', className='btn btn-secondary'),
            
            html.Button([
                html.I(className="fas fa-chart-bar", style={'marginRight': '8px'}),
                "Visão Completa"
            ], id='btn-complete', className='btn btn-dark'),
        ], style={'textAlign': 'center', 'marginBottom': '30px', 'display': 'flex', 
                 'justifyContent': 'center', 'flexWrap': 'wrap', 'gap': '10px'}),
        
        dcc.Graph(id='motor-graph', style={'height': '70vh', 'borderRadius': '8px', 
                                          'boxShadow': '0 4px 6px rgba(0,0,0,0.1)'}),
        
        dcc.Store(id='simulation-data')  # Para armazenar os dados da simulação
    ], style={'padding': '20px', 'backgroundColor': '#f8f9fa', 'minHeight': '100vh'})
])

# Funções para criar os gráficos (mantidas como antes)
def create_velocity_torque_plot(motor):
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Velocidade Mecânica", "Torque do Motor"),
        vertical_spacing=0.1
    )
    
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.velocidade, name='Velocidade', 
                            line=dict(color='#3498db')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.conjugado, name='Torque Elétrico', 
                            line=dict(color='#e74c3c')), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.conjcarga, name='Torque de Carga', 
                            line=dict(dash='dash', color='#7f8c8d')), row=2, col=1)
    
    fig.update_layout(height=700, title_text="Velocidade e Torque", 
                     template="plotly_white", showlegend=True)
    fig.update_xaxes(title_text="Tempo (s)", row=2, col=1)
    fig.update_yaxes(title_text="Velocidade (RPM)", row=1, col=1)
    fig.update_yaxes(title_text="Torque (Nm)", row=2, col=1)
    
    return fig

def create_currents_plot(motor):
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Correntes dq", "Correntes de Fase"),
        vertical_spacing=0.1
    )
    
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrented, name='Id', 
                            line=dict(color='#3498db')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.correnteq, name='Iq', 
                            line=dict(color='#e74c3c')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente1, name='Fase 1', 
                            line=dict(color='#3498db')), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente2, name='Fase 2', 
                            line=dict(color='#e74c3c')), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente3, name='Fase 3', 
                            line=dict(color='#2ecc71')), row=2, col=1)
    
    fig.update_layout(height=700, title_text="Correntes", 
                     template="plotly_white", showlegend=True)
    fig.update_xaxes(title_text="Tempo (s)", row=2, col=1)
    fig.update_yaxes(title_text="Corrente (A)", row=1, col=1)
    fig.update_yaxes(title_text="Corrente (A)", row=2, col=1)
    
    return fig

def create_voltages_plot(motor):
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Tensões de Controle dq", "Tensões de Fase"),
        vertical_spacing=0.1
    )
    
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensaosd, name='Vd', 
                            line=dict(color='#3498db')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensaosq, name='Vq', 
                            line=dict(color='#e74c3c')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao1, name='Fase 1', 
                            line=dict(color='#3498db')), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao2, name='Fase 2', 
                            line=dict(color='#e74c3c')), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao3, name='Fase 3', 
                            line=dict(color='#2ecc71')), row=2, col=1)
    
    fig.update_layout(height=700, title_text="Tensões", 
                     template="plotly_white", showlegend=True)
    fig.update_xaxes(title_text="Tempo (s)", row=2, col=1)
    fig.update_yaxes(title_text="Tensão (V)", row=1, col=1)
    fig.update_yaxes(title_text="Tensão (V)", row=2, col=1)
    
    return fig

def create_flux_temp_plot(motor):
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Fluxos Magnéticos", "Temperatura do Motor"),
        vertical_spacing=0.1
    )
    
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.fluxosd, name='Fluxo d', 
                            line=dict(color='#3498db')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.fluxosq, name='Fluxo q', 
                            line=dict(color='#e74c3c')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.temperatura, name='Temperatura', 
                            line=dict(color='#f39c12')), row=2, col=1)
    
    fig.update_layout(height=700, title_text="Fluxos e Temperatura", 
                     template="plotly_white", showlegend=True)
    fig.update_xaxes(title_text="Tempo (s)", row=2, col=1)
    fig.update_yaxes(title_text="Fluxo (Wb)", row=1, col=1)
    fig.update_yaxes(title_text="Temperatura (°C)", row=2, col=1)
    
    return fig

def create_control_plot(motor):
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Sinais de Controle FOC", "Erro de Velocidade"),
        vertical_spacing=0.1
    )
    
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vd_control, name='Vd control', 
                            line=dict(color='#3498db')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.vq_control, name='Vq control', 
                            line=dict(color='#e74c3c')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.speed_error, name='Erro de Velocidade', 
                            line=dict(color='#9b59b6')), row=2, col=1)
    
    fig.update_layout(height=700, title_text="Sinais de Controle", 
                     template="plotly_white", showlegend=True)
    fig.update_xaxes(title_text="Tempo (s)", row=2, col=1)
    fig.update_yaxes(title_text="Tensão (V)", row=1, col=1)
    fig.update_yaxes(title_text="Erro (rad/s)", row=2, col=1)
    
    return fig

def create_complete_plot(motor):
    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=(
            "Velocidade Mecânica",
            "Torque do Motor",
            "Correntes dq",
            "Correntes de Fase",
            "Tensões de Controle dq",
            "Tensões de Fase"
        ),
        vertical_spacing=0.08,
        horizontal_spacing=0.1
    )

    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.velocidade, name='Velocidade', 
                            line=dict(color='#3498db')), row=1, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.conjugado, name='Torque Elétrico', 
                            line=dict(color='#e74c3c')), row=1, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.conjcarga, name='Torque de Carga', 
                            line=dict(dash='dash', color='#7f8c8d')), row=1, col=2)

    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrented, name='Id', 
                            line=dict(color='#3498db')), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.correnteq, name='Iq', 
                            line=dict(color='#e74c3c')), row=2, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente1, name='Fase 1', 
                            line=dict(color='#3498db')), row=2, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente2, name='Fase 2', 
                            line=dict(color='#e74c3c')), row=2, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.corrente3, name='Fase 3', 
                            line=dict(color='#2ecc71')), row=2, col=2)

    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensaosd, name='Vd', 
                            line=dict(color='#3498db')), row=3, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensaosq, name='Vq', 
                            line=dict(color='#e74c3c')), row=3, col=1)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao1, name='Fase 1', 
                            line=dict(color='#3498db')), row=3, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao2, name='Fase 2', 
                            line=dict(color='#e74c3c')), row=3, col=2)
    fig.add_trace(go.Scatter(x=motor.tempo, y=motor.tensao3, name='Fase 3', 
                            line=dict(color='#2ecc71')), row=3, col=2)

    fig.update_layout(height=900, title_text="Visão Completa da Simulação", 
                     template="plotly_white", hovermode="x unified", showlegend=False)
    fig.update_xaxes(title_text="Tempo (s)", row=3, col=1)
    fig.update_xaxes(title_text="Tempo (s)", row=3, col=2)
    fig.update_yaxes(title_text="Velocidade (RPM)", row=1, col=1)
    fig.update_yaxes(title_text="Torque (Nm)", row=1, col=2)
    fig.update_yaxes(title_text="Corrente (A)", row=2, col=1)
    fig.update_yaxes(title_text="Corrente (A)", row=2, col=2)
    fig.update_yaxes(title_text="Tensão (V)", row=3, col=1)
    fig.update_yaxes(title_text="Tensão (V)", row=3, col=2)
    
    return fig

# Callback para executar a simulação e atualizar o gráfico
@app.callback(
    [Output('motor-graph', 'figure'),
     Output('simulation-data', 'data')],
    [Input('btn-velocity-torque', 'n_clicks'),
     Input('btn-currents', 'n_clicks'),
     Input('btn-voltages', 'n_clicks'),
     Input('btn-flux-temp', 'n_clicks'),
     Input('btn-control', 'n_clicks'),
     Input('btn-complete', 'n_clicks')]
)
def update_graph(btn_velocity, btn_currents, btn_voltages, btn_flux, btn_control, btn_complete):
    ctx = dash.callback_context
    if not ctx.triggered:
        # Executar simulação inicial
        motor = Motor(
            rs=0.04585,
            ld=0.00067,
            lq=0.00067,
            jm=0.05769,
            kf=0.1,
            lambda_m=0.13849,
            p=10,
            valor_mu=0.95,
            TL=False,
            torque=0.0
        )
        motor.speed_ref = 471.23
        motor.simulate()
        
        # Converter dados para JSON (apenas para demonstração)
        data = {
            'tempo': motor.tempo,
            'velocidade': motor.velocidade,
            'conjugado': motor.conjugado,
            # Adicione outras variáveis conforme necessário
        }
        
        return create_velocity_torque_plot(motor), data
    
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    # Executar simulação (ou usar dados armazenados)
    motor = Motor(
        rs=0.04585,
        ld=0.00067,
        lq=0.00067,
        jm=0.05769,
        kf=0.1,
        lambda_m=0.13849,
        p=10,
        valor_mu=0.95,
        TL=False,
        torque=0.0
    )
    motor.speed_ref = 471.23
    motor.simulate()
    
    # Converter dados para JSON (apenas para demonstração)
    data = {
        'tempo': motor.tempo,
        'velocidade': motor.velocidade,
        'conjugado': motor.conjugado,
        # Adicione outras variáveis conforme necessário
    }
    
    if button_id == 'btn-velocity-torque':
        return create_velocity_torque_plot(motor), data
    elif button_id == 'btn-currents':
        return create_currents_plot(motor), data
    elif button_id == 'btn-voltages':
        return create_voltages_plot(motor), data
    elif button_id == 'btn-flux-temp':
        return create_flux_temp_plot(motor), data
    elif button_id == 'btn-control':
        return create_control_plot(motor), data
    elif button_id == 'btn-complete':
        return create_complete_plot(motor), data

if __name__ == '__main__':
    app.run(debug=True)
