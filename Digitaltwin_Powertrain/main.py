import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import plotly.express as px
import dash_bootstrap_components as dbc
from models import BatteryPack, Inversor, Tire, Transmission, Vehicle, Motor
from simulation.Simulation import Simulation

def build_defaults():
    transmission = Transmission.Transmission(final_drive_ratio=4.0, efficiency=0.95)
    vehicle = Vehicle.Vehicle(mass=230.0, wheel_radius=0.275, drag_coeff=0.7789,
                      frontal_area=0.68, rolling_resistance=0.015, road_grade=0.0)
    battery = BatteryPack.BatteryPack(tipo_celula='Li-ion', n_serie=162, n_paralelo=1, soc_inicial=1.0)
    tire = Tire.Tire(pacejka_params=[0.333, 1.627, 1, 4.396, 931.4, 366.4], tire_friction_coef=1.45)
    inversor = Inversor.Inversor(eficiencia=0.95, freq_chaveamento=10000)
    motor = Motor(rs=0.04585, ld=0.00067, lq=0.00172, jm=0.05769, kf=0.1,
                  lambda_m=0.13849, p=10, valor_mu=0.99,
                  TL=False, torque=0.0,
                  speed_ref=471.23)

    sim = Simulation(motor=motor, vehicle=vehicle, transmission=transmission,
                     battery=battery, tire=tire, inversor=inversor, tmax=10, steps=10)
    return sim

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
        'vehicle': {'mass': 230, 'wheel_radius': 0.275, 'drag_coeff': 0.7789,
                    'frontal_area': 0.68, 'rolling_resistance': 0.015, 'road_grade': 0},
        'transmission': {'final_drive_ratio': 4.0, 'efficiency': 0.95},
        'battery': {'tipo_celula': 'Li-ion', 'n_serie': 162, 'n_paralelo': 1, 'soc_inicial': 1.0},
        'motor': {'rs': 0.04585, 'ld': 0.00067, 'lq': 0.00067, 'jm': 0.05769,
                  'kf': 0.1, 'lambda_m': 0.13849, 'p': 10, 'valor_mu': 0.95, 'velocidade_ref': 471.23},
        'inversor': {'eficiencia': 0.95, 'freq_chaveamento': 10000},
        'simulacao': {'tmax': 10},
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
                    dbc.Input(id='wheel-radius', type='number', value=0.275, step=0.01),
                    dbc.Label("Coef. de Arrasto"),
                    dbc.Input(id='drag-coeff', type='number', value=0.7789, step=0.01),
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
                    dbc.Input(id='motor-lq', type='number', value=0.00137, step=0.00001),
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
                    dbc.Input(id='sim-time', type='number', value=10, step=1),
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
def create_velocity_torque_plot(sim):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Velocidade Mecânica", "Torque do Motor"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.velocidade, name='Velocidade'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.conjugado, name='Torque Elétrico'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.conjcarga, name='Torque de Carga', line=dict(dash='dash')), row=2, col=1)
    fig.update_layout(height=900, title_text="Velocidade e Torque", template="plotly_white")
    fig.update_xaxes(title_text="Tempo (s)", row=2, col=1)
    fig.update_yaxes(title_text="Velocidade (RPM)", row=1, col=1)
    fig.update_yaxes(title_text="Torque (Nm)", row=2, col=1)
    return fig

def create_currents_plot(sim):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Correntes dq", "Correntes de Fase"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.corrented, name='Id'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.correnteq, name='Iq'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.corrente1, name='Fase 1'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.corrente2, name='Fase 2'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.corrente3, name='Fase 3'), row=2, col=1)
    fig.update_layout(height=900, title_text="Correntes", template="plotly_white")
    return fig

def create_voltages_plot(sim):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Tensões de Controle dq", "Tensões de Fase"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensaosd, name='Vd'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.vd_real, name='Vd real'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensaosq, name='Vq'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.vq_real, name='Vq real'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensao1, name='Fase 1'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensao2, name='Fase 2'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensao3, name='Fase 3'), row=2, col=1)
    fig.update_layout(height=900, title_text="Tensões", template="plotly_white")
    return fig

def create_flux_temp_plot(sim):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Fluxos Magnéticos", "Temperatura do sim"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.fluxosd, name='Fluxo d'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.fluxosq, name='Fluxo q'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.temperatura, name='Temperatura'), row=2, col=1)
    fig.update_layout(height=900, title_text="Fluxos e Temperatura", template="plotly_white")
    return fig

def create_control_plot(sim):
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Sinais de Controle FOC", "Erro de Velocidade"), vertical_spacing=0.1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.vd_control, name='Vd control'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.vq_control, name='Vq control'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.speed_error, name='Erro de Velocidade'), row=2, col=1)
    fig.update_layout(height=900, title_text="Sinais de Controle", template="plotly_white")
    return fig

def create_complete_plot(sim):
    fig = make_subplots(rows=3, cols=2, subplot_titles=(
        "Velocidade Mecânica", "Torque do sim", "Correntes", "Correntes de Fase", "Tensões de Controle", "Tensões de Fase"
    ), vertical_spacing=0.08, horizontal_spacing=0.1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.velocidade, name='Velocidade'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.conjugado, name='Torque Elétrico'), row=1, col=2)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.conjcarga, name='Torque de Carga', line=dict(dash='dash')), row=1, col=2)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.corrented, name='Id'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.correnteq, name='Iq'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.corrente1, name='Fase 1'), row=2, col=2)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.corrente2, name='Fase 2'), row=2, col=2)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.corrente3, name='Fase 3'), row=2, col=2)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensaosd, name='Vd'), row=3, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.vd_real, name='Vd real'), row=3, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensaosq, name='Vq'), row=3, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.vq_real, name='Vq real'), row=3, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensao1, name='Fase 1'), row=3, col=2)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensao2, name='Fase 2'), row=3, col=2)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tensao3, name='Fase 3'), row=3, col=2)
    fig.update_layout(height=900, title_text="Visão Completa da Simulação", template="plotly_white", hovermode="x unified", showlegend=False)
    return fig

def create_vehicle_plot(sim):
    fig = make_subplots(
        rows=2, cols=2, 
        subplot_titles=("Velocidade (km/h)", "Aceleração (m/s²)", "Forças (N)", "Torque na Roda (Nm)"), 
        vertical_spacing=0.15, horizontal_spacing=0.1
    )
    velocity_kmh = [v * 3.6 for v in sim.vehicle_velocity]
    fig.add_trace(go.Scatter(x=sim.tempo, y=velocity_kmh, name='Velocidade'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.vehicle_acceleration, name='Aceleração'), row=1, col=2)
    
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.tractive_force_hist, name='Força Trativa'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.resistive_force_hist, name='Forças Resistivas'), row=2, col=1)
    
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.wheel_torque, name='Torque na Roda'), row=2, col=2)
    
    fig.update_layout(height=900, title_text="Desempenho do Veículo", template="plotly_white", legend_tracegroupgap=180)
    fig.update_yaxes(title_text="Força (N)", row=2, col=1)
    return fig

def create_tire_plot(sim):
    fig = make_subplots(
        rows=1, cols=2, 
        subplot_titles=("Slip Ratio x Tempo", "Força Longitudinal x Slip Ratio")
    )

    fig.add_trace(
        go.Scatter(x=sim.tempo, y=sim.slip_ratio_hist, name='Slip Ratio', line=dict(color='blue')),
        row=1, col=1
    )
    fig.update_xaxes(title_text="Tempo [s]", row=1, col=1)
    fig.update_yaxes(title_text="Slip Ratio", row=1, col=1)

    if len(sim.slip_ratio_hist) > 0:
        df_tire = pd.DataFrame({
            'slip': sim.slip_ratio_hist,
            'fx': sim.longitudinal_force_hist
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

def create_battery_plot(sim):
    fig = make_subplots(rows=3, cols=1, subplot_titles=("Tensão do Banco (V)", "Corrente (A)", "SoC"), vertical_spacing=0.12)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.battery_voltage_hist, name='Tensão Banco'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.battery_current_hist, name='Corrente'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sim.tempo, y=sim.soc_hist, name='SoC'), row=3, col=1)
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
        'valor_mu': motor_modulation,
        'velocidade_ref': speed_ref
    }

    current_params['simulacao'] = {
        'tmax': sim_time
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
        torque=0.0 ,
        speed_ref= motor_params['velocidade_ref']
    )
    

 
    sim = Simulation(
        motor=motor, vehicle=vehicle, transmission=transmission,
        battery=battery, tire=tire, inversor=Inversor(eficiencia=0.95, freq_chaveamento=1000),
        tmax=simulacao_params['tmax'], steps=2000
    )

    sim.simulate()
    
    # Gerar todas as figuras
    figures = {
        'velocity_torque': create_velocity_torque_plot(sim),
        'currents': create_currents_plot(sim),
        'voltages': create_voltages_plot(sim),
        'flux_temp': create_flux_temp_plot(sim),
        'control': create_control_plot(sim),
        'complete': create_complete_plot(sim),
        'vehicle': create_vehicle_plot(sim),
        'tire': create_tire_plot(sim),
        'battery': create_battery_plot(sim)
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
    print("")
    app.run(debug=True)