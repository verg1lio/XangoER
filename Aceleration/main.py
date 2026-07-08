import dash
import numpy as np
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash_bootstrap_components as dbc
from Models import BatteryPack, Tire, Transmission, Vehicle, Motor, Pedal
from Simulation.Simulation import Simulation

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
    
    # 1. ATUALIZAÇÃO NO STORE: Adicionados L, h, dist_cg
    dcc.Store(id='parameters-store', data={
        'vehicle': {
            'mass': 240,
            'wheel_radius': 0.26,
            'wheel_mass': 5.0,
            'drag_coeff': 0.7789,
            'frontal_area': 0.68,
            'rolling_resistance': 0.015,
            'road_grade': 0,
            'L': 1.5,
            'h': 0.28,
            'dist_cg': 0.6,
        },
        'transmission': {
            'final_drive_ratio': 5.0,
            'Z_pinhao': None, 'Z_coroa': None,
            'chain_efficiency': 0.98, 'diff_efficiency': 0.97,
            'axle_inertia': 0.0015, 'diff_inertia': 0.012,
            'sprocket_inertia': 0.0002, 'coroa_inertia': 0.0015,
            'lsd_bias_ratio': 3.0,
        },
        'aero': {
            'has_wing': True,
            'lift_coeff_front': 1.0, 'area_front': 0.35,
            'lift_coeff_rear': 1.0,  'area_rear': 0.35,
        },
        'battery': {'tipo_celula': 'Li-ion', 'n_serie': 144, 'n_paralelo': 5, 'soc_inicial': 1.0},
        'motor': {'rs': 0.00706, 'ld': 0.0000965, 'lq': 0.0000965, 'jm': 0.02521,
                  'kf': 0.005, 'lambda_m': 0.04748, 'p': 10, 'valor_mu': 0.95, 'velocidade_ref': 575.96},
        'inversor': {'eficiencia': 0.95, 'freq_chaveamento': 10000},
        'simulacao': {'distancia_max': 75, 'pedal_t_pico': 3.0},
        'tire': {
            'pacejka_params': dict(Tire.Tire.HOOSIER_FSAE_LONG),
            'tire_friction_coef': 0.6,
            'Fz0': 654.0,
            'camber': 0.0,
        }
    }),
    dcc.Store(id='sidebar-state', data={'visible': True}),
    dcc.Store(id='stored-figures', data={}),
    dcc.Store(id='simulation-run-flag', data={'run': False}),
    dcc.Store(id='sim-history', data=[]),
    dcc.Store(id='last-params', data=None),

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
            html.Div([
                html.H4("Parâmetros", style={'margin': '0', 'color': '#2c3e50', 'fontWeight': '600'}),
                html.Small("Digital-Twin Xangô-E.Racing",
                           style={'color': '#6c757d', 'fontSize': '11px'}),
            ], style={'padding': '10px 4px 16px 32px', 'borderBottom': '1px solid #dee2e6',
                      'marginBottom': '12px'}),

            dbc.Accordion(
                [
                    dbc.AccordionItem(
                        [
                            dbc.Label("Massa (kg)"),
                            dbc.Input(id='vehicle-mass', type='number', value=240, step=1),
                            dbc.Label("Raio da Roda (m)"),
                            dbc.Input(id='wheel-radius', type='number', value=0.26, step=0.01),
                            dbc.Label("Massa da Roda (kg)"),
                            dbc.Input(id='wheel-mass', type='number', value=5.0, step=0.1),
                            dbc.Label("Coef. de Arrasto"),
                            dbc.Input(id='drag-coeff', type='number', value=0.7789, step=0.01),
                            dbc.Label("Área Frontal (m²)"),
                            dbc.Input(id='frontal-area', type='number', value=0.68, step=0.1),
                            dbc.Label("Resistência Rolamento"),
                            dbc.Input(id='rolling-resistance', type='number', value=0.015, step=0.001),
                            html.Hr(style={'margin': '12px 0'}),
                            html.Small("Geometria (CG / entre-eixos)",
                                       style={'color': '#6c757d', 'fontWeight': '600'}),
                            dbc.Label("Entre-eixos L (m)"),
                            dbc.Input(id='vehicle-L', type='number', value=1.5, step=0.01),
                            dbc.Label("Altura do CG h (m)"),
                            dbc.Input(id='vehicle-h', type='number', value=0.28, step=0.01),
                            dbc.Label("Dist.x CG (m)"),
                            dbc.Input(id='vehicle-dist-cg', type='number', value=0.6, step=0.01),
                        ],
                        title="🚗  Veículo",
                        item_id="acc-vehicle",
                    ),
                    dbc.AccordionItem(
                        [
                            html.Small("Relação de Transmissão",
                                       style={'color': '#6c757d', 'fontWeight': '600'}),
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Dentes Pinhão Z₁"),
                                    dbc.Input(id='trans-z-pinhao', type='number',
                                              value=None, step=1, min=1,
                                              placeholder='ex: 15'),
                                ], width=6),
                                dbc.Col([
                                    dbc.Label("Dentes Coroa Z₂"),
                                    dbc.Input(id='trans-z-coroa', type='number',
                                              value=None, step=1, min=1,
                                              placeholder='ex: 75'),
                                ], width=6),
                            ], style={'marginBottom': '4px'}),
                            html.Div(id='trans-ratio-display',
                                     style={'color': '#495057', 'fontSize': '12px',
                                            'marginBottom': '6px',
                                            'fontStyle': 'italic'}),
                            dbc.Label("Relação Manual (usado se Z₁/Z₂ vazios)"),
                            dbc.Input(id='final-drive-ratio', type='number',
                                      value=5.0, step=0.1),
                            html.Hr(style={'margin': '10px 0'}),
                            html.Small("Eficiências",
                                       style={'color': '#6c757d', 'fontWeight': '600'}),
                            dbc.Label("Eficiência Corrente"),
                            dbc.Input(id='trans-chain-efficiency', type='number',
                                      value=0.98, step=0.005, min=0.5, max=1.0),
                            dbc.Label("Eficiência Diferencial"),
                            dbc.Input(id='trans-diff-efficiency', type='number',
                                      value=0.97, step=0.005, min=0.5, max=1.0),
                            html.Hr(style={'margin': '10px 0'}),
                            html.Small("LSD",
                                       style={'color': '#6c757d', 'fontWeight': '600'}),
                            dbc.Label("Bias Ratio LSD (Torsen/Clutch-pack)"),
                            dbc.Input(id='trans-lsd-bias', type='number',
                                      value=3.0, step=0.5, min=1.0, max=10.0),
                            html.Hr(style={'margin': '10px 0'}),
                            html.Small("Inércias",
                                       style={'color': '#6c757d', 'fontWeight': '600'}),
                            dbc.Label("Pinhão — lado motor (kg·m²)"),
                            dbc.Input(id='sprocket_inertia', type='number',
                                      value=0.0002, step=0.00005, min=0),
                            dbc.Label("Coroa — lado roda (kg·m²)"),
                            dbc.Input(id='trans-coroa-inertia', type='number',
                                      value=0.0015, step=0.0001, min=0),
                            dbc.Label("Semi-eixo (kg·m²)"),
                            dbc.Input(id='axle_inertia', type='number',
                                      value=0.0015, step=0.0001),
                            dbc.Label("Diferencial (kg·m²)"),
                            dbc.Input(id='diff_inertia', type='number',
                                      value=0.012, step=0.0001),
                        ],
                        title="⚙️  Transmissão",
                        item_id="acc-trans",
                    ),
                    dbc.AccordionItem(
                        [
                            dbc.Row([
                                dbc.Col(dbc.Label("Possui Asa"), width=8,
                                        style={'display': 'flex', 'alignItems': 'center'}),
                                dbc.Col(
                                    dbc.Switch(id='wing-has-wing', value=True,
                                               className='mt-1'),
                                    width=4,
                                    style={'display': 'flex', 'justifyContent': 'flex-end'}),
                            ], className='mb-2'),
                            html.Div(id='wing-params-div', children=[
                                html.Small("Asa Dianteira",
                                           style={'color': '#6c757d', 'fontWeight': '600'}),
                                dbc.Row([
                                    dbc.Col([
                                        dbc.Label("Cl dianteiro"),
                                        dbc.Input(id='wing-cl-front', type='number',
                                                  value=1.0, step=0.1, min=0),
                                    ], width=6),
                                    dbc.Col([
                                        dbc.Label("Área (m²)"),
                                        dbc.Input(id='wing-area-front', type='number',
                                                  value=0.35, step=0.01, min=0),
                                    ], width=6),
                                ], style={'marginBottom': '8px'}),
                                html.Small("Asa Traseira",
                                           style={'color': '#6c757d', 'fontWeight': '600',
                                                  'display': 'block', 'marginBottom': '4px'}),
                                dbc.Row([
                                    dbc.Col([
                                        dbc.Label("Cl traseiro"),
                                        dbc.Input(id='wing-cl-rear', type='number',
                                                  value=1.0, step=0.1, min=0),
                                    ], width=6),
                                    dbc.Col([
                                        dbc.Label("Área (m²)"),
                                        dbc.Input(id='wing-area-rear', type='number',
                                                  value=0.35, step=0.01, min=0),
                                    ], width=6),
                                ]),
                                html.Div(id='wing-downforce-display',
                                         style={'color': '#495057', 'fontSize': '12px',
                                                'marginTop': '8px', 'fontStyle': 'italic'}),
                            ]),
                        ],
                        title="🪂  Aerodinâmica",
                        item_id="acc-aero",
                    ),
                    dbc.AccordionItem(
                        [
                            dbc.Label("Resistência (Ω)"),
                            dbc.Input(id='motor-rs', type='number', value=0.00706, step=0.0001),
                            dbc.Label("Indutância d (H)"),
                            dbc.Input(id='motor-ld', type='number', value=0.0000965, step=0.00001),
                            dbc.Label("Indutância q (H)"),
                            dbc.Input(id='motor-lq', type='number', value=0.0000965, step=0.00001),
                            dbc.Label("Inércia (kg·m²)"),
                            dbc.Input(id='motor-jm', type='number', value=0.02521, step=0.0001),
                            dbc.Label("Atrito (N·m·s)"),
                            dbc.Input(id='motor-kf', type='number', value=0.005, step=0.001),
                            dbc.Label("Fluxo λm (Wb)"),
                            dbc.Input(id='motor-lambda', type='number', value=0.04748, step=0.0001),
                            dbc.Label("Pares de Polos"),
                            dbc.Input(id='motor-poles', type='number', value=10, step=1),
                            dbc.Label("Modulação"),
                            dbc.Input(id='motor-modulation', type='number', value=0.95, step=0.01, min=0, max=1),
                        ],
                        title="⚡  Motor (PMSM)",
                        item_id="acc-motor",
                    ),
                    dbc.AccordionItem(
                        [
                            dbc.Label("Tipo de Célula"),
                            dbc.Select(
                                id='battery-type',
                                options=[{'label': 'Li-ion', 'value': 'Li-ion'},
                                         {'label': 'LiFePO4', 'value': 'LiFePO4'}],
                                value='Li-ion'
                            ),
                            dbc.Label("Células em Série"),
                            dbc.Input(id='n-serie', type='number', value=144, step=1),
                            dbc.Label("Células em Paralelo"),
                            dbc.Input(id='n-paralelo', type='number', value=5, step=1),
                            dbc.Label("SoC Inicial"),
                            dbc.Input(id='soc-inicial', type='number', value=1.0, step=0.01, min=0, max=1),
                        ],
                        title="🔋  Bateria",
                        item_id="acc-battery",
                    ),
                    dbc.AccordionItem(
                        [
                            html.Small("Macro (escalonamento)",
                                       style={'color': '#6c757d', 'fontWeight': '600'}),
                            dbc.Label("λμx — Calspan scrub (~0.55–0.70)"),
                            dbc.Input(id='tire-mu', type='number', value=0.6, step=0.05, min=0.1, max=3.0),
                            dbc.Label("Fz0 — carga nominal (N)"),
                            dbc.Input(id='tire-Fz0', type='number', value=654.0, step=10.0, min=50.0),
                            dbc.Label("Camber γ (rad)"),
                            dbc.Input(id='tire-camber', type='number', value=0.0, step=0.01),
                            html.Hr(style={'margin': '12px 0'}),
                            html.Small("Coeficientes PAC2002 — Fx",
                                       style={'color': '#6c757d', 'fontWeight': '600'}),
                            dbc.Label("PCX1 — shape Cfx"),
                            dbc.Input(id='tire-PCX1', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PCX1'], step=0.01),
                            dbc.Label("PDX1 — μx em Fznom"),
                            dbc.Input(id='tire-PDX1', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PDX1'], step=0.01),
                            dbc.Label("PDX2 — Δμx com carga"),
                            dbc.Input(id='tire-PDX2', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PDX2'], step=0.01),
                            dbc.Label("PDX3 — Δμx com camber"),
                            dbc.Input(id='tire-PDX3', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PDX3'], step=0.01),
                            dbc.Label("PEX1 — curvatura Efx em Fznom"),
                            dbc.Input(id='tire-PEX1', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PEX1'], step=0.01),
                            dbc.Label("PEX2 — ΔE com carga"),
                            dbc.Input(id='tire-PEX2', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PEX2'], step=1e-3),
                            dbc.Label("PEX3 — ΔE com carga²"),
                            dbc.Input(id='tire-PEX3', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PEX3'], step=1e-3),
                            dbc.Label("PEX4 — termo de E em tração"),
                            dbc.Input(id='tire-PEX4', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PEX4'], step=1e-3),
                            dbc.Label("PKX1 — Kfx/Fz em Fznom"),
                            dbc.Input(id='tire-PKX1', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PKX1'], step=0.1),
                            dbc.Label("PKX2 — Δ(Kfx/Fz) com carga"),
                            dbc.Input(id='tire-PKX2', type='number',
                                      value=Tire.Tire.HOOSIER_FSAE_LONG['PKX2'], step=0.1),
                        ],
                        title="🛞  Pneu (PAC2002)",
                        item_id="acc-tire",
                    ),
                    dbc.AccordionItem(
                        [
                            dbc.Label("Distância Máxima (m)"),
                            dbc.Input(id='dist-max', type='number', value=75, step=1, min=1),
                            dbc.Label("Velocidade Ref (rad/s)"),
                            dbc.Input(id='speed-ref', type='number', value=575.96, step=1),
                            dbc.Label("Tempo até o Pico do Pedal (s)"),
                            dbc.Input(id='pedal-peak-time', type='number', value=3.0,
                                      step=0.1, min=0, max=30),
                        ],
                        title="🎯  Simulação",
                        item_id="acc-sim",
                    ),
                ],
                start_collapsed=False,
                always_open=False,
                active_item="acc-vehicle",
                flush=True,
                id='param-accordion',
            ),

            html.Div([
                dbc.Button(
                    [html.I(className="fas fa-play me-2"), "Executar Simulação"],
                    id="run-simulation", n_clicks=0,
                    color="danger", size="lg",
                    style={'width': '100%', 'fontWeight': '600',
                           'backgroundColor': '#830B0B', 'borderColor': '#830B0B'},
                ),
                html.Div(id='sim-status',
                         style={'fontSize': '11px', 'color': '#6c757d',
                                'textAlign': 'center', 'marginTop': '8px', 'minHeight': '16px'}),
            ], style={'padding': '16px 12px', 'borderTop': '1px solid #dee2e6',
                      'backgroundColor': '#f8f9fa', 'position': 'sticky', 'bottom': 0}),

        ], id='sidebar', style={'width': '320px', 'backgroundColor': '#ffffff',
                               'height': '100vh', 'overflowY': 'auto', 'position': 'fixed',
                               'transition': 'left 0.3s ease', 'left': '0',
                               'boxShadow': '2px 0 8px rgba(0,0,0,0.06)',
                               'display': 'flex', 'flexDirection': 'column'}),

        # Área principal com tabs e gráfico
        html.Div([
            html.Div([
                html.H1("Digital-Twin de Powertrain — Xangô-E.Racing",
                        style={'color': '#2c3e50', 'margin': '0', 'fontSize': '22px',
                               'fontWeight': '600'}),
                html.Small("Simulação PMSM + Veículo + Bateria + Pneu",
                           style={'color': '#6c757d'}),
            ], style={'padding': '14px 24px', 'borderBottom': '1px solid #dee2e6',
                      'backgroundColor': '#fff', 'boxShadow': '0 1px 3px rgba(0,0,0,0.04)'}),

            html.Div([
                dbc.Tabs(
                    [
                        dbc.Tab(label="Velocidade & Torque", tab_id="velocity_torque"),
                        dbc.Tab(label="Correntes",           tab_id="currents"),
                        dbc.Tab(label="Tensões",             tab_id="voltages"),
                        dbc.Tab(label="Fluxo & Temp.",       tab_id="flux_temp"),
                        dbc.Tab(label="Controle FOC",        tab_id="control"),
                        dbc.Tab(label="Visão Completa",      tab_id="complete"),
                        dbc.Tab(label="Veículo",             tab_id="vehicle"),
                        dbc.Tab(label="Pneu",                tab_id="tire"),
                        dbc.Tab(label="Bateria",             tab_id="battery"),
                        dbc.Tab(label="📋 Histórico",        tab_id="history"),
                    ],
                    id="plot-tabs",
                    active_tab="velocity_torque",
                ),

                dcc.Loading(
                    id="loading-simulation",
                    type="dot",
                    color="#9B0909",
                    children=html.Div([
                        dcc.Graph(
                            id='motor-graph',
                            style={'height': '82vh', 'borderRadius': '8px',
                                   'boxShadow': '0 2px 6px rgba(0,0,0,0.08)',
                                   'backgroundColor': '#fff', 'marginTop': '12px'},
                            config={'displaylogo': False,
                                    'modeBarButtonsToRemove': ['lasso2d', 'select2d']},
                        ),
                        html.Div(id='history-area', style={'display': 'none'}),
                    ]),
                ),
            ], style={'padding': '16px 24px'}),

        ], id='main-content', style={'transition': 'margin-left 0.3s ease',
                                     'backgroundColor': '#f4f6f8', 'minHeight': '100vh'}),
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
    # Alterado para 2 linhas x 2 colunas para acomodar os novos gráficos
    fig = make_subplots(
        rows=2, cols=2, 
        subplot_titles=(
            "Slip Ratio x Tempo", 
            "Força Longitudinal (Fx) x Slip Ratio",
            "Força Normal Dinâmica (Fz) x Tempo",    # Novo Gráfico
            "Velocidade Angular da Roda x Tempo"     # Novo Gráfico
        ),
        vertical_spacing=0.15,
        horizontal_spacing=0.1
    )

    # --- Gráfico 1: Slip Ratio (Original) ---
    fig.add_trace(
        go.Scatter(x=sim.tempo, y=sim.slip_ratio_hist, name='Slip Ratio', line=dict(color='blue')),
        row=1, col=1
    )
    fig.update_xaxes(title_text="Tempo [s]", row=1, col=1)
    fig.update_yaxes(title_text="Slip Ratio", row=1, col=1)

    # --- Gráfico 2: Curva do Pneu (Fx vs slip) ---
    # Curva TEÓRICA do Pacejka avaliada num Fz fixo (média da operação).
    # Substitui o sort+drop_duplicates do código original, que costurava
    # pontos de instantes diferentes (cada um com seu Fz) numa linha só,
    # gerando o aspecto serrilhado pós-pico.
    if len(sim.slip_ratio_hist) > 0 and sim.tire is not None:
        slip_arr = np.asarray(sim.slip_ratio_hist)
        fz_arr   = np.asarray(sim.fz_hist)
        fz_ref   = float(np.mean(fz_arr))   # Fz médio da operação

        slip_curve = np.linspace(-0.05, max(float(slip_arr.max()), 0.5), 400)
        fx_curve   = [sim.tire.Tire_forces(fz_ref, s) for s in slip_curve]

        fig.add_trace(
            go.Scatter(
                x=slip_curve, y=fx_curve,
                name=f'Pacejka @ Fz={fz_ref:.0f} N (médio)',
                mode='lines', line=dict(color='green', width=2),
            ),
            row=1, col=2
        )
    fig.update_xaxes(title_text="Slip Ratio", row=1, col=2)
    fig.update_yaxes(title_text="Força Longitudinal (N)", row=1, col=2)

    # --- Gráfico 3: Fz Dinâmico x Tempo (Novo) ---
    # Usa o histórico de Fz calculado na simulação (sim.fz_hist)
    fig.add_trace(
        go.Scatter(x=sim.tempo, y=sim.fz_hist, name='Fz Dinâmico', line=dict(color='orange')),
        row=2, col=1
    )
    fig.update_xaxes(title_text="Tempo [s]", row=2, col=1)
    fig.update_yaxes(title_text="Força Normal (N)", row=2, col=1)

    # --- Gráfico 4: Velocidade Angular da Roda x Tempo (Novo) ---
    # Cálculo: Motor RPM -> Roda RPM -> Roda Rad/s
    if sim.transmission:
        ratio = sim.transmission.final_drive_ratio
        # sim.velocidade está em RPM (do motor)
        motor_rpm = np.array(sim.velocidade)
        wheel_rpm = motor_rpm / ratio
        wheel_rads = wheel_rpm * (2 * np.pi / 60) # Converte RPM para Rad/s
    else:
        wheel_rads = [0] * len(sim.tempo)

    fig.add_trace(
        go.Scatter(x=sim.tempo, y=wheel_rads, name='Vel. Angular Roda', line=dict(color='purple')),
        row=2, col=2
    )
    fig.update_xaxes(title_text="Tempo [s]", row=2, col=2)
    fig.update_yaxes(title_text="Velocidade Angular (rad/s)", row=2, col=2)

    # Ajuste de Layout
    fig.update_layout(height=900, title_text="Análise Detalhada do Pneu", template="plotly_white")
    
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
    base_sidebar = {'width': '320px', 'backgroundColor': '#ffffff',
                    'height': '100vh', 'overflowY': 'auto', 'position': 'fixed',
                    'transition': 'left 0.3s ease',
                    'boxShadow': '2px 0 8px rgba(0,0,0,0.06)',
                    'display': 'flex', 'flexDirection': 'column'}
    base_content = {'transition': 'margin-left 0.3s ease',
                    'backgroundColor': '#f4f6f8', 'minHeight': '100vh'}

    if n_clicks is None:
        visible = True
    else:
        visible = not sidebar_state.get('visible', True)

    if visible:
        sidebar_style = {**base_sidebar, 'left': '0'}
        content_style = {**base_content, 'marginLeft': '320px'}
        icon = html.I(className='fas fa-chevron-right')
    else:
        sidebar_style = {**base_sidebar, 'left': '-320px'}
        content_style = {**base_content, 'marginLeft': '0'}
        icon = html.I(className='fas fa-chevron-right', style={'transform': 'rotate(180deg)'})

    return sidebar_style, content_style, icon, {'visible': visible}

# ---------------------
# Helper: nomes amigáveis para o diff de parâmetros
# ---------------------
_PARAM_LABELS = {
    'vehicle.mass': 'Massa', 'vehicle.wheel_radius': 'Raio Roda',
    'vehicle.wheel_mass': 'Massa Roda', 'vehicle.drag_coeff': 'Cd',
    'vehicle.frontal_area': 'Área Frontal',
    'vehicle.rolling_resistance': 'Rolamento',
    'vehicle.L': 'Entre-eixos', 'vehicle.h': 'h CG', 'vehicle.dist_cg': 'Dist. CG',
    'transmission.final_drive_ratio': 'Relação Final',
    'transmission.Z_pinhao': 'Z Pinhão', 'transmission.Z_coroa': 'Z Coroa',
    'transmission.chain_efficiency': 'η Corrente',
    'transmission.diff_efficiency': 'η Diferencial',
    'transmission.axle_inertia': 'J semi-eixo',
    'transmission.diff_inertia': 'J diff',
    'transmission.sprocket_inertia': 'J pinhão',
    'transmission.coroa_inertia': 'J coroa',
    'transmission.lsd_bias_ratio': 'LSD bias',
    'aero.has_wing': 'Possui Asa',
    'aero.lift_coeff_front': 'Cl dianteiro', 'aero.area_front': 'Área front',
    'aero.lift_coeff_rear': 'Cl traseiro',   'aero.area_rear': 'Área rear',
    'battery.tipo_celula': 'Tipo Célula', 'battery.n_serie': 'Células Série',
    'battery.n_paralelo': 'Células Paral.', 'battery.soc_inicial': 'SoC Inicial',
    'motor.rs': 'Rs', 'motor.ld': 'Ld', 'motor.lq': 'Lq', 'motor.jm': 'Jm',
    'motor.kf': 'Kf', 'motor.lambda_m': 'λm', 'motor.p': 'Pares Polos',
    'motor.valor_mu': 'Modulação', 'motor.velocidade_ref': 'ω ref',
    'simulacao.distancia_max': 'Dist. Máx.', 'simulacao.pedal_t_pico': 't Pico Pedal (s)',
    'tire.tire_friction_coef': 'λμx', 'tire.Fz0': 'Fz0', 'tire.camber': 'Camber γ',
}

def _flatten(d, prefix=''):
    out = {}
    for k, v in (d or {}).items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(_flatten(v, key))
        else:
            out[key] = v
    return out

def _diff_params(old, new):
    """Retorna lista [(label, old_val, new_val)] de parâmetros alterados.

    Skips 'inversor' (não usado) e PAC2002 coefs internos (ruído visual)."""
    if old is None:
        return [('—', None, None)]
    old_flat, new_flat = _flatten(old), _flatten(new)
    changes = []
    for key, new_v in new_flat.items():
        if key.startswith('inversor.') or '.pacejka_params.' in key:
            continue
        old_v = old_flat.get(key)
        if old_v != new_v:
            label = _PARAM_LABELS.get(key, key)
            changes.append((label, old_v, new_v))
    return changes

def _fmt_val(v):
    if v is None:
        return '—'
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)

def _summarize_changes(changes):
    if not changes:
        return html.Span('— (sem alteração)', style={'color': '#6c757d'})
    if changes[0][0] == '—' and changes[0][1] is None:
        return html.Span('primeira simulação', style={'color': '#6c757d',
                                                      'fontStyle': 'italic'})
    first = changes[0]
    text = f"{first[0]}: {_fmt_val(first[1])} → {_fmt_val(first[2])}"
    if len(changes) > 1:
        text += f"  (+{len(changes)-1})"
    return text

# ---------------------
# Callback unificado: lê inputs, atualiza store, roda simulação, gera figuras
# ---------------------
@app.callback(
    [Output('stored-figures', 'data'),
     Output('simulation-run-flag', 'data'),
     Output('parameters-store', 'data'),
     Output('sim-status', 'children'),
     Output('sim-history', 'data'),
     Output('last-params', 'data')],
    [Input('run-simulation', 'n_clicks')],
    [State('vehicle-mass', 'value'),
     State('wheel-radius', 'value'),
     State('wheel-mass', 'value'),
     State('drag-coeff', 'value'),
     State('frontal-area', 'value'),
     State('rolling-resistance', 'value'),
     State('vehicle-L', 'value'),
     State('vehicle-h', 'value'),
     State('vehicle-dist-cg', 'value'),
     State('final-drive-ratio', 'value'),
     State('axle_inertia', 'value'),
     State('diff_inertia', 'value'),
     State('sprocket_inertia', 'value'),
     State('trans-z-pinhao', 'value'),
     State('trans-z-coroa', 'value'),
     State('trans-chain-efficiency', 'value'),
     State('trans-diff-efficiency', 'value'),
     State('trans-lsd-bias', 'value'),
     State('trans-coroa-inertia', 'value'),
     State('wing-has-wing', 'value'),
     State('wing-cl-front', 'value'),
     State('wing-area-front', 'value'),
     State('wing-cl-rear', 'value'),
     State('wing-area-rear', 'value'),
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
     State('dist-max', 'value'),
     State('speed-ref', 'value'),
     State('pedal-peak-time', 'value'),
     State('tire-mu', 'value'),
     State('tire-Fz0', 'value'),
     State('tire-camber', 'value'),
     State('tire-PCX1', 'value'),
     State('tire-PDX1', 'value'),
     State('tire-PDX2', 'value'),
     State('tire-PDX3', 'value'),
     State('tire-PEX1', 'value'),
     State('tire-PEX2', 'value'),
     State('tire-PEX3', 'value'),
     State('tire-PEX4', 'value'),
     State('tire-PKX1', 'value'),
     State('tire-PKX2', 'value'),
     State('parameters-store', 'data'),
     State('sim-history', 'data'),
     State('last-params', 'data')]
)
def run_simulation_once(n_clicks,
                        mass, wheel_radius, wheel_mass, drag_coeff, frontal_area, rolling_resistance,
                        L, h, dist_cg,
                        final_drive_ratio, axle_inertia, diff_inertia, sprocket_inertia,
                        trans_z_pinhao, trans_z_coroa,
                        trans_chain_eff, trans_diff_eff, trans_lsd_bias, trans_coroa_inertia,
                        wing_has_wing, wing_cl_front, wing_area_front, wing_cl_rear, wing_area_rear,
                        battery_type, n_serie, n_paralelo, soc_inicial,
                        motor_rs, motor_ld, motor_lq, motor_jm, motor_kf, motor_lambda, motor_poles, motor_modulation,
                        dist_max, speed_ref, pedal_peak_time,
                        tire_mu, tire_Fz0, tire_camber,
                        tire_PCX1, tire_PDX1, tire_PDX2, tire_PDX3,
                        tire_PEX1, tire_PEX2, tire_PEX3, tire_PEX4,
                        tire_PKX1, tire_PKX2,
                        current_params, history, last_params):
    if n_clicks is None or n_clicks == 0:
        return (dash.no_update, {'run': False}, current_params, "",
                history or [], last_params)

    vehicle_params = {
        'mass': mass, 'wheel_radius': wheel_radius, 'wheel_mass': wheel_mass,
        'drag_coeff': drag_coeff, 'frontal_area': frontal_area,
        'rolling_resistance': rolling_resistance, 'road_grade': 0,
        'L': L, 'h': h, 'dist_cg': dist_cg,
    }
    # Se Z₁ e Z₂ fornecidos, calculamos N a partir dos dentes;
    # caso contrário usamos final_drive_ratio manual.
    _use_teeth = (trans_z_pinhao and trans_z_coroa and
                  int(trans_z_pinhao) > 0 and int(trans_z_coroa) > 0)
    transmission_params = {
        'final_drive_ratio': final_drive_ratio,
        'Z_pinhao': int(trans_z_pinhao) if _use_teeth else None,
        'Z_coroa':  int(trans_z_coroa)  if _use_teeth else None,
        'chain_efficiency': trans_chain_eff or 0.98,
        'diff_efficiency':  trans_diff_eff  or 0.97,
        'axle_inertia': axle_inertia,
        'diff_inertia': diff_inertia,
        'sprocket_inertia': sprocket_inertia,
        'coroa_inertia': trans_coroa_inertia or 0.0015,
        'lsd_bias_ratio': trans_lsd_bias or 3.0,
    }
    battery_params = {
        'tipo_celula': battery_type, 'n_serie': n_serie,
        'n_paralelo': n_paralelo, 'soc_inicial': soc_inicial,
    }
    motor_params = {
        'rs': motor_rs, 'ld': motor_ld, 'lq': motor_lq, 'jm': motor_jm,
        'kf': motor_kf, 'lambda_m': motor_lambda, 'p': motor_poles,
        'valor_mu': motor_modulation, 'velocidade_ref': speed_ref,
    }
    simulacao_params = {'distancia_max': dist_max, 'pedal_t_pico': pedal_peak_time}
    tire_params = {
        'pacejka_params': {
            'PCX1': tire_PCX1, 'PDX1': tire_PDX1, 'PDX2': tire_PDX2, 'PDX3': tire_PDX3,
            'PEX1': tire_PEX1, 'PEX2': tire_PEX2, 'PEX3': tire_PEX3, 'PEX4': tire_PEX4,
            'PKX1': tire_PKX1, 'PKX2': tire_PKX2, 'PKX3': 0.0,
        },
        'tire_friction_coef': tire_mu,
        'Fz0': tire_Fz0,
        'camber': tire_camber,
    }

    aero_params = {
        'has_wing':        bool(wing_has_wing),
        'lift_coeff_front': wing_cl_front   or 1.0,
        'area_front':       wing_area_front or 0.35,
        'lift_coeff_rear':  wing_cl_rear    or 1.0,
        'area_rear':        wing_area_rear  or 0.35,
    }

    new_params = {
        'vehicle': vehicle_params,
        'transmission': transmission_params,
        'aero': aero_params,
        'battery': battery_params,
        'motor': motor_params,
        'inversor': current_params.get('inversor', {}),
        'simulacao': simulacao_params,
        'tire': tire_params,
    }

    transmission = Transmission.Transmission(**transmission_params)
    vehicle = Vehicle.Vehicle(**vehicle_params, **aero_params)
    
    battery = BatteryPack.BatteryPack(**battery_params)
    tire = Tire.Tire(**tire_params)

    motor = Motor.Motor(
        rs=motor_params['rs'],
        ld=motor_params['ld'],
        lq=motor_params['lq'],
        jm=motor_params['jm'],
        kf=motor_params['kf'],
        lambda_m=motor_params['lambda_m'],
        p=motor_params['p'],
        valor_mu=motor_params['valor_mu'],
        speed_ref=motor_params['velocidade_ref'],
        max_current=323.0,   # [A pico] → T_pico ≈ 230 N·m (datasheet EMRAX 228)
    )
    
    # ── Perfil do pedal — smoothstep 0 → 100 % no instante escolhido ──────
    # f(s) = 3s²−2s³: derivada zero no início e no fim, evita degrau de taxa
    # que excita oscilações no PI de velocidade ao fim de uma rampa linear.
    # "Tempo até o Pico do Pedal" (sidebar "🎯 Simulação") define o instante
    # em que o pedal atinge 100 %; daí em diante segura no máximo.
    # t_pico ≈ 0 → degrau (pedal a 100 % desde t=0).
    pedal = Pedal.Pedal(ganho=1.0)
    _t_pico = float(pedal_peak_time) if pedal_peak_time else 3.0
    if _t_pico > 1e-3:
        pedal_tempos, pedal_throttle = Pedal.Pedal.perfil_suave(0.0, _t_pico)
        pedal.set_profile(pedal_tempos, pedal_throttle)
    else:
        pedal.set_posicao(1.0)   # degrau: 100 % desde o início

    sim = Simulation(
        motor=motor, vehicle=vehicle, transmission=transmission,
        battery=battery, tire=tire,
        pedal=pedal,
        tmax=60.0, dmax=simulacao_params['distancia_max'], steps=12000,
        p_max_dc=80e3,   # FSAE EV.4.1 — limite de potência no barramento DC
    )

    sim.simulate()

    figures = {
        'velocity_torque': create_velocity_torque_plot(sim),
        'currents':        create_currents_plot(sim),
        'voltages':        create_voltages_plot(sim),
        'flux_temp':       create_flux_temp_plot(sim),
        'control':         create_control_plot(sim),
        'complete':        create_complete_plot(sim),
        'vehicle':         create_vehicle_plot(sim),
        'tire':            create_tire_plot(sim),
        'battery':         create_battery_plot(sim),
    }
    stored_figures = {key: fig.to_dict() for key, fig in figures.items()}

    t_final = float(sim.tempo[-1]) if len(sim.tempo) else 0.0
    d_final = float(sim.vehicle_position[-1]) if len(sim.vehicle_position) else 0.0
    status = f"✓ {t_final:.2f} s para {d_final:.1f} m  ·  run #{n_clicks}"

    changes = _diff_params(last_params, new_params)
    history = list(history or [])
    history.append({
        'n': n_clicks,
        'time': t_final,
        'dist': d_final,
        'changes': [(lbl, _fmt_val(ov), _fmt_val(nv)) for lbl, ov, nv in changes],
    })

    return (stored_figures, {'run': True}, new_params, status,
            history, new_params)

# ---------------------
# Render do histórico em tabela
# ---------------------
def _render_history_table(history):
    if not history:
        return html.Div(
            "Nenhuma simulação executada ainda.",
            style={'textAlign': 'center', 'padding': '60px 20px',
                   'color': '#6c757d', 'fontStyle': 'italic'},
        )

    header = html.Thead(html.Tr([
        html.Th("#",        style={'width': '60px'}),
        html.Th("Tempo (s)", style={'width': '110px'}),
        html.Th("Dist. (m)", style={'width': '110px'}),
        html.Th("Parâmetro alterado"),
    ]))

    rows = []
    for entry in reversed(history):
        changes = entry.get('changes') or []
        summary = _summarize_changes(
            [(lbl, ov, nv) for lbl, ov, nv in changes]
            if changes else []
        )
        rows.append(html.Tr([
            html.Td(f"#{entry['n']}",         style={'fontWeight': '600'}),
            html.Td(f"{entry['time']:.3f}"),
            html.Td(f"{entry['dist']:.1f}"),
            html.Td(summary),
        ]))

    return html.Div([
        html.Div([
            html.H5(f"Histórico de Simulações ({len(history)})",
                    style={'margin': '0', 'color': '#2c3e50'}),
            html.Small("Mais recentes no topo. Coluna 'Parâmetro alterado' "
                       "mostra a diferença vs. a run anterior.",
                       style={'color': '#6c757d'}),
        ], style={'padding': '14px 18px', 'borderBottom': '1px solid #dee2e6',
                  'backgroundColor': '#f8f9fa'}),
        dbc.Table([header, html.Tbody(rows)],
                  bordered=False, hover=True, striped=True,
                  responsive=True, size='sm',
                  style={'marginBottom': '0'}),
    ], style={'backgroundColor': '#fff', 'borderRadius': '8px',
              'boxShadow': '0 2px 6px rgba(0,0,0,0.08)',
              'marginTop': '12px', 'overflow': 'hidden'})

# ---------------------
# Callback de render: alterna gráfico ↔ histórico conforme a aba
# ---------------------
_GRAPH_STYLE_VISIBLE = {'height': '82vh', 'borderRadius': '8px',
                        'boxShadow': '0 2px 6px rgba(0,0,0,0.08)',
                        'backgroundColor': '#fff', 'marginTop': '12px'}
_GRAPH_STYLE_HIDDEN = {'display': 'none'}

@app.callback(
    [Output('motor-graph', 'figure'),
     Output('motor-graph', 'style'),
     Output('history-area', 'children'),
     Output('history-area', 'style')],
    [Input('plot-tabs', 'active_tab'),
     Input('stored-figures', 'data'),
     Input('sim-history', 'data')],
)
def render_active_tab(active_tab, stored_figures, history):
    if active_tab == 'history':
        return (go.Figure(), _GRAPH_STYLE_HIDDEN,
                _render_history_table(history or []),
                {'display': 'block'})

    if not stored_figures or not active_tab:
        fig = go.Figure(layout={
            'template': 'plotly_white',
            'annotations': [{
                'text': 'Clique em <b>Executar Simulação</b> para gerar os gráficos.',
                'xref': 'paper', 'yref': 'paper', 'x': 0.5, 'y': 0.5,
                'showarrow': False, 'font': {'size': 16, 'color': '#6c757d'},
            }],
            'xaxis': {'visible': False}, 'yaxis': {'visible': False},
        })
    else:
        fig = stored_figures.get(active_tab, {})

    return fig, _GRAPH_STYLE_VISIBLE, [], {'display': 'none'}

# ---------------------
# Callback: exibe relação calculada dos dentes
# ---------------------
@app.callback(
    Output('trans-ratio-display', 'children'),
    [Input('trans-z-pinhao', 'value'),
     Input('trans-z-coroa', 'value'),
     Input('final-drive-ratio', 'value')],
)
def update_ratio_display(z1, z2, ratio_manual):
    if z1 and z2 and int(z1) > 0:
        n = int(z2) / int(z1)
        return f"Relação calculada: {n:.3f} (Z₂/Z₁ = {int(z2)}/{int(z1)})"
    if ratio_manual:
        return f"Relação manual: {ratio_manual:.3f}"
    return "Preencha Z₁ e Z₂ ou a relação manual"


# ---------------------
# Callback: asa — oculta parâmetros e mostra downforce estimado
# ---------------------
@app.callback(
    [Output('wing-params-div', 'style'),
     Output('wing-downforce-display', 'children')],
    [Input('wing-has-wing', 'value'),
     Input('wing-cl-front', 'value'),
     Input('wing-area-front', 'value'),
     Input('wing-cl-rear', 'value'),
     Input('wing-area-rear', 'value')],
)
def update_wing_ui(has_wing, cl_front, area_front, cl_rear, area_rear):
    style = {} if has_wing else {'display': 'none'}
    if not has_wing or None in (cl_front, area_front, cl_rear, area_rear):
        return style, ""
    rho = 1.225
    v60  = 60  / 3.6
    v80  = 80  / 3.6
    v100 = 100 / 3.6
    def _f(v):
        return 0.5 * rho * ((cl_front * area_front) + (cl_rear * area_rear)) * v ** 2
    F_rear_60  = 0.5 * rho * cl_rear * area_rear * v60  ** 2
    F_rear_80  = 0.5 * rho * cl_rear * area_rear * v80  ** 2
    balance_80 = F_rear_80 / max(_f(v80), 1e-6) * 100
    msg = (f"@ 60 km/h: {_f(v60):.0f} N  |  "
           f"@ 80 km/h: {_f(v80):.0f} N ({balance_80:.0f}% traseiro)  |  "
           f"@ 100 km/h: {_f(v100):.0f} N")
    return style, msg


if __name__ == '__main__':
    print("")
    app.run(debug=True)