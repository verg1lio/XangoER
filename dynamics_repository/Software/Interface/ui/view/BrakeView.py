import flet as ft
from Software.Interface.ui.page.caixa_texto import caixa_texto
from Software.Interface.controllers import BrakeController
from Software.Interface.ui.page.select_model import head
from Software.Interface.ui.page.buttons import save_button, back_button, home_button, next_button, footer_buttons
from Software.Interface.controllers.envio_texto import envio
import time
import numpy as np

controller = BrakeController()

mensagem = envio()

def home_func(e, page):
    """
    Navega para a página inicial ("/").

    Args:
        e: Evento de clique.
        page (ft.Page): Página atual do Flet.
    """

    page.go("/")


def create_brake(page: ft.Page) -> ft.View:

    """
    Cria a view para cadastro de um novo modelo do sistema de freio.

    Exibe vários campos para entrada de parâmetros do modelo, com botões para salvar,
    voltar e ir para home.

    Args:
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com os campos e botões para criação do modelo de freio.
    """
    
    
    status_text = ft.Text("")  # Mensagem dinâmica de feedback

    # Caixas de texto (em ordem reversa porque vão sendo empilhadas no topo do Column)
    c_rr = caixa_texto("Coeficiente de resistência ao rolamento", page, mensagem)
    dist_eixos = caixa_texto("Distância entre eixos", page, mensagem, proximo=c_rr)
    m_tire = caixa_texto("Massa do pneu", page, mensagem, proximo=dist_eixos)
    m_wheel = caixa_texto("Massa da roda", page, mensagem, proximo=m_tire)
    massa_total = caixa_texto("Massa total do veículo", page, mensagem, proximo=m_wheel)
    raio_efetivo_disco = caixa_texto("Raio efetivo do disco", page, mensagem, proximo=massa_total)
    atrito_coeficiente = caixa_texto("Coeficiente de atrito da pastilha", page, mensagem, proximo=raio_efetivo_disco)
    Npast = caixa_texto("Número de pastilhas por disco", page, mensagem, proximo=atrito_coeficiente)
    Dwc = caixa_texto("Diâmetro do cilindro da roda em metros", page, mensagem, proximo=Npast)
    Dcm = caixa_texto("Diâmetro do Cilindro Mestre", page, mensagem, proximo=Dwc)
    raio_pneu = caixa_texto("Raio do Pneu", page, mensagem, proximo=Dcm)
    FzR = caixa_texto("Força de reação estática na traseira", page, mensagem, proximo=raio_pneu)
    FzF = caixa_texto("Força de reação estática na dianteira", page, mensagem, proximo=FzR)
    coef_atr_pneuSolo = caixa_texto("Coeficiente de atrito do pneu/solo", page, mensagem, proximo=FzF)
    centro_grav_y = caixa_texto("Altura do centro de gravidade", page, mensagem, proximo=coef_atr_pneuSolo)
    coef_at_pastdisc = caixa_texto("Coeficiente de atrito do contato pastilha/disco", page, mensagem, proximo=centro_grav_y)
    psi = caixa_texto("Distribuição de carga estática por eixo", page, mensagem, proximo=coef_at_pastdisc)
    desaceleracao = caixa_texto("Desaceleração", page, mensagem, proximo=psi)
    RedP = caixa_texto("Redução do Pedal", page, mensagem, proximo=desaceleracao)

    inputs_column = ft.Column(
        controls=[
            RedP, desaceleracao, psi, coef_at_pastdisc, centro_grav_y,
            coef_atr_pneuSolo, FzF, FzR, raio_pneu, Dcm, Dwc,
            Npast, atrito_coeficiente, raio_efetivo_disco, massa_total,
            m_wheel, m_tire, dist_eixos, c_rr
        ],
        height=page.height - 180,
        scroll=ft.ScrollMode.ADAPTIVE
    )

    # Função que atualiza texto de status e volta para a página inicial
    def modelo_criado():
        status_text.value = "Modelo Salvo!"
        view.controls = [ft.Container(content=status_text,alignment=ft.alignment.center,padding=20,expand=True)]

    
        page.update()
        time.sleep(2)
        page.go("/")



    def erro_criacao():
        status_text.value = f"Valor(es) Inválido(s)!"
        view.controls = [ft.Container(content=status_text,alignment=ft.alignment.center,padding=20,expand=True)]
        page.update()
        time.sleep(2)
        page.go("/create")



    def instance_brake(e):
        controller.instance_brake(
            dicionario=mensagem.dicionario,
            instanciado=modelo_criado,
            falha=erro_criacao
        )
    def back_create(e):
        page.go("/create")
    

    save = save_button(func=instance_brake)
    back = back_button(func=back_create)
    home = home_button(lambda e:home_func(e, page))


    # View retornada
    view = ft.View(
        route="/create_brake",
        controls=[
            ft.Text("Cadastro de Sistema de Freio", size=24, weight=ft.FontWeight.BOLD),
            inputs_column,
            ft.Row(controls=[back, home, save], alignment=ft.MainAxisAlignment.END)
        ],
        vertical_alignment=ft.MainAxisAlignment.START
    )

    return view


def run_brake(page:ft.Page) -> ft.View:

    """
    View para seleção do modelo de freio a ser executado.

    Args:
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com interface para escolha do modelo Brake a executar.
    """

    view = ft.View(
        route = "/run_brake",
        controls=[head("Brake", "run_brake_exec", run_brake_exec, page)]
    )

    return view


def run_brake_exec(model, page:ft.Page) -> ft.View:

    """
    View para entrada de dados para execução do modelo de freio.

    Apresenta campos para força do pedal, tempo e velocidade do veículo,
    além de botões para navegação.

    Args:
        model: Modelo de freio instanciado.
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com campos para execução do modelo e navegação.
    """

    pedal_forces = caixa_texto(hint_text="Força do pedal", view = page, variable=mensagem)
    
    tempo = caixa_texto(hint_text="Tempo", view = page, proximo=pedal_forces, variable=mensagem)
    
    velocidade = caixa_texto(hint_text="Velocidade do veículo", proximo=tempo, view = page, variable=mensagem)

    dicionario = mensagem.dicionario

    tela = ft.Row(
        controls = [
            velocidade, tempo, pedal_forces
        ], expand = 4
    )
    
    def back_func(e):

        page.go("/run_brake")

    def next_page(e):
        
        controller.result_data(dicionario=mensagem.dicionario, func=dados_results, page=page, brake_model=model)
    
    footer = footer_buttons(page=page,back_func=back_func, home_func=home_func,next_page=next_page)

    view = ft.View(route="/run_brake_exec",
                        controls=[tela, footer])

    return view



def dados_results(brake_model, pedal_forces, page:ft.Page) -> ft.View:

    """
    Exibe os resultados do cálculo do modelo de freio baseado na força do pedal.

    Inclui área para texto do resultado e botões para navegação entre telas.

    Args:
        brake_model: Instância do modelo de freio.
        pedal_forces (float): Valor da força do pedal.
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com os resultados e navegação.
    """
    
    text_result = brake_model.results(pedal_forces)
    
    content_area = ft.Container(
        content=ft.Column(
            controls=[
                ft.Text(value=text_result, no_wrap=False, selectable=True)
            ],
            alignment=ft.MainAxisAlignment.CENTER,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
            scroll=ft.ScrollMode.AUTO,
            expand=True,
        ),
        alignment=ft.alignment.center,
        padding=20,
        expand=True,
    )

    def back_func(e):
        page.go("/run_brake_exec")

    def home_func(e):
        page.go("/")
    
    def next_page(e):
        controller.graph(brake_model=brake_model, pedal_forces=pedal_forces, page=page,func=graph_brake)


    # Rodapé com botões distribuídos
    footer_buttons = ft.Row(
        controls=[
            back_button(func=back_func),
            home_button(home_func),
            next_button(func=next_page),
        ],
        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
        vertical_alignment=ft.CrossAxisAlignment.CENTER,
    )

    # Coluna principal com o conteúdo e o rodapé fixo
    return ft.View(
        route="/run_brake/result_data",
        controls=[
            ft.Column(
                controls=[
                    content_area,
                    footer_buttons
                ],
                expand=True
            )
        ]
    )

    
def graph_brake(brake_model, pedal_forces, page):

    """
    Exibe gráficos gerados a partir do modelo de freio e força do pedal.

    Mostra imagens base64 do gráfico e oferece navegação para tela anterior e home.

    Args:
        brake_model: Instância do modelo de freio.
        pedal_forces (float): Valor da força do pedal.
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com gráficos e botões de navegação.
    """
    
    time_intervals = np.arange(0, 4, 0.001)

    #recebendo os dados gráficos do modelo instanciado
    img_64 = brake_model.graph_2(1000, time_intervals, 10)

    #Criando uma função para voltar aos dados anteriores
    def voltar(e):
        page.go("/run_brake/result_data")
  
    #Cria uma row com paara receber os gráficos
    images = ft.Row(expand=2,
                    wrap=False,
                    scroll=ft.ScrollMode.ALWAYS,
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.START)
    

    tela = ft.Column(controls=[ft.Container(content=images, expand=True)])
    
    footer = footer_buttons(page=page, back_func=voltar, next_visibility=False, home_func=home_func)

    for img_data in img_64:
        images.controls.append(
            ft.Container(
                content=ft.Image(
                src_base64=img_data,
                fit=ft.ImageFit.CONTAIN,
                expand=True,
                repeat=ft.ImageRepeat.NO_REPEAT,
                border_radius=ft.border_radius.all(10),
            ),
            width=450,
            height=450
            )
            
        )
    return ft.View(route="/run_brake/graph",
                   controls = [
                       tela,
                       footer
                   ])
