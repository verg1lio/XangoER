import flet as ft
import time
from Software.Interface.ui.page.rolagem_tela import create_row_title, create_row_with_columns
from Software.Interface.ui.page.select_model import head
from Software.Interface.ui.page.caixa_texto import caixa_texto
from Software.Interface.controllers import DrivetrainController
from Software.Interface.controllers.envio_texto import envio
from Software.Interface.ui.page.buttons import save_button, home_button, back_button, footer_buttons



controller = DrivetrainController()

mensagem = envio()


def home_func(e, page):
    """
    Navega para a página inicial ("/").

    Args:
        e: Evento de clique.
        page (ft.Page): Página atual do Flet.
    """
    page.go("/")

def create_drivetrain(page: ft.Page) -> ft.View:

    """
    Cria a view para cadastro de um novo modelo de transmissão.

    Exibe campos encadeados para entrada dos parâmetros do modelo e botões para salvar,
    voltar e home. Exibe mensagens de status para sucesso ou erro.

    Args:
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com formulário para cadastro da transmissão.
    """

    status_text = ft.Text("")  # texto para feedback (erro ou sucesso)

    # Criação dos campos de entrada encadeados
    cp = caixa_texto(hint_text="Relação Coroa-Pinhão", variable=mensagem, view=page)
    torque = caixa_texto(hint_text="Torque", variable=mensagem, view=page, proximo=cp)
    rpm = caixa_texto(hint_text="RPM", proximo=torque, variable=mensagem, view=page)
    red_unica = caixa_texto(hint_text="Redução Única", proximo=rpm, variable=mensagem, view=page)
    red_primaria = caixa_texto(hint_text="Redução Primária", proximo=red_unica, variable=mensagem, view=page)
    aceleracao_ideal = caixa_texto(hint_text="Aceleração Ideal", proximo=red_primaria, variable=mensagem, view=page)
    raio_pneu = caixa_texto(hint_text="Raio do Pneu", proximo=aceleracao_ideal, variable=mensagem, view=page)
    coef_atrito = caixa_texto(hint_text="Coeficiente de Atrito", proximo=raio_pneu, variable=mensagem, view=page)
    distancia_eixos = caixa_texto(hint_text="Distância entre Eixos", proximo=coef_atrito, variable=mensagem, view=page)
    massa = caixa_texto(hint_text="Massa do Veículo", variable=mensagem, view=page, proximo=distancia_eixos)
    centro_gravidade_y = caixa_texto(hint_text="Centro de Gravidade em y", variable=mensagem, view=page, proximo=massa)
    centro_gravidade_x = caixa_texto(hint_text="Centro de Gravidade em X", variable=mensagem, view=page, proximo=centro_gravidade_y)

    inputs_column = ft.Column(
        controls=[
            centro_gravidade_x,
            centro_gravidade_y,
            massa,
            distancia_eixos,
            coef_atrito,
            raio_pneu,
            aceleracao_ideal,
            red_primaria,
            red_unica,
            rpm,
            torque,
            cp
        ],
        height=page.height - 180,
        scroll=ft.ScrollMode.ADAPTIVE
    )

    status_container = ft.Container(content=status_text, alignment=ft.alignment.center)

    def modelo_criado():
        status_text.value = "Modelo Salvo!"
        view.controls[:] = [status_container]
        page.update()
        time.sleep(3)
        page.go("/")

    def erro_criacao():
        status_text.value = "Valor(es) Inválido(s)!"
        view.controls[:] = [status_container]
        page.update()
        time.sleep(3)
        # volta para a mesma view
        page.go("/")


    def instance_dt(e):
        controller.instance_dt(
            dicionario=mensagem.dicionario,
            instanciado=modelo_criado,
            falha=erro_criacao
        )
    
    def back_create(e):
        page.go("/create")
    
    def home_func(e):
        page.go("/")

    save = save_button(func=instance_dt)
    back = back_button(func=back_create)
    home = home_button(home_func)



    # Criando a View com rota definida
    view = ft.View(
        route="/create_drivetrain",
        controls=[
            ft.Text("Cadastro de Transmissão", size=24, weight=ft.FontWeight.BOLD),
            inputs_column,
            ft.Row(controls=[back, home, save], alignment=ft.MainAxisAlignment.END)
        ],
        vertical_alignment=ft.MainAxisAlignment.START
    )

    return view



def run_dt(page:ft.Page) -> ft.View:

    """
    View para seleção do modelo de transmissão a ser executado.

    Args:
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com interface para escolha do modelo Drivetrain a executar.
    """

    view = ft.View(
        route = "/run_drivetrain",
        controls=[head("Drivetrain", "run_dt_exec", run_dt_exec, page)]
    )

    return view

def run_dt_exec(dt_model, page):

    """
    View para entrada de dados para execução do modelo de transmissão.

    Apresenta campos para torque final e RPM final, com botões de navegação.

    Args:
        dt_model: Modelo de transmissão instanciado.
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com campos para execução do modelo e navegação.
    """

    torque_final = caixa_texto(hint_text="Torque Final", view = page, variable=mensagem)
    rpm_final = caixa_texto(hint_text="RPM Final", view = page, variable=mensagem,proximo=torque_final)

    tela = ft.Row(controls=[rpm_final, torque_final])

    def voltar(e):
        page.go("/run")

    def info_tela(e):

        controller.info_func(dt_model, mensagem.dicionario, info, page)
        """dt_model.new_rpm = int(mensagem.dicionario["RPM Final"])

        page.go("")
        info(dt_model, page)"""

    footer = ft.Row(controls=[footer_buttons(page=page, back_func=voltar, home_func=home_func, next_page=info_tela)],
                    expand=4)

    view = ft.View(route="/run_dt_exec",
                   controls = [tela, footer])

    return view

def info(dt_model, page):

    """
    Exibe as informações/resultados do modelo de transmissão.

    Inclui texto com resultados e botões para navegar entre as telas.

    Args:
        dt_model: Modelo de transmissão instanciado.
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com informações do modelo e navegação.
    """

    info = dt_model.showResultsInfo()

    tela = ft.Column(controls=
                [ft.Text(value = info, no_wrap=False, max_lines=25)],
                alignment = ft.MainAxisAlignment.CENTER,
                horizontal_alignment = ft.CrossAxisAlignment.CENTER,
                scroll = ft.ScrollMode.ADAPTIVE
                )
    
    def back_func(e):
        controller.change_route(model=dt_model,route="/run_dt_exec", page=page, func=run_dt_exec)
    
    def next_page(e):

        controller.performance_run(dt_model, performance_run, page)

    footer = footer_buttons(back_func=back_func, home_func=home_func, next_page=next_page, page=page)

    return ft.View(route="/run_dt/result_data",
                   controls=[tela, footer])


def performance_run(dt_model, page):
    
    """
    Exibe a performance do veículo com base nos dados do modelo de transmissão.

    Mostra tabelas com forças e velocidades e permite navegação para próximas etapas.

    Args:
        dt_model: Modelo de transmissão instanciado.
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com os dados de performance e navegação.
    """

    rpm = dt_model.rpm
    new_rpm = dt_model.new_rpm

    dados, rpm_faixa = dt_model.CarPerformance()

    forca_trativa, velocidade_angular, velocidade_linear, forca_arrasto, resistencia_rolamento, forca_final = controller.dados(dados)


    linha1 = create_row_title("Força Trativa", "Velocidade Angular", "Velocidade Linear")
    linha2 = create_row_with_columns(forca_trativa, velocidade_angular, velocidade_linear)
    linha3 = create_row_title("Força de Arrasto", "Resistência ao Rolamento", "Força Final")
    linha4 = create_row_with_columns(forca_arrasto, resistencia_rolamento, forca_final)

    
    dt_model.new_rpm = new_rpm
    dt_model.rpm = rpm

    def info(e):
        page.go("/run_dt/result_data")

    def next_page(e):

        controller.slip_ratio(dt_model, slip_ratio, dados, rpm_faixa, page)



    footer = footer_buttons(page=page, back_func=info, home_func=home_func, next_page=next_page)
    
    layout = ft.Column(
        controls = [linha1, linha2, linha3, linha4, footer],
        expand=True,
        spacing=5
    )


    return ft.View(
        route="/run_dt/performance_run",
        controls=[ft.Container(
            content=layout,
            expand=True,
            padding=10,
        )]
    )

def slip_ratio(dt_model, dados, rpm_faixa, page:ft.Page):

    """
    Exibe a tela para entrada da velocidade do veículo e atualização do cálculo do slip ratio.

    Permite atualizar os dados e navegar para a próxima tela.

    Args:
        dt_model: Modelo de transmissão instanciado.
        dados: Dados necessários para o cálculo do slip ratio.
        rpm_faixa: Faixa de RPM para análise.
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com controles para velocidade, atualização e navegação.
    """

    velocidade_veiculo = caixa_texto("Velocidade do veiculo", view=page, variable=mensagem)

    def voltar(e):
        page.go("/run_dt/performance_run")

    def slip(e):
        
        controles = slip_controls(dados,int(velocidade_veiculo.value),dt_model.tire_radius, rpm_faixa)

        tela.controls = [velocidade_veiculo, atualizar_button, controles]
        
        page.update()

    def next_page(e):

        controller.halfshafts(dt_model,halfshafts,page)
    
    atualizar_button = ft.IconButton(icon=ft.Icons.NEXT_PLAN_ROUNDED, on_click=slip)
    
    tela = ft.Column(
        controls=[
            velocidade_veiculo,atualizar_button
        ],
        expand=4)

    
    footer = ft.Container(
        content=footer_buttons(page=page, back_func=voltar, home_func=home_func, next_page=next_page)
    )
    view = ft.View(
        route="/run_dt/slip_ratio",
        controls=[
            tela,
            footer
        ]
    )

    return view

def slip_controls(dados,velocidade_linear,raio_pneu, rpm_faixa):

    """
    Gera a interface com gráficos do slip ratio baseados nos dados fornecidos.

    Apresenta imagens com rolagem vertical.

    Args:
        dados: Dados utilizados para cálculo do slip ratio.
        velocidade_linear (float): Velocidade linear do veículo.
        raio_pneu (float): Raio do pneu.
        rpm_faixa: Faixa de RPM para análise.

    Returns:
        ft.Container: Container com imagens dos gráficos.
    """

    img_64 = controller.calcular_slip(dados, velocidade_linear, raio_pneu, rpm_faixa)

    # Coluna para exibir imagens com rolagem vertical
    images = ft.Column(
        expand=True,
        scroll=ft.ScrollMode.ALWAYS,
        alignment=ft.MainAxisAlignment.START,
        horizontal_alignment=ft.CrossAxisAlignment.CENTER  # Centraliza horizontalmente
    )

    # Container que ocupa toda a largura e altura disponível
    tela = ft.Container(
        content=images,
        expand=True
    )

    # Adiciona cada imagem, ocupando toda a largura disponível
    for img_data in img_64:
        images.controls.append(
            ft.Container(
                content=ft.Image(
                    src_base64=img_data,
                    fit=ft.ImageFit.FILL,
                    repeat=ft.ImageRepeat.NO_REPEAT,
                    border_radius=ft.border_radius.all(10),
                    expand=True
                ),
                width=float("inf"),   # ocupa toda a largura possível
                height=450,
                alignment=ft.alignment.center
            )
        )

    return tela


def halfshafts(dt_model, page:ft.Page) -> ft.View:

    """
    Exibe os resultados do dimensionamento do semi-eixo (Half Shafts).

    Inclui texto com os resultados e botões para navegação.

    Args:
        dt_model: Modelo de transmissão instanciado.
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com resultados do dimensionamento e navegação.
    """

    text = dt_model.resultados_HalfShaftsSizing()

    def voltar(e):
        page.go("/run_dt/slip_ratio")

    footer = ft.Row(
        controls= [
            footer_buttons(page=page, back_func=voltar, home_func=home_func, next_visibility=False),
        ],
        expand=1
    )
    tela  = ft.Container(
        content=ft.Column(
            controls=[
                ft.Text(value=text, no_wrap=False, selectable=True),
                footer
            ],
            alignment=ft.MainAxisAlignment.CENTER,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
            scroll=ft.ScrollMode.AUTO,
            expand=4,
        ),
        alignment=ft.alignment.center,
        padding=20,
        expand=True,
    )


    view = ft.View(
        route="/run_dt/halfshafts",
        controls=[
            tela
        ]
    )

    return view



