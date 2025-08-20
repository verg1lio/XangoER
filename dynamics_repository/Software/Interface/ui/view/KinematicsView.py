import flet as ft
import asyncio
import time
from Software.Interface.ui.page.caixa_texto import caixa_texto
from Software.Interface.controllers import KinematicsController
from Software.Interface.ui.page.buttons import save_button, back_button, home_button, footer_buttons, next_button
from Software.Interface.controllers.envio_texto import envio
from Software.Interface.ui.page.select_model  import head



controller = KinematicsController()
mensagem = envio()

def home_func(e, page):
    page.go("/")

def create_kinematics(page: ft.Page) -> ft.View:

    status_text = ft.Text("")  # Mensagem de feedback
    status_container = ft.Container(content=status_text, alignment=ft.alignment.center)

    # Campos de entrada
    spring_non_lin_coef = caixa_texto("Coeficiente não linear da mola: ", variable=mensagem, view=page, usar_none=True)
    spring_x = caixa_texto("Deslocamento da mola: ", proximo=spring_non_lin_coef, variable=mensagem, view=page)
    spring_k = caixa_texto("Rigidez da mola: ", proximo=spring_x, variable=mensagem, view=page)
    spring_type = caixa_texto("Tipo de mola: ", proximo=spring_k, variable=mensagem, view=page)

    L3 = caixa_texto("Comprimento da barra de saída: ", proximo=spring_type, variable=mensagem, view=page)
    L2 = caixa_texto("Comprimento da barra de acoplamento: ", proximo=L3, variable=mensagem, view=page)
    L1 = caixa_texto("Comprimento da barra de entrada: ", proximo=L2, variable=mensagem, view=page)
    L0 = caixa_texto("Comprimento da barra fixa: ", proximo=L1, variable=mensagem, view=page)

    damper_F_viscous = caixa_texto("Força viscosa do fluído", proximo=L0, variable=mensagem, view=page, usar_none=True)
    damper_K_friction = caixa_texto("Rigidez de fricção do Amortecedor: ", proximo=damper_F_viscous, variable=mensagem, view=page)
    damper_F_static = caixa_texto("Força Viscosa do fluído do Amortecedor: ", proximo=damper_K_friction, variable=mensagem, view=page)
    damper_type = caixa_texto("Tipo de Amortecedor: ", proximo=damper_F_static, variable=mensagem, view=page)

    inputs_column = ft.Column(
        controls=[
            damper_type,
            damper_F_static,
            damper_K_friction,
            damper_F_viscous,
            L0,
            L1,
            L2,
            L3,
            spring_type,
            spring_k,
            spring_x,
            spring_non_lin_coef
        ],
        height=page.height - 180,
        scroll=ft.ScrollMode.ADAPTIVE
    )

    

    # Funções de feedback com atraso
    def modelo_criado():
        status_text.value = "Modelo Salvo!"
        view.controls = [status_container]
        page.update()
        time.sleep(3)
        page.go("/")

    def erro_criacao():
        status_text.value = "Valor(es) Inválido(s)!"
        view.controls = [status_container]
        page.update()
        time.sleep(3)
        page.go("/")

    def instance_kmt(e):
        controller.instance_kmt(
            dicionario=mensagem.dicionario,
            instanciado=modelo_criado,
            falha=erro_criacao
        )

    def back_create(e):
        page.go("/create")
    
    def home_func(e):
        page.go("/")

    save = save_button(func=instance_kmt)
    back = back_button(func=back_create)
    home = home_button(home_func)

    # View será criada antes para podermos atualizar os controls dela
    view = ft.View(
        route="/create_kinematics",
        controls=[
            ft.Text("Cadastro de Cinemática", size=24, weight=ft.FontWeight.BOLD),
            inputs_column,
            ft.Row(controls=[back, home, save], alignment=ft.MainAxisAlignment.END)  # Placeholder para o botão
        ],
        vertical_alignment=ft.MainAxisAlignment.START
    )


    return view

def run_kmt(page:ft.Page) -> ft.View:
    
    view = ft.View(
        route = "/run_kinematics",
        controls=[head("Kinematics", "run_kinematics/cinematica", run_cinematica, page)]
    )

    return view

def run_cinematica(model, page):

    graphics, text = controller.run_kmt(model)  # lista de imagens base64

    img_64 = graphics[2]

    frame_time: float = 0.2
    playing = {"value": True}
    img = ft.Image(src_base64=img_64[0], expand=True)

    # Guarda a referência da task
    animation_task = {"task": None}

    async def animate():
        index = 0
        while True:
            if not playing["value"]:
                await asyncio.sleep(0.1)
                continue

            img.src_base64 = img_64[index]
            try:
                img.update()
            except AssertionError:
                # Controle já removido da página, interrompe animação
                break

            index = (index + 1) % len(img_64)
            await asyncio.sleep(frame_time)

    # Inicia a animação em segundo plano
    animation_task["task"] = page.run_task(animate)

    def play_click(e):
        playing["value"] = True

    def pause_click(e):
        playing["value"] = False

    def voltar(e):
        # Para a animação antes de mudar de página
        playing["value"] = False
        if animation_task["task"]:
            animation_task["task"].cancel()
        page.go("/run_kinematics")

    def next_page(e):

        playing["value"] = False
        if animation_task["task"]:
            animation_task["task"].cancel()

        controller.dados_cinematica(graphics=graphics, text = text, func=dados_cinematica, page=page)

    tela = ft.Container(content=img, alignment=ft.alignment.center, expand=True)
    animation_buttons = ft.Row(
        alignment=ft.MainAxisAlignment.CENTER,
        controls=[
            ft.ElevatedButton("Play", on_click=play_click),
            ft.ElevatedButton("Pause", on_click=pause_click),
        ],
    )
    footer = footer_buttons(page, back_func=voltar, next_page=next_page, home_func=home_func)
    
    return ft.View(
        route="/run_kinematics/cinematica",
        controls=[tela, animation_buttons, footer],
        vertical_alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
        horizontal_alignment=ft.CrossAxisAlignment.CENTER,
    )


def dados_cinematica(graphics, text, page):
    # Criando uma função para voltar aos dados anteriores
    def voltar(e):
        page.go("/run_kinematics/cinematica")

    def next_page(e):
        controller.run_graphs(graphics=graphics, run_plotagem=run_graph, page=page)

    # Texto centralizado
    tela = ft.Container(
        content=ft.Text(value=text, size=18, text_align=ft.TextAlign.CENTER),
        alignment=ft.alignment.center,
        expand=True,  # ocupa todo espaço disponível
    )

    footer = footer_buttons(page=page, back_func=voltar, next_page=next_page, home_func=home_func)

    return ft.View(
        route="/run_kinematics/dados_cinematica",
        controls=[
            ft.Column(
                controls=[tela, footer],
                expand=True,
                alignment=ft.MainAxisAlignment.SPACE_BETWEEN,  # texto no meio, footer embaixo
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
            )
        ]
    )




def run_graph(graphics: list[str], page: ft.Page ) -> ft.View:


    current_index = {"value": 0}  # dicionário para manter referência mutável

    # Imagem exibida
    image_display = ft.Image(src_base64=graphics[current_index["value"]], width=600)

    def update_image():
        image_display.src_base64 = graphics[current_index["value"]]
        image_display.update()

    def voltar(e):
        if current_index["value"] == 0:
            page.go("/run_kinematics/dados_cinematica")  # se for a primeira imagem, volta pra home
        else:
            current_index["value"] -= 1
            update_image()

    def next_page(e):
        if current_index["value"] < len(graphics) - 2:
            current_index["value"] += 1
            update_image()
        else:
            next.visible = False
            page.update()

    back = back_button(voltar)
    home = home_button(lambda e: home_func(e, page))
    next = next_button(func=next_page)
        
    footer_buttons = ft.Row(
        controls=[
            back,
            home,
            next
        ],
        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
        vertical_alignment=ft.CrossAxisAlignment.CENTER,
    )

    # Conteúdo centralizado
    tela = ft.Container(
        content=image_display,
        alignment=ft.alignment.center,
    )

    return ft.View(
        route="/run_kinematics/graph",
        controls=[tela, footer_buttons],
    )