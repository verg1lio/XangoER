import flet as ft
import time
import asyncio
from Software.Interface.ui.page.caixa_texto import caixa_texto
from Software.Interface.controllers import TireController
from Software.Interface.ui.page.buttons import save_button, back_button, home_button, footer_buttons
from Software.Interface.controllers.envio_texto import envio
from Software.Interface.ui.page.select_model  import head

controller = TireController()
mensagem = envio()

def home_func(e, page):
    page.go("/")

def create_tire(page: ft.Page) -> ft.View:

    status_text = ft.Text("")  # Mensagem de feedback
    status_container = ft.Container(content=status_text, alignment=ft.alignment.center)

    # Campos de entrada
    tire_friction_coef = caixa_texto("Coeficiente de fricção entre pneu e pista", variable=mensagem, view=page)
    tire_k = caixa_texto("Rigidez vertical do pneu", variable=mensagem, view=page, proximo=tire_friction_coef)
    track_y = caixa_texto("Deformação vertical do pneu", proximo=tire_k, variable=mensagem, view=page)
    rear_axle_length = caixa_texto("Comprimento do eixo fixo traseiro", proximo=track_y, variable=mensagem, view=page)
    WB = caixa_texto("Entre-eixos fixo", proximo=rear_axle_length, variable=mensagem, view=page)
    B3 = caixa_texto("Comprimento do segundo braço de direção", proximo=WB, variable=mensagem, view=page)
    B2 = caixa_texto("Comprimento da bitola", proximo=B3, variable=mensagem, view=page)
    B1 = caixa_texto("Comprimento do braço de direção", proximo=B2, variable=mensagem, view=page)
    B0 = caixa_texto("Comprimento da barra de direção", proximo=B1, variable=mensagem, view=page)
    tire_Ca = caixa_texto("Ângulo do camber do pneu", proximo=B0, variable=mensagem, view=page)
    tire_Fz = caixa_texto("Carga vertical no pneu", proximo=tire_Ca, variable=mensagem, view=page)

    inputs_column = ft.Column(
        controls=[
            tire_Fz, tire_Ca, B0, B1, B2, B3, WB, rear_axle_length, track_y, tire_k, tire_friction_coef
        ],
        height=page.height - 180,
        scroll=ft.ScrollMode.ADAPTIVE
    )

    # Feedbacks de sucesso/erro
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

    def instance_tire(e):
        print(mensagem.dicionario)
        controller.instance_tire(
            dicionario=mensagem.dicionario,
            instanciado=modelo_criado,
            falha=erro_criacao
        )

    def back_create(e):
        page.go("/create")
    
    def home_func(e):
        page.go("/")

    save = save_button(func=instance_tire)
    back = back_button(func=back_create)
    home = home_button(home_func)



    # View principal
    view = ft.View(
        route="/create_tire",
        controls=[
            ft.Text("Cadastro de Pneu", size=24, weight=ft.FontWeight.BOLD),
            inputs_column,
            ft.Row(controls=[back, home, save], alignment=ft.MainAxisAlignment.END)
        ],
        vertical_alignment=ft.MainAxisAlignment.START
    )

    return view


def run_tire(page:ft.Page) -> ft.View:
    
    """
    View para seleção do modelo de pneu a ser executado.

    Args:
        page (ft.Page): Página atual do Flet.

    Returns:
        ft.View: View com interface para escolha do modelo Brake a executar.
    """

    view = ft.View(
        route = "/run_tire",
        controls=[head("Tire", "run_tire/camber", run_camber, page)]
    )

    return view



def run_camber(model, page:ft.Page) -> ft.View:
    """
    View do Flet para exibir o GIF animado do mecanismo (desktop app).
    O arquivo deve estar dentro da pasta 'assets'.


    """

    #recebendo os dados gráficos do modelo instanciado
    img_64 = controller.run_data(model)

    #Criando uma função para voltar aos dados anteriores
    def voltar(e):
        page.go("/run_tire")

    def next_page(e):
        controller.run_deformacao(img_64, run_deformacao, page)

    #Cria uma row com paara receber os gráficos
    images = ft.Row(expand=2,
                    wrap=False,
                    scroll=ft.ScrollMode.ALWAYS,
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.START)
    

    tela = ft.Column(controls=[ft.Container(content=images, expand=True)])
    
    footer = footer_buttons(page=page, back_func=voltar, next_page=next_page, home_func=home_func)

    for img_data in img_64["camber"]:
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
    return ft.View(route="/run_tire/camber",
                controls = [
                    tela,
                    footer
                ])




def run_deformacao(graph, page: ft.Page) -> ft.View:

    # Criando uma função para voltar aos dados anteriores
    def voltar(e):
        page.go("/run_tire/camber")

    def next_page(e):
        controller.run_graph(graph, run_graph, page)

        # controller.run_deformacao(graph, run_deformacao, page)
    
    
    # Coluna que centraliza o conteúdo
    tela = ft.Container(
                    content=ft.Image(src_base64=graph["deformacao"]),
                    alignment=ft.alignment.center,
                )

    footer = footer_buttons(page=page, back_func=voltar, next_page=next_page, home_func=home_func)

    return ft.View(
        route="/run_tire/deformacao",
        controls=[
            tela,
            footer
        ]
    )




def run_graph(img_64, page:ft.Page) -> ft.View:

    #Criando uma função para voltar aos dados anteriores
    def voltar(e):
        page.go("/run_tire/deformacao")

    def next_page(e):
        
        controller.run_mechanism(img_64, run_mechanism, page)

    #Cria uma row com paara receber os gráficos
    images = ft.Row(expand=2,
                    wrap=False,
                    scroll=ft.ScrollMode.ALWAYS,
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.START)
    

    tela = ft.Column(controls=[ft.Container(content=images, expand=True)])
    
    footer = footer_buttons(page=page, back_func=voltar, next_page=next_page, home_func=home_func)

    for img_data in img_64["graph"]:
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
    return ft.View(route="run_tire/graph",
                controls = [
                    tela,
                    footer
                ])




def run_mechanism(img_64, page:ft.Page) -> ft.View:

    img_64 = img_64["mechanism"]

    frame_time: float = 0.2

    playing = {"value": True}
    
    img = ft.Image(src_base64=img_64, expand=True)

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
        page.go("/run_tire/graph")

   
    tela = ft.Container(content=img, alignment=ft.alignment.center, expand=True)
    animation_buttons = ft.Row(
        alignment=ft.MainAxisAlignment.CENTER,
        controls=[
            ft.ElevatedButton("Play", on_click=play_click),
            ft.ElevatedButton("Pause", on_click=pause_click),
        ],
    )
    footer = footer_buttons(page, back_func=voltar, next_visibility=False, home_func=home_func)
    
    return ft.View(
        route="/run_tire/run_mechanism",
        controls=[tela, animation_buttons, footer],
        vertical_alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
        horizontal_alignment=ft.CrossAxisAlignment.CENTER,
    )
