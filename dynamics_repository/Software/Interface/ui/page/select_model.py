import flet as ft
from Software.Interface.controllers.MainController import MainController
from routes import route_map
from Software.Interface.ui.page.buttons import back_button, home_button, footer_buttons

controller = MainController()

def home_func(e, page):
    """
    Função para navegar à página inicial ("/").

    Args:
        e: Evento de clique (passado automaticamente).
        page (ft.Page): Página atual do Flet.
    """
    page.go("/")

def head(system, route, func, page):

    """
    Cria a interface para seleção e execução de um modelo salvo para um sistema específico.

    Exibe o número de modelos salvos, lista-os em um dropdown e ao selecionar,
    executa a função associada no controlador principal, além de prover botões de navegação.

    Args:
        system (str): Nome do sistema (ex: "Freio", "Transmissão").
        route (str): Rota para a qual a aplicação deve navegar ao executar o modelo.
        func (callable): Função a ser executada com o modelo selecionado.
        page (ft.Page): Página atual do Flet, usada para navegação.

    Returns:
        ft.Column: Coluna contendo os componentes da interface.
    """
    
    def back_func(e):
        page.go("/run")

    def none_func(e):
        return

    size, lista = controller.models_saved(system)
    
    op = [ft.dropdown.Option(n) for n in lista]
    
    def select_model(e, system=system):
        
        num = int(e.control.value)

        if size:
            
            controller.run(system, num, route=route, func=func, page=page)
            
            
            
    text = ft.Column(controls= [ft.Text(f"Número de Modelos Salvos\n{size}")],
            alignment=ft.MainAxisAlignment.CENTER,
            horizontal_alignment = ft.CrossAxisAlignment.CENTER, expand=1)
    
    
    dropdown = ft.Column(
        controls = [ft.Dropdown(label="Selecione o modelo", options=op, on_change=select_model)],
        alignment=ft.MainAxisAlignment.CENTER,
        horizontal_alignment = ft.CrossAxisAlignment.CENTER,
        expand= 1)                
    
    
    top = ft.Row(controls=[text, dropdown], expand=1, alignment =  ft.alignment.center)
    footer = footer_buttons(page=page, back_func=back_func, home_func=home_func, next_page=none_func, next_visibility=False)
    

    tela = ft.Column(controls = [top, footer], expand = True)
    
    return tela