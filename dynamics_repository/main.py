import flet as ft
from routes import route_map
from Software.Interface.ui.view import system, data_base
from docs.docs import readme
from Software.Interface.controllers import MainController

controller = MainController()

def main_page(page: ft.Page):

    """
    Cria a view principal da aplicação com uma barra de navegação e appbar.

    Exibe mensagem de boas-vindas e oferece navegação entre as páginas "New",
    "Run" e "Show" através do NavigationBar.

    Args:
        page (ft.Page): Página atual do Flet para configuração e navegação.

    Returns:
        ft.View: View principal configurada com AppBar e NavigationBar.
    """

    def check_item_clicked(e):
        """
        Evento para navegar à página de informações (/info) ao clicar no ícone.

        Args:
            e: Evento de clique.
        """
        page.go("/info")

    def page_view(e):
        """
        Controla a navegação quando o usuário seleciona uma aba no NavigationBar.

        Args:
            e: Evento de mudança no NavigationBar.
        """

        page_selected = e.control.selected_index

        if page_selected == 0:
            page.go("/create")
        elif page_selected == 1:
            page.go("/run")
        elif page_selected == 2:
            page.go("/data_base")

    return ft.View(
        route="/",
        controls=[
            ft.Container(
                content=ft.Text(
                    "Bem-vindo(a) ao sistema de dinâmica!",
                    size=24
                ),
                alignment=ft.alignment.center,
                padding=20,
                expand=True
            )
        ], appbar=ft.AppBar(
        leading=ft.Icon(ft.Icons.ELECTRIC_CAR_ROUNDED),
        leading_width=10,
        title=ft.Text("Xangô eRacing"),
        center_title=True,
        bgcolor=ft.Colors.BLUE_GREY,
        actions=[
            ft.IconButton(ft.Icons.INFO_OUTLINE, on_click=check_item_clicked)
        ]
    ), navigation_bar=ft.NavigationBar(
                destinations=[
                    ft.NavigationBarDestination(icon=ft.Icons.ADD_CIRCLE_OUTLINE_OUTLINED, label="New"),
                    ft.NavigationBarDestination(icon=ft.Icons.PLAY_CIRCLE_OUTLINE_OUTLINED, label="Run"),
                    ft.NavigationBarDestination(icon=ft.Icons.TABLE_CHART_ROUNDED, label="Show")
                ],
                on_change=page_view
            )
    )


def main(page: ft.Page):

    """
    Função principal de inicialização da aplicação.

    Configura a página principal, tema, tamanho da janela, e gerencia
    o roteamento das views com base na URL atual.

    Args:
        page (ft.Page): Página principal do Flet.
    """

    page.title = "Dynamics App"
    page.window_width = 800
    page.window_height = 500
    page.theme_mode = ft.ThemeMode.DARK

    def route_change(e):
        page.views.clear()
        
        route_map.update({"/": lambda: main_page(page)})
        route_map.update({"/create": lambda: system(page, 0)})
        route_map.update({"/run": lambda: system(page, 1)})
        route_map.update({"/data_base": lambda: data_base(page)})
        route_map.update({"/info": lambda: readme(page)})
        view_fn = route_map.get(page.route)

        if view_fn:
            page.views.append(view_fn())
        else:
            
            page.views.append(ft.View(route=page.route, controls=[ft.Text("Página não encontrada")]))


        page.update()

    page.on_route_change = route_change
    page.go("/")
    
ft.app(target=main)

