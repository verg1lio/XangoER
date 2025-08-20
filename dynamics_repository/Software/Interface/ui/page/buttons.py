import flet as ft

def save_button(func, opt=True):
    """
    Cria um botão flutuante de salvar (ícone ADD).

    Args:
        func (callable): Função a ser chamada no clique.
        opt (bool, opcional): Define se o botão está visível. Padrão é True.

    Returns:
        ft.Container: Container contendo o FloatingActionButton configurado.
    """
    save_button = ft.Container(
        content=ft.FloatingActionButton(
            icon=ft.Icons.ADD,
            on_click=func,
            visible=opt
        ),
        alignment=ft.alignment.bottom_right,
        padding=10
    )
    return save_button

def next_button(func, opt=True):
    """
    Cria um botão flutuante para avançar (ícone NAVIGATE_NEXT).

    Args:
        func (callable): Função a ser chamada no clique.
        opt (bool, opcional): Define se o botão está visível. Padrão é True.

    Returns:
        ft.Container: Container contendo o FloatingActionButton configurado.
    """
    next_button = ft.Container(
        content=ft.FloatingActionButton(
            icon=ft.Icons.NAVIGATE_NEXT,
            on_click=func,
            visible=opt
        ),
        alignment=ft.alignment.bottom_right,
        padding=10
    )
    return next_button

def back_button(func, opt=True):
    """
    Cria um botão flutuante para voltar (ícone NAVIGATE_BEFORE).

    Args:
        func (callable): Função a ser chamada no clique.
        opt (bool, opcional): Define se o botão está visível. Padrão é True.

    Returns:
        ft.Container: Container contendo o FloatingActionButton configurado.
    """
    back_button = ft.Container(
        content=ft.FloatingActionButton(
            icon=ft.Icons.NAVIGATE_BEFORE,
            on_click=func,
            visible=opt
        ),
        alignment=ft.alignment.bottom_left,
        padding=10
    )
    return back_button

def home_button(func, opt=True):
    """
    Cria um botão flutuante para ir à página inicial (ícone HOME).

    Args:
        func (callable): Função a ser chamada no clique.
        opt (bool, opcional): Define se o botão está visível. Padrão é True.

    Returns:
        ft.Container: Container contendo o FloatingActionButton configurado.
    """
    back_button = ft.Container(
        content=ft.FloatingActionButton(
            icon=ft.Icons.HOME,
            on_click=func,
            visible=opt
        ),
        alignment=ft.alignment.bottom_left,
        padding=10
    )
    return back_button

def func_none(e):
    """
    Função vazia padrão que não realiza nenhuma ação.

    Args:
        e: Evento de clique (recebido automaticamente).

    Returns:
        None
    """
    return

def footer_buttons(page, back_func=func_none, home_func=func_none, next_page=func_none, back_visibility=True, home_visibility=True, next_visibility=True):
    """
    Cria uma linha com os botões de rodapé: voltar, home e avançar, configuráveis.

    Args:
        page (ft.Page): Página atual do Flet.
        back_func (callable, opcional): Função ao clicar no botão voltar.
        home_func (callable, opcional): Função ao clicar no botão home, recebe (evento, page).
        next_page (callable, opcional): Função ao clicar no botão avançar.
        back_visibility (bool, opcional): Visibilidade do botão voltar. Padrão True.
        home_visibility (bool, opcional): Visibilidade do botão home. Padrão True.
        next_visibility (bool, opcional): Visibilidade do botão avançar. Padrão True.

    Returns:
        ft.Row: Linha contendo os botões configurados.
    """
    footer_buttons = ft.Row(
        controls=[
            back_button(back_func, back_visibility),
            home_button(lambda e: home_func(e, page), home_visibility),
            next_button(func=next_page, opt=next_visibility)
        ],
        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
        vertical_alignment=ft.CrossAxisAlignment.CENTER,
    )
    return footer_buttons
