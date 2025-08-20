import flet as ft
from Software.Interface.ui.page.buttons import footer_buttons
from Software.Interface.ui.view import create_brake, create_drivetrain, create_kinematics, create_tire, run_brake, run_dt, run_tire, run_kmt
from routes import route_map



def home_func(e, page):
    page.go("/")


def none_func(e):
    return

def system(page: ft.Page, option: int) -> ft.View:

    
    route_map.update({"/create_brake": lambda: create_brake(page)})
    route_map.update({"/create_drivetrain": lambda: create_drivetrain(page)})
    route_map.update({"/create_kinematics": lambda: create_kinematics(page)})
    route_map.update({"/create_tire": lambda: create_tire(page)})
    route_map.update({"/run_brake": lambda: run_brake(page)})
    route_map.update({"/run_drivetrain": lambda: run_dt(page)})
    route_map.update({"/run_kinematics": lambda :run_kmt(page)})
    route_map.update({"/run_tire": lambda :run_tire(page)})

    def on_change(e):
        system_selected = e.control.value

        if option == 1:
            if system_selected == "Freio":
                page.go("/run_brake")
            elif system_selected == "Transmissão":
                page.go("/run_drivetrain")
            elif system_selected == "Cinemática":
                page.go("/run_kinematics")
            elif system_selected == "Pneu":
                page.go("/run_tire")

        if option == 0:
            if system_selected == "Freio":
                page.go("/create_brake")
            elif system_selected == "Transmissão":
                page.go("/create_drivetrain")
            elif system_selected == "Cinemática":
                page.go("/create_kinematics")
            elif system_selected == "Pneu":
                page.go("/create_tire")
        
        page.update()

    dropdown = ft.Dropdown(
        label="Selecione um subsistema",
        options=[
            ft.dropdown.Option("Cinemática"),
            ft.dropdown.Option("Freio"),
            ft.dropdown.Option("Pneu"),
            ft.dropdown.Option("Transmissão")
        ],
        on_change=on_change
    )

    footer = footer_buttons(back_func=none_func, next_page=none_func, home_func=home_func, back_visibility=False, next_visibility=False, page=page)

    return ft.View(
        route="/create",
        controls=[
            ft.Container(content=dropdown, alignment=ft.alignment.center, padding=20),
            ft.Row(controls=[footer], alignment=ft.MainAxisAlignment.CENTER, expand=4)
        ]
    )

def view_models():
    pass