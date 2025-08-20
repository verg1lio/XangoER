import flet as ft

def create_scrollable_column(label):
        return ft.Column(
            controls=[ft.Text(label, text_align=ft.TextAlign.CENTER)],
            scroll=ft.ScrollMode.AUTO,
            alignment=ft.MainAxisAlignment.CENTER,
            expand=True
        )

def create_column_title(label):
    return ft.Container(
        content=ft.Text(label, weight="bold"),
        alignment=ft.alignment.center,
        padding=5
    )

def create_row_title(c1,c2,c3):
    return ft.Row(
        controls=[
            ft.Container(create_column_title(f"{c1}"), expand=True),
            ft.Container(create_column_title(f"{c2}"), expand=True),
            ft.Container(create_column_title(f"{c3}"), expand=True),
        ],
        height=30 # altura fina
    )

def create_row_with_columns(c1,c2,c3):
    return ft.Row(
        controls=[
            ft.Container(create_scrollable_column(f"{c1}"), expand=True, padding=5),
            ft.Container(create_scrollable_column(f"{c2}"), expand=True, padding=5),
            ft.Container(create_scrollable_column(f"{c3}"), expand=True, padding=5),
        ],
        expand=True
    )


