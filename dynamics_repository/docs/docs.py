import flet as ft

def readme(page:ft.Page):


    conteudo_md = read_doc("README.md")

    

    def route(e):
        page.go("/")

    def tap_link(e):

        if e.data == "https://github.com/verg1lio/XangoER":
            page.launch_url(e.data)
        
        else:
            print(e.data[5:-3])
            nome = e.data
            text_markdown.value = read_doc(str(e.data))
            page.update()

    text_markdown = ft.Markdown(
                conteudo_md,
                on_tap_link=tap_link)

    view = ft.View(
        route="/info",
        controls=[ft.Column(
            controls=[
                text_markdown
            ],
            alignment=ft.MainAxisAlignment.CENTER,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
            scroll=ft.ScrollMode.AUTO,
            expand=True,
        )
        ], appbar=ft.AppBar(
            leading=ft.Icon(ft.Icons.ELECTRIC_CAR),
            leading_width=10,
            title=ft.Text("Project Info"),
            center_title=True,
            actions=[
                ft.IconButton(
                    ft.Icons.ARROW_BACK,
                    on_click=route)
            ]
        ))

    return view


def read_doc(nome="README.md"):
    with open(nome, "r", encoding="utf-8") as f:
        conteudo_md = f.read()

    return conteudo_md