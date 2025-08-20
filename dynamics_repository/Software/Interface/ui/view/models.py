import flet as ft
from Software.Interface.controllers.MainController import MainController

def data_base(page: ft.Page):

    controller = MainController()
    df_atual_nome = {"nome": None}
    tabela_area = ft.Column(expand=True, scroll=ft.ScrollMode.AUTO)

    # Estado de edição
    linha_editando = {"index": None, "campos": {}}

    def home_func(e):
        page.go("/")

    def salvar_edicao(linha_idx):
        # Salva os valores no controller
        for col, tf in linha_editando["campos"].items():
            valor = tf.value
            controller.atualizar_celula(df_atual_nome["nome"], linha_idx, col, valor)

        # Sai do modo de edição
        linha_editando["index"] = None
        linha_editando["campos"] = {}
        atualizar_tabela(df_atual_nome["nome"])

    def cancelar_edicao():
        linha_editando["index"] = None
        linha_editando["campos"] = {}
        atualizar_tabela(df_atual_nome["nome"])

    def editar_linha(linha_idx):
        # Entra no modo de edição
        linha_editando["index"] = linha_idx
        atualizar_tabela(df_atual_nome["nome"])

    def excluir_linha(linha_idx):
        controller.excluir_linha(df_atual_nome["nome"], linha_idx)
        atualizar_tabela(df_atual_nome["nome"])

    def atualizar_tabela(nome_df):
        df_atual_nome["nome"] = nome_df
        df = controller.get_dataframe(nome_df)
        tabela_area.controls.clear()

        if df.empty:
            tabela_area.controls.append(
                ft.Container(
                    ft.Text("Nenhum dado disponível para exibição.", italic=True),
                    alignment=ft.alignment.center,
                    expand=True
                )
            )
        else:
            tabela = ft.DataTable(
                columns=[ft.DataColumn(ft.Text(c)) for c in df.columns] + [ft.DataColumn(ft.Text("Ações"))],
                rows=[]
            )

            for idx, row in df.iterrows():
                cells = []

                if linha_editando["index"] == idx:
                    # Linha em edição
                    for col in df.columns:
                        # Cria apenas se ainda não existe (para manter os valores digitados)
                        if col not in linha_editando["campos"]:
                            tf = ft.TextField(
                                value=str(row[col]),
                                width=180,
                                dense=True,
                                text_align=ft.TextAlign.CENTER
                            )
                            linha_editando["campos"][col] = tf
                        cells.append(ft.DataCell(linha_editando["campos"][col]))

                    # Botões de salvar e cancelar
                    btn_salvar = ft.IconButton(
                        ft.Icons.SAVE,
                        icon_color="green",
                        tooltip="Salvar alterações",
                        on_click=lambda e, linha=idx: salvar_edicao(linha)
                    )
                    btn_cancelar = ft.IconButton(
                        ft.Icons.CANCEL,
                        icon_color="red",
                        tooltip="Cancelar edição",
                        on_click=lambda e: cancelar_edicao()
                    )
                    cells.append(ft.DataCell(ft.Row([btn_salvar, btn_cancelar], spacing=5)))
                    tabela.rows.append(ft.DataRow(cells=cells, color="lightyellow"))

                else:
                    # Linha normal
                    for col in df.columns:
                        cells.append(ft.DataCell(ft.Text(str(row[col]))))

                    btn_editar = ft.IconButton(
                        ft.Icons.EDIT,
                        icon_color="blue",
                        tooltip="Editar linha",
                        on_click=lambda e, linha=idx: editar_linha(linha)
                    )
                    btn_del = ft.IconButton(
                        ft.Icons.DELETE,
                        icon_color="red",
                        tooltip="Excluir linha",
                        on_click=lambda e, linha=idx: excluir_linha(linha)
                    )
                    cells.append(ft.DataCell(ft.Row([btn_editar, btn_del], spacing=5)))
                    tabela.rows.append(ft.DataRow(cells=cells))

            tabela_area.controls.append(
                ft.Container(
                    content=ft.Row([tabela], scroll=ft.ScrollMode.AUTO),
                    expand=True
                )
            )

        page.update()

    # Botões para trocar o dataframe exibido
    botoes = ft.Row(
        [
            ft.ElevatedButton(nome, on_click=lambda e, n=nome: atualizar_tabela(n))
            for nome in controller.get_all_dataframe_names()
        ],
        spacing=10
    )

    footer = ft.Row(
        controls=[
            ft.Container(
                ft.IconButton(
                    ft.Icons.NAVIGATE_BEFORE,
                    on_click=home_func
                )
            )
        ]
    )

    # Inicializa a tabela com o primeiro dataframe disponível
    nomes_dfs = controller.get_all_dataframe_names()
    if nomes_dfs:
        atualizar_tabela(nomes_dfs[0])

    return ft.View(
        route="/data_base",
        controls=[
            botoes,
            tabela_area,
            footer
        ]
    )
