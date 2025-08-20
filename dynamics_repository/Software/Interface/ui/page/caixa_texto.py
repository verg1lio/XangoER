import flet as ft

def caixa_texto(hint_text, view, variable, proximo=None, usar_none = False):
    # Campo de texto
    texto = ft.TextField(
        hint_text=hint_text,
    )

    # Checkbox para forçar None
    use_none = ft.Checkbox(label="Usar None", value=False, visible=usar_none)

    # Mantém sua lógica: vazio -> None
    def on_change(e):
        if e.control.value == "":
            e.control.value = None
        view.update()
    texto.on_change = on_change

    # Enviar com Enter
    def enviar(e):
        variable.value = None if use_none.value else texto.value
        variable.text = hint_text
        variable.atualizar_dict()

        if proximo:
            try:
                proximo.focus()
            except Exception:
                pass
        view.update()
    texto.on_submit = enviar

    # Alternar None
    def toggle_none(e):
        if use_none.value:
            texto.value = None
            texto.disabled = True
            # 🔹 força atualização do dicionário imediatamente
            variable.value = None
        else:
            texto.disabled = False
            texto.focus()
            # 🔹 restaura o valor atual do campo no dicionário
            variable.value = texto.value

        variable.text = hint_text
        variable.atualizar_dict()
        view.update()

    use_none.on_change = toggle_none

    # Empacota controles
    container = ft.Column(controls=[texto, use_none])

    # Encaminha o focus() da Column para o TextField
    def focus_forward():
        texto.focus()
    container.focus = focus_forward

    return container
