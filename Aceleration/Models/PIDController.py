# -*- coding: utf-8 -*-
"""
PIDController.py — Controlador PI(D) genérico com anti-windup
==============================================================

Usado pela Simulation como controlador de corrente dos eixos d e q do motor
(dois instâncias independentes, ganhos calculados por cancelamento de polo).
O termo derivativo existe mas fica desligado (kd=0) no projeto.

O ponto delicado deste controlador é o ANTI-WINDUP em dois estágios:

  Estágio (a) — interno, automático:
      A saída bruta kp·e + ki·∫e + kd·de é saturada em ±limit. Quando isso
      acontece, o excesso cortado é devolvido ao integrador (back-calculation)
      dentro do próprio ``update()``.

  Estágio (b) — externo, via ``back_calculate()``:
      Na Simulation, a saída do PI é SOMADA a um termo de feedforward
      (desacoplamento de back-EMF) e o total é cortado a ±Vdc — um limite
      dinâmico que o PID não conhece. O excesso desse corte externo é
      devolvido ao integrador chamando ``back_calculate(excesso)``.

Sem o estágio (b), o integrador acumularia comando que nunca chega à planta
("windup silencioso") e o controle apresentaria overshoot ao sair da
saturação de tensão.
"""

import numpy as np


class Controller:
    """Controlador PID com saturação de saída e anti-windup por back-calculation.

    Lei de controle (forma paralela):

        u_raw = kp·e + ki·∫e·dt + kd·de/dt
        u     = clip(u_raw, −limit, +limit)

    Quando u_raw ≠ u (saturou), o integrador é corrigido subtraindo
    (u_raw − u)/ki, de modo que o estado interno permaneça consistente com
    a saída realmente aplicada à planta.

    Parameters
    ----------
    kp : float
        Ganho proporcional.
    ki : float
        Ganho integral. Se 0, o anti-windup é ignorado (não há integrador
        a proteger).
    kd : float
        Ganho derivativo (0 no uso atual do projeto — derivada de erro de
        corrente amplificaria ruído numérico).
    limit : float
        Limite absoluto de saturação da saída (a Simulation usa a tensão
        nominal do pack como limite).
    Ts : float, optional
        Período de amostragem padrão [s], usado quando ``update()`` é
        chamado sem ``dt`` explícito. Default 0.01.

    Attributes
    ----------
    integral : float
        Acumulador do termo integral ∫e·dt (já corrigido pelo anti-windup).
    prev_error : float
        Erro do passo anterior — base do termo derivativo por diferença finita.
    prev_output : float
        Última saída aplicada (diagnóstico).
    anti_windup_enabled : bool
        Chave geral do anti-windup (True por default; desligar só para
        estudos de comparação com/sem windup).
    """

    def __init__(self, kp, ki, kd, limit, Ts=0.01):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.limit = limit
        self.Ts = Ts

        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_output = 0.0
        self.anti_windup_enabled = True

    def update(self, error, dt=None):
        """Calcula o sinal de controle para o erro atual.

        Sequência: integra provisoriamente → calcula derivada → soma os três
        termos → satura em ±limit → devolve ao integrador o excesso cortado
        (back-calculation). Integrar ANTES de saturar e corrigir depois é
        equivalente a integração condicional, mas converge de forma limpa e
        simétrica quando a saturação é liberada.

        Parameters
        ----------
        error : float
            Erro de controle (referência − medida).
        dt : float, optional
            Passo de tempo [s]. Se omitido, usa ``self.Ts``.

        Returns
        -------
        float
            Saída de controle saturada em [−limit, +limit].
        """
        if dt is None:
            dt = self.Ts

        # Integra provisoriamente; o back-calculation abaixo reverte
        # qualquer excesso que a saturação cortar.
        self.integral += error * dt

        # Derivada por diferença finita — só avaliada se kd ≠ 0 (evita
        # divisão desnecessária no caminho quente do loop de 10 kHz).
        derivative = 0.0
        if dt > 0 and self.kd != 0.0:
            derivative = self.kd * (error - self.prev_error) / dt

        output_raw = self.kp * error + self.ki * self.integral + derivative

        # Saturação simétrica em ±limit
        if output_raw > self.limit:
            output = self.limit
        elif output_raw < -self.limit:
            output = -self.limit
        else:
            output = output_raw

        # Anti-windup por back-calculation: quando a saída foi cortada,
        # subtrai do integrador exatamente o excesso (em unidades de
        # integral: excesso/ki). O próximo passo então parte de um estado
        # consistente com a saída REALMENTE aplicada — isso evita o
        # ciclo-limite assimétrico do anti-windup por congelamento
        # (clamping) e converge de forma limpa ao sair da saturação.
        if self.anti_windup_enabled and self.ki > 0.0 and output_raw != output:
            self.integral -= (output_raw - output) / self.ki

        self.prev_error = error
        self.prev_output = output

        return output

    def back_calculate(self, excess):
        """Gancho de anti-windup EXTERNO: corrige o integrador em ``excess/ki``.

        Usado quando a saída do controlador sofre uma saturação adicional
        FORA do PID — no projeto, após somar o feedforward de desacoplamento
        (±ωe·L·i) e cortar o total a ±Vdc, limite dinâmico que varia com a
        carga da bateria e que o PID não conhece.

        Parameters
        ----------
        excess : float
            Défice pós-corte: ``saida_bruta − saida_cortada``. Excesso
            positivo ENCOLHE o integrador (puxa a saída para baixo);
            negativo o expande. Zero é ignorado sem custo.
        """
        if self.anti_windup_enabled and self.ki > 0.0 and excess != 0.0:
            self.integral -= excess / self.ki

    def reset(self):
        """Zera os estados internos (integral, erro e saída anteriores).

        Chamado pela Simulation no início de cada ``simulate()`` para que
        execuções sucessivas partam do mesmo estado.
        """
        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_output = 0.0

    def set_parameters(self, kp, ki, kd):
        """Atualiza os ganhos do controlador em tempo de execução."""
        self.kp = kp
        self.ki = ki
        self.kd = kd
