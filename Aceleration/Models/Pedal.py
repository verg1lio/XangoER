# -*- coding: utf-8 -*-
"""
Pedal.py — Interface de comando do acelerador (throttle)
=========================================================

O pedal é a ÚNICA entrada de comando do piloto no digital twin. Ele produz
um sinal normalizado em [0, 1] que a Simulation converte diretamente em
referência de corrente de torque:

    iq_ref = max_current × posição_do_pedal   (controle de torque open-loop)

Não existe malha fechada de velocidade no caminho de controle — o carro
acelera com o torque que o pedal pede, limitado apenas pelos protetores
(limitador suave de rotação, teto de potência DC e derate térmico) dentro
de Simulation.py.

O modelo suporta dois modos:
  1. Posição estática  — valor fixo (ex.: degrau de 100 % desde t=0);
  2. Perfil temporal   — sequência (t, posição) interpolada linearmente,
     usada para modelar rampas de largada (launch control primitivo).
"""

import numpy as np


class Pedal:
    """Modelo do pedal do acelerador como modulador 0–1.

    Modos de operação
    -----------------
    • Posição estática — ``set_posicao(v)`` define um valor fixo que é
      devolvido em todas as chamadas subsequentes de ``update()``.
    • Perfil temporal — ``set_profile(t, p)`` define um perfil (instante,
      posição) interpolado linearmente por ``np.interp``. Antes do primeiro
      ponto repete-se o valor inicial; depois do último ponto repete-se o
      valor final (clamp nas bordas — o pedal "segura" a última posição).

    O gerador ``perfil_suave()`` cria rampas smoothstep (derivada nula nas
    extremidades), que são o perfil recomendado: uma rampa linear termina
    com descontinuidade de taxa e excitava oscilações nos controladores PI.

    Parameters
    ----------
    ganho : float, optional
        Fator multiplicativo aplicado à posição em ``get_referencia()``.
        Default 1.0 (sem escala).

    Attributes
    ----------
    posicao : float
        Posição corrente do pedal, sempre em [0.0, 1.0].
    ganho : float
        Ganho aplicado em get_referencia().
    """

    def __init__(self, ganho=1.0):
        self.posicao = 0.0
        self.ganho = ganho
        # Vetores do perfil temporal (None = modo estático)
        self._tprofile = None
        self._pprofile = None

    # ─── Modo estático ────────────────────────────────────────────────────
    def set_posicao(self, valor):
        """Define posição fixa do pedal, saturada em [0, 1].

        Usada para simular pedal-degrau (ex.: ``set_posicao(1.0)`` =
        acelerador no fundo desde o início da prova).
        """
        self.posicao = float(np.clip(valor, 0.0, 1.0))

    # ─── Modo perfil temporal ────────────────────────────────────────────
    def set_profile(self, tempos, posicoes):
        """Define um perfil temporal de posição do pedal.

        A partir daqui, ``update(t)`` passa a interpolar linearmente o
        perfil no instante pedido, ignorando a posição estática anterior.

        Parameters
        ----------
        tempos : array-like [s]
            Vetor de instantes (deve ser monotonicamente crescente —
            requisito do ``np.interp`` usado em ``update``).
        posicoes : array-like
            Posição do pedal em cada instante, em [0, 1] ou [0, 100] (%).
            Se o máximo do vetor exceder 1, assume-se que os valores estão
            em percentual e todos são divididos por 100 automaticamente.

        Raises
        ------
        ValueError
            Se os vetores tiverem tamanhos diferentes ou forem vazios.
        """
        t = np.asarray(tempos, dtype=float)
        p = np.asarray(posicoes, dtype=float)
        if t.size != p.size:
            raise ValueError("tempos e posicoes devem ter o mesmo tamanho")
        if t.size == 0:
            raise ValueError("perfil vazio")
        # Auto-detecta percentuais (0-100) e normaliza para 0-1
        if np.nanmax(p) > 1.0 + 1e-9:
            p = p / 100.0
        p = np.clip(p, 0.0, 1.0)
        self._tprofile = t
        self._pprofile = p
        # Inicializa posição corrente no primeiro ponto do perfil
        self.posicao = float(p[0])

    def update(self, t):
        """Atualiza e retorna a posição do pedal no instante ``t``.

        Chamada uma vez por passo de controle pela Simulation. Se nenhum
        perfil foi definido, apenas devolve a posição estática corrente.
        Fora do intervalo do perfil, aplica clamp: antes do primeiro ponto
        usa o valor inicial; depois do último ponto mantém o valor final
        (o pedal permanece "cravado" onde o perfil terminou).
        """
        if self._tprofile is None:
            return self.posicao
        self.posicao = float(np.interp(
            t, self._tprofile, self._pprofile,
            left=self._pprofile[0], right=self._pprofile[-1]))
        return self.posicao

    def has_profile(self):
        """True se um perfil temporal foi definido via ``set_profile``."""
        return self._tprofile is not None

    # ─── Perfil suavizado ────────────────────────────────────────────────
    @staticmethod
    def perfil_suave(t_inicio: float, t_fim: float, n: int = 100):
        """Gera arrays (tempos, posicoes) com rampa smoothstep cúbica 0 → 100 %.

        Usa o polinômio de Hermite f(s) = 3s² − 2s³ com s ∈ [0, 1], que
        garante f(0)=0, f(1)=1 e DERIVADA NULA em ambas as extremidades.

        Por que não uma rampa linear?
            Uma rampa linear termina com "quina": a taxa de variação do
            pedal salta de constante para zero instantaneamente. Esse degrau
            de taxa excitava oscilações nos controladores PI (resposta a
            entrada rampa-parada). O smoothstep entra e sai suavemente,
            eliminando o transiente.

        No otimizador, t_fim (= t_pico) funciona como launch control
        primitivo: atrasar o pico do pedal reduz wheelspin na largada.

        Parameters
        ----------
        t_inicio, t_fim : float  [s]
            Início e fim da rampa. Antes de t_inicio o perfil vale 0;
            depois de t_fim, o clamp do ``update`` mantém 100 %.
        n : int
            Número de pontos gerados (100 é mais que suficiente — o
            ``update`` interpola linearmente entre eles, e o passo da
            simulação é 1e-4 s).

        Returns
        -------
        tempos, posicoes : np.ndarray
            Vetores prontos para ``set_profile`` (posicoes em [0, 1]).
        """
        t = np.linspace(t_inicio, t_fim, n)
        s = (t - t_inicio) / (t_fim - t_inicio)
        posicoes = 3*s**2 - 2*s**3
        return t, posicoes

    # ─── Saída ────────────────────────────────────────────────────────────
    def get_referencia(self):
        """Retorna a referência final do pedal (posicao × ganho)."""
        return self.posicao * self.ganho
