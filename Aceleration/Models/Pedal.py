import numpy as np


class Pedal:
    """Modelo do pedal do acelerador como modulador 0–1.

    Suporta dois modos de operação:
      • Posição estática   — set_posicao(v) define um valor fixo
      • Perfil temporal    — set_profile(t, p) define um perfil (t, posição)
                             interpolado linearmente. Antes do primeiro ponto
                             repete-se o valor inicial; depois do último ponto
                             repete-se o valor final (clamping nas bordas).

    Parameters
    ----------
    ganho : float, optional
        Fator de ganho aplicado à posição. Default 1.0.

    Attributes
    ----------
    posicao : float
        Posição corrente do pedal [0.0, 1.0].
    ganho : float
        Ganho aplicado em get_referencia().
    """

    def __init__(self, ganho=1.0):
        self.posicao = 0.0
        self.ganho = ganho
        self._tprofile = None
        self._pprofile = None

    # ─── Modo estático ────────────────────────────────────────────────────
    def set_posicao(self, valor):
        """Define posição fixa, saturada em [0, 1]."""
        self.posicao = float(np.clip(valor, 0.0, 1.0))

    # ─── Modo perfil temporal ────────────────────────────────────────────
    def set_profile(self, tempos, posicoes):
        """Define um perfil temporal de posição do pedal.

        Parameters
        ----------
        tempos : array-like [s]
            Vetor de instantes (deve ser monotonicamente crescente).
        posicoes : array-like
            Posição do pedal em cada instante, em [0, 1] ou [0, 100] (%).
            Valores > 1 são automaticamente normalizados por 100.
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
        # Inicializa posição corrente no primeiro ponto
        self.posicao = float(p[0])

    def update(self, t):
        """Atualiza a posição do pedal a partir do perfil no instante t.

        Após o último ponto do perfil, mantém o valor final (clamp à direita).
        Retorna a posição atualizada.
        """
        if self._tprofile is None:
            return self.posicao
        self.posicao = float(np.interp(
            t, self._tprofile, self._pprofile,
            left=self._pprofile[0], right=self._pprofile[-1]))
        return self.posicao

    def has_profile(self):
        return self._tprofile is not None

    # ─── Perfil suavizado ────────────────────────────────────────────────
    @staticmethod
    def perfil_suave(t_inicio: float, t_fim: float, n: int = 100):
        """Gera arrays (tempos, posicoes) com rampa suavizada via smoothstep cúbico.

        Usa f(s) = 3s²−2s³, que garante derivada zero em t_inicio e t_fim.
        Isso evita o degrau de taxa que excita oscilações no controlador PI de
        velocidade quando uma rampa linear termina abruptamente.

        Parameters
        ----------
        t_inicio, t_fim : float  [s]
        n : int  número de pontos (100 é mais que suficiente para dt=1e-4 s)

        Returns
        -------
        tempos, posicoes : np.ndarray  (posicoes em [0, 1])
        """
        t = np.linspace(t_inicio, t_fim, n)
        s = (t - t_inicio) / (t_fim - t_inicio)
        posicoes = 3*s**2 - 2*s**3
        return t, posicoes

    # ─── Saída ────────────────────────────────────────────────────────────
    def get_referencia(self):
        """Retorna a referência (posicao × ganho)."""
        return self.posicao * self.ganho