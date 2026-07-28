# -*- coding: utf-8 -*-
"""
BatteryPack.py — Banco de baterias (modelo de Shepherd modificado)
===================================================================

Modela o acumulador do carro como um circuito equivalente: fonte de tensão
dependente do estado de carga (SoC) em série com uma resistência interna.
A equação de Shepherd modificada captura os três efeitos dominantes na
descarga de células de lítio:

  1. Queda ôhmica instantânea (R·i) — proporcional à corrente;
  2. Queda por polarização (K·Q/(Q−q)) — cresce à medida que a célula
     esvazia (q = carga já usada);
  3. Região exponencial (A·e^(−B·q)) — o "degrau" inicial de tensão logo
     no começo da descarga.

Papel na simulação:
  A Simulation calcula a corrente DC (I_bat) por balanço de potência e
  pergunta a este modelo: (a) qual a tensão terminal do pack sob essa
  corrente (→ limite de tensão ±Vdc dos comandos do inversor) e (b) quais
  as derivadas de SoC/Iast para integrar no RK4.

INVARIANTE CRÍTICO (não regredir):
  O modelo de Shepherd é avaliado POR CÉLULA e depois escalado por n_serie.
  A corrente e o Iast que chegam em nível de PACK são divididos por
  n_paralelo antes de entrar na equação — sem isso, o sag de tensão fica
  n_paralelo× superestimado e o barramento DC colapsa (ver docstring de
  ``calcular_tensao``).
"""

import numpy as np
class BatteryPack:
    """Banco de baterias baseado no modelo de Shepherd modificado simplificado.

    Representa um pack composto por células idênticas em arranjo
    série/paralelo (ex.: 144S5P = 144 células em série × 5 fileiras em
    paralelo). Modela a tensão terminal, a dinâmica do estado de carga
    (SoC), a corrente acumulada (efeito de polarização) e propriedades
    físicas agregadas (massa e volume totais).

    Arranjo série/paralelo — regras de escala:
      • Tensão   → multiplica por n_serie   (células em série somam tensão)
      • Corrente → divide por n_paralelo    (fileiras dividem a corrente)
      • Capacidade → multiplica por n_paralelo [Ah]

    Parameters
    ----------
    tipo_celula : str
        Tipo eletroquímico da célula. Opções implementadas no catálogo:
        - 'Li-ion'   (E0=3.7 V, Q=3.0 Ah, R=20 mΩ)
        - 'LiFePO4'  (E0=3.2 V, Q=2.8 Ah, R=8 mΩ)
    n_serie : int
        Número de células em série (define a tensão do pack).
    n_paralelo : int
        Número de fileiras em paralelo (define capacidade e corrente máxima).
    soc_inicial : float, optional
        Estado de carga inicial do pack (0.0–1.0). Default 1.0 (cheio).

    Attributes
    ----------
    E0 : float
        Tensão de circuito aberto nominal de UMA célula [V].
    K : float
        Constante de polarização [V] — escala a queda K·Q/(Q−q).
    Q : float
        Capacidade de UMA célula [Ah].
    A : float
        Amplitude da região exponencial [V].
    B : float
        Constante de decaimento exponencial [Ah⁻¹].
    R : float
        Resistência interna de UMA célula [Ω].
    peso : float
        Massa de uma célula [kg].
    volume : float
        Volume de uma célula [L].
    n_serie, n_paralelo : int
        Configuração do arranjo.
    soc : float
        Estado de carga corrente (0–1).
    carga_total : float
        Capacidade total do pack [Ah] (= Q · n_paralelo).
    inv_carga_total : float
        1/(capacidade em Coulombs) pré-computado — converte corrente [A]
        diretamente em dSoC/dt no loop de 10 kHz.
    Iast : float
        Corrente acumulada ∫i·dt [C] (alimenta o termo de polarização).
    tempo_acumulado : float
        Tempo decorrido [s] (usado no ramo de carga do modelo).
    tensao_hist, corrente_hist, soc_hist, tempo : list
        Históricos (mantidos por compatibilidade; a Simulation faz o
        logging principal em seus próprios arrays).
    """

    def __init__(self, tipo_celula, n_serie, n_paralelo, soc_inicial=1.0):
        self.parametros = self.definir_celula(tipo_celula)
        if self.parametros is None:
            raise ValueError(f"Unknown cell type: {tipo_celula}")

        # Desempacota os parâmetros da célula escolhida em atributos diretos
        params = self.parametros
        self.E0 = params['E0']
        self.K = params['K']
        self.Q = params['Q']
        self.A = params['A']
        self.B = params['B']
        self.R = params['R']
        self.peso = params['peso']
        self.volume = params['volume']

        self.n_serie = n_serie
        self.n_paralelo = n_paralelo

        self.soc = soc_inicial
        # Capacidade total do pack [Ah] e seu inverso em Coulombs
        # (3600 converte Ah → C); pré-computado para o caminho quente.
        self.carga_total = self.Q * n_paralelo  # Ah
        self.inv_carga_total = 1.0 / (self.carga_total * 3600.0)

        # Estados de integração (avançados pelo RK4 da Simulation)
        self.Iast = 0.0  # corrente acumulada ∫i·dt [Coulombs]
        self.tempo_acumulado = 0.0  # tempo decorrido [s]

        # Históricos (compatibilidade)
        self.tensao_hist = []
        self.corrente_hist = []
        self.soc_hist = []
        self.tempo = []

    def definir_celula(self, tipo):
        """Retorna os parâmetros de uma célula eletroquímica pré-definida.

        Catálogo interno com parâmetros de Shepherd típicos por química.
        Para adicionar uma química nova (ex.: a célula real do acumulador
        quando homologada), basta acrescentar uma entrada ao dicionário.

        Parameters
        ----------
        tipo : str
            Identificador da química ('Li-ion' ou 'LiFePO4').

        Returns
        -------
        dict or None
            Dicionário de parâmetros se o tipo existir, senão None
            (o __init__ converte None em ValueError).
        """
        catalogo = {
            'Li-ion': {
                'E0': 3.7,       # tensão nominal de circuito aberto [V]
                'K': 0.005,      # constante de polarização [V]
                'Q': 3.0,        # capacidade [Ah]
                'A': 0.1,        # amplitude exponencial [V]
                'B': 0.01,       # decaimento exponencial [Ah⁻¹]
                'R': 0.02,       # resistência interna [Ω]
                'peso': 0.045,   # massa [kg]
                'volume': 0.025  # volume [L]
            },
            'LiFePO4': {
                'E0': 3.2,
                'K': 0.003,
                'Q': 2.8,
                'A': 0.05,
                'B': 0.015,
                'R': 0.008,
                'peso': 0.07,
                'volume': 0.03
            }
        }
        return catalogo.get(tipo)

    def calcular_derivadas(self, corrente):
        """Derivadas dos estados da bateria para a integração numérica (RK4).

        Três estados evoluem no vetor de estado da Simulation:
          • SoC:  dSoC/dt = −i / (capacidade em C)   → descarga linear
          • Iast: dIast/dt = i                       → integral da corrente,
            usada pelo termo de polarização de calcular_tensao()
          • tempo_acumulado: dt/dt = 1               → relógio interno do
            ramo de carga do modelo

        Parameters
        ----------
        corrente : float
            Corrente DC aplicada [A] (positiva = descarga, negativa = carga).
            IMPORTANTE: deve ser a corrente DC do barramento (balanço de
            potência), nunca isq — ver invariante no cabeçalho do módulo.

        Returns
        -------
        tuple (dsoc_dt, dIast_dt, dtime_dt)
        """
        dsoc_dt = -corrente * self.inv_carga_total
        dIast_dt = corrente
        dtime_dt = 1.0
        return dsoc_dt, dIast_dt, dtime_dt

    def calcular_tensao(self, corrente, soc, Iast, tempo_acumulado):
        """Tensão terminal do pack [V] sob a corrente e o estado dados.

        Equação de Shepherd modificada (ramo de descarga), por célula:

            Vt = E0 − R·i − K·(Q/(Q−q))·(Iast/3600) − K·(Q/(Q−q))·q + A·e^(−B·q)
                 └──┬──┘   └──────────┬───────────┘   └─────┬─────┘   └───┬───┘
                 queda      polarização dinâmica       polarização     região
                 ôhmica      (histórico de corrente)    estática      exponencial

        onde q = (1 − SoC)·Q é a carga já retirada da célula [Ah].

        Parameters
        ----------
        corrente : float
            Corrente aplicada [A] em nível de PACK (positiva = descarga).
        soc : float
            Estado de carga corrente (0.0–1.0).
        Iast : float
            Corrente acumulada [C] em nível de PACK.
        tempo_acumulado : float
            Tempo decorrido [s] (só usado no ramo de carga).

        Returns
        -------
        float
            Tensão terminal do PACK [V] (= tensão da célula × n_serie).

        Notes
        -----
        O modelo de Shepherd é avaliado POR CÉLULA e o resultado escalado
        por n_serie.  Os estados chegam em nível de PACK (corrente total,
        Iast total, capacidade total), então são convertidos para célula:
            i_cell    = corrente / n_paralelo
            q_cell    = (1 − soc) · Q          [Ah por célula]
            Iast_cell = Iast / n_paralelo
        Sem essa divisão, a corrente do pack inteiro era aplicada à
        resistência de UMA célula — com n_paralelo > 1 o sag de tensão
        ficava n_paralelo× superestimado e o barramento colapsava ao piso.
        (Correção crítica — NÃO regredir.)
        """
        # Conversão pack → célula (ver Notes)
        i_cell = corrente / self.n_paralelo
        q_cell = (1 - soc) * self.Q                # carga usada por célula [Ah]
        Iast_cell = Iast / self.n_paralelo

        # Termos da equação de Shepherd (denom1 protegido contra célula vazia)
        denom1 = max(self.Q - q_cell, 1e-6)
        term_exp = self.A * np.exp(-self.B * q_cell)     # região exponencial
        term_K = self.K * self.Q / denom1                # ganho de polarização
        term_carga = term_K * q_cell                     # polarização estática
        term_Iast = term_K * (Iast_cell / 3600)          # polarização dinâmica

        if corrente >= 0:  # Descarga (caso normal da prova de aceleração)
            Vt = self.E0 - self.R * i_cell - term_Iast - term_carga + term_exp
        else:  # Carga (regeneração — não usada no modelo atual, mantida por completude)
            i_cell_abs = abs(i_cell)
            denom2 = max(i_cell_abs * tempo_acumulado - 0.1 * self.Q, 1e-6)
            term_K_charge = self.K * self.Q / denom2 * (Iast_cell / 3600)
            Vt = self.E0 - self.R * i_cell_abs - term_K_charge - term_carga + term_exp

        # Escala célula → pack: células em série somam tensão
        return Vt * self.n_serie

    def calcular_peso_total(self):
        """Massa total do pack [kg] (n_serie × n_paralelo × massa da célula)."""
        return self.n_serie * self.n_paralelo * self.peso

    def calcular_volume_total(self):
        """Volume total do pack [L] (n_serie × n_paralelo × volume da célula)."""
        return self.n_serie * self.n_paralelo * self.volume

    def calcular_tensao_nominal(self):
        """Tensão nominal do pack [V] (E0 × n_serie, sem carga).

        Usada pela Simulation como denominador do balanço de potência
        (I_bat = P_ac / V_nominal) e como base do piso de tensão Vdc_floor.
        Ex.: 144S Li-ion → 144 × 3.7 = 532.8 V (abaixo do teto FSAE de 600 V).
        """
        return self.E0 * self.n_serie
