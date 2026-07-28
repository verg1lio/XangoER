# -*- coding: utf-8 -*-
"""
constants.py — Constantes matemáticas pré-computadas do projeto
================================================================

Este módulo centraliza constantes numéricas usadas repetidamente pelas
transformações de coordenadas do motor (Park/Clarke) e por outros cálculos
trigonométricos do digital twin.

Por que pré-computar?
---------------------
O loop de simulação roda a 10 kHz (passo de 1e-4 s) — cada operação evitada
dentro do loop importa. Avaliar `np.sqrt(3)/2` ou `2*np.pi/3` uma única vez
na importação do módulo elimina milhares de recomputações idênticas.

Constantes disponíveis
----------------------
SQRT3      : √3            — aparece nas transformações trifásicas (Clarke).
SQRT3_2    : √3/2 ≈ 0.866  — projeção das fases b/c no plano αβ; usada na
                             transformada inversa de Park (Motor.py).
PI         : π             — alias de np.pi por conveniência.
PI23       : 2π/3 ≈ 2.094  — defasagem angular de 120° entre as três fases
                             do motor; usada para gerar ia, ib, ic a partir
                             do referencial dq (Motor.abc_currents_from_dq).
RQ23       : √(2/3)        — fator de escala da transformada de Park
                             INVARIANTE EM POTÊNCIA. Mantido como referência;
                             o projeto usa a forma invariante em AMPLITUDE
                             (sem esse fator), consistente com Kt = 1.5·p·λm.
TWO_THIRDS : 2/3           — fator da transformada de Clarke invariante em
                             amplitude (abc → αβ).
"""

import numpy as np


SQRT3 = np.sqrt(3)
SQRT3_2 = SQRT3 / 2
PI = np.pi
PI23 = 2 * PI / 3
RQ23 = np.sqrt(2 / 3)
TWO_THIRDS = 2 / 3
