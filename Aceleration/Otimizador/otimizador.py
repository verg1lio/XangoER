"""
otimizador.py — Evolução Diferencial para relação de transmissão + pedal

Minimiza o tempo para cruzar ``dmax`` metros variando:
  x[0] = N      — relação de transmissão final (Z_coroa / Z_pinhao)
  x[1] = t_pico — instante [s] em que o pedal atinge 100 % (rampa
                  smoothstep 3s²−2s³; t_pico=0 → degrau de pedal)

Mudanças v2 (2026-06):
  • SIM_CONFIG centraliza todos os parâmetros — edite aqui para sincronizar
    com os valores do dashboard (main.py).
  • Transmissão atualizada: chain_efficiency, diff_efficiency, sprocket_inertia,
    coroa_inertia, lsd_bias_ratio (parâmetros da arquitetura atual).
  • Veículo atualizado: asa dianteira/traseira (lift_coeff_front/rear, area_front/rear)
    e resistência ao rolamento velocidade-dependente Cr(v).
  • Penalidade de slip removida do fitness.
    Motivo: o Pacejka calibrado (PCX1=1.5, PEX1=-0.5, κ_peak=0.204) já penaliza
    patinagem via queda de Fx — adicionar penalidade no fitness criaria viés duplo
    contra relações que naturalmente operam com mais slip no lançamento.
    Adicionalmente, a penalidade era silenciosamente nula na v1 por bug:
    results.get('slip_ratio_hist') → chave correta é 'slip_ratio'.
  • find_integer_ratio(): pós-otimização, converte N* contínuo no par (Z₁, Z₂)
    inteiro mais próximo, com tabela de candidatos ordenados por erro.
  • sweep_diagnostico(): varredura ASCII do landscape N × t_cruzamento antes
    da DE — confirma que o ótimo não está em platô ou descontinuidade.
  • Progresso da DE exibe t_75m, slip_max e SoC por iteração.

Mudanças v3 (2026-06) — revisão física do envelope do powertrain:
  • speed_ref 1200 → 575.96 rad/s (5500 RPM, limite real do EMRAX 228).
    Com 1200 rad/s a DE convergia para N* ≈ 12, pois o ótimo de aceleração
    é N* ≈ ω_max·r/v_final e o motor "girava" o dobro do fisicamente possível.
  • p_max_dc = 80 kW (FSAE EV.4.1) — antes não havia limite de potência e o
    modelo entregava torque constante até o limitador (>300 kW).
  • max_current = 323 A pico → torque de pico 230 N·m (datasheet), antes 302.
  • kf 0.1 → 0.005 N·m·s/rad (0.1 dissipava ~33 kW em atrito viscoso).
  • Bateria 264S1P (976 V!) → 144S5P (532.8 V, dentro do teto FSAE de 600 V).
  • wheel_radius 0.275 → 0.26 m (Hoosier 20.5x7.0-13).
  • Tração do eixo corrigida em Simulation: teto = 2 × Fx(Fz_por_roda).

Mudanças v4 (2026-06) — otimização conjunta transmissão + pedal:
  • Espaço de busca 2D: [N, t_pico_pedal].  A rampa do pedal funciona como
    launch control primitivo — suaviza o "penhasco" de wheelspin que antes
    forçava N* a ficar logo abaixo do limiar de patinagem com pedal-degrau.
  • DE/rand/1/bin vetorial com crossover binomial por dimensão (j_rand
    garante ao menos um gene do mutante).
  • sweep_diagnostico agora imprime grade 2D N × t_pico.
  • simulation() aceita escalar N (retrocompatível, t_pico=0) ou [N, t_pico].
"""

import numpy as np
import sys
import os
import io
import contextlib

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Models import BatteryPack, Tire, Transmission, Vehicle, Motor, Pedal
from Simulation.Simulation import Simulation


# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURAÇÃO — sincronize com os valores do dashboard (main.py)
# ─────────────────────────────────────────────────────────────────────────────
SIM_CONFIG = {
    # ── Transmissão ───────────────────────────────────────────────────────────
    'chain_efficiency':   0.98,     # corrente lubrificada FSAE
    'diff_efficiency':    0.97,     # diferencial LSD
    'axle_inertia':       0.0015,   # [kg·m²] por semi-eixo
    'diff_inertia':       0.012,    # [kg·m²] corpo do diferencial
    'sprocket_inertia':   0.0002,   # [kg·m²] pinhão (lado motor, sem reflexão)
    'coroa_inertia':      0.0015,   # [kg·m²] coroa (lado diferencial, refletida)
    'lsd_bias_ratio':     3.0,      # razão de torque bias do LSD
    # ── Veículo ───────────────────────────────────────────────────────────────
    'mass':               240.0,    # [kg] massa total (piloto incluso)
    'wheel_radius':       0.26,     # [m] raio efetivo — Hoosier 20.5x7.0-13
    'wheel_mass':         5.0,      # [kg] por conjunto roda+pneu
    'drag_coeff':         0.7789,   # Cd aerodinâmico
    'frontal_area':       0.68,     # [m²] área frontal
    'rolling_resistance': 0.015,    # Cr0 base; Cr(v) = Cr0·(1 + v/150)
    'L':                  1.5,      # [m] entre-eixos
    'h':                  0.28,     # [m] altura do CG
    'dist_cg':            0.6,      # [m] distância eixo dianteiro → CG
    'has_wing':           True,
    'lift_coeff_front':   1.0,      # Cl asa dianteira
    'area_front':         0.35,     # [m²] área asa dianteira
    'lift_coeff_rear':    1.0,      # Cl asa traseira
    'area_rear':          0.35,     # [m²] área asa traseira
    # ── Bateria ───────────────────────────────────────────────────────────────
    # 144S → 532.8 V nominal (abaixo do teto de 600 V da FSAE EV).
    # 5P → R_pack = 144×0.02/5 = 0.576 Ω e ~30 A por célula a 80 kW.
    # SUBSTITUA pela configuração real do acumulador quando definida.
    'tipo_celula':        'Li-ion',
    'n_serie':            144,
    'n_paralelo':         5,
    'soc_inicial':        1.0,
    # ── Motor — EMRAX 228 ─────────────────────────────────────────────────────
    'rs':                 0.00706,
    'ld':                 0.0000965,
    'lq':                 0.0000965,
    'jm':                 0.02521,
    'kf':                 0.005,    # [N·m·s/rad] atrito viscoso (rolamentos/ventilação).
                                    # 0.1 era irreal: 57 N·m (~33 kW) a 5500 RPM, e as
                                    # perdas no ferro já são modeladas separadamente.
    'lambda_m':           0.04748,
    'p':                  10,
    'valor_mu':           0.99,
    'speed_ref':          575.96,   # [rad/s] = 5500 RPM — limite do EMRAX 228.
                                    # (1200 rad/s = 11459 RPM era 2× o limite real e
                                    #  inflava N* para ~12: N* ≈ ω_max·r/v_final.)
    'max_current':        323.0,    # [A pico] → T_pico = 1.5·p·λm·323 ≈ 230 N·m
                                    # (datasheet EMRAX 228; 300√2=424 A dava 302 N·m)
    # ── Powertrain / regulamento ──────────────────────────────────────────────
    'p_max_dc':           80e3,     # [W] FSAE EV.4.1 — limite no barramento DC
    # ── Simulação ─────────────────────────────────────────────────────────────
    'tmax':               30.0,     # [s] limite de segurança
    'dmax':               75.0,     # [m] distância alvo (Acceleration FSAE)
}


# ─────────────────────────────────────────────────────────────────────────────
# Funções auxiliares internas
# ─────────────────────────────────────────────────────────────────────────────

def _build_simulation(i_ratio: float, cfg: dict,
                      t_pico: float = 0.0) -> Simulation:
    """Instancia todos os modelos com N=i_ratio e retorna um objeto Simulation.

    t_pico > 0 → pedal com rampa smoothstep 0→100 % em t_pico segundos.
    t_pico = 0 → sem pedal (Simulation usa throttle 100 % constante).
    """
    trans = Transmission.Transmission(
        final_drive_ratio=float(i_ratio),
        chain_efficiency=cfg['chain_efficiency'],
        diff_efficiency=cfg['diff_efficiency'],
        axle_inertia=cfg['axle_inertia'],
        diff_inertia=cfg['diff_inertia'],
        sprocket_inertia=cfg['sprocket_inertia'],
        coroa_inertia=cfg['coroa_inertia'],
        lsd_bias_ratio=cfg['lsd_bias_ratio'],
    )
    veh = Vehicle.Vehicle(
        mass=cfg['mass'],
        wheel_radius=cfg['wheel_radius'],
        wheel_mass=cfg['wheel_mass'],
        drag_coeff=cfg['drag_coeff'],
        frontal_area=cfg['frontal_area'],
        rolling_resistance=cfg['rolling_resistance'],
        L=cfg['L'], h=cfg['h'], dist_cg=cfg['dist_cg'],
        has_wing=cfg['has_wing'],
        lift_coeff_front=cfg['lift_coeff_front'], area_front=cfg['area_front'],
        lift_coeff_rear=cfg['lift_coeff_rear'],   area_rear=cfg['area_rear'],
    )
    bat = BatteryPack.BatteryPack(
        tipo_celula=cfg['tipo_celula'],
        n_serie=cfg['n_serie'],
        n_paralelo=cfg['n_paralelo'],
        soc_inicial=cfg['soc_inicial'],
    )
    tire = Tire.Tire(
        pacejka_params=Tire.Tire.HOOSIER_FSAE_LONG,
        tire_friction_coef=0.6,
    )
    mot = Motor.Motor(
        rs=cfg['rs'], ld=cfg['ld'], lq=cfg['lq'],
        jm=cfg['jm'], kf=cfg['kf'], lambda_m=cfg['lambda_m'],
        p=cfg['p'], valor_mu=cfg['valor_mu'], speed_ref=cfg['speed_ref'],
        max_current=cfg.get('max_current'),
    )
    ped = None
    if t_pico is not None and float(t_pico) > 1e-3:
        ped = Pedal.Pedal(ganho=1.0)
        tp, pp = Pedal.Pedal.perfil_suave(0.0, float(t_pico))
        ped.set_profile(tp, pp)
    return Simulation(
        motor=mot, vehicle=veh, transmission=trans,
        battery=bat, tire=tire, pedal=ped,
        tmax=cfg['tmax'], dmax=cfg['dmax'],
        p_max_dc=cfg.get('p_max_dc', 80e3),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Função de aptidão
# ─────────────────────────────────────────────────────────────────────────────

def simulation(x, cfg: dict = None):
    """Avalia a aptidão de um candidato.

    Parameters
    ----------
    x : float ou sequência
        Escalar N (retrocompatível, pedal-degrau) ou vetor [N, t_pico]:
          x[0] = relação de transmissão final N
          x[1] = instante [s] em que o pedal atinge 100 % (0 → degrau)

    Fitness = t_75m [s] sem penalidade de slip.

    O modelo de Pacejka calibrado (PCX1=1.5, PEX1=-0.5) penaliza
    patinagem fisicamente — candidatos que geram slip excessivo produzem
    menos tração e, portanto, maior t_75m de forma natural.

    Returns
    -------
    fitness : float
        Tempo de cruzamento [s], ou constante grande em falha/não-chegada.
    info : dict
        reached, t_75m, slip_max, slip_mean, soc_final.
    """
    if cfg is None:
        cfg = SIM_CONFIG

    if np.isscalar(x):
        i_ratio, t_pico = float(x), 0.0
    else:
        x = np.asarray(x, dtype=float).ravel()
        i_ratio = float(x[0])
        t_pico  = float(x[1]) if x.size > 1 else 0.0

    dmax      = cfg['dmax']
    _no_reach = dict(reached=False, t_75m=None,
                     slip_max=None, slip_mean=None, soc_final=None)

    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            sim     = _build_simulation(i_ratio, cfg, t_pico=t_pico)
            results = sim.simulate(t0=0, tf=cfg['tmax'])
    except Exception as e:
        sys.stderr.write(f"[EXCECAO N={i_ratio:.3f} t_pico={t_pico:.2f}]: {e}\n")
        return 1e9, _no_reach

    if results is None:
        return 1e9, _no_reach

    pos = results.get('vehicle_position', np.array([]))
    t   = results.get('t',                np.array([]))

    if len(pos) == 0 or float(pos[-1]) < 1.0:
        return 1e5, _no_reach

    if float(pos[-1]) < dmax:
        return 100.0 + (dmax - float(pos[-1])), _no_reach

    idx      = int(np.argmax(pos >= dmax))
    t_cruzou = float(t[idx])
    if np.isnan(t_cruzou):
        return 1e9, _no_reach

    slip_arr = np.abs(np.asarray(results.get('slip_ratio', []), dtype=float)[:idx + 1])
    soc_arr  = np.asarray(results.get('soc', []), dtype=float)

    info = dict(
        reached   = True,
        t_75m     = t_cruzou,
        slip_max  = float(np.max(slip_arr))  if slip_arr.size > 0 else None,
        slip_mean = float(np.mean(slip_arr)) if slip_arr.size > 0 else None,
        soc_final = float(soc_arr[idx])      if soc_arr.size > idx else None,
    )
    return t_cruzou, info


# ─────────────────────────────────────────────────────────────────────────────
# Snap para par inteiro de dentes
# ─────────────────────────────────────────────────────────────────────────────

def find_integer_ratio(n_opt: float,
                       z1_range: tuple = (10, 20),
                       z2_range: tuple = (20, 100),
                       n_candidates: int = 5):
    """Converte relação contínua N* no par (Z₁, Z₂) de dentes inteiros mais próximo.

    Percorre Z_pinhao ∈ z1_range, calcula Z_coroa ideal e testa floor/ceil.
    Ordena candidatos por erro percentual |N_inteiro − N*| / N*.

    Parameters
    ----------
    n_opt : float
        Relação de transmissão contínua ótima encontrada pela DE.
    z1_range : tuple
        (min, max) de dentes do pinhão. Padrão (10, 20) — sprockets FSAE típicos.
    z2_range : tuple
        (min, max) de dentes da coroa. Padrão (20, 100).
    n_candidates : int
        Número máximo de candidatos retornados.

    Returns
    -------
    list of dict
        Candidatos com chaves: Z1, Z2, N, erro_pct.
    """
    candidatos = []
    for z1 in range(z1_range[0], z1_range[1] + 1):
        z2_ideal = n_opt * z1
        for z2 in (int(np.floor(z2_ideal)), int(np.ceil(z2_ideal))):
            if z2_range[0] <= z2 <= z2_range[1]:
                N_real = z2 / z1
                erro   = abs(N_real - n_opt)
                candidatos.append(dict(Z1=z1, Z2=z2, N=N_real,
                                       erro_pct=100.0 * erro / max(n_opt, 1e-9)))

    candidatos.sort(key=lambda c: c['erro_pct'])
    vistos, resultado = set(), []
    for c in candidatos:
        key = (c['Z1'], c['Z2'])
        if key not in vistos:
            vistos.add(key)
            resultado.append(c)
            if len(resultado) >= n_candidates:
                break
    return resultado


# ─────────────────────────────────────────────────────────────────────────────
# Varredura de diagnóstico (landscape)
# ─────────────────────────────────────────────────────────────────────────────

def sweep_diagnostico(distancia: float,
                      n_min: float = 3.0,
                      n_max: float = 8.0,
                      n_points: int = 8,
                      t_picos: tuple = (0.0, 0.5, 1.0, 2.0),
                      cfg: dict = None):
    """Varre a grade N × t_pico e exibe o landscape t_75m em ASCII.

    Permite verificar se o ótimo está em uma região suave ou se há
    descontinuidades no espaço de busca antes de rodar a DE.

    Returns
    -------
    list of (N, t_pico, fitness, info)
    """
    if cfg is None:
        cfg = SIM_CONFIG

    valores_N = np.linspace(n_min, n_max, n_points)
    registros = []
    total     = n_points * len(t_picos)

    print(f"\n  Varredura landscape: grade {n_points}x{len(t_picos)}  "
          f"N=[{n_min:.1f}, {n_max:.1f}]  t_pico={list(t_picos)} s  "
          f"dmax={distancia} m")
    k = 0
    grade = {}
    for N in valores_N:
        for tp in t_picos:
            k += 1
            sys.stdout.write(f"\r  [{k:02d}/{total}] N={N:5.2f} t_pico={tp:.2f} ...")
            sys.stdout.flush()
            f, info = simulation([N, tp], cfg)
            registros.append((N, tp, f, info))
            grade[(N, tp)] = (f, info)
    print()

    alcancou = [(N, tp, f) for (N, tp, f, info) in registros if info['reached']]
    if not alcancou:
        print("  Nenhuma configuracao atingiu a distancia alvo nesta faixa.")
        return registros

    # Tabela: linhas = N, colunas = t_pico, célula = t_75m [s]
    header = "  ".join(f"t={tp:4.2f}s" for tp in t_picos)
    print(f"\n  {'N':>6}  {header}")
    print(f"  {'-'*6}  " + "  ".join('-' * 7 for _ in t_picos))
    for N in valores_N:
        celulas = []
        for tp in t_picos:
            f, info = grade[(N, tp)]
            celulas.append(f"{f:7.4f}" if info['reached'] else '    ---')
        print(f"  {N:6.2f}  " + "  ".join(celulas))

    N_b, tp_b, f_b = min(alcancou, key=lambda r: r[2])
    print(f"\n  Melhor ponto da grade: N={N_b:.2f}  t_pico={tp_b:.2f} s  "
          f"t_75m={f_b:.4f} s\n")
    return registros


# ─────────────────────────────────────────────────────────────────────────────
# Evolução Diferencial — DE/rand/1/bin
# ─────────────────────────────────────────────────────────────────────────────

def diferentialEvo(
        distancia_desejada: float = 75.0,
        i_min: float = 3.0,
        i_max: float = 20.0,
        t_pico_min: float = 0.0,
        t_pico_max: float = 3.0,
        pop_size: int = 20,
        max_iter: int = 30,
        CR: float = 0.9,
        F_scale: float = 0.8,
        tol: float = 1e-4,
        cfg: dict = None,
        executar_sweep: bool = False):
    """Otimiza [N, t_pico_pedal] via Evolução Diferencial DE/rand/1/bin.

    Variáveis de decisão:
      x[0] = N      ∈ [i_min, i_max]          relação de transmissão final
      x[1] = t_pico ∈ [t_pico_min, t_pico_max] instante [s] em que o pedal
             atinge 100 % (smoothstep; 0 → degrau).  A rampa atua como
             launch control: sem ela, o pedal-degrau força N* a ficar logo
             abaixo do limiar de wheelspin.

    Parameters
    ----------
    executar_sweep : bool
        Se True, roda sweep_diagnostico (grade 2D) antes da DE.
        Recomendado na primeira execução para validar os limites de busca.

    Returns
    -------
    best_x : np.ndarray shape (2,)
        [N*, t_pico*] ótimo encontrado.
    best_fitness : float
        t_75m correspondente [s].
    fitness_history : list
        Melhor fitness por iteração (para analisar convergência).
    x_history : list of np.ndarray
        Melhor [N, t_pico] por iteração.
    """
    if cfg is None:
        cfg = SIM_CONFIG

    lo  = np.array([i_min, t_pico_min], dtype=float)
    hi  = np.array([i_max, t_pico_max], dtype=float)
    dim = 2

    print(f"\n{'='*62}")
    print(f"  Otimizador de Transmissao + Pedal - Xango E-Racing")
    print(f"  Distancia alvo : {distancia_desejada} m")
    print(f"  Espaco de busca: N in [{i_min:.1f}, {i_max:.1f}]   "
          f"t_pico in [{t_pico_min:.1f}, {t_pico_max:.1f}] s")
    print(f"  Populacao={pop_size}   max_iter={max_iter}   CR={CR}   F={F_scale}   tol={tol}")
    print(f"{'='*62}\n")

    if executar_sweep:
        sweep_diagnostico(distancia_desejada,
                          n_min=i_min, n_max=min(i_max, 8.0),
                          cfg=cfg)

    # ── Inicializa população (pop_size × 2) ───────────────────────────────────
    pop     = np.random.uniform(lo, hi, size=(pop_size, dim))
    fit     = np.full(pop_size, np.inf)
    inf_pop = [None] * pop_size

    print(f"  Avaliando populacao inicial ({pop_size} individuos)...")
    for k in range(pop_size):
        sys.stdout.write(f"\r  [{k+1:02d}/{pop_size}] "
                         f"N={pop[k,0]:.3f} t_pico={pop[k,1]:.2f}  ")
        sys.stdout.flush()
        fit[k], inf_pop[k] = simulation(pop[k], cfg)
    print()

    best_idx     = int(np.argmin(fit))
    best_fitness = float(fit[best_idx])
    best_x       = pop[best_idx].copy()
    best_info    = inf_pop[best_idx]

    fitness_history = []
    x_history       = []

    # ── Loop DE principal ─────────────────────────────────────────────────────
    def _progress(iteration, i_indiv, n_indiv, trial_x, result_info):
        """Atualiza a linha de progresso após cada avaliação individual."""
        t75_str  = (f"{result_info['t_75m']:.4f}s"
                    if result_info and result_info.get('t_75m') else '    --- ')
        slip_str = (f"{result_info['slip_max']:.3f}"
                    if result_info and result_info.get('slip_max') is not None else '---')
        sys.stdout.write(
            f"\r  iter {iteration+1:02d}/{max_iter} "
            f"[{i_indiv+1:02d}/{n_indiv}] "
            f"trial=({trial_x[0]:.3f}, {trial_x[1]:.2f}s)  "
            f"best: N={best_x[0]:.4f} t_pico={best_x[1]:.3f}s "
            f"t={t75_str} slip={slip_str}  ")
        sys.stdout.flush()

    for iteration in range(max_iter):
        for i in range(pop_size):
            idxs = [j for j in range(pop_size) if j != i]
            if len(idxs) < 3:
                continue
            a, b, c = np.random.choice(idxs, 3, replace=False)

            # Mutação: DE/rand/1
            mutant = pop[a] + F_scale * (pop[b] - pop[c])
            # Cruzamento binomial por dimensão; j_rand garante ao menos
            # um gene vindo do mutante (DE/rand/1/bin canônico)
            cross = np.random.rand(dim) < CR
            cross[np.random.randint(dim)] = True
            trial = np.where(cross, mutant, pop[i])
            trial = np.clip(trial, lo, hi)

            f_trial, info_trial = simulation(trial, cfg)
            if f_trial < fit[i]:
                pop[i]     = trial
                fit[i]     = f_trial
                inf_pop[i] = info_trial
                if f_trial < best_fitness:
                    best_fitness = f_trial
                    best_x       = trial.copy()
                    best_info    = info_trial

            # Atualiza progresso a cada avaliação — sem silêncio entre iterações
            _progress(iteration, i, pop_size, trial, best_info)

        fitness_history.append(best_fitness)
        x_history.append(best_x.copy())
        print()   # quebra de linha ao final de cada iteração completa

        # Convergência: std das últimas 10 melhores < tol
        if len(fitness_history) >= 10:
            if np.std(fitness_history[-10:]) < tol:
                print(f"  Convergencia na iteracao {iteration + 1}.")
                break

    print(f"\n")

    # ── Relatório final ───────────────────────────────────────────────────────
    best_N, best_tp = float(best_x[0]), float(best_x[1])
    print(f"{'='*62}")
    print(f"  RESULTADO FINAL")
    print(f"{'-'*62}")
    print(f"  Relacao continua otima  : N* = {best_N:.4f}")
    print(f"  Pico do pedal otimo     : t_pico* = {best_tp:.3f} s "
          f"{'(degrau)' if best_tp < 1e-3 else '(rampa smoothstep)'}")
    if best_info and best_info.get('t_75m'):
        print(f"  Tempo t_75m             : {best_info['t_75m']:.4f} s")
        sm = best_info.get('slip_max')
        sv = best_info.get('slip_mean')
        sc = best_info.get('soc_final')
        print(f"  Slip maximo (ate 75 m)  : {sm:.3f}" if sm is not None else
              f"  Slip maximo (ate 75 m)  : ---")
        print(f"  Slip medio  (ate 75 m)  : {sv:.3f}" if sv is not None else
              f"  Slip medio  (ate 75 m)  : ---")
        print(f"  SoC ao cruzar 75 m      : {sc:.4f}" if sc is not None else
              f"  SoC ao cruzar 75 m      : ---")
    else:
        print(f"  Fitness bruto           : {best_fitness:.4f}")
    print(f"{'-'*62}")

    # ── Snap para dentes inteiros (avaliado com t_pico*) ─────────────────────
    print(f"  Snap para pares inteiros  (Z_pinhao in [9, 18], t_pico={best_tp:.3f} s):")
    print(f"\n  {'Z1':>4}  {'Z2':>4}  {'N real':>8}  {'Erro':>7}  Avaliado")
    print(f"  {'-'*4}  {'-'*4}  {'-'*8}  {'-'*7}  {'-'*10}")
    candidatos = find_integer_ratio(best_N,
                                    z1_range=(9, 18),
                                    z2_range=(20, 180))
    avaliados  = []
    for c in candidatos:
        f_int, info_int = simulation([c['N'], best_tp], cfg)
        c['t_75m'] = info_int.get('t_75m')
        avaliados.append(c)
        t_str = f"{info_int['t_75m']:.4f} s" if info_int.get('t_75m') else '  ---'
        print(f"  {c['Z1']:>4}  {c['Z2']:>4}  {c['N']:8.4f}  "
              f"{c['erro_pct']:6.2f}%  {t_str}")

    # Melhor candidato inteiro
    valid = [c for c in avaliados if c.get('t_75m') is not None]
    if valid:
        best_int = min(valid, key=lambda c: c['t_75m'])
        delta    = best_int['t_75m'] - (best_info.get('t_75m') or best_fitness)
        print(f"\n  Melhor par inteiro: Z1={best_int['Z1']}  Z2={best_int['Z2']}"
              f"  N={best_int['N']:.4f}  t={best_int['t_75m']:.4f} s"
              f"  (delta={delta:+.4f} s vs. N* continuo)")

    print(f"{'='*62}\n")

    return best_x, best_fitness, fitness_history, x_history


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    best_x, best_fitness, fit_hist, x_hist = diferentialEvo(
        distancia_desejada = SIM_CONFIG['dmax'],
        i_min              = 3.0,
        i_max              = 20.0,
        t_pico_min         = 0.0,
        t_pico_max         = 3.0,
        pop_size           = 20,
        max_iter           = 30,
        CR                 = 0.9,
        F_scale            = 0.8,
        tol                = 1e-4,
        executar_sweep     = True,
    )
