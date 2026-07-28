# DigitalTwin Aceleration — Guia Didático Completo
### Xangô E-Racing · Fórmula SAE Elétrico · UFBA
**Versão didática v1.0 — baseada na documentação técnica v0.1.5**
Autor do projeto e da documentação: Marco Affonso de Carvalho Santos

---

## Para quem é este documento

Para **você que acabou de entrar na equipe** e nunca viu nada disso na vida.

Não pressupomos que você saiba o que é um PMSM, uma transformada de Park, um Runge-Kutta ou a Magic Formula de Pacejka. Não pressupomos nem que você já tenha feito Dinâmica ou Circuitos Elétricos. O que pressupomos é: você sabe que força é massa vezes aceleração, sabe o que é uma derivada (mesmo que só de nome) e consegue ler um pedaço de código Python sem entrar em pânico.

Tudo o mais está explicado aqui, do zero, na ordem certa.

> **Aviso honesto sobre o tamanho.** Este documento é longo de propósito. Ele não foi feito para ser lido de uma sentada — foi feito para ser lido em partes e depois consultado. As Partes I e II você lê uma vez. As Partes III e IV você vai reabrir toda vez que precisar mexer num módulo específico.

---

## Como usar este guia

O documento tem seis partes, em ordem crescente de dificuldade:

| Parte | Conteúdo | Para quem |
|---|---|---|
| **0** | O que é o projeto, qual problema ele resolve | Todo mundo |
| **I** | Fundamentos: física, matemática, eletricidade, controle | Quem nunca viu o assunto |
| **II** | Arquitetura do software — como as peças se conectam | Todo mundo, antes de mexer no código |
| **III** | Os 7 módulos de modelo, um por um | Quem vai mexer em algum modelo |
| **IV** | O motor de simulação (`Simulation.py`) em profundidade | Quem vai mexer no núcleo |
| **V** | Dashboard e otimizador | Quem vai *usar* a ferramenta |
| **VI** | Como ler resultados, limitações, exercícios, glossário | Todo mundo |

**Três trilhas sugeridas:**

- 🟢 **"Só quero usar a ferramenta"** → Parte 0 → Parte II (§8) → Parte V → Parte VI (§28)
- 🟡 **"Quero entender a física"** → Parte 0 → Parte I inteira → Parte III → Parte VI
- 🔴 **"Vou mexer no código"** → tudo, na ordem

---

## Convenções deste documento

- **Blocos "🧠 Intuição"** dão a analogia antes da matemática. Se a equação te assustar, leia primeiro a intuição.
- **Blocos "⚠️ Não regredir"** marcam decisões de projeto que já causaram bug grave uma vez. Se você for mexer nelas, leia o histórico antes.
- **Blocos "🔍 No código"** mostram exatamente onde, em qual arquivo e qual linha, aquela física vira Python.
- **Blocos "✏️ Exercício"** são tarefas curtas para você fazer rodando o código. É a melhor forma de fixar.
- Símbolos matemáticos aparecem em Unicode (ω, λ, κ, η, ρ) porque é assim que estão nos comentários do código também.

---

## Índice

**Parte 0 — O problema e a ideia**
- 0.1 O que é um "gêmeo digital"? · 0.2 A prova Acceleration · 0.3 O carro simulado · 0.4 O que o programa faz · 0.5 Um resultado de referência

**Parte I — Fundamentos**
- [1. Física do movimento](#1-física-do-movimento-as-duas-leis-de-newton-que-usamos) — Newton translacional e rotacional, inércia, potência
- [2. Derivadas, ODEs e integração numérica](#2-derivadas-odes-e-integração-numérica) — vetor de estado, Euler, RK4, por que 10 kHz
- [3. Eletricidade e magnetismo](#3-eletricidade-e-magnetismo-o-mínimo-necessário) — tensão, corrente, indutância, fluxo, back-EMF
- [4. Como um PMSM produz torque](#4-como-um-pmsm-produz-torque) — campo girante, pares de polos, Kt
- [5. O referencial dq (Park)](#5-o-referencial-dq-transformada-de-park--a-ideia-mais-importante-do-controle) — a analogia do carrossel, FOC, convenção amplitude-invariante
- [6. Controle: PI, saturação e windup](#6-controle-pi-saturação-e-windup) — malha aberta vs fechada, cancelamento de polo, anti-windup
- [7. Unidades e sinais](#7-nota-sobre-unidades-e-sinais)

**Parte II — Arquitetura do software**
- [8. O mapa do território](#8-o-mapa-do-território)
- [9. O princípio de design central](#9-o-princípio-de-design-central) — por que os modelos não têm estado
- [10. O ciclo de vida de uma simulação](#10-o-ciclo-de-vida-de-uma-simulação)
- [11. Quem chama quem](#11-quem-chama-quem-o-diagrama-de-dependências)

**Parte III — Os módulos, um por um**
- [12. `Motor.py`](#12-motorpy--o-motor-pmsm-emrax-228) — parâmetros, perdas no ferro, Rs(T), térmica
- [13. `BatteryPack.py`](#13-batterypackpy--a-bateria-modelo-de-shepherd-modificado) — Shepherd, ⚠️ escala pack↔célula
- [14. `Transmission.py`](#14-transmissionpy--o-trem-de-força) — relação N, reflexão de inércia por 1/N²
- [15. `Vehicle.py`](#15-vehiclepy--dinâmica-longitudinal-do-carro) — resistências, transferência de carga, ⚠️ dupla contagem de massa
- [16. `Tire.py`](#16-tirepy--o-pneu-pacejka-pac2002--magic-formula) — slip ratio, Magic Formula, ⚠️ calibração, wheelspin emergente
- [17. `Pedal.py`](#17-pedalpy--a-entrada-do-piloto) — smoothstep como launch control
- [18. `PIDController.py`](#18-pidcontrollerpy--o-controlador-pi-genérico) — anti-windup em dois estágios

**Parte IV — O motor de simulação**
- [19. `__init__` — a montagem do sistema](#19-__init__--a-montagem-do-sistema)
- [20. O loop principal — os 8 passos](#20-o-loop-principal--os-8-passos)
- [21. A ODE linha a linha](#21-a-ode-_physics_ode-linha-a-linha)
- [22. Os três limitadores de `iq_ref`](#22-os-três-limitadores-de-iq_ref)
- [23. Field weakening](#23-field-weakening-enfraquecimento-de-campo)
- [24. **Um passo inteiro, com números reais**](#24-um-passo-inteiro-com-números-reais) ← *comece por aqui se quiser ver tudo junto*
- [25. Os bugs históricos e o que NÃO regredir](#25-os-bugs-históricos-e-por-que-não-regredir)

**Parte V — As ferramentas**
- [26. O dashboard web](#26-o-dashboard-web-mainpy)
- [27. O otimizador de transmissão](#27-o-otimizador-de-transmissão-otimizadorotimizadorpy)

**Parte VI — Usando, criticando e estendendo**
- [28. Como ler os resultados](#28-como-ler-os-resultados-e-desconfiar-deles) — testes de sanidade, sintomas e causas
- [29. Limitações conhecidas](#29-limitações-conhecidas--o-que-este-modelo-não-faz)
- [30. Exercícios e respostas](#30-exercícios--e-as-respostas)
- [31. Glossário](#31-glossário)
- [32. FAQ e primeiros passos](#32-faq-e-primeiros-passos)

---

# Parte 0 — O problema e a ideia

## 0.1 O que é um "gêmeo digital"?

Um **gêmeo digital** (*digital twin*) é um modelo computacional de um sistema físico real, detalhado o bastante para que você possa fazer perguntas ao computador em vez de fazer ao carro.

Perguntas do tipo:

> *"Se eu trocar o pinhão de 13 para 12 dentes, o carro faz os 75 metros em menos tempo?"*

Você tem três formas de responder isso:

| Forma | Custo | Confiabilidade |
|---|---|---|
| Achismo / intuição de garagem | zero | baixa |
| Fabricar o pinhão e testar na pista | R$ + dias + risco de quebrar | alta, mas só para *aquele* pinhão |
| **Simular** | minutos de CPU | boa, se o modelo for honesto |

O gêmeo digital não substitui o teste em pista. Ele **reduz o espaço de busca**: em vez de testar 8 relações de transmissão na pista, você simula 600 no computador, escolhe as 2 melhores e testa essas.

> 🧠 **Intuição**
> Pense no gêmeo digital como um "carro de mentira que obedece às mesmas leis da física do carro de verdade". Quanto mais leis você ensinar a ele — atrito do pneu, queda de tensão da bateria, aquecimento do motor — mais o carro de mentira se comporta como o de verdade. E cada lei que você *não* ensina é uma fonte de erro conhecida (listamos todas na §29).

## 0.2 A prova que estamos simulando: Acceleration

A Fórmula SAE tem várias provas dinâmicas. A mais simples de todas — e por isso a primeira a ser modelada — é a **Acceleration**:

```
   PARTIDA                                                  CHEGADA
   (carro parado)                                            
      │                                                         │
      ├─────────────────── 75 metros em linha reta ─────────────┤
      │                                                         │
   t = 0                                                    t = t₇₅
```

Regras que importam para o modelo:
- Largada **parada** (velocidade zero, motor parado);
- Pista **reta** e **plana**;
- A pontuação depende do **tempo** para cruzar os 75 m;
- Limite regulamentar (FSAE EV.4.1): a potência no barramento DC não pode passar de **80 kW**.

Por que essa prova é boa para começar:
- **Sem curvas** → só precisamos modelar a dinâmica *longitudinal* (para frente/para trás). Nada de forças laterais, nada de esterçamento, nada de transferência lateral de carga.
- **Sem freios** → só aceleração.
- **Curta** → a bateria descarrega quase nada (SoC cai ~0,6 %), então erros no modelo de bateria têm pouco impacto.
- **Objetivo único e claro**: minimizar t₇₅.

E, mesmo sendo "simples", ela já força o modelo a acertar todos os subsistemas ao mesmo tempo: se o pneu estiver errado, o carro patina de mentira; se a bateria estiver errada, o barramento colapsa; se a transmissão estiver errada, o motor não chega na rotação certa.

## 0.3 O carro simulado

| Subsistema | O que é | Números |
|---|---|---|
| **Motor** | PMSM EMRAX 228 (motor síncrono de ímãs permanentes) | 230 N·m de pico, 5500 RPM máx, 10 pares de polos |
| **Bateria** | Acumulador de lítio, arranjo 144S5P | 532,8 V nominal, 15 Ah, 32,4 kg |
| **Transmissão** | Pinhão → corrente → coroa → diferencial LSD → semi-eixos | Relação única N (~4 a 5), η ≈ 0,95 |
| **Pneus** | Hoosier FSAE slick 20.5×7.0-13 | raio 0,26 m, modelo Pacejka PAC2002 |
| **Chassi** | Monoposto FSAE com asas | 240 kg (com piloto), entre-eixos 1,5 m |
| **Inversor** | Ideal (sem PWM), controle FOC | η agregada 0,90 |

## 0.4 O que o programa faz, em uma frase

> Ele pega o estado do carro num instante (correntes no motor, rotação, velocidade, temperatura, carga da bateria), calcula que forças e torques agem sobre ele, avança tudo **0,0001 segundo** no tempo, e repete ~50.000 vezes até o carro cruzar os 75 m.

É isso. Todo o resto do documento é detalhar **quais** forças, **como** calculá-las e **por que** avançar 0,0001 s por vez.

## 0.5 Um resultado, para você saber onde vamos chegar

Configuração padrão (N = 5,0; pedal subindo suavemente até 100 % em 3 s):

```
✅ Concluído — 51036 pontos  t_final=5.104 s  dist_final=75.00 m
   SoC final=0.994   Vdc final=529.9 V
🏁 t_75m = 5.104 s  |  v_final = 104.1 km/h
```

Cinco segundos e um décimo. Com o pedal otimizado (rampa em 1 s em vez de 3 s), cai para **4,44 s**. Essa diferença de 0,66 s — obtida **sem trocar uma única peça do carro**, só mudando como o piloto pisa no acelerador — é exatamente o tipo de descoberta que justifica o projeto existir.

---

# Parte I — Fundamentos

*Se você já domina algum destes tópicos, pule. Mas leia pelo menos os blocos "🧠 Intuição" — eles fixam o vocabulário que o resto do documento usa.*

---

## 1. Física do movimento: as duas leis de Newton que usamos

Todo este simulador é, no fundo, **duas equações de Newton acopladas**. Se você entender essas duas, entendeu a espinha dorsal do projeto.

### 1.1 Movimento de translação (o carro andando para frente)

A famosa segunda lei de Newton:

```
F = m · a
```

- **F** = soma de todas as forças [N — newtons]
- **m** = massa [kg]
- **a** = aceleração [m/s²]

Reescrita do jeito que o simulador usa (aceleração é a derivada da velocidade):

```
m · dv/dt = F_tração − F_resistência
```

Ou seja: *a massa do carro vezes a taxa de variação da velocidade é igual à força que empurra menos a força que segura.*

> 🧠 **Intuição**
> Empurrar um carrinho de supermercado cheio: se você empurra com 100 N e o atrito puxa para trás com 30 N, sobram 70 N para acelerar. Se o carrinho tem 35 kg, ele ganha 2 m/s de velocidade a cada segundo. É exatamente isso.

**Nossas forças:**
- **F_tração** — o que o pneu empurra contra o chão (§16). Sempre positiva na prova de aceleração.
- **F_resistência** — arrasto do ar + resistência ao rolamento dos pneus + componente da gravidade se a pista fosse inclinada (§15).

### 1.2 Movimento de rotação (o rotor do motor girando)

Existe uma versão "girante" da mesma lei:

```
T = J · α        ou       J · dω/dt = ΣT
```

- **T** = torque [N·m — newton-metro]
- **J** = momento de inércia [kg·m²]
- **α = dω/dt** = aceleração angular [rad/s²]
- **ω** (ômega) = velocidade angular [rad/s]

| Translação | Rotação | Significado |
|---|---|---|
| força F [N] | torque T [N·m] | o que causa a mudança |
| massa m [kg] | inércia J [kg·m²] | a "preguiça" de mudar |
| velocidade v [m/s] | velocidade angular ω [rad/s] | o quão rápido |
| posição x [m] | ângulo θ [rad] | onde está |

> 🧠 **Intuição sobre momento de inércia**
> Massa é "quanto custa mudar a velocidade". Momento de inércia é "quanto custa mudar a *rotação*" — e depende não só da massa, mas de **quão longe do eixo** ela está. Um patinador no gelo gira mais rápido quando encolhe os braços: a massa é a mesma, mas J diminuiu porque a massa ficou mais perto do eixo. Uma roda de bicicleta com aro pesado é mais difícil de acelerar que uma com raio de mesma massa concentrada no centro.

**Sobre radianos:** ω em rad/s é a unidade "natural" da física. Mas todo mundo na equipe fala em RPM. A conversão:

```
RPM = ω · 60 / (2π)          ω [rad/s] = RPM · 2π / 60
```

Exemplos do nosso carro: 5500 RPM = 575,96 rad/s (é de onde vem o `speed_ref=575.96` no código).

### 1.3 Os dois corpos e o que os liga

O simulador trata **motor** e **carro** como dois corpos separados:

```
Corpo 1 (rotacional — o rotor do motor):
   J_eff · dω/dt = T_em − T_carga − k_f·ω − T_ferro

Corpo 2 (translacional — o carro inteiro):
   m · dv/dt = F_tração − F_resistência
```

E quem faz a ponte entre os dois é o **pneu**. Isso é uma escolha de modelagem importante, e vamos voltar nela várias vezes.

> 🧠 **Por que não um corpo só?**
> Um modelo mais simples assumiria que a roda **nunca patina**: nesse caso ω e v estão travados um ao outro (v = ω·r/N) e você pode escrever uma equação só. Mais simples, mas incapaz de representar wheelspin — e wheelspin é exatamente o fenômeno que decide uma largada de FSAE. Ao separar os dois corpos, deixamos que ω e v evoluam com alguma liberdade, e a diferença entre eles é o **escorregamento** (slip), que alimenta o modelo de pneu.

### 1.4 Potência

Potência é energia por unidade de tempo [W — watts]:

```
Translação:  P = F · v
Rotação:     P = T · ω
```

Isso vai ser usado o tempo todo. Por exemplo, o limite de 80 kW da FSAE: se o motor está a 300 rad/s, o torque máximo que o regulamento permite (ignorando eficiências) é 80.000/300 ≈ 267 N·m.

> ✏️ **Exercício 1.** Nosso motor tem torque de pico 230 N·m. A que rotação (em RPM) ele atinge 80 kW mecânicos? *(Resposta ao final, §30.)*

---

## 2. Derivadas, ODEs e integração numérica

Esta seção responde a pergunta: **por que o programa avança de 0,0001 s em 0,0001 s?**

### 2.1 O que é uma derivada, para o nosso uso

Derivada = **taxa de variação**. `dv/dt` é "quanto a velocidade muda por segundo". Se seu carro está a 10 m/s e `dv/dt = 4 m/s²`, daqui a um décimo de segundo ele estará (aproximadamente) a 10,4 m/s.

Nós escrevemos "d alguma coisa / dt" o tempo todo. Sempre significa a mesma coisa: **a velocidade com que aquela grandeza está mudando agora**.

### 2.2 O que é uma EDO (ODE)

Uma **Equação Diferencial Ordinária** é uma equação que te diz a *derivada* de uma grandeza em função do *valor atual* dela e de outras coisas. Exemplo do nosso simulador:

```
dv/dt = (F_tração(v, ω, ...) − F_resistência(v)) / m
```

Repare: para saber a *variação* da velocidade, preciso saber a velocidade *atual* (a resistência do ar depende de v²). Isso é o que caracteriza uma ODE — e é o que torna impossível resolver "de uma vez": você precisa avançar aos poucos.

> 🧠 **Intuição**
> Imagine dirigir de olhos vendados, mas com alguém do lado dizendo a cada instante "você está acelerando a 3 m/s² agora". Se você anotar a velocidade inicial e for somando `3 × (tempo decorrido)` em passinhos pequenos, você consegue reconstruir a viagem inteira. É exatamente o que o integrador numérico faz.

### 2.3 O vetor de estado

Nosso sistema tem **11 grandezas** que evoluem no tempo. Juntamos todas num vetor:

```
x = [ isd, isq, iso,   ← correntes do motor no referencial dq0 [A]
      ωm, θm,          ← velocidade [rad/s] e ângulo [rad] do rotor
      T,               ← temperatura do enrolamento [°C]
      v, x_pos,        ← velocidade [m/s] e posição [m] do carro
      SoC, I*, t_acc ] ← estado de carga, corrente acumulada, relógio da bateria
```

E o simulador sabe calcular `dx/dt` — as 11 derivadas — a partir de `x`. Essa função é o coração do programa e se chama `_physics_ode()`.

> 🧠 **O que é "estado"?**
> Estado é o **conjunto mínimo de números que você precisa saber agora para prever o futuro**. Se eu te der esses 11 números, você não precisa saber nada sobre o passado do carro — pode calcular todo o resto. Grandezas que *não* são estado (torque, força de tração, slip) são calculadas *a partir* do estado a cada instante.

### 2.4 Método de Euler: o jeito ingênuo

O jeito mais simples de avançar no tempo:

```
x[k+1] = x[k] + h · f(x[k])
```

onde `h` é o passo de tempo e `f` é a função que devolve as derivadas.

Ou seja: "assume que a derivada fica constante durante o passo inteiro". Funciona, mas o erro se acumula rápido — para uma trajetória curva, é como aproximá-la por uma sequência de retinhas tangentes, sempre errando para fora da curva.

### 2.5 Runge-Kutta de 4ª ordem (RK4): o que usamos

O RK4 avalia a derivada **quatro vezes** dentro de cada passo e faz uma média ponderada:

```
k1 = f( x                )       ← inclinação no início do passo
k2 = f( x + (h/2)·k1     )       ← inclinação no meio, usando k1
k3 = f( x + (h/2)·k2     )       ← inclinação no meio, usando k2 (refinada)
k4 = f( x + h·k3         )       ← inclinação no fim

x[k+1] = x[k] + (h/6)·(k1 + 2·k2 + 2·k3 + k4)
```

> 🧠 **Intuição**
> Você quer estimar quanto vai andar na próxima hora. Euler pergunta a velocidade **agora** e multiplica por uma hora. RK4 pergunta a velocidade agora, estima a velocidade daqui a meia hora (duas vezes, refinando), estima a do fim da hora, e faz uma média com peso maior para as estimativas do meio. Muito mais preciso, ao custo de 4 avaliações em vez de 1.

**Quanto mais preciso?** O erro do Euler por passo cresce com h²; o do RK4, com h⁵. Na prática: para a mesma precisão, o RK4 permite passos ~100× maiores. Como o custo é só 4×, ele sai muito mais barato.

🔍 **No código:** `Simulation._rk4()`, linhas 561-581.

```python
k1 = self._physics_ode(t,        x,             vd, vq, I_bat, telemetry)
k2 = self._physics_ode(t + dt/2, x + (dt/2)*k1, vd, vq, I_bat)
k3 = self._physics_ode(t + dt/2, x + (dt/2)*k2, vd, vq, I_bat)
k4 = self._physics_ode(t + dt,   x + dt*k3,     vd, vq, I_bat)
return x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
```

### 2.6 Por que o passo é fixo em 1e-4 s (10 kHz)?

Duas razões, e as duas importam:

**(a) Rigidez do sistema.** As correntes elétricas do motor mudam MUITO mais rápido que a velocidade do carro. A constante de tempo elétrica é τ_e = L/Rs = 9,65e-5 / 0,00706 ≈ **13,7 ms**... mas a malha de controle de corrente foi sintonizada para 500 Hz de banda, o que significa dinâmica de ~2 ms. E o back-EMF acopla as equações de forma que transientes bem mais rápidos aparecem. Um passo de 0,1 ms dá ~20 pontos por ciclo de 500 Hz — o mínimo confortável.

**(b) O controlador é digital.** Num carro real, o inversor calcula os comandos de tensão a uma taxa fixa (tipicamente 8–20 kHz) e mantém esses comandos constantes até o próximo cálculo. Isso se chama **ZOH** (*zero-order hold*, retentor de ordem zero). Nosso simulador reproduz isso: o controle roda 1× por passo, e as saídas (vd, vq) ficam **congeladas** durante as 4 avaliações do RK4.

> ⚠️ **Isso não é um detalhe.** Se os PIs rodassem *dentro* da ODE, cada uma das 4 avaliações do RK4 veria um comando diferente e o RK4 deixaria de ser matematicamente válido (ele pressupõe que `f` é uma função pura do estado). Foi um bug real da versão v1 do simulador, corrigido na v2. Veja §25.3.

**O custo:** 75 m em ~5 s a 10 kHz = ~51.000 passos × 4 avaliações = ~204.000 chamadas da ODE. Leva de 30 s a 2 min por simulação.

---

## 3. Eletricidade e magnetismo: o mínimo necessário

### 3.1 Tensão, corrente, resistência

- **Corrente (i)** [A — ampères]: fluxo de carga elétrica. Analogia hidráulica: vazão de água num cano.
- **Tensão (V)** [V — volts]: "pressão" elétrica, a diferença de potencial que empurra a corrente. Analogia: pressão da água.
- **Resistência (R)** [Ω — ohms]: o quanto o material se opõe à passagem. Analogia: cano estreito.

**Lei de Ohm:** `V = R · i`

**Potência dissipada em calor numa resistência:** `P = R · i²`

> ⚠️ Note o **quadrado**. Dobrar a corrente **quadruplica** o calor. É por isso que o modelo térmico do motor importa, e por isso que motores elétricos sofrem em aceleração máxima prolongada.

### 3.2 Indutância

Um **indutor** é uma bobina. A grandeza é a **indutância L** [H — henry], e a lei é:

```
V = L · di/dt
```

Ou seja: a tensão sobre um indutor é proporcional à *taxa de variação* da corrente, não à corrente em si.

> 🧠 **Intuição**
> Indutância é a "inércia" da corrente elétrica. Assim como você não consegue mudar instantaneamente a velocidade de um carro pesado, você não consegue mudar instantaneamente a corrente num indutor — teria que aplicar tensão infinita. Em termos hidráulicos: é como um cano longo cheio de água em movimento, que "não quer" parar nem acelerar de repente.

Isso é essencial no motor: os enrolamentos do estator **são** indutores (Ld e Lq no nosso modelo, ambos ≈ 96,5 μH). Por isso o controle de corrente precisa de um controlador — a corrente não vai instantaneamente para onde você quer.

Reescrevendo para uso na ODE (o que queremos é a derivada):

```
di/dt = V / L
```

Com L pequeno (96,5 μH), `1/L ≈ 10.363` — números grandes, dinâmica rápida. É a razão principal do passo de 0,1 ms.

### 3.3 Fluxo magnético e a lei de Faraday

- **Fluxo concatenado (λ)** [Wb — weber]: quantidade de campo magnético "abraçado" pela bobina.
- **Lei de Faraday:** um fluxo magnético **variando** induz tensão numa bobina: `V = dλ/dt`.

No nosso motor, os ímãs permanentes do rotor produzem um fluxo fixo **λm = 0,04748 Wb**. Quando o rotor gira, esse fluxo, do ponto de vista das bobinas do estator, **varia** — e isso induz uma tensão chamada **força contra-eletromotriz** (back-EMF):

```
V_backEMF = ωe · λm
```

onde ωe é a velocidade **elétrica** (§4.2).

> ⚠️ **Consequência prática crucial.** O back-EMF cresce com a rotação e se opõe à tensão que você aplica. A 5500 RPM: ωe = 575,96 × 10 = 5759,6 rad/s elétricos, logo V_backEMF = 5759,6 × 0,04748 ≈ **273 V**. Sua bateria dá 532,8 V. Quanto sobra para forçar corrente? Cada vez menos, conforme o motor acelera. Chega um ponto em que não sobra nada — a **velocidade-base**. Acima dela, é preciso a trapaça chamada *field weakening* (§23).

---

## 4. Como um motor PMSM produz torque

**PMSM** = *Permanent Magnet Synchronous Motor* — motor síncrono de ímãs permanentes.

### 4.1 A ideia básica

```
       ESTATOR (parado, tem as bobinas)
      ╔═══════════════════════════════╗
      ║   ┌───┐                ┌───┐  ║
      ║   │ A │    ROTOR       │ A'│  ║
      ║   └───┘   ┌───────┐    └───┘  ║
      ║           │ N ─ S │           ║   ← ímãs permanentes no rotor
      ║   ┌───┐   └───────┘    ┌───┐  ║
      ║   │ B │                │ B'│  ║
      ║   └───┘                └───┘  ║
      ╚═══════════════════════════════╝
```

1. O **rotor** (que gira) tem ímãs permanentes fixos nele — daí o "PM". Ele tem um campo magnético próprio, sempre presente, sem gastar energia.
2. O **estator** (parado, é a carcaça) tem três conjuntos de bobinas (fases a, b, c). Ao passar corrente por elas em sequência, cria-se um **campo magnético girante**.
3. Os ímãs do rotor tentam se alinhar com esse campo girante. Se o campo do estator gira, o rotor "persegue" — e essa perseguição é o **torque**.
4. "Síncrono" significa que o rotor gira exatamente na mesma frequência do campo do estator (dividida pelos pares de polos). Não há escorregamento como num motor de indução.

> 🧠 **Intuição**
> Imagine uma bússola (o rotor) e você segurando um ímã por fora, girando devagar em volta (o campo do estator). A agulha da bússola persegue seu ímã. Se você acelerar a rotação do seu ímã, a agulha acelera junto — enquanto o torque magnético for suficiente. Se você girar rápido demais, a agulha "perde o passo" e o torque desmorona. Um motor síncrono real também perde a sincronia se sobrecarregado — mas nosso controle FOC (§5) mede a posição do rotor e ajusta o campo em tempo real, então isso nunca acontece.

### 4.2 Pares de polos

Nosso EMRAX 228 tem **p = 10 pares de polos** (20 polos magnéticos). Isso significa que, para o rotor dar **uma volta mecânica**, o campo magnético precisa dar **10 voltas elétricas**.

```
ωe = p · ωm          θe = p · θm
```

- **ωm** = velocidade **mecânica** [rad/s] — o que você mede com um tacômetro.
- **ωe** = velocidade **elétrica** [rad/s] — a frequência das correntes nas bobinas.

Exemplo: 5500 RPM mecânicos = 575,96 rad/s mecânicos = 5759,6 rad/s elétricos = 916,7 Hz. As correntes nas fases oscilam a quase 1 kHz.

🔍 **No código:** aparece como `we = self.p * wm` em `_physics_ode`, linha 425.

### 4.3 A equação de torque

Para o nosso motor (ímãs superficiais, Ld ≈ Lq, sem torque de relutância):

```
T_em = 1,5 · p · λm · isq
```

Ou, definindo a **constante de torque** Kt:

```
Kt = 1,5 · p · λm = 1,5 × 10 × 0,04748 = 0,7122 N·m/A
T_em = Kt · isq
```

**Isso é lindo e é o ponto central do controle vetorial:** o torque é *diretamente proporcional* a uma única corrente, `isq`. Se você quer 230 N·m, você pede isq = 230/0,7122 = 323 A. É por isso que `max_current = 323.0` no código.

Mas... o que é `isq`?

---

## 5. O referencial dq (transformada de Park) — a ideia mais importante do controle

Esta é a parte que costuma travar quem chega. Vamos com calma.

### 5.1 O problema

Nas três fases do motor, as correntes são **senóides defasadas de 120°**:

```
ia(t) = I·cos(ωe·t)
ib(t) = I·cos(ωe·t − 120°)
ic(t) = I·cos(ωe·t + 120°)
```

A 5500 RPM, essas senóides oscilam a ~917 Hz. Tentar controlar diretamente três senóides de alta frequência com controladores PI é ruim: um PI só zera erro em regime permanente para sinais **constantes**. Perseguir uma senóide sempre deixa erro de amplitude e atraso de fase.

### 5.2 A solução: mudar de ponto de vista

> 🧠 **A analogia do carrossel**
> Você está no chão vendo um carrossel girar. Um cavalinho passa por você subindo e descendo, indo e voltando — do seu ponto de vista, tudo oscila. Agora **suba no carrossel** e gire junto com ele. De repente o cavalinho está parado em relação a você: está sempre "ali, na frente e à direita". A física não mudou nada; **só o seu referencial mudou** — e no novo referencial tudo virou constante.
>
> A transformada de Park faz exatamente isso: ela coloca o observador girando junto com o rotor. As três correntes senoidais viram **duas correntes contínuas**.

### 5.3 Os dois eixos

No referencial girante, definimos:

- **Eixo d** (*direct*) — alinhado com o fluxo dos ímãs do rotor. Corrente `isd` neste eixo **reforça ou enfraquece** o campo magnético, mas **não produz torque**.
- **Eixo q** (*quadrature*) — perpendicular (90° elétricos) ao eixo d. Corrente `isq` neste eixo produz **todo o torque**.

```
              q  (produz torque)
              ↑
              │
              │
    ──────────┼──────────→ d  (alinhado com os ímãs — não produz torque)
              │
```

> 🧠 **Por que 90°?**
> Torque magnético máximo acontece quando os dois campos estão perpendiculares. Se você empurra uma porta pela maçaneta, empurra **perpendicular** à porta — empurrar na direção da dobradiça não faz nada. Eixo d = direção da dobradiça. Eixo q = perpendicular, onde a força vira torque.

### 5.4 O que a estratégia FOC faz

**FOC** = *Field-Oriented Control* (controle orientado por campo). A estratégia é:

1. Meça a posição do rotor (θm) → calcule θe = p·θm.
2. Transforme as correntes medidas (ia, ib, ic) para o referencial dq usando θe.
3. Controle `isq` para dar o torque que você quer (`isq_ref = T_desejado / Kt`).
4. Controle `isd` = 0 normalmente (não desperdice corrente sem produzir torque) — exceto em field weakening, quando `isd < 0` (§23).
5. Transforme os comandos de tensão (vd, vq) de volta para as três fases e mande para o inversor.

**Resultado:** o motor AC de três fases passa a ser controlado como se fosse um motor DC de escovas — que é a coisa mais fácil do mundo de controlar. Essa é a mágica do FOC, e é por isso que ele domina a indústria de tração elétrica.

### 5.5 Convenção "invariante em amplitude" (⚠️ ler antes de mexer)

Existem duas convenções para escrever a transformada de Park, e elas diferem por um fator √(2/3). Nosso projeto usa a **invariante em amplitude**, o que implica:

| Consequência | Valor |
|---|---|
| `|i_dq|` é a **amplitude de pico** da corrente de fase | isq = 323 A ⟹ pico de fase = 323 A ⟹ I_rms = 323/√2 ≈ 228 A |
| Torque tem fator **1,5** | `T = 1,5 · p · λm · isq` |
| Perdas no cobre têm fator **1,5** | `P_cu = 1,5 · Rs · (isd² + isq²)` |
| `max_current` é **pico**, não rms | 323 A pico |

> ⚠️ **Não regredir.** Esses três 1,5 vêm todos da mesma escolha. Se alguém "corrigir" um deles isoladamente (por exemplo, achar que perdas trifásicas deveriam ser `3·Rs·i²`), o modelo fica internamente inconsistente: o torque e as perdas passariam a falar de correntes diferentes. Se você mudar a convenção, mude **os três juntos** e revalide o torque de pico contra o datasheet.

### 5.6 As equações do PMSM em dq

Com tudo isso, as equações elétricas do motor ficam:

```
Ld · d(isd)/dt = vd − Rs·isd + ωe·Lq·isq
Lq · d(isq)/dt = vq − Rs·isq − ωe·(Ld·isd + λm)
```

Leitura termo a termo da segunda equação (a que importa para torque):

| Termo | O que é |
|---|---|
| `vq` | tensão que o inversor aplica no eixo q — é o que **você comanda** |
| `−Rs·isq` | queda ôhmica no cobre do enrolamento |
| `−ωe·Ld·isd` | **acoplamento cruzado**: o eixo d "vaza" para o eixo q |
| `−ωe·λm` | **back-EMF** — cresce com a rotação e rouba tensão |

> ⚠️ Repare que as duas equações são **acopladas**: isd aparece na equação de isq e vice-versa. Isso complica o controle. A solução é o **desacoplamento por feedforward** (§20.5): calculamos os termos de acoplamento e os somamos ao comando, cancelando-os. Cada eixo vira então uma planta de 1ª ordem independente e simples.

🔍 **No código:** `_physics_ode`, linhas 438-440.

```python
d_isd = (vd - rs_eff*isd + we*self.lq*isq) * self.inv_ld
d_isq = (vq - rs_eff*isq - we*(self.ld*isd + self.lambda_m)) * self.inv_lq
d_iso = (-rs_eff * iso) / self.L0
```

(`iso` é a componente de "sequência zero" — existe por completude matemática, decai sozinha para zero e não afeta torque. Ignore.)

---

## 6. Controle: PI, saturação e windup

### 6.1 Malha aberta vs. malha fechada

- **Malha aberta**: você manda o comando e torce. "Pedal a 50 % → corrente de 161 A." Sem verificação.
- **Malha fechada**: você mede o resultado, compara com o desejado e corrige. "Quero 161 A, estou medindo 140 A, erro = 21 A, aumento a tensão."

Nosso simulador usa **os dois, em níveis diferentes**:

```
   PEDAL  ──(malha ABERTA)──▶  iq_ref  ──(malha FECHADA: PI)──▶  vq  ──▶  motor
```

O pedal define quanta corrente queremos, **sem** medir a velocidade do carro (malha aberta — assim como no carro real: o piloto é a realimentação). Mas a corrente em si é controlada em malha fechada por um PI, porque a indutância impede que ela apareça instantaneamente.

> 🧠 **Por que não um PI de velocidade?**
> Seria tentador: "quero 5500 RPM, PI de velocidade gera a corrente". Mas na largada o erro é gigante (5500 RPM de erro) durante segundos. O integrador do PI acumula um valor enorme, e quando o carro finalmente chega perto da referência, esse acúmulo tem que ser "descarregado" — gerando um overshoot violento e oscilação. Esse fenômeno se chama **windup** e é explicado abaixo. A escolha de projeto foi eliminar a causa: sem PI de velocidade, sem windup de velocidade. A proteção de sobre-rotação virou um limitador suave e proporcional (§22.2).

### 6.2 O controlador PID

```
u(t) = kp·e(t)  +  ki·∫e(t)dt  +  kd·de(t)/dt
        ────┬───    ─────┬────      ─────┬────
       Proporcional   Integral      Derivativo
```

onde `e = referência − medida` (o erro).

| Termo | O que faz | Analogia |
|---|---|---|
| **P** — proporcional | reage ao erro **atual**. Grande erro → resposta grande. Sozinho, sempre deixa erro residual. | Você acelera mais quanto mais longe está do carro da frente |
| **I** — integral | acumula o erro ao longo do tempo. Elimina o erro residual, mas adiciona atraso e pode oscilar. | "Estou atrasado há 10 minutos, vou acelerar um pouco mais para compensar" |
| **D** — derivativo | reage à **velocidade de mudança** do erro. Antecipa e amortece. Amplifica ruído. | Você tira o pé antes de chegar perto, porque vê que está fechando rápido |

**No nosso projeto, kd = 0.** Só usamos PI. Motivo: derivar um erro de corrente amplificaria ruído numérico sem benefício real na malha rápida.

### 6.3 Como escolhemos os ganhos: cancelamento de polo

Não chutamos os ganhos. Usamos um método analítico.

A planta elétrica de cada eixo é um sistema de **primeira ordem**: `L·di/dt = v − Rs·i`. Sua constante de tempo é `τe = L/Rs`.

Escolhendo:

```
kp = ωc · L
ki = ωc · Rs
```

o **zero** do controlador (que fica em `ki/kp = Rs/L = 1/τe`) cai exatamente em cima do **polo** da planta e o cancela. Resultado: a malha fechada vira um sistema de primeira ordem com banda `ωc` — comportamento previsível e sem oscilação.

Com ωc = 2π·500 Hz = 3141,6 rad/s:
```
kp = 3141,6 × 9,65e-5 = 0,303
ki = 3141,6 × 0,00706 = 22,18
```

Confira: são exatamente os valores que a Simulation imprime ao iniciar.

🔍 **No código:** `Simulation.__init__`, linhas 266-275.

> 🧠 **Por que 500 Hz?** É uma escolha de engenharia. Alta o bastante para a corrente responder muito mais rápido que a mecânica (o carro acelera em segundos; a corrente precisa se estabelecer em milissegundos). Baixa o bastante para não bater no limite do passo de simulação (10 kHz é 20× a banda — regra de bolso confortável).

### 6.4 Saturação e windup — o ponto delicado

**Saturação**: a saída do controlador não pode ser qualquer coisa. A tensão que o inversor pode aplicar é limitada pela bateria (±Vdc). Pedir 800 V numa bateria de 532 V é impossível — o comando é **cortado**.

**Windup** é o que acontece quando você satura e o integrador não sabe:

```
Situação:  o PI pede 800 V, mas só 532 V chegam ao motor.
           O erro continua grande (a corrente não sobe o esperado).
           O integrador continua somando erro... 850... 900... 1200 V...
           
Depois:    o erro finalmente zera. Mas o integrador está em 1200 V!
           Ele precisa "descarregar" somando erro NEGATIVO por um tempão.
           → overshoot violento e oscilação
```

> 🧠 **Analogia**
> Você está empurrando um carro atolado com força máxima e não sai do lugar. Você pensa "preciso empurrar mais forte, mais forte, MAIS FORTE" (integrador subindo), mas seus músculos já estão no limite (saturado). Quando o carro finalmente sai, você está com toda essa "vontade acumulada" e sai empurrando muito mais do que precisa, derrubando o carro ladeira abaixo.

**A solução: back-calculation.** Quando a saída é cortada, você devolve o excesso ao integrador:

```
integral ← integral − (saída_bruta − saída_cortada) / ki
```

Assim o estado interno fica sempre consistente com o que **de fato** chegou à planta.

### 6.5 Anti-windup em dois estágios (nosso caso específico)

Nosso projeto tem uma sutileza: existem **duas** saturações.

```
       PI  ──▶ vq_pid ──┐
                        ├──▶ vq_total ──▶ [corta a ±Vdc] ──▶ vq (aplicado)
   feedforward ──▶ dec_q┘
       (desacoplamento)
```

- **Estágio (a) — interno ao PI:** o próprio `update()` corta a saída bruta a ±limit e faz back-calculation. Está em `PIDController.update()`.
- **Estágio (b) — externo:** depois de somar o feedforward, o total é cortado a ±Vdc — um limite **dinâmico** (a tensão da bateria muda com a carga!) que o PI não conhece. O excesso desse segundo corte é devolvido chamando `back_calculate()`.

🔍 **No código:** `Simulation.simulate()`, linhas 789-798:

```python
vd_total = vd_pid + dec_d
vq_total = vq_pid + dec_q

vd = float(np.clip(vd_total, -Vlim, Vlim))
vq = float(np.clip(vq_total, -Vlim, Vlim))

if vd_total != vd:
    self.id_ctrl.back_calculate(vd_total - vd)
if vq_total != vq:
    self.iq_ctrl.back_calculate(vq_total - vq)
```

> ⚠️ **Sem o estágio (b)** o integrador acumularia comando que nunca chega ao motor — um "windup silencioso", difícil de diagnosticar porque não aparece em nenhum log óbvio. Foi por isso que `back_calculate()` existe como método público do controlador.

---

## 7. Nota sobre unidades e sinais

Antes de seguir, fixe estas convenções — 90 % dos bugs de simulação são erro de unidade ou de sinal.

| Grandeza | Unidade no código | Cuidado comum |
|---|---|---|
| Velocidade angular | **rad/s** | RPM só aparece no log e no dashboard |
| Ângulo | **rad** | — |
| Corrente do motor | **A pico** (não rms) | ver §5.5 |
| Capacidade da bateria | **Ah** no catálogo, **C** (coulomb) nas integrais | 1 Ah = 3600 C — daí os `/3600` no código |
| Carga vertical Fz | **N, positiva (magnitude)** | O ajuste PAC2002 original usa Fz negativo (ISO 8855); o código usa `abs()` e `|μx|` — equivalente |
| Slip ratio κ | **adimensional**, ∈ [−1, +1] | positivo = patinando; negativo = travando |
| Temperatura | **°C** | as leis usam ΔT, então °C e K dão o mesmo resultado |
| Eficiência η | **fração** (0,95), não % | — |

---

# Parte II — Arquitetura do software

## 8. O mapa do território

```
┌──────────────────────────────────────────────────────────────────────┐
│  CAMADA DE APLICAÇÃO                                                 │
│                                                                      │
│   main.py                          Otimizador/otimizador.py          │
│   (dashboard web Dash)             (Evolução Diferencial)            │
│        │                                   │                         │
│        │  instancia os modelos             │  instancia os modelos   │
│        │  e chama simulate()               │  ~620 vezes             │
│        └───────────────┬───────────────────┘                         │
├────────────────────────▼─────────────────────────────────────────────┤
│  CAMADA DE SIMULAÇÃO                                                 │
│                                                                      │
│   Simulation/Simulation.py                                           │
│   • guarda TODO o estado dinâmico (vetor de 11 variáveis)            │
│   • roda o loop RK4 a 10 kHz                                         │
│   • executa o controle FOC                                           │
│        │                                                             │
│        │  lê parâmetros de                                           │
├────────▼─────────────────────────────────────────────────────────────┤
│  CAMADA DE MODELOS  (Models/)                                        │
│                                                                      │
│   Motor.py   BatteryPack.py   Transmission.py   Vehicle.py           │
│   Tire.py    Pedal.py         PIDController.py                       │
│                                                                      │
│   → contêineres de PARÂMETROS + funções puras. Sem estado dinâmico.  │
├──────────────────────────────────────────────────────────────────────┤
│  Constants/constants.py — constantes matemáticas pré-computadas      │
└──────────────────────────────────────────────────────────────────────┘
```

## 9. O princípio de design central

> **Os modelos em `Models/` são contêineres de parâmetros sem estado dinâmico próprio. TODO o estado vive na `Simulation`, que o avança no tempo.**

Isso parece uma escolha arbitrária de organização, mas **não é** — é uma exigência matemática.

### 9.1 Por que os modelos não podem ter estado

Lembre-se: o RK4 chama a função de derivadas **4 vezes por passo**, com valores de estado **diferentes e intermediários** (`x`, `x + h/2·k1`, `x + h/2·k2`, `x + h·k3`). Três desses quatro pontos são "de mentira" — são sondagens matemáticas, não estados reais pelos quais o sistema passou.

Se o modelo de bateria guardasse `self.soc` internamente e o atualizasse a cada chamada, o RK4 estaria corrompendo o estado com valores intermediários. O resultado seria numericamente errado de um jeito silencioso e dificílimo de diagnosticar.

**A solução:** a ODE é uma **função pura**.

> 🧠 **Função pura** = mesma entrada sempre dá a mesma saída, e ela não modifica nada fora dela mesma. `_physics_ode(t, x, vd, vq, I_bat)` só lê `x` e devolve `dx/dt`. Não escreve em lugar nenhum.

Compare as duas assinaturas do BatteryPack:

```python
# ❌ Como NÃO é (teria estado interno):
def calcular_tensao(self):
    return f(self.soc, self.Iast)      # lê estado interno — impuro

# ✅ Como É:
def calcular_tensao(self, corrente, soc, Iast, tempo_acumulado):
    return f(corrente, soc, Iast, ...)  # recebe tudo por parâmetro — puro
```

O `BatteryPack` **tem** um atributo `self.soc`, mas ele só é usado como **condição inicial** — a Simulation o lê uma vez em `_reset_ic()` e daí em diante o SoC vive no vetor de estado.

### 9.2 Onde ficam os controladores

Os PIs **têm** estado (o integrador!). Por isso eles rodam **fora** da ODE, uma vez por passo, e suas saídas ficam congeladas. Isso mantém a ODE pura e ainda por cima é fisicamente mais realista (é o ZOH de um controlador digital real — §2.6).

```
     ┌─── 1 vez por passo ──────────────┐   ┌──── 4 vezes por passo ────┐
     │  controle (tem estado: PIs)      │   │  ODE (pura, sem estado)   │
     │  → produz vd, vq, I_bat          │──▶│  → produz dx/dt           │
     └──────────────────────────────────┘   └───────────────────────────┘
```

## 10. O ciclo de vida de uma simulação

```
1. O chamador (dashboard ou otimizador) instancia os 6 modelos
   com os parâmetros do carro.
        │
        ▼
2. Simulation.__init__(...)
   • lê os parâmetros dos modelos com getattr
   • monta J_eff (a inércia efetiva) — §19.2
   • sintoniza os PIs de corrente por cancelamento de polo
   • pré-aloca os arrays de log
        │
        ▼
3. sim.simulate()
   ┌──── repete ~51.000 vezes ────────────────────────┐
   │  1. desempacota x[k]                             │
   │  2. calcula I_bat (balanço de potência)          │
   │  3. calcula Vdc da bateria                       │
   │  4. monta iq_ref (pedal + 3 limitadores)         │
   │     e id_ref (field weakening)                   │
   │  5. roda os PIs → vd, vq                         │
   │  6. RK4: x[k+1] = f(x[k], vd, vq, I_bat)         │
   │  7. loga x[k] e a telemetria                     │
   │  8. para se posição ≥ 75 m                       │
   └──────────────────────────────────────────────────┘
        │
        ▼
4. Retorna:
   • um dict de 28 séries temporais (chaves públicas — 't', 'isq',
     'vehicle_velocity', ...) → para o otimizador e scripts externos
   • os atributos sim.* (nomes internos) → para o dashboard
```

> ⚠️ **Pegadinha de nomenclatura.** Os nomes internos e as chaves públicas **não são iguais**: `self.corrented` ↔ `'isd'`, `self.correnteq` ↔ `'isq'`, `self.conjugado` ↔ `'torque_electromagnetic'`. Isso é herança histórica (o código nasceu em português). Regra prática: **se você está escrevendo código fora da Simulation, use o dict.** A lista completa de campos internos está na tupla `_LOG_FIELDS` (linhas 324-338) e a de chaves públicas no `return` do `simulate()` (linhas 898-927).

## 11. Quem chama quem: o diagrama de dependências

```
Simulation
   ├── Motor            → lê rs, ld, lq, jm, kf, lambda_m, p,
   │                      max_current, speed_ref, coef. térmicos,
   │                      coef. de perdas no ferro
   │                    → chama inverse_park_transform() e
   │                      abc_currents_from_dq() (só para LOG)
   │
   ├── Vehicle          → chama calculate_resistance_forces(v)
   │                      calculate_load_transfer(a, v)
   │                      calculate_reflected_inertia(trans)
   │
   ├── Transmission     → chama motor_to_wheel_torque()
   │                            motor_to_wheel_speed()
   │                    → lê final_drive_ratio, efficiency, inércias
   │
   ├── Tire             → chama SlipRatio() e Tire_forces()
   │
   ├── BatteryPack      → chama calcular_tensao()
   │                            calcular_derivadas()
   │                            calcular_tensao_nominal()
   │
   ├── Pedal            → chama update(t)
   │
   └── PIDController    → 2 instâncias (eixo d e eixo q)
                        → chama update(), back_calculate(), reset()
```

Repare que **os modelos nunca chamam uns aos outros**. Tudo passa pela Simulation. Isso torna cada modelo testável isoladamente — e é o que permite as células de demonstração da documentação original rodarem cada modelo sozinho.

Duas exceções aparentes, que na verdade confirmam a regra:
- `Vehicle.calculate_reflected_inertia(transmission)` **recebe** a transmissão como argumento — não a guarda, só lê parâmetros dela naquele instante.
- `Vehicle.calculate_load_torque(v, transmission)` idem — e essa função sequer é usada pela Simulation (é conveniência para análises isoladas).

---

# Parte III — Os módulos, um por um

*Cada seção segue a mesma estrutura: **o que é fisicamente** → **a matemática** → **como está no código** → **os parâmetros** → **armadilhas**.*

---

## 12. `Motor.py` — o motor PMSM (EMRAX 228)

### 12.1 O que este módulo faz (e o que NÃO faz)

**Faz:** guarda os parâmetros elétricos, mecânicos e térmicos do motor, e fornece duas funções de transformação de coordenadas.

**Não faz:** nenhuma dinâmica. As equações do PMSM, a integração, o controle — tudo isso está na `Simulation`. Este módulo é uma **ficha técnica executável**.

> 🧠 Se você abrir `Motor.py` esperando encontrar "o motor girando", vai se decepcionar. Pense nele como o datasheet do EMRAX transcrito para Python, com duas calculadoras de brinde.

### 12.2 O motor real

O EMRAX 228 é um motor axial de ímãs permanentes muito usado em Formula Student. Características que entram no modelo:

| Parâmetro | Símbolo | Valor no projeto | O que significa |
|---|---|---|---|
| Resistência de estator | Rs | 0,00706 Ω | resistência de UMA fase, medida a 25 °C |
| Indutância eixo d | Ld | 96,5 μH | "inércia" da corrente no eixo d |
| Indutância eixo q | Lq | 96,5 μH | idem, eixo q. **Ld ≈ Lq** ⟹ sem torque de relutância |
| Inércia do rotor | Jm | 0,02521 kg·m² | quão difícil é acelerar o rotor |
| Atrito viscoso | kf | 0,005 N·m·s | rolamentos + ventilação |
| Fluxo dos ímãs | λm | 0,04748 Wb | define Kt e o back-EMF |
| Pares de polos | p | 10 | 20 polos magnéticos |
| Corrente máx. | max_current | 323 A (pico) | ⟹ 230 N·m de pico (datasheet) |
| Rotação máx. | speed_ref | 575,96 rad/s | = 5500 RPM, limite mecânico |

### 12.3 Ld ≈ Lq: por que isso simplifica tudo

Em motores com ímãs **enterrados** no rotor (IPM), Ld ≠ Lq, e surge um torque extra chamado **torque de relutância**:

```
T = 1,5·p·[ λm·isq + (Ld − Lq)·isd·isq ]
                     └──────┬───────┘
                    torque de relutância
```

O EMRAX 228 tem ímãs **superficiais** (SPM): Ld ≈ Lq, o segundo termo some, e ficamos com a fórmula limpa `T = Kt·isq`. Menos física, mas fielmente representa este motor.

> ⚠️ Se um dia a equipe trocar por um motor IPM, **essa simplificação quebra** — e o ponto de operação de máximo torque por ampère (MTPA) deixaria de ser `isd = 0`. Seria necessário rever a estratégia de referência de corrente.

### 12.4 Efeitos secundários modelados

O modelo não para na equação ideal. Três efeitos de segunda ordem entram:

**(a) Perdas no ferro (Steinmetz)**

O campo magnético alternado no núcleo de ferro do estator dissipa energia por dois mecanismos:

```
P_ferro = k_h·|f_e|  +  k_e·f_e²
          ────┬────     ───┬────
         histerese    Foucault
                    (correntes parasitas)
```

- **Histerese**: o ferro "resiste" a ser magnetizado e desmagnetizado. Perde energia por ciclo → proporcional à frequência.
- **Foucault (eddy currents)**: correntes induzidas no próprio ferro. Proporcional ao quadrado da frequência.

Calibração: `k_h = 0.9 W/Hz`, `k_e = 9e-4 W/Hz²`. A 4000 RPM (f_e = 666,7 Hz): P_ferro = 0,9×666,7 + 9e-4×666,7² ≈ 600 + 400 = **1000 W**. Bate com a divisão 60/40 histerese/Foucault do datasheet.

**(b) Rs varia com a temperatura**

O cobre fica mais resistivo quando esquenta, ~0,393 % por kelvin:

```
Rs(T) = Rs0 · [ 1 + α·(T − T_ref) ]        α = 0,00393 /K,  T_ref = 25 °C
```

Consequência em cascata: mais temperatura → mais Rs → mais perdas → mais temperatura. E também: mais Rs muda a planta que o PI de corrente está controlando (embora a sintonia continue estável, com margem sobrando).

> 🔍 **Detalhe importante:** o código guarda **dois** atributos: `rs0` (o valor frio, imutável, referência da lei) e `rs` (o efetivo). Se você precisar da resistência nominal para alguma conta, use `rs0`, nunca `rs`.

**(c) Modelo térmico de massa concentrada**

```
dT/dt = ( P_cobre + P_ferro − P_resfriamento ) / (m·C)

P_cobre        = 1,5 · Rs(T) · (isd² + isq²)
P_resfriamento = h_cool · A_cool · (T − T_ambiente)
```

`m·C = 13,5 kg × 385 J/(kg·K) = 5197 J/K` — ou seja, 5197 joules para subir 1 °C.

> 🧠 **"Massa concentrada" (lumped)** significa que tratamos motor inteiro como **um único bloco à mesma temperatura**. Na realidade o enrolamento esquenta muito mais rápido que a carcaça. Nosso modelo subestima a temperatura de pico do cobre e superestima a da carcaça. Para uma prova de 5 s isso é aceitável; para o Endurance de 22 km, não seria (§29).

**Defaults de resfriamento:** `h_cool = 10 W/(m²·K)` (convecção natural), `A_cool = 0,15 m²`. Referência: ar forçado ≈ 50–200; líquido forçado ≈ 500–2000.

### 12.5 As duas funções de transformação

Ambas são **exclusivamente para logging** — elas não realimentam a dinâmica, que roda inteiramente em dq.

**`inverse_park_transform(vd, vq, theta_e)`** — dq → abc (tensões de fase que o inversor aplicaria):

```python
valpha = vd * cos(θe) - vq * sin(θe)      # Park inversa: dq → αβ
vbeta  = vd * sin(θe) + vq * cos(θe)
vs1 = valpha                               # Clarke inversa: αβ → abc
vs2 = -0.5*valpha + (√3/2)*vbeta
vs3 = -0.5*valpha - (√3/2)*vbeta
```

**`abc_currents_from_dq(isd, isq, theta_e, flux_d)`** — dq → abc (correntes de fase), forma direta:

```python
is1 = isd·cos(θe)        − isq·sin(θe)
is2 = isd·cos(θe − 2π/3) − isq·sin(θe − 2π/3)
is3 = isd·cos(θe + 2π/3) − isq·sin(θe + 2π/3)
```

> 🧠 **O que você vê no dashboard** quando abre a aba "Correntes": três senóides de amplitude igual a `|i_dq|`, defasadas de 120°, com frequência crescendo conforme o motor acelera. É a reconstrução do que realmente circularia nos cabos.

### 12.6 Armadilhas

| Armadilha | Consequência |
|---|---|
| Aumentar `kf` "para representar perdas" | As perdas no ferro **já** são modeladas separadamente. `kf = 0,1` dissiparia ~33 kW e dupla-contaria. Calibrado em **0,005**. |
| Confundir `max_current` com corrente rms | É **pico**. 323 A pico = 228 A rms. Se você colocar o valor rms do datasheet aqui, perde √2 de torque. |
| Mudar λm sem revisar max_current | Kt = 1,5·p·λm. Mudar λm muda o torque de pico. Sempre revalide: `Kt × max_current` deve dar o torque do datasheet. |
| Usar `rs` onde deveria usar `rs0` | Você pega a resistência quente do último passo, não a nominal. |

> ✏️ **Exercício 2.** Instancie o motor no Python e calcule: (a) o torque de pico; (b) a velocidade-base (onde ωe·λm = 0,9·532,8 V), em RPM. Compare com o que a documentação original imprime.

---

## 13. `BatteryPack.py` — a bateria (modelo de Shepherd modificado)

### 13.1 O que uma bateria "faz" do ponto de vista do simulador

Duas coisas, e só:

1. **Fornece uma tensão** que cai conforme você puxa corrente e conforme ela descarrega. Essa tensão é o **teto do que o inversor pode aplicar** no motor.
2. **Descarrega**: o SoC (estado de carga) cai proporcionalmente à corrente integrada no tempo.

### 13.2 O circuito equivalente

```
        ┌──── R_interna ────┐
   ●────┤                   ├──── + terminal
        │  E(SoC)  (fonte)  │
   ●────┴───────────────────┴──── − terminal
```

Uma **fonte de tensão que depende do estado de carga**, em série com uma **resistência interna**. Quando você puxa corrente i, perde `R·i` de tensão nos terminais — o famoso **sag** (afundamento) de tensão.

### 13.3 A equação de Shepherd modificada

Ela captura três efeitos que a bateria real apresenta ao descarregar:

```
Vt = E0 − R·i − K·(Q/(Q−q))·(I*/3600) − K·(Q/(Q−q))·q + A·e^(−B·q)
     └┬─┘ └┬─┘  └──────────┬─────────┘  └───────┬─────┘  └───┬───┘
   tensão  queda      polarização           polarização   região
   nominal ôhmica     DINÂMICA               ESTÁTICA   exponencial
                   (histórico de i)         (quanto já
                                             descarregou)
```

onde `q = (1 − SoC)·Q` é a carga já retirada da célula [Ah], e `I*` é a corrente acumulada ∫i·dt [C].

**Traduzindo cada termo:**

| Termo | Efeito físico | Quando domina |
|---|---|---|
| `E0` | tensão de circuito aberto nominal | sempre (é a base) |
| `−R·i` | queda ôhmica instantânea | **sempre que há corrente** — é o sag que você vê no gráfico |
| `−K·(Q/(Q−q))·q` | polarização estática: fica mais difícil tirar carga conforme a célula esvazia | no fim da descarga (a "queda de joelho" do gráfico) |
| `−K·(Q/(Q−q))·(I*/3600)` | polarização dinâmica: histórico de corrente também penaliza | descargas longas e pesadas |
| `+A·e^(−B·q)` | região exponencial: o "degrau" característico logo no início | primeiros % de descarga |

> 🧠 **Por que o denominador `Q − q`?** Conforme q → Q (célula vazia), esse termo explode. É a modelagem matemática do "joelho" — o momento em que a tensão despenca. O código protege isso com `denom1 = max(self.Q - q_cell, 1e-6)`.

### 13.4 ⚠️ O invariante crítico: escala pack ↔ célula

**Este é o ponto mais importante do módulo. Leia com atenção.**

O nosso pack é **144S5P**: 144 células em série, 5 fileiras em paralelo. Isso escala assim:

```
Tensão      →  × n_série     (células em série SOMAM tensão)
Corrente    →  ÷ n_paralelo  (fileiras DIVIDEM a corrente)
Capacidade  →  × n_paralelo  [Ah]
```

A equação de Shepherd é escrita **para uma célula**. Mas a Simulation entrega os valores em **nível de pack**. Então o código converte antes de aplicar a equação:

```python
i_cell    = corrente / self.n_paralelo      # 200 A no pack → 40 A por célula
q_cell    = (1 - soc) * self.Q
Iast_cell = Iast / self.n_paralelo
# ... aplica Shepherd com valores de CÉLULA ...
return Vt * self.n_serie                    # e escala o resultado de volta
```

> ⚠️ **NÃO REGREDIR.** Sem a divisão por `n_paralelo`, a corrente do pack inteiro cairia sobre a resistência de **uma** célula. Com n_paralelo = 5, o sag de tensão sairia **5× superestimado** e o barramento DC colapsaria para o piso. Este foi um bug real, corrigido, e está documentado no cabeçalho do módulo.

**Verificação numérica rápida:** pack a 200 A, SoC = 1,0.
- Errado: 200 A × 0,02 Ω = 4 V de queda por célula → 4 × 144 = **576 V de sag** (!)
- Certo: 40 A × 0,02 Ω = 0,8 V por célula → 0,8 × 144 = **115 V de sag** ✓

O gráfico da documentação original confirma: a 200 A, o pack sai de 532,8 V nominal para ~430 V. Coerente.

### 13.5 O catálogo de células

```python
'Li-ion':  E0=3.7 V, K=0.005, Q=3.0 Ah, A=0.1,  B=0.01,  R=0.020 Ω, 45 g, 25 mL
'LiFePO4': E0=3.2 V, K=0.003, Q=2.8 Ah, A=0.05, B=0.015, R=0.008 Ω, 70 g, 30 mL
```

Para adicionar a célula real quando ela for homologada, basta acrescentar uma entrada ao dicionário em `definir_celula()`.

### 13.6 O pack 144S5P do projeto

| Grandeza | Cálculo | Resultado |
|---|---|---|
| Tensão nominal | 3,7 × 144 | **532,8 V** (teto FSAE: 600 V ✓) |
| Capacidade | 3,0 × 5 | **15 Ah** |
| Energia | 532,8 × 15 | ≈ 8 kWh |
| Massa | 0,045 × 144 × 5 | **32,4 kg** |
| Nº de células | 144 × 5 | 720 |

> 🧠 **Por que 144 em série?** Porque 144 × 3,7 = 532,8 V, que está confortavelmente abaixo do teto de 600 V da FSAE, e é tensão alta o suficiente para reduzir a corrente (menos perdas nos cabos, menos calor). Regra geral: em EV, tensão alta = corrente baixa = menos perdas ôhmicas.

### 13.7 Os três estados integrados

```python
def calcular_derivadas(self, corrente):
    dsoc_dt  = -corrente * self.inv_carga_total   # SoC cai proporcional à corrente
    dIast_dt = corrente                            # integral de corrente [C]
    dtime_dt = 1.0                                 # relógio interno
    return dsoc_dt, dIast_dt, dtime_dt
```

`inv_carga_total = 1/(Q_total × 3600)` é pré-computado no `__init__` — converte corrente [A] diretamente em dSoC/dt sem divisão no loop de 10 kHz.

> ⚠️ **`corrente` aqui DEVE ser a corrente DC do barramento**, calculada por balanço de potência — **nunca** `isq`. Este é o "bug raiz" da §25.1.

### 13.8 Armadilhas

| Armadilha | Consequência |
|---|---|
| Passar `isq` como corrente | Colapso do barramento DC. Ver §25.1. |
| Esquecer a conversão pack→célula ao editar | Sag n_paralelo× superestimado. |
| Achar que o SoC importa muito nesta prova | Em 5 s, o SoC cai de 1,000 para 0,994. Praticamente nada. O que **importa** é o **sag de tensão**, que limita a margem do inversor. |
| Confundir Ah e C nas integrais | 1 Ah = 3600 C. Daí os `/3600` espalhados pela equação. |

> ✏️ **Exercício 3.** Calcule à mão a tensão terminal do pack com SoC = 1,0, I = 150 A, I* = 0. Depois rode `bat.calcular_tensao(150, 1.0, 0.0, 0.0)` e compare. Dica: com SoC = 1, q = 0, então vários termos zeram.

---

## 14. `Transmission.py` — o trem de força

### 14.1 A cadeia física

```
Motor ──▶ Pinhão (Z₁) ──corrente──▶ Coroa (Z₂) ──▶ Diferencial LSD ──▶ Semi-eixos ──▶ Rodas
        (gira rápido)                          (gira devagar)
```

Relação única (não há caixa de marchas):

```
N = Z_coroa / Z_pinhão
```

Exemplo do projeto: Z₁ = 13, Z₂ = 52 → **N = 4,0**.

### 14.2 As quatro conversões

| Grandeza | Motor → Roda | Roda → Motor |
|---|---|---|
| **Torque** | `T_roda = T_motor · N · η` | `T_motor = T_roda / (N·η)` |
| **Velocidade** | `ω_roda = ω_motor / N` | `ω_motor = ω_roda · N` |

> 🧠 **A regra de ouro das transmissões:** você **multiplica torque** e **divide velocidade** pelo mesmo fator. Potência se conserva (menos as perdas). É a alavanca: mais força, menos deslocamento.

**Sobre a eficiência η:** ela **sempre reduz o torque útil**, não importa o sentido. Por isso multiplica num sentido e divide no outro. Nosso η = η_corrente × η_diferencial = 0,98 × 0,97 = **0,9506**.

### 14.3 As inércias, separadas por lado

Este é o detalhe sutil do módulo. Cada peça gira a uma velocidade diferente, e isso muda como ela é "sentida" pelo motor.

```
   LADO RÁPIDO (velocidade do motor)      LADO LENTO (velocidade da roda)
   ────────────────────────────────      ────────────────────────────────
   • pinhão (sprocket_inertia)            • coroa (coroa_inertia)
                                          • diferencial (diff_inertia)
     ⟹ soma DIRETO em J_eff               • semi-eixos (axle_inertia)
                                          • rodas
                                          
                                            ⟹ refletidas por 1/N²
```

**Por que 1/N²?** Por conservação de energia cinética:

```
½·J_lento·ω_roda²  =  ½·J_equivalente·ω_motor²

com ω_roda = ω_motor/N:

½·J_lento·(ω_motor/N)² = ½·J_equiv·ω_motor²
⟹  J_equiv = J_lento / N²
```

> 🧠 **Intuição.** Com N = 4, a roda gira 4× mais devagar que o motor. Uma peça girando na roda "custa" muito menos ao motor, porque para o motor girar 4 voltas ela só gira 1. E como energia cinética vai com ω², a redução é 4² = **16×**.

### 14.4 O LSD

`lsd_bias_ratio = 3.0` documenta o diferencial autoblocante. **Não afeta a simulação atual** — em linha reta, com tração simétrica nas duas rodas, o LSD se comporta como um diferencial aberto. O parâmetro está lá registrado para quando o modelo evoluir para curvas.

### 14.5 Prioridade de definição da relação

```python
if Z_pinhao is not None and Z_coroa is not None and int(Z_pinhao) > 0:
    self.final_drive_ratio = float(Z_coroa) / float(Z_pinhao)   # 1ª opção
else:
    self.final_drive_ratio = float(final_drive_ratio) or 5.0     # 2ª opção
```

Par de dentes tem **precedência** sobre valor manual, porque dentes são inteiros e representam hardware que existe. É esse o formato que o otimizador entrega ao final (§27).

### 14.6 Armadilhas

| Armadilha | Consequência |
|---|---|
| Somar a inércia da coroa direto em J_eff | Superestima J_eff em N² vezes — o carro fica "pesado" de mentira |
| Passar tanto `efficiency` quanto `chain_efficiency` | O `efficiency` explícito **sobrescreve** o produto. Comportamento intencional (retrocompatibilidade), mas confunde. |
| Z_pinhao = 0 | Protegido pelo teste `int(Z_pinhao) > 0` (evita divisão por zero vinda do dashboard) |

> ✏️ **Exercício 4.** Com N = 4,0, r = 0,26 m e o motor a 5500 RPM, qual a velocidade do carro em km/h assumindo que não há escorregamento? (Confira com o valor impresso na documentação original: 135 km/h.) Agora repita com N = 5,0. O que isso te diz sobre o compromisso entre aceleração e velocidade máxima?

---

## 15. `Vehicle.py` — dinâmica longitudinal do carro

Este módulo responde **três perguntas** para a Simulation, a cada passo.

### 15.1 Pergunta 1: quais forças resistem ao movimento?

```
F_resist = ½·ρ·Cd·A·v²  +  Cr(v)·m·g·cos(θ)  +  m·g·sin(θ)
           └─────┬─────┘    └──────┬───────┘    └────┬────┘
              arrasto           rolamento          rampa
```

**Arrasto aerodinâmico.** Proporcional a **v²** — dobra a velocidade, quadruplica o arrasto.
- ρ = 1,225 kg/m³ (ar ao nível do mar, 15 °C)
- Cd = 0,7789, A = 0,68 m² → Cd·A = 0,53 m². Alto para um carro, normal para um FSAE com asas (asa gera downforce *e* arrasto).

**Resistência ao rolamento.** O pneu se deforma ao rolar e dissipa energia. Modelada com um coeficiente que cresce com a velocidade:

```
Cr(v) = Cr0 · (1 + v/v_ref)        Cr0 = 0,015,  v_ref = 150 m/s
```

Efeito suave: +6,7 % a 10 m/s, +13 % a 20 m/s, +20 % a 30 m/s.

**Rampa.** `m·g·sin(θ)`. Com `road_grade = 0` (pista plana), some. Está lá para simular subidas.

**Ordem de grandeza:** a 104 km/h (28,9 m/s, velocidade final da prova), o arrasto vale ½×1,225×0,7789×0,68×28,9² ≈ **271 N** e o rolamento ≈ 42 N. Contra uma tração da ordem de 2000 N. Ou seja: **na aceleração, a resistência é secundária** — o que limita é o pneu e a potência.

### 15.2 Pergunta 2: quanta carga vertical há na roda motriz?

Esta é a pergunta que alimenta o pneu. Três contribuições:

```
Fz_por_roda = ½ · [ m·g·d_cg/L  +  m·a·h/L  +  ½ρ·Cl_tras·A_tras·v² ]
                    └────┬────┘    └───┬──┘    └────────┬──────────┘
                      estático      transferência    downforce
                                     de carga         traseiro
```

**(a) Estático.** Momento em torno do eixo dianteiro: `m·g·d_cg/L`. Com d_cg = 0,6 m e L = 1,5 m → 40 % do peso no eixo traseiro. (240 kg × 9,81 × 0,4 = 942 N no eixo, 471 N por roda.)

**(b) Transferência longitudinal.** `m·a·h/L`.

> 🧠 **Por que acelerar "senta" o carro na traseira?** O motor empurra o carro pelo contato pneu-solo (na altura do chão), mas a massa está concentrada no CG, a h = 0,28 m de altura. Isso cria um **binário** que tende a girar o carro para trás — descarrega o eixo dianteiro e carrega o traseiro. É por isso que motos empinam e carros de arrancada levantam a frente.

Com a = 10 m/s²: 240 × 10 × 0,28/1,5 = **448 N** transferidos para o eixo traseiro. Comparado aos 942 N estáticos, é quase 50 % a mais de carga na roda motriz.

**(c) Downforce traseiro.** `½ρ·Cl_tras·A_tras·v²`. Cresce com v², então é irrelevante na largada e importante no fim. A 28,9 m/s: ½×1,225×1,0×0,35×28,9² ≈ **179 N** no eixo traseiro.

> ⚠️ **Só a asa TRASEIRA entra.** O carro é RWD (tração traseira). O downforce da asa dianteira carrega o eixo dianteiro, que não tem tração — não ajuda em nada na aceleração (e ainda gera arrasto). Isso está explícito no código: `calculate_load_transfer` usa apenas `lift_coeff_rear` e `area_rear`.

### 15.3 O acoplamento realimentado (e como resolvemos)

Repare no problema:

```
   Fz depende da aceleração  a   (transferência de carga)
                ↕
   a depende de Fz   (porque Fz limita a tração via pneu)
```

Isso é um **loop algébrico**: as duas grandezas dependem uma da outra *no mesmo instante*. Não dá para calcular uma sem a outra.

**Solução adotada: 2 iterações de ponto fixo.**

```python
# Iteração 0: chuta a = 0
Fz0 = self.vehicle.calculate_load_transfer(0.0, vv)
Fx0 = self.tire.Tire_forces(Fz0, slip) * n_drv
Ft0 = np.sign(Ftração_ideal) * min(abs(Ftração_ideal), abs(Fx0))
a0  = (Ft0 - Fresist) / m

# Iteração 1: refina Fz com a aceleração estimada
Fz     = self.vehicle.calculate_load_transfer(a0, vv)
Fx_max = self.tire.Tire_forces(Fz, slip) * n_drv
Ft     = np.sign(Ftração_ideal) * min(abs(Ftração_ideal), abs(Fx_max))
```

> 🧠 **Por que 2 iterações bastam?** Porque a correção converge rápido: a segunda iteração já traz a estimativa para dentro de ~1 % do ponto fixo. Uma terceira iteração custaria 33 % mais tempo de CPU para ganhar quase nada — e lembre que isso roda 204.000 vezes por simulação. É um compromisso deliberado.

🔍 **No código:** `Simulation._physics_ode`, linhas 477-486.

### 15.4 Pergunta 3: qual a inércia refletida ao motor?

```
J_refletida = J_pinhão + (J_rodas + J_eixos + J_diff + J_coroa)/N² + m·(r/N)²
              └───┬───┘   └───────────────┬─────────────────────┘   └───┬───┘
              lado motor              lado lento / N²            massa translacional
```

**Inércia das rodas.** `J_roda = k·m_roda·r²` com **k = 0,80**.

> 🧠 **O fator de forma k.** Para um anel fino (toda massa no aro), k = 1,0. Para um disco maciço uniforme, k = 0,5. Uma roda de FSAE (aro + pneu, com raio interno moderado) fica em k ≈ 0,8. Todas as **4** rodas entram, porque todas giram, mesmo as não-motrizes.

**O termo translacional `m·(r/N)²`.**

```
Energia cinética do carro:  ½·m·v²
Com v = ω_motor·r/N:        ½·m·(r/N)²·ω_motor²
⟹ inércia equivalente:      m·(r/N)²
```

> ⚠️ **Erro clássico:** escrever `m·r²` e depois dividir tudo por N². Isso dá o mesmo resultado *por acidente* se você dividir o total, mas a função devolve valor errado se chamada isoladamente. Foi corrigido na v2 do módulo (está no histórico do cabeçalho).

### 15.5 ⚠️ O modelo de dois corpos e a dupla contagem

Aqui está a armadilha mais fina do projeto inteiro.

`calculate_reflected_inertia()` devolve o **total**, incluindo `m·(r/N)²`. Isso é correto para um **modelo de corpo único** (acoplamento rígido, sem escorregamento).

Mas a nossa `Simulation` usa um **modelo de dois corpos**:

```
J_eff · dω/dt = T_em − T_carga − ...      ← motor (rotacional)
m     · dv/dt = F_tração − F_resist       ← carro (translacional)
```

Neste arranjo, a massa do carro **já aparece** na segunda equação. Se ela também entrasse em J_eff, estaríamos contando a massa **duas vezes**. Então a Simulation subtrai:

```python
j_all         = vehicle.calculate_reflected_inertia(transmission)
j_translation = vehicle.mass * (r / N) ** 2
self.j_eff    = max(self.jm + j_all - j_translation, self.jm)
```

**Os números do projeto — e uma lição de leitura de documentação:**

| Parcela | com N = 4,0 | com N = 5,0 |
|---|---|---|
| J das 4 rodas (0,8·5·0,26² × 4) | 1,0816 | 1,0816 |
| + semi-eixos + diff + coroa | 1,0981 | 1,0981 |
| → refletido: J_lento / N² | 0,0686 | 0,0439 |
| + pinhão (lado motor) | 0,0688 | 0,0441 |
| Massa translacional m·(r/N)² | 1,0140 | 0,6490 |
| **Total** (`calculate_reflected_inertia`) | **1,0828** | 0,6931 |
| **J_eff = Jm + rotacional** | 0,0940 | **0,0693** |

> 🧠 A documentação original mostra `J refletida total = 1,0828` (exemplo com N = 4,0, da seção de transmissão) e `J_eff = 0,0693` (saída da simulação, que usa N = 5,0). **São configurações diferentes.** Ao comparar números entre seções de qualquer documentação de simulação, sempre confira **qual configuração estava em uso**. Este tipo de confusão já custou horas de depuração a mais de uma equipe.

O módulo também deixa as parcelas disponíveis separadamente, nos atributos `j_rot_reflected` e `j_translational`, para quem preferir montar J_eff pelo caminho direto.

### 15.6 Armadilhas

| Armadilha | Consequência |
|---|---|
| Usar o retorno completo de `calculate_reflected_inertia` como J_eff | Dupla contagem da massa: o carro fica ~15× mais "inerte" no eixo do motor |
| Passar downforce dianteiro para a roda motriz | Tração fictícia — o carro acelera mais do que poderia |
| Esquecer que Fz é **por roda** e Fx é **de um pneu** | Fator 2 no teto de tração |
| Não passar `L`, `h`, `dist_cg` | O código cai num fallback de 50 % do peso por eixo, **sem transferência de carga** |

> ✏️ **Exercício 5.** Rode `veh.calculate_load_transfer(a, v)` para a ∈ {0, 5, 10, 15} m/s² a v = 0. Quanto Fz cresce da largada estática até 15 m/s²? Que efeito isso tem sobre a tração máxima disponível?

---

## 16. `Tire.py` — o pneu (Pacejka PAC2002 / Magic Formula)

Este é o módulo que decide se o carro **anda** ou **patina**. É provavelmente o mais importante para o tempo final.

### 16.1 A pergunta central

> *Dada a carga vertical Fz sobre o pneu e o escorregamento κ, quanta força longitudinal Fx ele consegue transmitir ao chão?*

### 16.2 Primeiro: o que é slip ratio (κ)?

Um pneu **não** transmite força sem escorregar um pouco. Isso não é falha — é o mecanismo. A borracha na região de contato se deforma; parte dela adere e parte desliza. A medida disso é o **slip ratio**:

```
κ = (velocidade periférica da roda − velocidade do carro) / referência
```

| κ | Significado |
|---|---|
| 0 | roda gira exatamente na velocidade do carro — **força zero** |
| 0,1 | roda gira 10 % mais rápido — transmitindo boa força |
| ~0,2 | **pico de força** para pneu slick FSAE |
| 0,5 | patinando bastante — força já **caiu** |
| 1,0 | roda girando no vazio (carro parado) — força bem baixa |

> 🧠 **Contraintuitivo mas essencial:** patinar mais **não** dá mais tração. Existe um ponto ótimo (κ ≈ 0,2) e, além dele, a força **cai**. É por isso que existe controle de tração, e é por isso que um piloto bom não "queima" o pneu na largada.

### 16.3 A fórmula robusta do slip (e por que não é a clássica)

A definição clássica é `κ = (ωr − v)/v`. Problema: **na largada v = 0** e isso explode para infinito.

Nossa versão:

```
κ = (ω·r − v) / max( |ω·r|, |v|, 0,1 )
```

```python
v_wheel = omega * raio_pneu
denom   = np.maximum.reduce([np.abs(v_wheel), np.abs(v), np.full_like(v, eps)])
slip    = (v_wheel - v) / denom
return np.clip(slip, -slip_max, slip_max)
```

Dividindo pelo **maior** dos termos, o resultado fica naturalmente limitado a [−1, +1]:
- Roda girando no vazio (v = 0, ω > 0): denominador = |ω·r| → **κ = +1** exatamente.
- Carro acompanhando a roda: κ → 0.
- Ambos parados: o piso `eps = 0,1` evita 0/0.

> ⚠️ **Bug histórico corrigido.** Com a fórmula clássica, o slip explodia no primeiro instante da largada, jogando a Magic Formula na cauda descendente da curva. O carro "perdia" tração fictícia justo na largada — o momento que mais importa.

### 16.4 A Magic Formula

```
Fx = Dx · sin( Cx · arctan( Bx·κ − Ex·(Bx·κ − arctan(Bx·κ)) ) )
```

Ela é chamada "mágica" porque **não vem de nenhuma dedução física** — é uma forma funcional empírica que, com 4 coeficientes bem escolhidos, reproduz curvas de pneu reais com precisão notável. Foi proposta por Hans Pacejka nos anos 1980 e virou padrão da indústria.

**O papel de cada coeficiente:**

```
    Fx │
       │        ╭─── Dx (altura do pico)
       │      ╱   ╲___________
       │    ╱                 ← Cx controla a queda pós-pico
       │  ╱ ← Kxκ (inclinação na origem)
       │╱
       └─────┬─────────────────── κ
           κ_peak (Ex ajusta onde fica)
```

| Coeficiente | Nome | Efeito |
|---|---|---|
| **Dx** | valor de pico | altura máxima da curva = `\|μx\|·λμx·Fz` |
| **Cx** | fator de forma | controla a assíntota: `Dx·sin(Cx·π/2)`. Quanto maior Cx, mais a curva cai depois do pico |
| **Bx** | fator de rigidez | `Kxκ/(Cx·Dx)` — junto com Cx e Dx define onde fica o pico |
| **Ex** | fator de curvatura | ajusta a forma perto do pico (assimetria tração/frenagem) |
| **Kxκ** | rigidez longitudinal | inclinação da curva na origem |

E os parâmetros dependem da carga através de `dfz = (Fz − Fz0)/Fz0`:

```
Cx  = PCX1
μx  = (PDX1 + PDX2·dfz)·(1 − PDX3·γ²)          γ = camber
Dx  = |μx| · λμx · Fz
Kxκ = Fz·(PKX1 + PKX2·dfz)·exp(PKX3·dfz)
Bx  = Kxκ/(Cx·Dx)
Ex  = (PEX1 + PEX2·dfz + PEX3·dfz²)·(1 − PEX4·sgn(κ))     e Ex ≤ 1
```

### 16.5 ⚠️ A calibração: dois dicionários e por quê

O módulo tem **dois** conjuntos de parâmetros. Isso confunde muita gente. Entenda a diferença:

**`HOOSIER_FSAE_LONG_CALSPAN`** — o ajuste **bruto de bancada** (Tire Test Consortium / Calspan). Está lá **só como registro histórico da origem dos dados**. **Não use.**

Por que não? Porque, aplicado tal e qual, ele produz uma curva **sem pico**:

- `PEX1 = 1,346` → a Magic Formula clampa Ex a 1,0 (requisito matemático, Pacejka 2002 eq. 4.E12) → curva monotonicamente crescente → **"patinar mais sempre dá mais tração"**. Fisicamente falso.
- `PCX1 = 1,279` → assíntota em 0,906·Dx → queda máxima de apenas 9,4 %.
- `PKX1 = 58,52` → Bx = 21,7 → pico teórico em κ ≈ 0,13, irrealisticamente baixo.

**`HOOSIER_FSAE_LONG`** — o conjunto **operacional**, corrigido. É este que o dashboard e o otimizador usam.

| Parâmetro | Calspan | Operacional | Motivo |
|---|---|---|---|
| PCX1 | 1,279 | **1,50** | assíntota cai para 0,707·Dx → queda máxima de 29,3 %, criando pico real |
| PKX1 | 58,52 | **23,0** | Bx = 8,53 → κ_peak ≈ 0,20 (slick FSAE típico) |
| PKX2 | 5,48 | **2,15** | escalonado proporcionalmente (mantém a mesma sensibilidade à carga) |
| PEX1 | 1,346 | **−0,50** | de platô plano para queda progressiva. PEX1 < 0 é fisicamente correto para pneus de competição (Pacejka 2002, §4.3.2) |

**⚠️ Repare no que NÃO foi mexido:** PDX1, PDX2, PDX3 — os coeficientes de **nível de atrito**. As correções mudaram a **forma** da curva, não o **quanto** de atrito o pneu tem. Isso é metodologicamente importante: não estamos inventando aderência, estamos consertando a forma.

**Resultado verificado numericamente:**
```
κ_peak = 0,204   |   queda a κ=0,5: −11,2 %   |   queda a κ=1,0: −19,7 %
```

### 16.6 O fator λμx (`tire_friction_coef`)

`tire_friction_coef = 0.6`. Ele multiplica direto o pico Dx.

> 🧠 **Por que 0,6 e não 1,0?** Os dados de bancada TTC/Calspan são medidos numa **esteira de lixa** (*sandpaper belt*), que gera aderência bem maior que asfalto real. O pico |PDX1| ≈ 3,5 é ~50 % inflado. A correção padrão da comunidade é aplicar um fator de escala (*scrub factor*) de 0,55–0,70.
>
> ⚠️ **Nunca use λμx > 1,0 com o dataset `HOOSIER_FSAE_LONG`** — seria dizer que o pneu adere melhor no asfalto do que na lixa.

### 16.7 A convenção de sinal de Fz

O ajuste PAC2002 original da Hoosier usa **Fz negativo** (apontando para baixo, norma ISO 8855). Com PDX1 ≈ −3,51 e Fz < 0, o produto D = μx·Fz sai positivo naturalmente.

Nosso código usa **Fz positivo (magnitude)** e toma `|μx|`. Numericamente idêntico:

```python
Fz_safe = max(abs(float(Fz)), 1e-6)
mux = (p['PDX1'] + p['PDX2'] * dfz) * (1.0 - p['PDX3'] * self.camber ** 2)
Dx = abs(mux) * self.tire_friction_coef * Fz_safe
```

> Isso **não é um workaround** — é a outra metade da mesma convenção, aplicada consistentemente.

### 16.8 Fz0 e a extrapolação em baixa carga

`Fz0 = 654 N` (~150 lbf) é a carga nominal de referência do ajuste — padrão TTC para slick FSAE de 13".

> ⚠️ **Cuidado conhecido.** O dataset Hoosier responde mal em Fz baixo. `PDX2 = +0,633` faz μx **subir** quando dfz < 0 (carga abaixo da nominal) — contra a tendência física real (pneus perdem coeficiente de atrito em carga baixa? Não — na verdade *ganham*; o problema é a magnitude da extrapolação fora da faixa medida).
>
> **Mitigação:** escolha Fz0 próximo da carga estática operacional real do carro (m·g/n_rodas). Com 240 kg e 40 % atrás: 471 N por roda traseira estática. Fz0 = 654 N é razoável considerando que sob aceleração a carga sobe para ~700–900 N.

### 16.9 A consequência: wheelspin emergente

**Este é o resultado mais elegante do projeto.**

A Simulation usa Fx como **teto**, não como valor:

```python
Ft = sign(F_ideal) · min( |F_ideal|,  n_rodas · Fx_Pacejka(Fz, κ) )
```

Quando o motor pede mais força do que o pneu transmite:
1. O excesso **não acelera o carro**;
2. Mas **acelera a roda** (ela é um corpo separado! §1.3);
3. A roda acelerando faz κ subir;
4. κ subindo além de 0,204 joga a curva de Pacejka na **cauda descendente**;
5. Fx **cai** → menos tração ainda.

**Resultado: wheelspin emerge da física, sem uma única linha de código dizendo "se patinar, faça X".** Não há lógica especial, não há flag `is_spinning`. É o modelo se comportando corretamente.

> 🧠 **Por que isso é importante para o otimizador.** Porque significa que **patinar já custa tempo naturalmente**. Não precisamos adicionar uma penalidade artificial de slip na função objetivo — a física cobra a conta sozinha. (E, aliás, a penalidade antiga que existia era silenciosamente nula por um bug de chave de dicionário — descobriu-se isso justamente ao validar essa propriedade.)

### 16.10 Armadilhas

| Armadilha | Consequência |
|---|---|
| Usar `HOOSIER_FSAE_LONG_CALSPAN` | Curva sem pico → patinar dá mais tração → otimizador escolhe relações absurdas |
| λμx > 1,0 | Aderência fisicamente impossível |
| "Consertar" PDX1 em vez de usar λμx | Você mistura nível de atrito com forma de curva; perde a rastreabilidade dos dados originais |
| Esquecer o `× n_driven_wheels` | Metade da tração real (Fz é por roda, Fx é de um pneu) |
| Limitar κ artificialmente antes da MF | Impede a operação na cauda descendente, que é fisicamente real |

> ✏️ **Exercício 6.** Plote `Tire_forces(654, κ)` para κ de 0 a 1. Encontre o pico numericamente. Depois troque para o dataset CALSPAN e plote de novo. Explique em uma frase por que o segundo gráfico é fisicamente impossível.

---

## 17. `Pedal.py` — a entrada do piloto

### 17.1 O papel

O pedal é a **única entrada de comando** do "piloto" no gêmeo digital. Ele produz um número entre 0 e 1, que vira referência de corrente:

```
iq_ref = max_current × posição_do_pedal
```

Controle de **torque em malha aberta** — exatamente como um carro elétrico real (o piloto é a realimentação).

### 17.2 Dois modos

**Modo estático:** `set_posicao(1.0)` → pedal cravado no fundo desde t = 0. Simula uma largada "brutal".

**Modo perfil temporal:** `set_profile(tempos, posicoes)` → uma sequência (t, posição) interpolada linearmente por `np.interp`. Fora do intervalo, aplica **clamp**: antes do primeiro ponto usa o valor inicial; depois do último, mantém o valor final (o pedal "segura" onde parou).

Detalhe de conveniência: se o vetor de posições tiver máximo > 1, o código assume percentual e divide tudo por 100 automaticamente.

### 17.3 A rampa smoothstep — e por que não uma rampa linear

```python
@staticmethod
def perfil_suave(t_inicio, t_fim, n=100):
    t = np.linspace(t_inicio, t_fim, n)
    s = (t - t_inicio) / (t_fim - t_inicio)
    posicoes = 3*s**2 - 2*s**3        # polinômio de Hermite
    return t, posicoes
```

O polinômio `f(s) = 3s² − 2s³` tem três propriedades:
- `f(0) = 0`, `f(1) = 1` (vai de 0 a 100 %);
- `f'(0) = 0` e `f'(1) = 0` — **derivada nula nas duas extremidades**.

```
 100% ┤                        ╭─────  smoothstep
      │                     ╭──╯
      │                  ╭──╯   ╱  rampa linear
      │              ╭───╯    ╱
      │          ╭───╯      ╱
      │      ╭───╯        ╱
   0% ┼──────╯──────────╱────────────
      0                        t_fim
                            ↑ a linear tem "quina" aqui
```

> ⚠️ **Por que isso importa.** Uma rampa linear termina com **descontinuidade de taxa**: a velocidade de variação do pedal salta de constante para zero instantaneamente. Isso é uma entrada tipo "rampa-parada", que excita oscilações nos controladores PI. Foi um problema real observado no simulador. O smoothstep entra e sai suavemente e elimina o transiente.

### 17.4 O pedal como launch control

No otimizador, `t_fim` (chamado `t_pico`) é uma **variável de decisão**. Atrasar o pico do pedal é um **launch control primitivo**:

- Pedal-degrau (t_pico = 0): torque máximo com o carro parado → wheelspin violento → κ dispara → tração cai → tempo ruim.
- Pedal rampeado (t_pico = 1 s): o torque cresce junto com a velocidade → o slip fica perto do ótimo por mais tempo.

**Resultado medido:** t_pico = 3,0 s → **5,104 s**; t_pico = 1,0 s → **4,444 s**. Só mudando como o piloto pisa.

> 🧠 **Isso vale como conclusão de engenharia real.** Antes de gastar dinheiro em hardware, olhe para a largada. O ganho de 0,66 s aqui é maior do que muitos ganhos que você conseguiria trocando componentes.

---

## 18. `PIDController.py` — o controlador PI genérico

### 18.1 O papel

Duas instâncias independentes: uma controla `isd`, outra controla `isq`. Ambas com os mesmos ganhos (a planta é a mesma, já que Ld = Lq).

`kd = 0` em ambas — só PI, sem derivativo (§6.2).

### 18.2 A sequência do `update()`

```python
self.integral += error * dt                     # 1. integra provisoriamente

derivative = 0.0                                # 2. derivada (só se kd ≠ 0)
if dt > 0 and self.kd != 0.0:
    derivative = self.kd * (error - self.prev_error) / dt

output_raw = self.kp*error + self.ki*self.integral + derivative   # 3. soma

if output_raw > self.limit:      output = self.limit              # 4. satura
elif output_raw < -self.limit:   output = -self.limit
else:                            output = output_raw

if self.anti_windup_enabled and self.ki > 0.0 and output_raw != output:
    self.integral -= (output_raw - output) / self.ki              # 5. devolve
```

> 🧠 **Por que integrar ANTES de saturar e corrigir depois?**
> A alternativa seria "integração condicional" (congelar o integrador quando saturado). As duas abordagens são equivalentes em regime, mas a back-calculation **converge de forma mais limpa e simétrica** ao sair da saturação. O congelamento tende a produzir um ciclo-limite assimétrico: ele para de integrar num limite mas não no outro.

### 18.3 `back_calculate()` — o gancho externo

```python
def back_calculate(self, excess):
    if self.anti_windup_enabled and self.ki > 0.0 and excess != 0.0:
        self.integral -= excess / self.ki
```

Existe porque a Simulation aplica uma **segunda** saturação, a ±Vdc, depois de somar o feedforward de desacoplamento (§6.5). Esse limite é **dinâmico** (a tensão da bateria muda com a carga) e o PI não tem como saber dele.

**Convenção de sinal:** excesso **positivo encolhe** o integrador (puxa a saída para baixo); negativo o expande. Zero é ignorado.

### 18.4 `reset()`

Chamado pela Simulation no início de **cada** `simulate()`. Garante que execuções sucessivas — o otimizador roda ~620 delas! — partam do mesmo estado. Sem isso, a segunda simulação começaria com o integrador contaminado da primeira.

### 18.5 Armadilhas

| Armadilha | Consequência |
|---|---|
| Desligar `anti_windup_enabled` | Overshoot violento ao sair de saturação de tensão |
| Esquecer de chamar `back_calculate` após o corte externo | "Windup silencioso" — não aparece em nenhum log óbvio |
| Chamar `update()` mais de uma vez por passo | Integrador anda o dobro; corrompe a sintonia |
| Esquecer `reset()` entre simulações | Contaminação entre execuções do otimizador |

---

# Parte IV — O motor de simulação (`Simulation.py`)

*Este é o arquivo mais longo do projeto (~930 linhas) e o único que tem estado dinâmico. Vamos por partes.*

---

## 19. `__init__` — a montagem do sistema

Nada de física acontece aqui. O construtor **prepara o terreno** para o loop rodar rápido e correto.

### 19.1 Leitura dos parâmetros por `getattr`

```python
self.p        = getattr(motor, 'p',        1)
self.ld       = getattr(motor, 'ld',       1.0)
self.rs0      = getattr(motor, 'rs0', getattr(motor, 'rs', 0.0))
...
```

Todos os modelos são **opcionais**. Se você passar `battery=None`, a simulação simplesmente não integra os estados de bateria e usa `Vdc = motor.Vdc` (600 V de fallback). Isso permite testar subsistemas isoladamente.

O `getattr` com default também protege contra modelos incompletos vindos do dashboard.

### 19.2 Montagem de `J_eff` (a subtração crítica)

```python
if vehicle is not None and transmission is not None:
    j_all         = vehicle.calculate_reflected_inertia(transmission)
    r, N          = vehicle.wheel_radius, transmission.final_drive_ratio
    j_translation = vehicle.mass * (r / N) ** 2
    self.j_eff    = max(self.jm + j_all - j_translation, self.jm)
else:
    self.j_eff = self.jm
```

Já explicado em §15.5. O `max(..., self.jm)` é uma rede de segurança: garante J_eff ≥ J_rotor mesmo com parâmetros degenerados que alguém digite no dashboard.

### 19.3 Pré-computação para o caminho quente

```python
self.inv_ld  = 1.0 / self.ld
self.inv_lq  = 1.0 / self.lq
self.inv_jm  = 1.0 / self.j_eff
```

> 🧠 **Por que isso importa.** A ODE é chamada **4 vezes por passo × 51.000 passos ≈ 204.000 vezes**. Divisão em ponto flutuante é ~4× mais cara que multiplicação. Trocar 3 divisões por 3 multiplicações em cada chamada economiza tempo real de CPU. É uma otimização legítima porque L e J não mudam durante a simulação.

Mesmo raciocínio em `BatteryPack.inv_carga_total` e `Vehicle.half_rho`.

### 19.4 Sintonia automática dos PIs

```python
wc_curr = 2.0 * np.pi * 500.0      # 500 Hz de banda
kp_curr = wc_curr * self.ld        # = 0,303
ki_curr = wc_curr * self.rs        # = 22,18

self.id_ctrl = PID.Controller(kp_curr, ki_curr, 0.0, limit=self.Vdc_nominal, Ts=self.hp)
self.iq_ctrl = PID.Controller(kp_curr, ki_curr, 0.0, limit=self.Vdc_nominal, Ts=self.hp)
```

Os ganhos são calculados **dos parâmetros reais do motor** (§6.3). Consequência prática valiosa: se alguém trocar o motor no dashboard, o tuning se ajusta sozinho. Nada de ganhos hardcoded que ficam errados na primeira mudança de configuração.

### 19.5 O piso de tensão

```python
self.Vdc_nominal = battery.calcular_tensao_nominal()   # 532,8 V
self.Vdc_floor   = 0.3 * self.Vdc_nominal              # 159,8 V
```

Rede de proteção contra o colapso descrito em §25.1. Se por algum motivo a bateria retornar tensão absurdamente baixa, o simulador não deixa cair abaixo de 30 % da nominal — o que evita a divergência numérica em cascata, ainda que o resultado daquela simulação deva ser considerado suspeito.

### 19.6 Constante de tempo do filtro de velocidade

```python
self.tau_wm = 5e-3     # τ = 5 ms  →  corner ≈ 31,8 Hz
```

Um filtro passa-baixa de primeira ordem aplicado a `wm`, usado **somente** pelos limitadores de rotação e potência. Explicado em §22.4.

### 19.7 O print de diagnóstico

```
[Simulation] Kt=0.7122 N·m/A  J_eff=0.0693 kg·m²  Vdc_nom=532.8 V
             PI corrente kp=0.303 ki=22.180
             (ref. PI velocidade kp=3.058 ki=0.918 ζ≈6.06)
```

> 🧠 **Leia sempre esse print.** Ele é a sua primeira verificação de sanidade. Se Kt sair errado, seu λm ou p estão errados. Se J_eff sair muito grande, provavelmente há dupla contagem de massa. Se Vdc_nom sair muito baixo, o arranjo da bateria está errado.
>
> Os ganhos do "PI de velocidade" são só informativos — **não existe** PI de velocidade no laço de controle. Ficaram no print como referência de diagnóstico.

---

## 20. O loop principal — os 8 passos

Este é o coração. Vamos passo a passo, e depois faremos um exemplo numérico completo.

```
┌─── PARA CADA PASSO (dt = 1e-4 s) ────────────────────────────────┐
│                                                                   │
│  1. Desempacota x[k]                                              │
│  2. I_bat = P_mec / (η_drive · V_nominal)     ← balanço de potência│
│  3. Vdc = bateria(I_bat, SoC, I*)             ← tensão deste passo │
│  4. iq_ref = pedal → [3 limitadores]                              │
│     id_ref = field weakening (se necessário)                      │
│  5. vd, vq = PIs de corrente + desacoplamento, cortados a ±Vdc    │
│  6. x[k+1] = RK4(x[k], vd, vq, I_bat)   ← entradas CONGELADAS      │
│  7. Loga x[k] e a telemetria                                      │
│  8. Para se posição ≥ 75 m                                        │
└───────────────────────────────────────────────────────────────────┘
```

> ⚠️ **A ORDEM IMPORTA.** Não é arbitrária. Repare na cadeia de dependências: a corrente depende do estado; a tensão depende da corrente; as referências dependem da tensão (field weakening!); os comandos dependem das referências e da tensão; a física depende dos comandos. Trocar dois passos de lugar quebra alguma dessas dependências.

### 20.1 Passo 1 — desempacotar

```python
(isd, isq, iso, wm, theta_m, temp, vv, vp, soc, Iast, tacc) = self._unpack(x)
theta_e = self.p * theta_m
we      = self.p * wm
```

`_unpack()` sempre devolve os **11** valores, mesmo que o vetor real tenha menos (subsistemas ausentes recebem defaults neutros). Isso mantém o resto do código independente da composição do vetor.

### 20.2 Passo 2 — a corrente da bateria (⚠️ o invariante mais importante)

```python
Ce_now = self.Kt * isq                                # torque atual [N·m]
P_mec  = Ce_now * wm                                  # potência mecânica [W]
P_ac   = max(P_mec / self.eta_drive, 0.0)             # potência elétrica [W]
I_bat  = P_ac / max(self.Vdc_nominal, 1.0)            # corrente DC [A]
```

**A cadeia lógica, em palavras:**
1. O torque vezes a rotação é a potência mecânica que sai pelo eixo.
2. Para produzir isso, o conjunto inversor+motor consome mais, porque tem perdas: divide por η_drive = 0,90.
3. Essa potência elétrica sai da bateria como corrente DC: divide pela tensão.

**Três detalhes deliberados:**

- **Usa `V_nominal` no denominador, não a tensão terminal.** Se usasse a terminal, teríamos referência circular (a tensão depende da corrente que depende da tensão...). A aproximação é boa: o erro é da ordem do sag relativo (~5 %).
- **`max(..., 0.0)`** — apenas potência **positiva** é extraída. Não há regeneração neste modelo.
- **Na primeira iteração** (isd = isq = 0) temos I_bat = 0 e Vdc = tensão de circuito aberto. Início limpo.

> ⚠️ **NÃO REGREDIR.** Isto é a correção do "bug raiz". Detalhes completos em §25.1.

### 20.3 Passo 3 — a tensão da bateria

```python
if self.battery is not None and 0.0 < soc <= 1.0:
    try:
        Vdc = self.battery.calcular_tensao(I_bat, soc, Iast, tacc)
        Vdc = max(Vdc, self.Vdc_floor)
    except Exception:
        Vdc = self.Vdc_nominal
else:
    Vdc = self.Vdc_nominal
Vlim = Vdc
```

`Vlim` é o **teto dos comandos do inversor neste passo**. Sob carga pesada a bateria afunda e o FOC perde margem de tensão — o que motiva o field weakening (§23).

### 20.4 Passo 4 — as referências de corrente

Detalhado em §22 (limitadores de iq) e §23 (field weakening).

### 20.5 Passo 5 — os PIs e o desacoplamento

```python
dec_d = -we * self.lq * isq                    # feedforward eixo d
dec_q =  we * (self.ld * isd + self.lambda_m)  # feedforward eixo q (inclui back-EMF)

vd_pid = self.id_ctrl.update(id_ref - isd, dt=dt)
vq_pid = self.iq_ctrl.update(iq_ref - isq, dt=dt)

vd_total = vd_pid + dec_d
vq_total = vq_pid + dec_q

vd = float(np.clip(vd_total, -Vlim, Vlim))
vq = float(np.clip(vq_total, -Vlim, Vlim))

if vd_total != vd: self.id_ctrl.back_calculate(vd_total - vd)
if vq_total != vq: self.iq_ctrl.back_calculate(vq_total - vq)
```

**Por que o feedforward?** Volte a §5.6: as equações dq são acopladas. Os termos `+ωe·Lq·isq` e `−ωe·(Ld·isd + λm)` fazem cada eixo perturbar o outro.

A estratégia: **calcule esses termos e some-os ao comando**, com sinal invertido em relação a como aparecem na planta. Eles se cancelam, e cada PI passa a controlar uma planta de 1ª ordem simples e independente.

> 🧠 **Analogia.** Você está mirando com um rifle num dia de vento lateral. Você **sabe** que o vento vai empurrar a bala 10 cm para a direita. Duas estratégias: (a) atirar e corrigir depois de ver onde acertou (só realimentação — lento); (b) **já mirar 10 cm à esquerda** (feedforward) e usar a realimentação só para o erro residual. O feedforward mata a perturbação conhecida antes que ela apareça.
>
> É por isso que `dec_q` inclui `ωe·λm` — o back-EMF é uma perturbação **perfeitamente conhecida**. Não faz sentido esperar o PI descobrir sozinho.

### 20.6 Passo 6 — o RK4 com entradas congeladas

```python
telem = {}
x = self._rk4(t, x, vd, vq, I_bat, dt, telemetry=telem)
```

`vd`, `vq` e `I_bat` ficam fixos nas 4 avaliações (§2.6). O dicionário `telem` é uma otimização inteligente:

> 🔍 A ODE calcula slip, Fz, Fx, Ft, resistência — tudo isso o logger também quer. Em vez de recomputar (rodando a Magic Formula duas vezes por passo!), a ODE escreve esses intermediários em `telem` — **mas só na primeira avaliação (k1)**, que é a que corresponde exatamente a x[k]. As avaliações k2, k3, k4 recebem `telemetry=None` e não escrevem nada.
>
> Isso garante que o log é consistente com o estado **antes** do passo, e economiza ~25 % do custo de pneu.

### 20.7 Passo 7 — logging em arrays pré-alocados

```python
k = self.k_log
self.tempo[k]      = t
self.corrented[k]  = isd
self.correnteq[k]  = isq
...
self.k_log = k + 1
```

**Sem `.append()`.** Os arrays são pré-alocados com `np.zeros(N)` em `_init_storage(N)` e escritos por índice.

> 🧠 **Por quê?** `list.append()` faz realocação amortizada — a cada tanto ele copia a lista inteira para um bloco maior. Com ~40 variáveis × ~51.000 passos, isso viraria overhead perceptível. Pré-alocar e escrever por índice é O(1) garantido. No final, os arrays são truncados para o tamanho realmente usado.

**Note também:** loga-se `x[k]` — o estado **antes** do passo — junto com os comandos de controle daquele passo. Isso mantém alinhamento temporal correto entre estado e comando.

### 20.8 Passo 8 — o critério de parada

```python
if self.dmax is not None and vp >= self.dmax:
    break
```

Para ao cruzar 75 m. O `tmax` (30 s por padrão) é apenas uma rede de segurança para configurações que nunca chegariam lá.

---

## 21. A ODE (`_physics_ode`) linha a linha

Sete blocos, em ordem.

### 21.1 Rs(T)

```python
rs_eff = self.rs0 * (1.0 + self.alpha_cu * (temp - self.T_ref))
```

Avaliado **dentro** da ODE, e não fora, para que cada estágio do RK4 use o Rs consistente com o `temp` **daquele estágio**. Detalhe pequeno, mas é o tipo de coisa que separa uma integração correta de uma aproximada.

### 21.2 Dinâmica elétrica dq

```python
d_isd = (vd - rs_eff*isd + we*self.lq*isq) * self.inv_ld
d_isq = (vq - rs_eff*isq - we*(self.ld*isd + self.lambda_m)) * self.inv_lq
d_iso = (-rs_eff * iso) / self.L0
```

Já detalhado em §5.6. `iso` (sequência zero) decai exponencialmente e não afeta torque; existe por completude. `L0` é aproximado como 10 % da média de Ld e Lq, porque o datasheet não fornece esse valor.

### 21.3 Acoplamento pneu–veículo

O bloco mais rico. Sequência:

```python
# 1. Força que o motor PEDE no contato (ainda sem limite de aderência)
Tmotor_roda   = self.transmission.motor_to_wheel_torque(self.Kt * isq)
Ftração_ideal = Tmotor_roda / r

# 2. Slip entre a roda e o carro
omega_roda = self.transmission.motor_to_wheel_speed(wm)
slip       = Tire.SlipRatio(omega_roda, r, vv)
Fresist    = self.vehicle.calculate_resistance_forces(vv)

# 3. Loop algébrico de 2 iterações (Fz ↔ aceleração) — §15.3
Fz0 = self.vehicle.calculate_load_transfer(0.0, vv)
Fx0 = self.tire.Tire_forces(Fz0, slip) * n_drv
Ft0 = np.sign(Ftração_ideal) * min(abs(Ftração_ideal), abs(Fx0))
a0  = (Ft0 - Fresist) / m

Fz     = self.vehicle.calculate_load_transfer(a0, vv)
Fx_max = self.tire.Tire_forces(Fz, slip) * n_drv
Ft     = np.sign(Ftração_ideal) * min(abs(Ftração_ideal), abs(Fx_max))

# 4. Reação no eixo do motor
Tcarga = Ft * r / (N * eta)

# 5. Equação translacional do carro
Fr   = Ft - Fresist
d_vv = 0.0 if (vv < 0.01 and Fr < 0) else Fr / m
```

**Comentários sobre pontos específicos:**

- **`min(|F_ideal|, |Fx_max|)`** — é AQUI que o wheelspin nasce (§16.9).
- **`* n_drv`** — Fz é por roda, Fx é de um pneu, o eixo tem `n_driven_wheels` = 2. Em linha reta com tração simétrica as duas contribuem igualmente.
- **`vv < 0.01 and Fr < 0`** → `d_vv = 0`. Isso impede o carro de "andar para trás" parado no grid, caso a resistência supere a tração inicial. É uma restrição unilateral simples, não um modelo de freio.
- **`try/except`** envolvendo tudo — se algo falhar no pneu (parâmetros absurdos vindos do dashboard), imprime e segue com forças zeradas em vez de derrubar a simulação inteira.

### 21.4 Perdas no ferro → torque de freio

```python
f_elec = abs(we) * (1.0 / (2.0 * np.pi))
P_iron = self.k_h_iron * f_elec + self.k_e_iron * f_elec * f_elec
wm_safe = max(abs(wm), self.wm_floor)
T_iron_brake = (P_iron / wm_safe) * np.sign(wm) if wm != 0.0 else 0.0
```

> 🧠 **Como uma perda em watts vira torque?** Por definição, P = T·ω. Se o ferro dissipa P_iron watts, o equivalente mecânico é um torque `T = P_iron/ω` **opondo-se** ao movimento. O `wm_floor = 1,0 rad/s` evita a singularidade quando ω → 0 (na largada, P_iron também é ~0, então isso é seguro).

Essa energia perdida vai **integralmente** para o aquecimento — ela reaparece somada às perdas do cobre na equação térmica logo abaixo.

### 21.5 Mecânica do eixo

```python
Ce      = self.Kt * isq
d_wm    = (Ce - Tcarga - self.kf * wm - T_iron_brake) * self.inv_jm
d_theta = wm
```

A equação de Newton rotacional, com quatro termos:
- `+Ce` — torque eletromagnético, o que empurra;
- `−Tcarga` — reação da força de tração, o que o carro "cobra";
- `−kf·wm` — atrito viscoso (rolamentos, ventilação);
- `−T_iron_brake` — freio equivalente das perdas no ferro.

### 21.6 Térmica

```python
P_cu   = 1.5 * rs_eff * (isd**2 + isq**2)
P_cool = self.h_cool * self.A_cool * max(temp - self.T_ambient, 0.0)
d_temp = (P_cu + P_iron - P_cool) * self.inv_mC
```

O `1,5` vem da convenção amplitude-invariante — **o mesmo 1,5 do Kt** (§5.5).

O `max(temp − T_ambient, 0)` impede resfriamento "negativo" (o ar não aquece o motor se este estiver mais frio que o ambiente — simplificação razoável).

### 21.7 Bateria

```python
if self.battery:
    dsoc, dIast, dtacc = self.battery.calcular_derivadas(I_bat)
```

Com `I_bat` **congelado** — a mesma corrente que foi usada para calcular a tensão. Isso garante consistência entre tensão terminal e evolução do SoC dentro do passo.

### 21.8 Montagem do vetor de derivadas

```python
dx = [d_isd, d_isq, d_iso, d_wm, d_theta, d_temp]
if self.vehicle and self.transmission:
    dx += [d_vv, vv]        # ← note: a derivada da POSIÇÃO é a VELOCIDADE
if self.battery:
    dx += [dsoc, dIast, dtacc]
return np.array(dx, dtype=float)
```

> 🧠 O par `[d_vv, vv]` costuma confundir. O estado é `[velocidade, posição]`; portanto as derivadas são `[aceleração, velocidade]`. É o truque padrão de reduzir uma equação de segunda ordem a duas de primeira ordem.

---

## 22. Os três limitadores de `iq_ref`

O pedal pede uma corrente. Três protetores podem reduzi-la, e vale **o mais restritivo dos três**.

```python
iq_ref = float(np.clip(min(iq_pedal, iq_from_speed, iq_from_power), 0.0, self.max_current))
```

### 22.1 Limitador térmico

```python
if temp >= self.T_max:        thermal_derate = 0.0
elif temp > self.T_alarm:     thermal_derate = (self.T_max - temp) / (self.T_max - self.T_alarm)
else:                         thermal_derate = 1.0

iq_pedal = self.max_current * pedal_pos * thermal_derate
```

Rampa linear de 1 → 0 entre T_alarm = 130 °C e T_max = 160 °C (limites típicos de isolamento classe H). Acima de T_max, corrente zero.

Numa prova de 5 s isso **nunca dispara** — mas está lá para quando o modelo for estendido ao Endurance.

### 22.2 Limitador de sobre-rotação (smoothstep)

```python
_lo = 0.90 * self.speed_ref            # 518,4 rad/s
_hi = self.speed_ref                   # 575,96 rad/s
_prog  = clip((self.wm_filt - _lo) / (_hi - _lo), 0.0, 1.0)
_s     = 1.0 - _prog                   # 1 em _lo, 0 em _hi
_fator = _s * _s * (3.0 - 2.0 * _s)    # smoothstep
iq_from_speed = self.max_current * _fator
```

Fade suave de 100 % → 0 % da corrente entre 90 % e 100 % de `speed_ref`.

> 🧠 **Por que smoothstep e não uma rampa linear?**
> Uma rampa linear tem **"joelho"**: em 0,9·speed_ref o ganho salta descontinuamente de 0 para um valor finito. Existe um laço de realimentação implícito aqui — `wm → iq_ref → T_em → wm` — e um ganho descontínuo dentro de um laço realimentado é receita para ciclo-limite (oscilação sustentada). O smoothstep tem derivada nula nos dois extremos e não excita esse modo.
>
> Este é o **mesmo raciocínio** do perfil de pedal (§17.3). Não é coincidência: descontinuidade de derivada dentro de laço realimentado sempre dá problema.

### 22.3 Limitador de potência DC (o regulamento)

```python
if self.p_max_dc is not None and self.Kt > 0.0:
    wm_pwr        = max(self.wm_filt, self.wm_floor)
    iq_from_power = (self.p_max_dc * self.eta_drive / (self.Kt * wm_pwr))
```

Vem de inverter o balanço de potência:

```
P_dc ≈ P_mec/η = Kt·iq·ω/η  ⟹  iq_max = P_max_dc · η / (Kt · ω)
```

Isso produz **o perfil clássico de tração elétrica**:

```
  Torque │
    230  ┤━━━━━━━━━━━━┓
         │            ┃╲                  ← região de potência constante
         │            ┃  ╲___             (torque ∝ 1/ω)
         │            ┃      ╲______
       0 ┼────────────┸─────────────── ω
         0         ω_base
         └──── torque constante ────┘
```

Com os nossos números: ω_base ≈ 80.000 × 0,90 / (0,7122 × 323) ≈ **313 rad/s** ≈ 2990 RPM. Abaixo disso o pedal manda; acima, o teto de 80 kW manda.

> 🧠 **É por isso que a curva de `iq_ref` no gráfico tem uma "corcova" seguida de queda hiperbólica.** A corcova é o pedal subindo; a queda ∝ 1/ω é o limitador de potência entrando. Aquele formato é a assinatura visual do regulamento FSAE EV.4.1 no seu gráfico.

### 22.4 O filtro passa-baixa em `wm`

```python
_alpha_wm = dt / (self.tau_wm + dt)
self.wm_filt += _alpha_wm * (wm - self.wm_filt)
```

Filtro de primeira ordem discreto, τ = 5 ms (corner ≈ 31,8 Hz).

> ⚠️ **O filtro é usado APENAS pelos limitadores 22.2 e 22.3.** A dinâmica do motor continua usando o `wm` exato, sem filtro. Isso é essencial: filtrar a velocidade na física introduziria um atraso fictício e mudaria a resposta do sistema.
>
> **Por que existe:** havia ripple numérico de torque que realimentava em `iq_ref` e disparava ciclo-limite justamente no "joelho" do limitador. O filtro atenua o ripple sem afetar a dinâmica real (31,8 Hz está bem acima da banda mecânica de ~5 Hz e bem abaixo da frequência do ripple).

---

## 23. Field weakening (enfraquecimento de campo)

### 23.1 O problema

Volte a §3.3. O back-EMF `ωe·λm` cresce com a rotação e come a margem de tensão disponível. A **velocidade-base** é onde ele iguala o limite:

```
ωe_base = Vlim / λm     ⟹    ωm_base = Vlim / (p · λm)
```

Com Vlim = 532,8 V (e o fator de segurança 0,9): ωm_base ≈ 0,9 × 532,8 / (10 × 0,04748) ≈ **1010 rad/s** = 9644 RPM.

> 🧠 **Boa notícia para o nosso carro:** 9644 RPM está **muito acima** do limite mecânico de 5500 RPM do EMRAX 228. Ou seja: **na prova de aceleração o field weakening nunca é acionado.** No gráfico de correntes, `id_ref` fica cravado em zero a prova inteira.
>
> Então por que ele existe no código? (a) Porque o dashboard permite trocar bateria e motor, e outras combinações **podem** cruzar a velocidade-base; (b) porque o `Vlim` é dinâmico — se a bateria afundar muito, a velocidade-base cai; (c) porque é fisicamente correto ter.

### 23.2 A ideia

Se o fluxo dos ímãs é o problema, **reduza o fluxo efetivo**. Injetando corrente **negativa** no eixo d, criamos um campo que se opõe parcialmente aos ímãs:

```
λ_efetivo = λm + Ld·isd        com isd < 0 ⟹ λ_efetivo < λm
```

Menos fluxo → menos back-EMF → mais margem de tensão → dá para continuar controlando iq e mantendo torque em rotação alta.

**O preço:** a corrente de eixo d não produz torque nenhum. Você gasta ampères (e calor!) só para comprar margem de tensão. É um trade-off consciente, não um almoço grátis.

### 23.3 A implementação

```python
we_abs = abs(we)
id_ref = 0.0
if (we_abs > self.wm_floor * self.p and self.lambda_m > 0.0 and self.ld > 0.0):
    fw_limit = 0.90 * Vlim
    if we_abs * self.lambda_m > fw_limit:
        id_fw  = (fw_limit / we_abs - self.lambda_m) / self.ld
        id_ref = float(np.clip(id_fw, -self.max_current * 0.5, 0.0))
```

Três proteções:
- **`fw_limit = 0,90·Vlim`** — operamos dentro de 90 % do limite, deixando margem para o PI de corrente trabalhar.
- **`clip(..., −0,5·max_current, 0)`** — id_fw nunca passa de metade da corrente máxima. Proteção contra **desmagnetização dos ímãs** (campo reverso forte pode danificar permanentemente ímãs de neodímio) e contra sobrecorrente.
- **`we_abs > wm_floor·p`** — evita divisão por ~zero em rotação baixa.

---

## 24. Um passo inteiro, com números reais

Nada amarra melhor a teoria do que seguir **um único passo** do começo ao fim. Vamos pegar o instante **t ≈ 2,3 s** da simulação padrão (N = 5,0; pedal smoothstep com pico em 3 s) — o momento de pico de torque e de pico de wheelspin.

**O estado nesse instante** (lido dos gráficos da simulação):

```
isq  ≈ 280 A          isd = 0 A            wm  ≈ 316 rad/s (3018 RPM)
v    ≈ 10,7 m/s       κ   ≈ 0,35           temp ≈ 30 °C
SoC  ≈ 0,997          Iast ≈ pequeno
```

---

### Passo 1 — desempacotar

```
we = p · wm = 10 × 316 = 3160 rad/s elétricos     (f_e = 503 Hz)
θe = 10 × θm
```

### Passo 2 — corrente da bateria

```
Ce    = Kt · isq   = 0,7122 × 280   = 199,4 N·m
P_mec = Ce · wm    = 199,4 × 316    = 63,0 kW
P_ac  = P_mec/η    = 63,0 / 0,90    = 70,0 kW
I_bat = P_ac/V_nom = 70.000 / 532,8 = 131 A
```

> ✅ **Confira no gráfico:** a curva de potência DC mostra ~70 kW nesse instante, e a de corrente de bateria mostra ~130–150 A. Bate.

### Passo 3 — tensão da bateria

```
i_cell  = 131 / 5           = 26,2 A            ← divisão por n_paralelo!
queda   = 26,2 × 0,02       = 0,524 V/célula
sag     = 0,524 × 144       = 75,5 V
Vdc     ≈ 532,8 − 75,5      ≈ 457 V
```

> ✅ **Confira no gráfico:** a tensão da bateria está em ~460 V nesse trecho. Bate.
>
> ⚠️ E veja o que aconteceria **sem** a divisão por n_paralelo: 131 × 0,02 × 144 = 377 V de sag → Vdc = 156 V → colapso. É a §13.4 em ação.

### Passo 4 — as referências de corrente

**Pedal** (smoothstep, t_pico = 3 s, em t = 2,3 s → s = 0,767):
```
f(s) = 3s² − 2s³ = 3(0,588) − 2(0,451) = 0,862
iq_pedal = 323 × 0,862 × 1,0 (sem derate térmico) = 278 A
```

**Limitador de rotação:** wm_filt ≈ 316 ≪ 0,90 × 575,96 = 518 rad/s
```
iq_from_speed = 323 A       (não limita)
```

**Limitador de potência:**
```
iq_from_power = 80.000 × 0,90 / (0,7122 × 316) = 320 A
```

**Resultado:**
```
iq_ref = min(278, 323, 320) = 278 A     ← o PEDAL ainda manda, por pouco
```

> 🧠 Estamos exatamente no **cotovelo**. Mais alguns décimos de segundo e wm cresce, `iq_from_power` cai abaixo de 278 e o **regulamento** passa a mandar. É a "corcova" seguida de queda hiperbólica que você vê no gráfico de `iq_ref`.

**Field weakening:**
```
we · λm  = 3160 × 0,04748 = 150 V
0,90 · Vlim = 0,90 × 457  = 411 V
150 < 411  ⟹  id_ref = 0        (não ativa — como esperado, §23.1)
```

### Passo 5 — comandos de tensão

```
dec_d = −we·Lq·isq        = −3160 × 9,65e-5 × 280 = −85,4 V
dec_q =  we·(Ld·isd + λm) =  3160 × (0 + 0,04748) = +150,0 V   ← back-EMF

erro_q = iq_ref − isq = 278 − 280 = −2 A      (o PI está rastreando bem)
vq_pid ≈ pequeno; o integrador carrega o grosso do comando
vq_total = vq_pid + 150,0        ⟹  bem abaixo de ±457 V, não satura
```

> 🧠 Note que o **feedforward faz quase todo o trabalho**: dos comandos totais, a maior parcela é o desacoplamento (que cancela o back-EMF conhecido) e o PI cuida só do erro residual de 2 A. É assim que deve ser numa malha bem projetada.

### Passo 6 — a física (dentro do RK4)

**O que o motor pede:**
```
T_roda   = 199,4 × 5,0 × 0,9506 = 947,7 N·m
F_ideal  = 947,7 / 0,26         = 3.645 N
```

**O que o pneu permite:**
```
ω_roda   = 316 / 5,0            = 63,2 rad/s
v_roda   = 63,2 × 0,26          = 16,4 m/s
κ        = (16,4 − 10,7)/16,4   = 0,35      ← a roda gira 53% mais rápido que o carro

Fz (por roda traseira):
   estático     = 240 × 9,81 × 0,6/1,5 / 2  = 471 N
   transferência (a ≈ 12) = 240 × 12 × 0,28/1,5 / 2 = 269 N
   downforce traseiro     = ½×1,225×1,0×0,35×10,7² / 2 = 12 N
   Fz ≈ 752 N por roda

Fx(752 N, κ = 0,35) ≈ 1.500 N por pneu       ← já na CAUDA DESCENDENTE
Fx_max do eixo = 2 × 1.500 = 3.000 N
```

**A saturação:**
```
Ft = sign(3645) · min(3645, 3000) = 3.000 N     ← SATURADO PELO PNEU
```

> 🔴 **Isto é wheelspin, em números.** O motor pede 3.645 N. O pneu entrega 3.000 N. Os 645 N de excesso **não empurram o carro**: eles aceleram a roda, aumentam κ, e — como já estamos além de κ_peak = 0,204 — fazem Fx cair ainda mais. Nenhuma linha de código diz "detecte patinagem". A física fez sozinha.

**As duas equações de Newton:**
```
Resistência: aero = ½×1,225×0,7789×0,68×10,7² = 37 N
             rolamento ≈ 38 N        →  F_resist ≈ 75 N

Perdas no ferro: f_e = 503 Hz
   P_iron = 0,9×503 + 9e-4×503² = 453 + 228 = 681 W
   T_ferro = 681/316 = 2,15 N·m

ROTACIONAL:  d_wm = (199,4 − 164,1 − 0,005×316 − 2,15) / 0,0693
                   = 31,6 / 0,0693 = +456 rad/s²
             onde Tcarga = 3000 × 0,26/(5,0 × 0,9506) = 164,1 N·m

TRANSLACIONAL: d_vv = (3.000 − 75) / 240 = +12,2 m/s²
```

**E aqui está o ponto:**

```
A roda acelera:   d(v_roda)/dt = 456/5,0 × 0,26 = 23,7 m/s²
O carro acelera:  d(v)/dt      =                  12,2 m/s²
                                                  ─────────
                  a roda ganha 11,5 m/s² sobre o carro
                  ⟹ κ CONTINUA CRESCENDO
```

É exatamente o que o gráfico de slip mostra: κ subindo até o pico de 0,35 por volta de t = 2,4 s, e só depois caindo quando o limitador de potência corta `iq_ref` e o carro "alcança" a roda.

### Passos 7 e 8

Loga tudo, avança 0,0001 s, verifica se `vp ≥ 75`. Ainda não (estamos em ~14 m). **Repete mais 28.000 vezes.**

---

> 🧠 **A moral desta seção.** Cada número acima é uma conta simples — multiplicação, divisão, uma senóide aqui e ali. Nenhuma delas é difícil. A dificuldade está em fazer as **onze** grandezas evoluírem juntas, de forma consistente, cinquenta mil vezes. É por isso que existe o gêmeo digital: nenhuma conta de guardanapo consegue rastrear esse acoplamento por 5 segundos.

## 25. Os bugs históricos (e por que NÃO regredir)

Esta seção é a memória institucional do projeto. Se você for mexer no código, **leia antes**.

### 25.1 ⚠️ BUG RAIZ — a corrente da bateria

**O que era (v1 e v2):** passar `abs(isq)` como corrente para `calcular_tensao()` e `calcular_derivadas()`.

**Por que está errado:** `isq` é a corrente no eixo q do **referencial dq rotativo do motor**. Ela **não é** a corrente DC do banco de baterias. São grandezas de circuitos diferentes, separadas pelo inversor.

**A cascata de destruição:**

```
1. R_interna (calculada errado, sem escala paralela) = 264 × 0,02 = 5,28 Ω
2. isq ≈ 700 A
3. Queda = 5,28 × 700 = 3.696 V   ≫  V_nominal
4. Vdc  <  0
5. Clamp: Vdc = max(Vdc, 1.0) = 1 V
6. Com 1 V não se cancela o desacoplamento de eixo d (≈ ωe·Lq·isq ≈ 473 V):
      d_isd/dt = (±1 − 0 + 473)/9,65e-5   →   isd → −∞
7. isd → −1000 A  →  torque oscilando ±600 N·m
8. P_cobre = 1,5 × Rs × (isd² + isq²) ≈ 15.800 W  →  temperatura → 4000 °C
9. Tensão da bateria oscilando entre −4000 V e +100 V
```

Uma simulação que "explode" desse jeito é fácil de detectar. O perigoso é o caso intermediário, em que a explosão é parcial e os resultados parecem plausíveis mas estão errados.

**A correção:** balanço de potência (§20.2). `I_bat` é calculada **uma vez** por passo e usada **tanto** para a tensão **quanto** para a ODE — mantendo consistência.

> ⚠️ **Sintomas de regressão:** picos de torque sem explicação, temperatura absurda, tensão de bateria negativa, `isd` grande com `id_ref = 0`. Se você vir qualquer um desses, verifique primeiro o que está sendo passado como corrente para a bateria.

### 25.2 ⚠️ Escala pack → célula no Shepherd

Detalhado em §13.4. Sem a divisão por `n_paralelo`, o sag de tensão sai 5× superestimado e o barramento colapsa.

### 25.3 ⚠️ PIDs dentro da ODE (v1)

Os controladores rodavam dentro da função de derivadas. Cada avaliação do RK4 atualizava o integrador — ou seja, o integrador avançava **4 vezes por passo**, com estados intermediários fictícios. O RK4 deixava de ser matematicamente válido.

**Corrigido na v2:** controladores fora da ODE, saídas congeladas (§2.6, §9.2).

### 25.4 ⚠️ Dupla contagem da massa em J_eff

Detalhado em §15.5. Usar o retorno completo de `calculate_reflected_inertia()` como J_eff conta a massa do carro duas vezes.

### 25.5 ⚠️ Slip explodindo na largada

Detalhado em §16.3. A fórmula clássica `κ = (ωr − v)/v` diverge em v = 0.

### 25.6 ⚠️ Pacejka sem pico

Detalhado em §16.5. O dataset Calspan bruto produz curva monotonicamente crescente — "patinar sempre dá mais tração". O otimizador, sob esse modelo, escolheria relações de transmissão absurdas.

### 25.7 ⚠️ Rampa linear no pedal e no limitador

Detalhado em §17.3 e §22.2. Descontinuidade de derivada dentro de laço realimentado gera ciclo-limite.

### 25.8 O quadro-resumo: decisões que NÃO devem ser revertidas

| # | Decisão | Motivo |
|---|---|---|
| 1 | `I_bat` por balanço de potência (nunca `abs(isq)`) | colapso do barramento DC, picos fantasmas |
| 2 | Conversão pack→célula no Shepherd | sag de tensão 5× superestimado |
| 3 | `speed_ref = 575,96 rad/s` (5500 RPM) | limite real do EMRAX 228; usar 1200 inflava N* para ~12 |
| 4 | `p_max_dc = 80 kW` | regulamento FSAE EV.4.1 |
| 5 | Teto de tração = n_rodas × Fx(Fz_por_roda) | Fz é por roda, Fx é de um pneu |
| 6 | `kf = 0,005` | 0,1 dissipava ~33 kW e dupla-contava as perdas no ferro |
| 7 | Bateria 144S5P (532,8 V) | teto FSAE de 600 V |
| 8 | `max_current = 323 A` | reproduz os 230 N·m do datasheet |
| 9 | Pacejka `HOOSIER_FSAE_LONG` calibrado | dataset Calspan bruto não tinha pico |
| 10 | J_eff sem massa translacional | dupla contagem no modelo de dois corpos |
| 11 | PIs fora da ODE, saídas congeladas | validade matemática do RK4 |
| 12 | Smoothstep no pedal e no limitador de rotação | evita ciclo-limite por descontinuidade de derivada |

---

# Parte V — As ferramentas construídas em volta

## 26. O dashboard web (`main.py`)

```bash
python main.py       # abre em http://127.0.0.1:8050/
```

### 26.1 O que é

Uma aplicação web feita com **Dash** (framework Python da Plotly) que permite mexer em qualquer parâmetro do carro pela interface, rodar a simulação e ver 10 abas de gráficos — **sem escrever uma linha de código**.

É a porta de entrada para quem chega na equipe. Você não precisa entender `Simulation.py` para usar o dashboard.

### 26.2 A estrutura

**Sidebar recolhível** — 7 cartões accordion:

| Cartão | O que você ajusta |
|---|---|
| 🚗 Veículo | massa, raio e massa da roda, Cd, área frontal, Cr, L, h, dist_cg |
| ⚙️ Transmissão | Z_pinhão, Z_coroa (ou relação manual), eficiências, inércias |
| 🪂 Aerodinâmica | tem asa?, Cl e área dianteira e traseira |
| ⚡ Motor | Rs, Ld, Lq, Jm, kf, λm, p, max_current, speed_ref |
| 🔋 Bateria | tipo de célula, n_série, n_paralelo, SoC inicial |
| 🛞 Pneu (PAC2002) | coeficientes Pacejka, λμx, Fz0, camber |
| 🎯 Simulação | tmax, dmax, p_max_dc, perfil de pedal |

E o botão **Executar Simulação**.

**10 abas de resultado:**

| Aba | O que mostra | Quando usar |
|---|---|---|
| Velocidade & Torque | v(t), T_em, T_carga | primeira olhada, sempre |
| Correntes | isd, isq, ia/ib/ic | verificar se o FOC está rastreando |
| Tensões | vd, vq, va/vb/vc | verificar saturação de tensão |
| Fluxo & Temp. | λd, λq, temperatura | diagnóstico térmico |
| Controle FOC | iq_ref vs isq, id_ref | ver qual limitador está mandando |
| Visão Completa | painel resumo | apresentações |
| Veículo | posição, velocidade, aceleração, forças | dinâmica longitudinal |
| Pneu | slip, Fz, Fx | diagnosticar wheelspin |
| Bateria | Vdc, I_bat, SoC, potência DC | verificar o limite de 80 kW |
| 📋 Histórico | comparação entre execuções | estudos de sensibilidade |

### 26.3 O fluxo de callbacks

```
[Executar Simulação]  ──▶  run_simulation_once()
                            1. lê TODOS os campos da sidebar (via State)
                            2. instancia os 6 modelos
                            3. sim = Simulation(...); sim.simulate()
                            4. gera as 9 figuras Plotly
                            5. serializa no store 'stored-figures'
                            6. acrescenta entrada em 'sim-history'
                               com o DIFF de parâmetros

[troca de aba]        ──▶  render_active_tab()
                            apenas busca a figura pronta no store
                            → NÃO re-simula
```

> 🧠 **Por que essa separação importa.** Simular leva 30 s a 2 min. Se cada troca de aba disparasse uma nova simulação, o dashboard seria inutilizável. A estratégia de "simular uma vez, guardar as figuras, servir sob demanda" é o padrão certo para qualquer aplicação Dash com cálculo pesado.

### 26.4 A aba Histórico

Cada execução vira uma linha da tabela, com:
- tempo final e distância;
- **quais parâmetros mudaram** em relação à execução anterior (o diff é calculado por `_diff_params()`).

> 🧠 Isso transforma o dashboard num caderno de laboratório automático. Você muda um parâmetro, roda, e a tabela registra "N: 5,0 → 4,5 · t₇₅: 5,104 → 4,98". É perfeito para estudos de sensibilidade manuais, e evita o clássico "não lembro qual configuração deu aquele resultado bom".

### 26.5 ⚠️ Sincronização manual

> **Os defaults da sidebar do dashboard e o `SIM_CONFIG` do otimizador são arquivos separados e precisam ser mantidos em sincronia manualmente.**
>
> Se você mudar a massa do carro no dashboard mas esquecer de mudar no otimizador, as duas ferramentas simularão carros diferentes — e você vai perder um bom tempo tentando entender por que os números não batem. **Esta é uma dívida técnica conhecida do projeto.**

---

## 27. O otimizador de transmissão (`Otimizador/otimizador.py`)

```bash
python Otimizador/otimizador.py       # ⏱️ horas — são ~620 simulações
```

É **standalone**: não usa o dashboard, instancia tudo por conta própria.

### 27.1 O problema formulado

```
    min      t₇₅ₘ(N, t_pico)
  N, t_pico

  sujeito a:  N ∈ [3, 20]
              t_pico ∈ [0, 3] s
```

Duas variáveis de decisão:
- **N** — a relação final de transmissão (o hardware);
- **t_pico** — o instante em que o pedal atinge 100 % (a técnica de largada).

### 27.2 Por que otimizar o pedal junto com a transmissão?

Esta é a decisão de projeto mais interessante do otimizador.

> 🧠 Com pedal-degrau (t_pico = 0), o "penhasco" de wheelspin prende o ótimo de N **logo abaixo do limiar de patinagem**. Qualquer relação mais agressiva faz o carro patinar na largada e perder tempo. O otimizador conclui "N = 4,2 é o melhor" e para por aí.
>
> Mas isso é uma conclusão sobre **o par (hardware, técnica)**, não sobre o hardware. Com a rampa smoothstep agindo como **launch control**, relações mais agressivas voltam a ser viáveis: o torque cresce junto com a velocidade, o slip fica perto do ótimo por mais tempo, e o penhasco desaparece.
>
> **Lição de engenharia:** otimizar um subsistema isoladamente encontra o ótimo local do subsistema, não do sistema. Se você fixar a técnica de largada, encontra a melhor transmissão *para aquela largada ruim*.

### 27.3 O algoritmo: Evolução Diferencial (DE/rand/1/bin)

**Evolução Diferencial** é um algoritmo de otimização **estocástico** e **sem derivadas**. Isso importa: nossa função objetivo é uma simulação inteira — não temos a derivada de t₇₅ em relação a N, e a superfície tem descontinuidades (o penhasco de wheelspin).

População: **20 indivíduos**, até **30 gerações**. Cada indivíduo é um par `x = [N, t_pico]`.

Para cada indivíduo `xᵢ`:

**1. Mutação**
```
v = xₐ + F · (x_b − x_c)          com a, b, c distintos e ≠ i;  F = 0,8
```
> 🧠 Pega um indivíduo qualquer e o desloca na direção da diferença entre outros dois. É a sacada da DE: a "escala" da mutação se **auto-ajusta** — quando a população está espalhada, as diferenças são grandes e a busca é ampla; quando ela converge, as diferenças encolhem e a busca vira refino local. Não há parâmetro de "tamanho de passo" para você calibrar.

**2. Crossover binomial**
```
para cada gene j:  u[j] = v[j] se rand() < CR, senão x[j]      CR = 0,9
(j_rand garante que pelo menos um gene venha do mutante)
```

**3. Seleção gulosa**
```
se t₇₅(u) < t₇₅(xᵢ):  xᵢ ← u          senão:  xᵢ permanece
```
> A população **nunca piora**. É uma propriedade forte: o melhor indivíduo é monotonicamente não-crescente ao longo das gerações.

**4. Convergência**
```
desvio-padrão dos últimos 10 melhores < 1e-4  ⟹  para
```

**Custo:** cada avaliação de fitness = **uma simulação completa** do gêmeo digital. 20 × 30 = 600, mais o pós-processamento ≈ **620 simulações**. Daí as horas.

### 27.4 A escada de penalidades

O otimizador precisa de um valor de fitness mesmo quando o candidato é inviável. A escada:

| Situação | Fitness |
|---|---|
| exceção na simulação | 10⁹ |
| carro andou < 1 m | 10⁵ |
| não cruzou os 75 m | 100 + metros faltantes |
| cruzou | **t₇₅ real [s]** |

> 🧠 **Por que "100 + metros faltantes" e não um valor fixo?** Porque assim a penalidade é **informativa**: um candidato que andou 70 m tem fitness melhor que um que andou 20 m. Isso dá um **gradiente** para a DE seguir mesmo na região inviável, guiando a busca de volta para a região boa. Uma penalidade constante seria um platô — o algoritmo ficaria vagando cego.

### 27.5 Sem penalidade de slip (decisão documentada)

O otimizador **não** penaliza patinagem explicitamente. Motivo: o Pacejka calibrado já cobra o preço da patinação em tempo (§16.9) — patinar joga a curva na cauda descendente, a tração cai e o carro leva mais tempo. Penalizar de novo seria contar duas vezes.

> 🧠 **Um detalhe honesto do histórico:** a penalidade de slip que existia antes era **silenciosamente nula** por causa de um bug de chave de dicionário — ela nunca chegou a fazer efeito. Isso foi descoberto justamente ao validar a propriedade acima. É um bom exemplo de por que vale a pena entender *por que* algo funciona, e não só constatar que funciona.

### 27.6 Pós-processamento

**`find_integer_ratio(N*)`** — a DE devolve um N **contínuo** (ex.: 4,83). Mas você não pode fabricar uma coroa com 58,3 dentes. Esta função converte no par de dentes realizável mais próximo:

```
Pares de dentes mais próximos de N* = 4,83:
  Z1=12  Z2= 58  → N=4,8333  (erro 0,07 %)
  Z1=11  Z2= 53  → N=4,8182  (erro 0,24 %)
  Z1=13  Z2= 63  → N=4,8462  (erro 0,33 %)
  Z1=14  Z2= 68  → N=4,8571  (erro 0,56 %)
```

E cada candidato é **reavaliado** com o t_pico ótimo — porque nada garante que o melhor t_pico para N = 4,83 continue sendo o melhor para N = 4,8182.

**`sweep_diagnostico()`** — imprime uma grade ASCII de N × t_pico. Serve para **inspecionar o landscape antes** de rodar a DE.

> 🧠 **Faça isso.** Sempre. Antes de gastar horas numa otimização estocástica, olhe a superfície. Você vai ver se ela tem um mínimo claro, se tem vários mínimos locais, ou se tem um penhasco. Isso te diz se o resultado da DE é confiável e onde procurar.

### 27.7 Reprodutibilidade

- A **simulação** é determinística (sem sementes aleatórias). Mesma entrada, mesma saída, sempre.
- A **DE** usa `np.random` **sem seed**. Duas execuções dão resultados ligeiramente diferentes.

> ⚠️ Se você precisa reproduzir exatamente uma otimização (para um relatório, para debug), fixe `np.random.seed(...)` antes de rodar. Anote a seed junto com o resultado.

### 27.8 Uma nota de calibração histórica

> Um resultado anterior do otimizador devolvia **N ≈ 11,95**, o que era incompatível com a arquitetura RWD de motor central (times comparáveis usam 3–5:1; relações de 13–15:1 aparecem só em arquiteturas 4WD com motores de cubo e redutores planetários). A investigação revelou que `speed_ref` estava em 1200 rad/s em vez dos 575,96 rad/s reais do EMRAX 228 — ou seja, o otimizador estava explorando uma faixa de rotação que o motor físico não alcança, e "comprava" relação longa que na prática seria inatingível.
>
> **A lição:** quando um otimizador devolve um resultado que destoa do benchmark do setor, a hipótese mais provável **não** é que você descobriu algo genial. É que um limite físico está mal representado no modelo. Sempre valide o resultado contra o que times comparáveis fazem.

---

# Parte VI — Usando, criticando e estendendo o modelo

## 28. Como ler os resultados (e desconfiar deles)

### 28.1 O painel de seis gráficos

A simulação padrão produz este conjunto. Aprenda a ler cada um:

**① Velocidade do carro [km/h]**
> Cresce quase linearmente enquanto a tração está **saturada no pneu**, depois a taxa cai quando o teto de 80 kW passa a mandar. Se você vir um degrau ou uma quina, algo está errado.

**② Torques no eixo do motor [N·m]**
> `T_em = Kt·isq` sobe com a rampa do pedal até ~205 N·m. A **diferença** entre T_em e T_carga é o que acelera o rotor. Se as duas curvas se encostam, o rotor parou de acelerar (regime).

**③ Slip ratio κ**
> Pico na largada (a roda dispara com o carro parado) e convergência para valores pequenos conforme v cresce. A linha em κ = 0,204 marca o pico de tração. **Tempo acima dessa linha = tempo desperdiçado.**

**④ Correntes de referência [A]**
> A "corcova" inicial é o pedal; a queda ∝ 1/ω depois é o limitador de potência. `id_ref` fica em 0 até a velocidade-base (que neste carro nunca chega).

**⑤ Bateria — tensão e corrente DC**
> Vdc afunda quando I_bat sobe. Se Vdc encostar no `Vdc_floor` (159,8 V), **desconfie**: algo está errado (§25.1).

**⑥ Potência no barramento DC [kW]**
> Deve saturar **exatamente** em 80 kW, formando um platô. Se ultrapassar, o limitador não está funcionando. Se ficar bem abaixo, você não está usando todo o regulamento — provavelmente o pneu ou o pedal estão limitando antes.

### 28.2 Testes de sanidade — faça sempre

Antes de acreditar em qualquer resultado, verifique:

| ✓ | Verificação | Valor esperado |
|---|---|---|
| ☐ | Kt no print inicial | 0,7122 N·m/A |
| ☐ | J_eff no print inicial | ~0,07 kg·m² (ordem de grandeza) |
| ☐ | Vdc_nominal | 532,8 V |
| ☐ | Potência DC máxima | ≤ 80 kW, com platô |
| ☐ | Vdc mínimo | > 400 V (nunca perto do floor) |
| ☐ | Temperatura final | < 60 °C numa prova de 5 s |
| ☐ | SoC final | ~0,99 (a prova é curta) |
| ☐ | Torque de pico | ≤ 230 N·m |
| ☐ | Slip máximo | < 1,0, e não ficando cravado em 1,0 |
| ☐ | Velocidade final | 100–120 km/h |
| ☐ | Distância final | exatamente 75,00 m |

### 28.3 Sintomas e causas prováveis

| Sintoma | Causa provável |
|---|---|
| Temperatura absurda (centenas ou milhares de °C) | Corrente errada na bateria (§25.1) ou isd divergindo |
| Tensão de bateria negativa ou no floor | Escala pack→célula (§13.4) ou corrente errada |
| Torque oscilando violentamente | isd divergindo por falta de margem de tensão |
| Carro não sai do lugar | Tração < resistência; verifique λμx, Fz, relação N |
| Slip cravado em 1,0 a prova toda | Relação N muito longa, ou pedal-degrau, ou λμx muito baixo |
| Otimizador devolvendo N absurdo | `speed_ref` errado (§27.8) ou Pacejka sem pico (§16.5) |
| Simulação lenta demais | Normal: 30 s a 2 min. Se demorar muito mais, `tmax` está alto e o carro não está cruzando 75 m |
| Resultado muda entre execuções | A simulação é determinística. Se mudar, é o otimizador (sem seed) ou você mudou algo |

### 28.4 O ceticismo saudável

> 🧠 **Um modelo que concorda com você não está te ajudando.** O valor do gêmeo digital está em ele **discordar** de vez em quando, te forçando a descobrir quem está errado — ele ou sua intuição. Quando o resultado bate perfeitamente com o que você esperava, desconfie: talvez você tenha calibrado o modelo até ele dizer o que você queria ouvir.
>
> **Regra prática da equipe:** todo resultado do simulador que for usado numa decisão de projeto deve ser (a) validado contra um benchmark externo (outras equipes, datasheets, literatura) ou (b) validado em pista assim que possível.

---

## 29. Limitações conhecidas — o que este modelo NÃO faz

Ser explícito sobre limitações é o que separa um modelo de engenharia de um chute com gráficos bonitos.

### 29.1 Física ausente

| Limitação | Impacto | Quando importaria |
|---|---|---|
| **Sem regeneração** | só potência positiva sai da bateria | Endurance, Autocross (frenagens) |
| **Inversor ideal** | sem PWM, harmônicos, perdas de comutação; η = 0,90 agregada | análise de eficiência fina, EMI, ripple de torque |
| **Pneu só longitudinal** | sem forças laterais, sem elipse de atrito | qualquer prova com curva |
| **Pneu sem relaxamento** | a força aparece instantaneamente com o slip | transientes muito rápidos (< 20 ms) |
| **Sem dinâmica de suspensão** | transferência de carga é instantânea e rígida | pitch dinâmico na largada, pista irregular |
| **LSD inativo** | parâmetro documentado mas sem efeito em reta | curvas, tração assimétrica |
| **Térmica de massa única** | não separa cobre/ferro/carcaça | Endurance, gestão térmica |
| **Sem modelo de BMS** | sem limites de corrente/temperatura por célula | dimensionamento de acumulador |
| **Sem escorregamento de embreagem/corrente** | trem de força rígido | análise de vibração torcional |
| **Pista perfeita** | plana, aderência uniforme, sem vento | condições reais de competição |

### 29.2 Simplificações numéricas

- **2 iterações de ponto fixo** no loop Fz ↔ aceleração (§15.3) — erro residual ~1 %;
- **`V_nominal` no denominador** do balanço de potência em vez da tensão terminal (§20.2) — erro da ordem do sag relativo (~5 %);
- **`L0` estimado** como 10 % da média de Ld/Lq (o datasheet não fornece) — irrelevante, já que iso decai a zero.

### 29.3 Detalhes de implementação a saber

- Os `print` da Simulation usam caracteres Unicode (🚀, ζ, ≈). Em consoles Windows com codepage cp1252 é preciso `PYTHONIOENCODING=utf-8` para não levantar `UnicodeEncodeError`.
- `Models/`, `Simulation/` e `Constants/` funcionam como **namespace packages** (sem `__init__.py`) — as importações são feitas módulo a módulo.
- Nomes internos ≠ chaves públicas do dict de retorno (§10).
- Defaults do dashboard e `SIM_CONFIG` do otimizador precisam de sincronização **manual** (§26.5).

### 29.4 Onde o modelo pode ser estendido (ideias de projeto)

Em ordem crescente de esforço:

1. **Fácil** — adicionar a química real da célula ao catálogo do `BatteryPack`; adicionar perfil de pedal de dois estágios; adicionar Skidpad como critério de parada alternativo.
2. **Médio** — modelo térmico de dois nós (cobre + carcaça); modelo de inversor com perdas de comutação; validação contra dados de telemetria de pista.
3. **Difícil** — pneu combinado (elipse de atrito) para simular Autocross; dinâmica de suspensão; simulação de Endurance com gestão de energia.

---

## 30. Exercícios — e as respostas

Faça estes com o código aberto. Levam de 10 a 40 minutos cada.

### 30.1 Os exercícios do texto

**Exercício 1.** Nosso motor tem torque de pico 230 N·m. A que rotação ele atinge 80 kW mecânicos?

<details><summary>Resposta</summary>

`ω = P/T = 80.000/230 = 347,8 rad/s = 3.321 RPM`.
Acima disso, o regulamento (não o motor) passa a limitar o torque. É a velocidade-base **de potência** — não confundir com a velocidade-base de **tensão** do field weakening (9.644 RPM), que é outro fenômeno.
</details>

**Exercício 2.** Calcule (a) o torque de pico do motor e (b) a velocidade-base de tensão em RPM.

<details><summary>Resposta</summary>

(a) `Kt = 1,5 × 10 × 0,04748 = 0,7122 N·m/A`; `T = 0,7122 × 323 = 230,0 N·m` ✓ (datasheet)
(b) `ωe_base = 0,9 × 532,8/0,04748 = 10.099 rad/s elétricos`; `ωm = 10.099/10 = 1.010 rad/s = 9.644 RPM`.
Como o limite mecânico é 5.500 RPM, o field weakening nunca ativa nesta configuração.
</details>

**Exercício 3.** Tensão terminal do pack com SoC = 1,0, I = 150 A, I* = 0.

<details><summary>Resposta</summary>

Com SoC = 1 → q = 0, então `term_carga = 0`, `term_Iast = 0` e `A·e⁰ = A = 0,1`.
`i_cell = 150/5 = 30 A`
`Vt_célula = 3,7 − 0,02×30 − 0 − 0 + 0,1 = 3,2 V`
`Vt_pack = 3,2 × 144 = 460,8 V`
Confira no gráfico da documentação original: a curva de 150 A começa em ~460 V. ✓
</details>

**Exercício 4.** Velocidade do carro a 5500 RPM com N = 4,0 e com N = 5,0.

<details><summary>Resposta</summary>

N = 4: `ω_roda = 575,96/4 = 144,0 rad/s`; `v = 144,0 × 0,26 = 37,4 m/s = 134,8 km/h`
N = 5: `ω_roda = 115,2 rad/s`; `v = 29,9 m/s = 107,8 km/h`

**A lição:** N maior = mais torque na roda (aceleração melhor) mas menor velocidade máxima. Para uma prova de 75 m, você nunca chega perto de 135 km/h — então N = 5 "desperdiça" menos velocidade máxima do que parece. É por isso que o ótimo da prova de aceleração fica em N maior do que a intuição de "carro rápido" sugere.
</details>

**Exercício 5.** Fz por roda traseira para a ∈ {0, 5, 10, 15} m/s², a v = 0.

<details><summary>Resposta</summary>

`Fz = [m·g·d_cg/L + m·a·h/L]/2 = [942 + 44,8·a]/2`

| a [m/s²] | Fz por roda [N] | ganho |
|---|---|---|
| 0 | 471 | — |
| 5 | 583 | +24 % |
| 10 | 695 | +48 % |
| 15 | 807 | +71 % |

Como Fx é aproximadamente proporcional a Fz, **acelerar forte gera mais tração disponível** — um efeito de realimentação positiva que ajuda na largada. É por isso que o loop de ponto fixo (§15.3) é necessário.
</details>

**Exercício 6.** Plote `Tire_forces(654, κ)` para κ ∈ [0, 1] com os dois datasets.

<details><summary>Resposta</summary>

Com `HOOSIER_FSAE_LONG`: pico em κ ≈ 0,204, caindo ~11 % em κ = 0,5 e ~20 % em κ = 1,0.
Com `HOOSIER_FSAE_LONG_CALSPAN`: curva **monotonicamente crescente** — Fx máximo em κ = 1,0.

Por que o segundo é impossível: significaria que **um pneu girando no vazio transmite mais força que um pneu em contato ótimo**. Se fosse verdade, controle de tração seria uma tecnologia inútil e todo piloto deveria queimar pneu na largada.
</details>

### 30.2 Exercícios adicionais

**Exercício 7 (fácil).** Rode a simulação padrão. Depois mude só `t_pico` de 3,0 para 1,0 s. Quanto o tempo melhora? Olhe o gráfico de slip antes e depois — o que mudou?

**Exercício 8 (fácil).** Rode com `p_max_dc = None` (sem limite de potência). Quanto o carro ganharia se o regulamento não existisse? Qual a potência de pico que ele usaria?

**Exercício 9 (médio).** Varra `tire_friction_coef` de 0,4 a 0,8 em passos de 0,1 e plote t₇₅ em função de λμx. A relação é linear? O que isso te diz sobre a importância de escolher bem os pneus?

**Exercício 10 (médio).** Coloque `has_wing=False`. Quanto tempo o carro perde? Agora aumente só `lift_coeff_rear` — em que ponto o ganho de tração começa a ser compensado pelo arrasto? *(Dica: o Cd atual já embute o arrasto induzido das asas, então esse exercício é aproximado. Pense em como você melhoraria o modelo para responder isso direito.)*

**Exercício 11 (médio).** Rode `sweep_diagnostico()` do otimizador. Descreva o formato do landscape N × t_pico. Onde estão os penhascos? Faz sentido rodar uma DE nessa superfície, ou uma busca em grade bastaria?

**Exercício 12 (difícil).** Comente a linha do filtro passa-baixa (`self.wm_filt += ...`) e use `wm` direto nos limitadores. Rode e olhe o gráfico de `iq_ref` de perto na região do joelho do limitador. Consegue ver o ciclo-limite que motivou o filtro?

**Exercício 13 (difícil).** Reintroduza o bug: passe `abs(isq)` para o modelo de bateria em vez de `I_bat`. Rode e documente a cascata de destruição da §25.1. **Depois desfaça.** *(Fazer o bug acontecer de propósito, uma vez, é a melhor forma de nunca mais reintroduzi-lo por acidente.)*

**Exercício 14 (projeto).** Implemente uma penalidade de temperatura no otimizador e veja se ela muda o ótimo. Depois argumente se ela deveria existir, dado o que você sabe sobre a duração da prova.

---

## 31. Glossário

### Símbolos

| Símbolo | Nome | Unidade | Onde aparece |
|---|---|---|---|
| ω, ωm | velocidade angular mecânica | rad/s | mecânica do motor |
| ωe | velocidade angular elétrica | rad/s | `ωe = p·ωm` |
| θm, θe | ângulo mecânico / elétrico | rad | transformadas |
| p | pares de polos | — | 10 no EMRAX 228 |
| λm | fluxo dos ímãs | Wb | 0,04748 |
| Kt | constante de torque | N·m/A | `1,5·p·λm = 0,7122` |
| isd, isq | correntes de eixo d e q | A (pico) | FOC |
| vd, vq | tensões de eixo d e q | V | comandos do inversor |
| Ld, Lq | indutâncias de eixo d e q | H | 96,5 μH |
| Rs | resistência de estator | Ω | 0,00706 (a 25 °C) |
| J, Jm, J_eff | inércias | kg·m² | rotor / efetiva |
| kf | atrito viscoso | N·m·s | 0,005 |
| N | relação final de transmissão | — | Z_coroa/Z_pinhão |
| η | eficiência | — | 0,9506 (transmissão) |
| κ (kappa) | slip ratio | — | ∈ [−1, +1] |
| Fz | carga vertical no pneu | N | por roda |
| Fx | força longitudinal do pneu | N | de um pneu |
| μx, λμx | atrito de pico / fator de escala | — | Pacejka |
| ρ (rho) | densidade do ar | kg/m³ | 1,225 |
| Cd, Cl | coef. de arrasto / sustentação | — | aerodinâmica |
| SoC | estado de carga | 0–1 | bateria |
| I* (Iast) | corrente acumulada ∫i·dt | C | polarização |
| Q | capacidade de célula | Ah | 3,0 |
| E0 | tensão de circuito aberto | V | 3,7 por célula |

### Termos

| Termo | Significado |
|---|---|
| **Anti-windup** | Técnica para impedir que o integrador de um PI acumule comando que não chega à planta |
| **Back-calculation** | Forma de anti-windup: devolve ao integrador o excesso cortado pela saturação |
| **Back-EMF** | Força contra-eletromotriz — tensão induzida pela rotação do rotor, `ωe·λm` |
| **Barramento DC** | O circuito de corrente contínua entre a bateria e o inversor |
| **Clarke (transformada)** | abc → αβ (três eixos para dois, ainda estacionários) |
| **Derate** | Redução deliberada de um limite para proteger um componente (aqui, térmica) |
| **DE** | Evolução Diferencial — algoritmo de otimização estocástico sem derivadas |
| **Digital twin** | Gêmeo digital — modelo computacional fiel de um sistema físico |
| **Downforce** | Força aerodinâmica para baixo, gerada por asas — aumenta a aderência |
| **EDO / ODE** | Equação diferencial ordinária |
| **Feedforward** | Compensação antecipada de uma perturbação conhecida, sem esperar o erro aparecer |
| **Field weakening** | Injeção de corrente negativa no eixo d para reduzir o fluxo efetivo e ganhar margem de tensão |
| **FOC** | *Field-Oriented Control* — controle vetorial no referencial dq |
| **FSAE / Formula Student** | Competição universitária de projeto e construção de carros de corrida |
| **Ganho de malha** | Quanto o controlador amplifica o erro |
| **Landscape** | A "paisagem" da função objetivo num problema de otimização |
| **Magic Formula** | Modelo empírico de pneu de Hans Pacejka |
| **Massa concentrada (lumped)** | Simplificação térmica: tratar tudo como um bloco à mesma temperatura |
| **MTPA** | *Maximum Torque Per Ampere* — estratégia de referência de corrente ótima |
| **Park (transformada)** | αβ → dq (referencial estacionário para referencial girante) |
| **PMSM** | *Permanent Magnet Synchronous Motor* — motor síncrono de ímãs permanentes |
| **Ponto fixo (iteração de)** | Método de resolver um loop algébrico chutando e refinando |
| **Relutância (torque de)** | Torque extra em motores com Ld ≠ Lq. Ausente no nosso motor |
| **RK4** | Runge-Kutta de 4ª ordem — integrador numérico |
| **RWD** | *Rear-Wheel Drive* — tração traseira |
| **Sag** | Afundamento de tensão da bateria sob carga |
| **Shepherd (modelo de)** | Modelo empírico de tensão de bateria em função do SoC |
| **Slip ratio** | Escorregamento relativo entre a roda e o solo |
| **Smoothstep** | Curva `3s²−2s³` com derivada nula nas extremidades |
| **SoC** | *State of Charge* — estado de carga da bateria |
| **Steinmetz** | Forma empírica das perdas no ferro: histerese + Foucault |
| **TTC / Calspan** | *Tire Test Consortium* — bancada de ensaio de pneus |
| **Velocidade-base** | Rotação acima da qual o back-EMF consome toda a margem de tensão |
| **Wheelspin** | Patinagem da roda motriz |
| **Windup** | Saturação do integrador de um PI durante saturação de saída |
| **ZOH** | *Zero-Order Hold* — retentor de ordem zero, comando congelado entre amostras |

---

## 32. FAQ e primeiros passos

### 32.1 Perguntas frequentes

**"Preciso saber tudo isso para contribuir?"**
Não. Se você vai usar o dashboard, leia a Parte 0, §8 e a Parte V. Se você vai mexer em **um** módulo, leia a Parte I, §9 e a seção daquele módulo. A Parte IV só é obrigatória se você for mexer no núcleo.

**"Por que tanto código em português misturado com inglês?"**
Herança. O projeto nasceu em português (`corrented`, `conjugado`, `calcular_tensao`) e foi ganhando nomes em inglês nas partes novas. **Não "padronize" isso sem combinar com a equipe** — renomear atributos quebra o dashboard, o otimizador e qualquer script externo. Se for fazer, faça de uma vez, com testes.

**"A simulação está lenta. Posso aumentar o passo?"**
Cuidado. Você pode testar `HP = 2e-4` (5 kHz) e comparar t₇₅ com o resultado de referência. Se a diferença for < 0,5 %, provavelmente é seguro para varreduras exploratórias. Mas **valide sempre com o passo original** antes de usar um resultado numa decisão. E lembre: o passo também define a taxa do controlador digital, então mudá-lo muda o sistema simulado, não só a precisão numérica.

**"Meu resultado não bate com o da documentação."**
Cheque, nesta ordem: (1) a relação N — muitos exemplos usam N = 4,0 e a simulação padrão usa N = 5,0 (§15.5); (2) o perfil de pedal (t_pico = 3 s ou 1 s?); (3) o dataset de pneu; (4) `speed_ref`.

**"Posso simular Autocross/Endurance com isso?"**
Não como está. O pneu só tem dinâmica longitudinal (§29.1). Estender para forças combinadas é um bom projeto de extensão, mas é trabalho sério.

**"Onde eu começo se quero contribuir?"**
Faça os Exercícios 7 a 11. Eles te dão familiaridade com a ferramenta sem risco de quebrar nada. Depois escolha uma das extensões "fáceis" da §29.4.

### 32.2 Checklist de primeiro dia

```
☐ Clonar o repositório
☐ Instalar dependências (numpy, dash, plotly, dash-bootstrap-components)
☐ Rodar `python main.py` e abrir http://127.0.0.1:8050/
☐ Clicar em "Executar Simulação" sem mudar nada
☐ Conferir os testes de sanidade da §28.2
☐ Passear pelas 10 abas e identificar cada gráfico usando a §28.1
☐ Mudar UM parâmetro, rodar de novo, olhar a aba Histórico
☐ Fazer os Exercícios 7 e 8
```

### 32.3 Antes de mandar um Pull Request

```
☐ Li a §25 (bugs históricos) e confirmei que não estou regredindo nenhum
☐ Rodei a simulação padrão e comparei com a referência (t₇₅ = 5,104 s)
☐ Conferi os testes de sanidade da §28.2
☐ Se mexi num modelo, atualizei a docstring do módulo
☐ Se mudei um default do dashboard, sincronizei o SIM_CONFIG do otimizador (§26.5)
☐ Se mudei uma decisão da tabela §25.8, documentei o porquê
```

---

## Agradecimentos

Do autor do projeto:

- **Leonardo Gonçalves** — pela parceria no desenvolvimento e pelas contribuições às primeiras versões do simulador do powertrain;
- **Igor Maia** — pelo apoio, pelas discussões técnicas e pelas contribuições ao programa de powertrain da equipe;
- E a toda a equipe **Xangô E-Racing**, cujo trabalho torna este projeto possível. 🏁

---

*Documento didático v1.0 · baseado na Documentação Técnica DigitalTwin Aceleration v0.1.5*
*Xangô E-Racing · Escola Politécnica · UFBA*
