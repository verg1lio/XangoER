# -*- coding: utf-8 -*-
"""
V2_Dimensionamento_Freio_Estilo_Exemplo_v5.py

Este script é focado no DIMENSIONAMENTO (cálculo e dimensionamento) 
do sistema de freios.

ATUALIZAÇÃO:
- Corrigidos os erros de digitação (f/r trocados) nas funções de 
  hidráulica, torque e verificação.
- Adicionado um "PASSO 3" dedicado a printar apenas os resultados ESTÁTICOS (a=0).
- O "PASSO 5" (Verificação) agora calcula a desaceleração alvo usando 
  a fórmula explícita do usuário: a = (Ff + Fr) / mv, 
  onde Ff/Fr = W_static * mu.
"""

import numpy as np
import sys

# Importa a classe BrakeSystem do V1 para "ligar" os códigos,
# conforme solicitado, e evitar repetições.
try:
    from V1 import BrakeSystem
except ImportError:
    print("Aviso: Não foi possível encontrar o arquivo V1.PY.")
    print("Algumas funcionalidades podem não estar disponíveis ou podem gerar erros.")
    # Define uma classe placeholder caso V1.PY não esteja acessível
    class BrakeSystem:
        def __init__(self, params):
            print("Classe BrakeSystem (Placeholder) inicializada.")
            self.params = params

# Constante
G = 9.81  # Aceleração da gravidade (m/s²)

class BrakeSizing:
    """
    Classe dedicada ao dimensionamento e verificação dos componentes 
    mecânicos e hidráulicos do sistema de freio.
    
    A inicialização segue o padrão exato do exemplo fornecido.
    """
    
    def __init__(self, params):
        """
        Inicializa a classe desempacotando todos os parâmetros 
        do dicionário 'params' para atributos 'self'.
        """
        self.params = params
        
        # --- Parâmetros do Veículo (Vehicle) ---
        self.massa_total = params['massa_total']
        self.HCG = params['HCG']
        self.entre_eixos = params['entre_eixos']
        self.dist_cg_dianteiro = params['dist_cg_dianteiro']
        self.dist_cg_traseiro = params['dist_cg_traseiro']
        self.raio_pneu = params['raio_pneu']
        self.mu_pneu = params['mu_pneu']

        # --- Parâmetros de Freio (Brake Sizing) ---
        self.L_pedal = params['L_pedal']
        self.L_pivot_cm = params['L_pivot_cm']
        self.dist_bb_dianteiro = params['dist_bb_dianteiro']
        self.dist_bb_traseiro = params['dist_bb_traseiro']
        self.D_cm = params['D_cm']
        self.D_pin_f = params['D_pin_f']
        self.N_pistoes_f_um_lado = params['N_pistoes_f_um_lado']
        self.D_pin_r = params['D_pin_r']
        self.N_pistoes_r_um_lado = params['N_pistoes_r_um_lado']
        self.tipo_pinca_f = params['tipo_pinca_f']
        self.tipo_pinca_r = params['tipo_pinca_r']
        self.mu_pastilha = params['mu_pastilha']
        self.raio_efetivo_disco_f = params['raio_efetivo_disco_f']
        self.raio_efetivo_disco_r = params['raio_efetivo_disco_r']
        
        # Dicionário para armazenar os resultados dos cálculos
        self.results = {}


    def calcular_relacao_pedal(self, F_piloto_alvo_N):
        """
        Calcula o ganho mecânico do pedal (Relação de Pedal).
        """
        if self.L_pivot_cm == 0:
            raise ValueError("Distância do pivot ao cilindro mestre (L_pivot_cm) não pode ser zero.")
            
        pedal_ratio = self.L_pedal / self.L_pivot_cm
        F_saida_pedal = pedal_ratio * F_piloto_alvo_N
        
        self.results['pedal_ratio'] = pedal_ratio
        self.results['F_piloto_aplicada_N'] = F_piloto_alvo_N
        self.results['F_saida_pedal_N'] = F_saida_pedal
        
        print(f"[Pedal] Relação de Pedal (Ganho): {pedal_ratio:.2f}:1")
        print(f"[Pedal] Força do Piloto Aplicada: {F_piloto_alvo_N:.1f} N")
        print(f"[Pedal] Força de Saída (na Balance Bar): {F_saida_pedal:.1f} N")
        
        return F_saida_pedal

    def calcular_distribuicao_balance_bar(self, F_saida_pedal_N):
        """
        Calcula a distribuição de força entre os cilindros mestres.
        """
        L_total_bb = self.dist_bb_dianteiro + self.dist_bb_traseiro
        if L_total_bb == 0:
            raise ValueError("Distância total da Balance Bar não pode ser zero.")

        # Força no cilindro mestre TRASEIRO (pivot no dianteiro)
        F_cm_traseiro = F_saida_pedal_N * (self.dist_bb_dianteiro / L_total_bb)
        # Força no cilindro mestre DIANTEIRO (pivot no traseiro)
        F_cm_dianteiro = F_saida_pedal_N * (self.dist_bb_traseiro / L_total_bb)
        
        # Evita divisão por zero se a força for 0
        soma_forcas = F_cm_dianteiro + F_cm_traseiro
        distribuicao_frente_pct = (F_cm_dianteiro / soma_forcas) * 100 if soma_forcas > 0 else 0

        self.results['F_cm_dianteiro_N'] = F_cm_dianteiro
        self.results['F_cm_traseiro_N'] = F_cm_traseiro
        self.results['distribuicao_frenagem_pct'] = distribuicao_frente_pct
        
        # Print de debug que você tinha colocado, agora corrigido
        print(f"Força CM Traseiro: {F_cm_traseiro:.1f} N, Força CM Dianteiro: {F_cm_dianteiro:.1f} N")

        print(f"[Balance Bar] Distribuição de Frenagem: {distribuicao_frente_pct:.1f}% Dianteira / {100-distribuicao_frente_pct:.1f}% Traseira")
        
    def calcular_transferencia_carga(self, desaceleracao_ms2):
        """
        Calcula a transferência de carga longitudinal durante a frenagem.
        """
        if self.entre_eixos == 0:
            raise ValueError("Entre-eixos (L) não pode ser zero.")
        if (self.dist_cg_dianteiro + self.dist_cg_traseiro) != self.entre_eixos:
            print(f"Aviso: Soma das distâncias do CG ({self.dist_cg_dianteiro + self.dist_cg_traseiro}m) não bate com entre-eixos ({self.entre_eixos}m).")

        peso_total_N = self.massa_total * G
        W_f_static_N = peso_total_N * (self.dist_cg_traseiro / self.entre_eixos)
        W_r_static_N = peso_total_N * (self.dist_cg_dianteiro / self.entre_eixos)
        
        delta_W_N = (self.massa_total * desaceleracao_ms2 * self.HCG) / self.entre_eixos
        
        W_f_dynamic_N = W_f_static_N + delta_W_N
        W_r_dynamic_N = W_r_static_N - delta_W_N
        
        if W_r_dynamic_N < 0:
            print("AVISO: Carga traseira negativa. Risco de capotamento (stoppie)!")
            W_r_dynamic_N = 0 
            
        self.results['W_f_static_N'] = W_f_static_N
        self.results['W_r_static_N'] = W_r_static_N
        self.results['W_f_dynamic_N'] = W_f_dynamic_N
        self.results['W_r_dynamic_N'] = W_r_dynamic_N
        
        if desaceleracao_ms2 > 0: # Não printar para o cálculo estático (a=0)
            print(f"[Carga] Desaceleração Alvo: {desaceleracao_ms2/G:.2f} G")
            print(f"[Carga] Transferência de Carga (Delta W): {delta_W_N:.1f} N")
            # --- CORRIGIDO (Prints estavam trocados no seu código base) ---
            print(f"[Carga] Dinâmica Dianteira: {W_f_dynamic_N:.1f} N (Estática: {W_f_static_N:.1f} N)")
            print(f"[Carga] Dinâmica Traseira: {W_r_dynamic_N:.1f} N (Estática: {W_r_static_N:.1f} N)")
        
        return W_f_dynamic_N, W_r_dynamic_N

    def calcular_pressao_hidraulica(self):
        """
        Calcula as pressões de linha com base na força da balance bar 
        e no diâmetro do cilindro mestre.
        """
        F_cm_dianteiro = self.results['F_cm_dianteiro_N']
        F_cm_traseiro = self.results['F_cm_traseiro_N']
        
        A_cm = np.pi * (self.D_cm / 2)**2
        
        # --- CORRIGIDO (Variáveis trocadas no seu código base) ---
        # Pressão Dianteira usa Força Dianteira
        P_linha_f_bar = (F_cm_dianteiro / A_cm) * 1e-5
        # Pressão Traseira usa Força Traseira
        P_linha_r_bar = (F_cm_traseiro / A_cm) * 1e-5
        
        self.results['A_cm_m2'] = A_cm
        self.results['P_linha_f_bar'] = P_linha_f_bar
        self.results['P_linha_r_bar'] = P_linha_r_bar
        
        print(f"[Hidráulica] Área Cilindro Mestre: {A_cm * 1e6:.1f} mm² (D={self.D_cm*1000:.1f} mm)")
        # --- CORRIGIDO (Prints trocados no seu código base) ---
        print(f"[Hidráulica] Pressão Linha Dianteira: {P_linha_f_bar:.2f} bar")
        print(f"[Hidráulica] Pressão Linha Traseira: {P_linha_r_bar:.2f} bar")

    def _calcular_forca_fechamento(self, P_linha_Pa, D_pin, N_pistoes_um_lado, tipo_pinca):
        """Helper para calcular força de fechamento (clamping force)."""
        A_pin_um_lado = N_pistoes_um_lado * (np.pi * (D_pin / 2)**2)
        F_hidraulica = P_linha_Pa * A_pin_um_lado
        
        if tipo_pinca == 'flutuante':
            F_fec_N = F_hidraulica * 2
        elif tipo_pinca == 'fixa':
            F_fec_N = F_hidraulica
        else:
            raise ValueError(f"Tipo de pinça '{tipo_pinca}' desconhecido. Use 'flutuante' or 'fixa'.")
            
        return F_fec_N, A_pin_um_lado

    def calcular_torque_frenagem(self):
        """
        Calcula o torque de frenagem gerado em cada disco.
        """
        P_linha_f_Pa = self.results['P_linha_f_bar'] * 1e5
        P_linha_r_Pa = self.results['P_linha_r_bar'] * 1e5
        
        # Dianteira
        F_fec_f, A_pin_f = self._calcular_forca_fechamento(
            P_linha_f_Pa,
            self.D_pin_f,
            self.N_pistoes_f_um_lado,
            self.tipo_pinca_f
        )
        
        # Traseira
        F_fec_r, A_pin_r = self._calcular_forca_fechamento(
            P_linha_r_Pa,
            self.D_pin_r,
            self.N_pistoes_r_um_lado,
            self.tipo_pinca_r
        )
        
        # Torque no Disco = 2 * F_fec * mu_pastilha * R_efetivo
        T_disco_f = 2 * F_fec_f * self.mu_pastilha * self.raio_efetivo_disco_f
        T_disco_r = 2 * F_fec_r * self.mu_pastilha * self.raio_efetivo_disco_r
        
        self.results['F_fec_f_N'] = F_fec_f
        self.results['F_fec_r_N'] = F_fec_r
        self.results['A_pin_f_total_um_lado_m2'] = A_pin_f
        self.results['A_pin_r_total_um_lado_m2'] = A_pin_r
        self.results['T_disco_f_Nm'] = T_disco_f
        self.results['T_disco_r_Nm'] = T_disco_r

        # --- CORRIGIDO (Prints trocados no seu código base) ---
        print(f"[Torque] Força Fechamento Diant: {F_fec_f:.1f} N")
        print(f"[Torque] Força Fechamento Tras: {F_fec_r:.1f} N")
        print(f"[Torque] Torque Disco Diant (por roda): {T_disco_f:.1f} Nm")
        print(f"[Torque] Torque Disco Tras (por roda): {T_disco_r:.1f} Nm")

    def calcular_desaceleracao_eficiencia_max(self):
        """
        Calcula a desaceleração máxima teórica e atingível,
        e a eficiência de frenagem.
        """
        B = self.results['distribuicao_frenagem_pct'] / 100.0
        W_f_static_N = self.results['W_f_static_N']
        W_r_static_N = self.results['W_r_static_N']

        a_max_teorica_ms2 = self.mu_pneu * G
        
        den_f = self.massa_total * (B - (self.mu_pneu * self.HCG / self.entre_eixos))
        den_r = self.massa_total * ((1 - B) + (self.mu_pneu * self.HCG / self.entre_eixos))

        a_lock_f_ms2 = (self.mu_pneu * W_f_static_N) / den_f if den_f > 0 else float('inf')
        a_lock_r_ms2 = (self.mu_pneu * W_r_static_N) / den_r if den_r > 0 else float('inf')
            
        a_max_atingivel_ms2 = min(a_max_teorica_ms2, a_lock_f_ms2, a_lock_r_ms2)
        
        if a_max_atingivel_ms2 == a_lock_f_ms2:
            limite = "Travamento Dianteiro"
        elif a_max_atingivel_ms2 == a_lock_r_ms2:
            limite = "Travamento Traseiro"
        else:
            limite = "Aderência Total (Ideal)"
            
        eficiencia_pct = (a_max_atingivel_ms2 / a_max_teorica_ms2) * 100

        self.results['a_max_teorica_ms2'] = a_max_teorica_ms2
        self.results['a_max_atingivel_ms2'] = a_max_atingivel_ms2
        self.results['eficiencia_frenagem_pct'] = eficiencia_pct
        self.results['limite_frenagem'] = limite

        print(f"[Desaceleração] Máx. Teórica (Ideal): {a_max_teorica_ms2/G:.3f} G")
        print(f"[Desaceleração] Limite Travamento Diant: {a_lock_f_ms2/G:.3f} G")
        print(f"[Desaceleração] Limite Travamento Tras: {a_lock_r_ms2/G:.3f} G")
        print(f"[Desaceleração] Máx. Atingível (Real): {a_max_atingivel_ms2/G:.3f} G (Limite: {limite})")
        print(f"[Eficiência] Eficiência de Frenagem: {eficiencia_pct:.1f}%")

    def verificar_forca_piloto_para_desaceleracao(self, desaceleracao_alvo_g):
        """
        Calcula a força que o piloto precisa aplicar para atingir 
        uma desaceleração alvo. (Cálculo inverso)
        
        :param desaceleracao_alvo_g: Desaceleração desejada (em G's).
        """
        a_alvo_ms2 = desaceleracao_alvo_g * G
        
        # 1. Carga dinâmica
        W_f_dyn, W_r_dyn = self.calcular_transferencia_carga(a_alvo_ms2)
        
        # 2. Torque máximo (limite) no pneu
        T_pneu_max_f = W_f_dyn * self.mu_pneu * self.raio_pneu
        T_pneu_max_r = W_r_dyn * self.mu_pneu * self.raio_pneu
        
        # 3. Torque de freio necessário
        F_total_necessaria = self.massa_total * a_alvo_ms2
        
        B_default = self.params.get('distribuicao_frenagem_pct_fixa', 60.0)
        B = self.results.get('distribuicao_frenagem_pct', B_default) / 100.0
        
        F_frenagem_f_necessaria = F_total_necessaria * B
        F_frenagem_r_necessaria = F_total_necessaria * (1.0 - B)
        
        T_disco_f_necessario = (F_frenagem_f_necessaria / 2) * self.raio_pneu
        T_disco_r_necessario = (F_frenagem_r_necessaria / 2) * self.raio_pneu

        if T_disco_f_necessario > T_pneu_max_f:
            print(f"AVISO: Eixo dianteiro travará antes de atingir {desaceleracao_alvo_g:.3f} G.")
            T_disco_f_necessario = T_pneu_max_f
        if T_disco_r_necessario > T_pneu_max_r:
            print(f"AVISO: Eixo traseiro travará antes de atingir {desaceleracao_alvo_g:.3f} G.")
            T_disco_r_necessario = T_pneu_max_r

        # 4. Força de fechamento (F_fec) necessária
        F_fec_f_nec = T_disco_f_necessario / (2 * self.mu_pastilha * self.raio_efetivo_disco_f) if self.raio_efetivo_disco_f > 0 else 0
        F_fec_r_nec = T_disco_r_necessario / (2 * self.mu_pastilha * self.raio_efetivo_disco_r) if self.raio_efetivo_disco_r > 0 else 0

        # 5. Pressão de linha necessária
        # (Areas A_pin_f/r são calculadas no Passo 2, 'calcular_torque_frenagem')
        A_pin_f = self.results['A_pin_f_total_um_lado_m2']
        A_pin_r = self.results['A_pin_r_total_um_lado_m2']
        
        F_hidraulica_f_nec = F_fec_f_nec / 2 if self.tipo_pinca_f == 'flutuante' else F_fec_f_nec
        F_hidraulica_r_nec = F_fec_r_nec / 2 if self.tipo_pinca_r == 'flutuante' else F_fec_r_nec
        
        P_linha_f_nec_Pa = F_hidraulica_f_nec / A_pin_f if A_pin_f > 0 else 0
        P_linha_r_nec_Pa = F_hidraulica_r_nec / A_pin_r if A_pin_r > 0 else 0

        # 6. Força no Cilindro Mestre necessária
        A_cm = self.results['A_cm_m2']
        F_cm_f_nec = P_linha_f_nec_Pa * A_cm
        F_cm_r_nec = P_linha_r_nec_Pa * A_cm
        
        # 7. Força no pedal necessária
        # (Usa a força dianteira como base para o cálculo)
        F_saida_pedal_nec = F_cm_f_nec / (B) if B > 0 else 0
        
        # 8. Força do piloto necessária
        pedal_ratio = self.results['pedal_ratio']
        F_piloto_necessaria_N = F_saida_pedal_nec / pedal_ratio if pedal_ratio > 0 else 0

        print(f"--- VERIFICAÇÃO PARA {desaceleracao_alvo_g:.3f} G (Meta) ---")
        # --- CORRIGIDO (Prints trocados no seu código base) ---
        print(f"Torque Diant. Necessário (por roda): {T_disco_f_necessario:.1f} Nm (Max Aderência: {T_pneu_max_f:.1f} Nm)")
        print(f"Torque Tras. Necessário (por roda): {T_disco_r_necessario:.1f} Nm (Max Aderência: {T_pneu_max_r:.1f} Nm)")
        print(f"Pressão Linha Diant. Necessária: {P_linha_f_nec_Pa*1e-5:.2f} bar")
        print(f"Pressão Linha Tras. Necessária: {P_linha_r_nec_Pa*1e-5:.2f} bar")
        print(f"Força CM Diant. Necessária: {F_cm_f_nec:.1f} N")
        print(f"Força CM Tras. Necessária: {F_cm_r_nec:.1f} N")
        print(f"FORÇA DO PILOTO NECESSÁRIA (para {desaceleracao_alvo_g:.3f} G): {F_piloto_necessaria_N:.1f} N")
        print(f"Relação força necessária por desaceleração {F_piloto_necessaria_N / desaceleracao_alvo_g:.3f} N/G")
        
        if desaceleracao_alvo_g == 1.0: # Mantém a checagem de 1G se o alvo for 1G
            if 445 <= F_piloto_necessaria_N <= 489:
                print("Resultado (1G): ÓTIMO (Dentro da faixa ideal de 445-489 N para 1G)")
            elif F_piloto_necessaria_N < 445:
                print("Resultado (1G): ACEITÁVEL (Força abaixo do ideal, pode ser sensível demais)")
            else:
                print("Resultado (1G): RUIM (Força excessiva, acima de 489 N para 1G)")
        if desaceleracao_alvo_g == desaceleracao_alvo_g:
            if F_piloto_necessaria_N / desaceleracao_alvo_g > 263 and F_piloto_necessaria_N / desaceleracao_alvo_g < 445:
                print("PedaL MUITO BOM")
            elif F_piloto_necessaria_N / desaceleracao_alvo_g > 445 and F_piloto_necessaria_N / desaceleracao_alvo_g < 668:
                print("Pedal ACEITÁVEL")
            elif F_piloto_necessaria_N / desaceleracao_alvo_g > 668:
                print("Pedal RUIM")
        
        self.results['F_piloto_para_meta_g'] = F_piloto_necessaria_N


# --- Exemplo de Uso do Script de Dimensionamento ---
if __name__ == "__main__":
    
    print("Iniciando Script de Dimensionamento de Freios (V5 - Lógica Corrigida)...")
    print("="*60)

    # 1. PARÂMETROS DE ENTRADA ORGANIZADOS
    params = {
        # Parâmetros do Veículo
        'massa_total': 300.0,  # kg (Veículo + Piloto)
        'HCG': 0.35,           # m (Altura do Centro de Gravidade)
        'entre_eixos': 1.6,    # m (Wheelbase)
        'dist_cg_dianteiro': 1.1, # m (Distância do CG ao eixo Dianteiro)
        'dist_cg_traseiro': 0.5,  # m (Distância do CG ao eixo Traseiro)
        'raio_pneu': 0.3,      # m (Raio do pneu, ex: 10 polegadas)
        'mu_pneu': 1.5,        # Coeficiente de atrito pneu/solo (pista seca)
        
        # Parâmetros de Freio (Dimensionamento)
        'L_pedal': 220 / 1000.0,         # m (Comprimento total do pedal, L2)
        'L_pivot_cm': 50 / 1000.0,       # m (Distância pivot à balance bar, L1)
        'dist_bb_dianteiro': 11.2 / 1000.0, # m (Distância pivot BB ao CM Diant)
        'dist_bb_traseiro': 20 / 1000.0,  # m (Distância pivot BB ao CM Tras)
        'D_cm': 0.015875, # m (Cilindro Mestre 5/8 polegada = 15.875mm)
        'D_pin_f': 0.030, # m (Diâmetro pistão pinça diant, ex: 30mm)
        'N_pistoes_f_um_lado': 1, # (ex: 2 pistões de um lado)
        'D_pin_r': 0.028, # m (Diâmetro pistão pinça tras, ex: 28mm)
        'N_pistoes_r_um_lado': 1, # (ex: 1 pistão de um lado)
        'tipo_pinca_f': 'flutuante', # 'flutuante' ou 'fixa'
        'tipo_pinca_r': 'flutuante', # 'flutuante' ou 'fixa'
        'mu_pastilha': 0.45,           # Coef. atrito pastilha/disco
        'raio_efetivo_disco_f': 0.1, # m (100mm)
        'raio_efetivo_disco_r': 0.1, # m (100mm)
    }

    # 2. PARÂMETROS DA SIMULAÇÃO
    F_PILOTO_ALVO = 873.0 # Newtons (Força alvo para os cálculos diretos)

    # --- EXECUÇÃO DOS CÁLCULOS ---
    
    try:
        # Inicializa a classe de dimensionamento com o dicionário 'params' plano
        dimensionamento = BrakeSizing(params)
        
        print("\n--- PASSO 1: CÁLCULO DE PEDAL E BALANCE BAR ---")
        F_saida = dimensionamento.calcular_relacao_pedal(F_PILOTO_ALVO)
        dimensionamento.calcular_distribuicao_balance_bar(F_saida)

        print("\n--- PASSO 2: CÁLCULO HIDRÁULICO E TORQUE (BASEADO NA FORÇA DO PILOTO) ---")
        dimensionamento.calcular_pressao_hidraulica()
        dimensionamento.calcular_torque_frenagem()

        print("\n--- PASSO 3: ANÁLISE ESTÁTICA (a=0) --- (SOLICITADO)")
        # Calcula e printa apenas as cargas estáticas
        dimensionamento.calcular_transferencia_carga(0) # Calcula e armazena W_static
        print(f"Carga Estática Dianteira: {dimensionamento.results['W_f_static_N']:.1f} N")
        print(f"Carga Estática Traseira: {dimensionamento.results['W_r_static_N']:.1f} N")
        
        print("\n--- PASSO 4: ANÁLISE DE DESACELERAÇÃO MÁXIMA E EFICIÊNCIA ---")
        # Esta função usa os valores estáticos calculados no Passo 3
        dimensionamento.calcular_desaceleracao_eficiencia_max()
        
        print("\n--- PASSO 5: VERIFICAÇÃO (COM DESACELERAÇÃO MÁXIMA ALVO) --- (SOLICITADO)")
        # Cálculo da desaceleração alvo conforme sua fórmula: a = (Ff + Fr) / mv
        print("Calculando desaceleração alvo (a = (W_f_static*mu + W_r_static*mu) / m)...")
        mu_pneu = params['mu_pneu']
        Wf_stat = dimensionamento.results['W_f_static_N']
        Wr_stat = dimensionamento.results['W_r_static_N']
        massa_total = params['massa_total']
        
        Ff = Wf_stat * mu_pneu
        Fr = Wr_stat * mu_pneu
        
        a_alvo_ms2 = (Ff + Fr) / massa_total
        a_alvo_g = a_alvo_ms2 / G # Converte para G's para passar para a função
        print(a_alvo_ms2)
        print(f"Desaceleração Alvo (baseada em F_estatica * mu): {a_alvo_ms2:.3f} G")
        
        # Qual a força do piloto necessária para atingir a desaceleração alvo?
        dimensionamento.verificar_forca_piloto_para_desaceleracao(a_alvo_g)

        print("Dimensionamento concluído com sucesso.")

    except ValueError as e:
        print(f"\nERRO NO DIMENSIONAMENTO: {e}", file=sys.stderr)
    except KeyError as e:
        print(f"\nERRO: Parâmetro ausente no dicionário 'params': {e}", file=sys.stderr)
        print("Verifique se o dicionário 'params' está completo.")
    except Exception as e:
        print(f"\nOcorreu um erro inesperado: {e}", file=sys.stderr)
