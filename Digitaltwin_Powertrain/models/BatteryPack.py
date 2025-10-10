import numpy as np
class BatteryPack:
    """A Battery Pack object based on a simplified Modified Shepherd model.

    This class represents a battery pack composed of cells arranged in series 
    and parallel. It models the terminal voltage, state-of-charge (SOC) dynamics, 
    accumulated current, and provides physical properties such as total mass 
    and volume.

    Parameters
    ----------
    tipo_celula : str
        Type of electrochemical cell. Implemented options:
        - 'Li-ion'
        - 'LiFePO4'
    n_serie : int
        Number of cells in series.
    n_paralelo : int
        Number of cells in parallel.
    soc_inicial : float, optional
        Initial state of charge of the pack (0.0–1.0). Default is 1.0.

    Attributes
    ----------
    E0 : float
        Nominal open-circuit voltage of a single cell [V].
    K : float
        Polarization constant [V].
    Q : float
        Capacity of a single cell [Ah].
    A : float
        Exponential voltage coefficient [V].
    B : float
        Exponential capacity coefficient [Ah^-1].
    R : float
        Internal resistance [Ω].
    peso : float
        Mass of a single cell [kg].
    volume : float
        Volume of a single cell [L].
    n_serie : int
        Number of series-connected cells.
    n_paralelo : int
        Number of parallel-connected cells.
    soc : float
        Current state of charge (SOC).
    carga_total : float
        Total pack capacity [Ah].
    inv_carga_total : float
        Precomputed inverse of capacity in [1/Coulombs].
    Iast : float
        Accumulated current [Coulombs].
    tempo_acumulado : float
        Elapsed time [s].
    tensao_hist : list
        History of terminal voltage values.
    corrente_hist : list
        History of current values.
    soc_hist : list
        History of SOC values.
    tempo : list
        History of time values.

    Methods
    -------
    definir_celula(tipo):
        Returns electrochemical parameters for a given cell type.
    calcular_derivadas(corrente):
        Computes time derivatives of SOC, accumulated current, and time.
    calcular_tensao(corrente, soc, Iast, tempo_acumulado):
        Computes pack terminal voltage given state and load current.
    calcular_peso_total():
        Returns total pack mass [kg].
    calcular_volume_total():
        Returns total pack volume [L].
    calcular_tensao_nominal():
        Returns nominal pack voltage [V].
    """

    def __init__(self, tipo_celula, n_serie, n_paralelo, soc_inicial=1.0):
        self.parametros = self.definir_celula(tipo_celula)
        if self.parametros is None:
            raise ValueError(f"Unknown cell type: {tipo_celula}")

        # Extract parameters
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
        self.carga_total = self.Q * n_paralelo  # Ah
        self.inv_carga_total = 1.0 / (self.carga_total * 3600.0)

        # Integration states
        self.Iast = 0.0  # Accumulated current [Coulombs]
        self.tempo_acumulado = 0.0  # Time [s]

        # History
        self.tensao_hist = []
        self.corrente_hist = []
        self.soc_hist = []
        self.tempo = []

    def definir_celula(self, tipo):
        """Return parameters of a predefined electrochemical cell.

        Parameters
        ----------
        tipo : str
            Cell type identifier.

        Returns
        -------
        dict or None
            Dictionary with cell parameters if type exists, otherwise None.
        """
        catalogo = {
            'Li-ion': {
                'E0': 3.7,
                'K': 0.005,
                'Q': 3.0,
                'A': 0.1,
                'B': 0.01,
                'R': 0.02,
                'peso': 0.045,
                'volume': 0.025
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
        """Compute state derivatives for numerical integration.

        Parameters
        ----------
        corrente : float
            Applied current [A] (positive for discharge, negative for charge).

        Returns
        -------
        tuple
            dsoc_dt : float
                Derivative of state of charge.
            dIast_dt : float
                Derivative of accumulated current.
            dtime_dt : float
                Derivative of elapsed time (constant = 1.0).
        """
        dsoc_dt = -corrente * self.inv_carga_total
        dIast_dt = corrente
        dtime_dt = 1.0
        return dsoc_dt, dIast_dt, dtime_dt

    def calcular_tensao(self, corrente, soc, Iast, tempo_acumulado):
        """Compute the terminal voltage of the battery pack.

        Parameters
        ----------
        corrente : float
            Applied current [A] (positive for discharge, negative for charge).
        soc : float
            Current state of charge (0.0–1.0).
        Iast : float
            Accumulated current [Coulombs].
        tempo_acumulado : float
            Elapsed time [s].

        Returns
        -------
        float
            Terminal voltage of the pack [V].
        """
        carga_utilizada = (1 - soc) * self.carga_total
        denom1 = max(self.carga_total - carga_utilizada, 1e-6)
        inv_denom1 = 1.0 / denom1

        term_exp = self.A * np.exp(-self.B * carga_utilizada)
        term_K = self.K * self.carga_total * inv_denom1
        term_carga = term_K * carga_utilizada
        term_Iast = term_K * (Iast / 3600)

        if corrente >= 0:  # Discharge
            if corrente < 100 and soc >= 0.0:
                Vt = self.E0 - self.R * corrente - term_Iast - term_carga + term_exp
            else:  # Limit max discharge current
                Vt = self.E0 - self.R * corrente * 0.2 - term_Iast - term_carga + term_exp
        else:  # Charge
            corrente_abs = abs(corrente)
            denom2 = max(corrente_abs * tempo_acumulado - 0.1 * self.carga_total, 1e-6)
            term_K_charge = self.K * self.carga_total / denom2 * (Iast / 3600)
            Vt = self.E0 - self.R * corrente_abs - term_K_charge - term_carga + term_exp

        return Vt * self.n_serie

    def calcular_peso_total(self):
        """Return the total pack mass.

        Returns
        -------
        float
            Total mass [kg].
        """
        return self.n_serie * self.n_paralelo * self.peso

    def calcular_volume_total(self):
        """Return the total pack volume.

        Returns
        -------
        float
            Total volume [L].
        """
        return self.n_serie * self.n_paralelo * self.volume

    def calcular_tensao_nominal(self):
        """Return the nominal pack voltage.

        Returns
        -------
        float
            Nominal voltage [V].
        """
        return self.E0 * self.n_serie
