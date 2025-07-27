import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import time

class SimulacaoBateria:
    def __init__(self):
        # --- Parâmetros da Bateria ---
        self.rho_bateria = 2500
        self.cp_bateria = 780
        self.k_bateria = 1.5
        self.alpha = self.k_bateria / (self.rho_bateria * self.cp_bateria)

        # --- Dimensões ---
        self.L, self.H, self.D = 0.254, 0.1524, 0.1016
        self.nx, self.ny, self.nz = 30, 15, 10
        self.x = np.linspace(0, self.L, self.nx)
        self.y = np.linspace(0, self.H, self.ny)
        self.z = np.linspace(0, self.D, self.nz)
        self.dx = np.diff(self.x)
        self.dy = np.diff(self.y)
        self.dz = np.diff(self.z)
        self.X, self.Y, self.Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # --- Condições Iniciais ---
        self.T_initial = 323
        self.T_inlet = 300
        self.T_outlet = 302
        self.T_amb = 300
        self.h_coolant = 1000
        self.h_air_top = 15
        self.h_air_sides = 25

        # --- Tempo e Convergência ---
        self.dt = 3
        self.max_iter = 3000
        self.convergence_limit = 1e-4

        # --- Inicializa Temperatura ---
        self.T = np.full((self.nx, self.ny, self.nz), self.T_initial, dtype=float)

        # --- Armazenamento ---
        self.temp_max_time = []
        self.temp_avg_time = []
        self.central_temp_over_time = []
        self.time_steps_to_plot = [0, 25, 50, 100, 500, 1000, 3000]
        self.temp_profiles = []

    def q_dot(self, t):
        """
        Simula a geração volumétrica de calor da bateria (W/m³) ao longo do tempo,
        durante uma volta de 2 minutos (120 segundos), com aceleração, frenagens etc.
        """
        t = t % 120  # repete a cada volta

        # Picos de calor típicos por situação (em W/m³)
        base = 3000  # base constante mínima
        pico_aceleracao = 15000
        medio = 5000
        recuperacao = 2000
        parado = 500

        if t < 10:
            # Largada — pico rápido
            return base + (pico_aceleracao - base) * np.exp(-((t - 5) ** 2) / 4)
        elif 10 <= t < 30 or 60 <= t < 80 or 100 <= t < 115:
            # Retas — geração média-alta com oscilações rápidas
            return medio + 1500 * np.sin(2 * np.pi * t / 5)
        elif 30 <= t < 60 or 80 <= t < 100:
            # Curvas / frenagens — geração baixa, com variação leve
            return recuperacao + 500 * np.cos(2 * np.pi * t / 15)
        elif 115 <= t < 120:
            # Resfriamento — quase nada gerado
            return parado
        else:
            return base

    def simular(self):
        start_time = time.time()

        for iter in range(self.max_iter):
            T_old = self.T.copy()

            # --- Cond. de Contorno ---
            self.T[:, :, 0] = self.T_inlet + (self.T_outlet - self.T_inlet) * self.Y[:, :, 0] / self.H
            self.T[:, :, -1] = (self.h_coolant * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dz[-1] * T_old[:, :, -1]) / (self.rho_bateria * self.cp_bateria * self.dz[-1] + self.h_coolant * self.dt)
            self.T[0, :, :] = (self.h_air_top * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dy[0] * T_old[0, :, :]) / (self.rho_bateria * self.cp_bateria * self.dy[0] + self.h_air_top * self.dt)
            self.T[:, 0, :] = (self.h_air_sides * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dx[0] * T_old[:, 0, :]) / (self.rho_bateria * self.cp_bateria * self.dx[0] + self.h_air_sides * self.dt)
            self.T[:, -1, :] = (self.h_air_sides * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dx[-1] * T_old[:, -1, :]) / (self.rho_bateria * self.cp_bateria * self.dx[-1] + self.h_air_sides * self.dt)

            # --- Diferenças Finitas ---
            for i in range(1, self.nx - 1):
                for j in range(1, self.ny - 1):
                    for k in range(1, self.nz - 1):
                        self.T[i, j, k] += self.alpha * self.dt * (
                            (T_old[i+1, j, k] - 2*T_old[i, j, k] + T_old[i-1, j, k]) / self.dx[i]**2 +
                            (T_old[i, j+1, k] - 2*T_old[i, j, k] + T_old[i, j-1, k]) / self.dy[j]**2 +
                            (T_old[i, j, k+1] - 2*T_old[i, j, k] + T_old[i, j, k-1]) / self.dz[k]**2
                        ) + self.q_dot(iter * self.dt) * self.dt / (self.rho_bateria * self.cp_bateria)

            # --- Registro ---
            self.temp_max_time.append(np.max(self.T))
            self.temp_avg_time.append(np.mean(self.T))
            ci, cj, ck = self.nx//2, self.ny//2, self.nz//2
            self.central_temp_over_time.append(self.T[ci, cj, ck])

            if iter in self.time_steps_to_plot:
                self.temp_profiles.append(self.T.copy())

            if np.max(np.abs(self.T - T_old)) < self.convergence_limit:
                print(f"Convergência atingida em {iter} iterações.")
                break

            if iter % 10 == 0:
                elapsed = time.time() - start_time
                print(f"Iteração: {iter}, Temp máxima: {np.max(self.T):.2f} K, Tempo decorrido: {elapsed:.2f} s")

        print(f"Simulação finalizada em {time.time() - start_time:.2f} segundos.")

    def plot1(self):
        print(self.temp_profiles)
        plt.figure(figsize=(10, 6))
        tempo = np.arange(len(self.temp_max_time)) * self.dt
        plt.plot(tempo, self.temp_max_time, label="Temperatura Máxima (K)", color='r')
        plt.plot(tempo, self.temp_avg_time, label="Temperatura Média (K)", color='b')
        plt.plot(tempo, self.central_temp_over_time, label="Temperatura Central (K)", color='g')
        plt.xlabel("Tempo (s)")
        plt.ylabel("Temperatura (K)")
        plt.title("Evolução da Temperatura da Bateria")
        plt.legend()
        plt.grid(True)
        plt.show()

    def plot3D(self):

        for idx, T_profile in enumerate(self.temp_profiles):
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection='3d')

            temp_celsius = T_profile - 273.15
            norm = mcolors.Normalize(vmin=np.min(temp_celsius), vmax=np.max(temp_celsius))
            cmap = plt.cm.hot_r  # versão invertida para cores mais frias claras

            voxels = np.ones(T_profile.shape, dtype=bool)
            colors = np.empty(T_profile.shape + (4,), dtype=float)

            it = np.nditer(temp_celsius, flags=['multi_index'])
            while not it.finished:
                idx_xyz = it.multi_index
                colors[idx_xyz] = cmap(norm(it[0]))
                it.iternext()

            ax.voxels(voxels, facecolors=colors, edgecolor='k')

            ax.set_xlabel('Comprimento (x) [cm]')
            ax.set_ylabel('Altura (y) [cm]')
            ax.set_zlabel('Profundidade (z) [cm]')
            plt.title(f'Distribuição de Temperatura por Célula (Iteração {self.time_steps_to_plot[idx]})')

            mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            mappable.set_array(temp_celsius)
            fig.colorbar(mappable, ax=ax, label="Temperatura (°C)")

            plt.show()

    def plot_plano_central(self, step_idx=-1):

        T_profile = self.temp_profiles[step_idx]
        temp_celsius = T_profile - 273.15

        k = self.nz // 2
        plano = temp_celsius[:, :, k]

        plt.figure(figsize=(8, 6))
        cp = plt.contourf(self.x, self.y, plano.T, 20, cmap='hot_r')
        plt.colorbar(cp, label='Temperatura (°C)')
        plt.xlabel('Comprimento (cm)')
        plt.ylabel('Altura (cm)')
        plt.title(f'Distribuição de Temperatura - Plano central (z) - Iteração {self.time_steps_to_plot[step_idx]}')
        plt.gca().invert_yaxis()
        plt.show()

simulacao = SimulacaoBateria()
simulacao.simular()
simulacao.plot1()
simulacao.plot3D()        
simulacao.plot_plano_central() 
