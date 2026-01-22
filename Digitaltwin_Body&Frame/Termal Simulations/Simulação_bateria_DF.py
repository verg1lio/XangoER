import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import time
import pandas

class SimulacaoBateria:
    def __init__(self,velocidade,correntes,battery_n_serie,battery_n_paralelo):

        # --- Parâmetros do Carro ---

        self.v = velocidade
        self.lenv = len(velocidade)

        # --- Parâmetros da Bateria ---
        self.rho_bateria = 2500
        self.cp_bateria = 780
        self.k_bateria = 1.5
        self.alpha = self.k_bateria / (self.rho_bateria * self.cp_bateria)
        self.battery_R = 0.2
        self.battery_n_serie = battery_n_serie
        self.battery_n_paralelo = battery_n_paralelo
        self.correntes = correntes

        # --- Parâmetros do Ar ---
        self.rho = 1.184
        self.mi = 1.849e-5
        self.ni = self.mi / self.rho
        self.k_ar = 0.02551
        self.c = 1007
        self.Pr = self.mi * self.c / self.k_ar


        # --- Dimensões ---
        self.tipo_geometria = "placa"
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
        self.T_initial = 300
        self.T_inlet = 300
        self.T_outlet = 302
        self.T_amb = 300
        self.h_coolant = 2000

        # --- Tempo e Convergência ---
        self.dt = 0.1
        self.max_iter = 6000
        self.convergence_limit = 1e-6

        # --- Inicializa Temperatura ---
        self.T = np.full((self.nx, self.ny, self.nz), self.T_initial, dtype=float)

        # --- Armazenamento ---
        self.temp_max_time = []
        self.temp_avg_time = []
        self.central_temp_over_time = []
        self.time_steps_to_plot = [0, 25, 50, 100, 500, 1000, 3000]
        self.temp_profiles = []

    def reynolds(self,idx):
        return self.rho * self.v[idx] * self.L / self.mi

    def calc_Nu(self,idx):
        Rey = self.reynolds(idx)
        if self.tipo_geometria == 'placa':
            if Rey < 200000:
                Nu = 0.332 * Rey**0.5 * self.Pr**(1/3)
            else:
                X = 200000 * self.ni / self.v[idx]
                Rey_X = X * self.v[idx] / self.ni
                A = 0.037 * Rey_X**0.8 - 0.664 * Rey_X**0.5
                Nu = 0.037 * (Rey**0.8 - A) * self.Pr**(1/3)
        elif self.tipo_geometria == 'cilindro':
            Nu = 0.3 + (0.62 * Rey**0.5 * self.Pr**(1/3)) / (1 + (0.4/self.Pr)**(2/3))**0.25
        else:
            raise ValueError("Tipo de geometria não suportada.")
        return Nu
    
    def calc_h(self, t):
        t = t%60
        idx = int(t*(self.lenv/60))
        h_air_sides = (self.calc_Nu(idx) * self.k_ar) / self.H
        h_air_top = (self.calc_Nu(idx) * self.k_ar) / self.D
        return h_air_sides, h_air_top


    def calcular_qdot(self, t):

        I = self.correntes[int(t)]

        R_pack = self.battery_R * self.battery_n_serie / self.battery_n_paralelo

        qdot = I**2 * R_pack * self.battery_n_paralelo * self.battery_n_serie

        return qdot



    def simular(self):
        start_time = time.time()

        for iter in range(self.max_iter):
            T_old = self.T.copy()

            h_air_sides, h_air_top = self.calc_h(iter*self.dt)

            # --- Cond. de Contorno ---
            self.T[:, :, 0] = self.T_inlet + (self.T_outlet - self.T_inlet) * self.Y[:, :, 0] / self.H
            self.T[:, :, -1] = (self.h_coolant * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dz[-1] * T_old[:, :, -1]) / (self.rho_bateria * self.cp_bateria * self.dz[-1] + self.h_coolant * self.dt)
            self.T[0, :, :] = (h_air_top * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dy[0] * T_old[0, :, :]) / (self.rho_bateria * self.cp_bateria * self.dy[0] + h_air_top * self.dt)
            self.T[:, 0, :] = (h_air_sides * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dx[0] * T_old[:, 0, :]) / (self.rho_bateria * self.cp_bateria * self.dx[0] + h_air_sides * self.dt)
            self.T[:, -1, :] = (h_air_sides * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dx[-1] * T_old[:, -1, :]) / (self.rho_bateria * self.cp_bateria * self.dx[-1] + h_air_sides * self.dt)
            self.T[-1, :, :] = (self.h_coolant * self.dt * self.T_amb + self.rho_bateria * self.cp_bateria * self.dz[-1] * T_old[-1, :, :]) / (self.rho_bateria * self.cp_bateria * self.dz[-1] + self.h_coolant * self.dt)

            # --- Diferenças Finitas ---
            for i in range(1, self.nx - 1):
                for j in range(1, self.ny - 1):
                    for k in range(1, self.nz - 1):
                        self.T[i, j, k] += self.alpha * self.dt * (
                            (T_old[i+1, j, k] - 2*T_old[i, j, k] + T_old[i-1, j, k]) / self.dx[i]**2 +
                            (T_old[i, j+1, k] - 2*T_old[i, j, k] + T_old[i, j-1, k]) / self.dy[j]**2 +
                            (T_old[i, j, k+1] - 2*T_old[i, j, k] + T_old[i, j, k-1]) / self.dz[k]**2
                        ) + self.calcular_qdot(iter * self.dt) * self.dt / (self.rho_bateria * self.cp_bateria)

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

    def plot_planos_centrais(self, step_idx=-1):

        T_profile = self.temp_profiles[step_idx]
        temp_celsius = T_profile - 273.15
        i = self.nx // 2
        j = self.ny // 2
        k = self.nz // 2
        plano_x = temp_celsius[i, :, :]
        plano_y = temp_celsius[:, j, :]
        plano_z = temp_celsius[:, :, k]

        plt.figure(figsize=(8, 6))
        cp = plt.contourf(self.y, self.z, plano_x.T, 20, cmap='hot_r')
        plt.colorbar(cp, label='Temperatura (°C)')
        plt.ylabel('Altura (m)')
        plt.xlabel('Profundidade (m)')
        plt.title(f'Distribuição de Temperatura - Plano central (x fixo) - Iteração {self.time_steps_to_plot[step_idx]}')
        plt.show()

        plt.figure(figsize=(8, 6))
        cp = plt.contourf(self.x, self.z, plano_y.T, 20, cmap='hot_r')
        plt.colorbar(cp, label='Temperatura (°C)')
        plt.xlabel('Largura (m)')
        plt.ylabel('Profundidade (m)')
        plt.title(f'Distribuição de Temperatura - Plano central (y fixo) - Iteração {self.time_steps_to_plot[step_idx]}')
        plt.show()

        plt.figure(figsize=(8, 6))
        cp = plt.contourf(self.x, self.y, plano_z.T, 20, cmap='hot_r')
        plt.colorbar(cp, label='Temperatura (°C)')
        plt.xlabel('Largura (m)')
        plt.ylabel('Altura (m)')
        plt.title(f'Distribuição de Temperatura - Plano central (z fixo) - Iteração {self.time_steps_to_plot[step_idx]}')
        plt.show()

'''
velocidade = [22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 
         20.38337567903568, 
         21.545306504458342, 
         18.44870403992181, 
         14.532464758634132, 
         10.132314142090229, 
         13.645807317358926, 
         10.98521599685414, 
         11.871116846023876, 
         11.931902433677735, 
         9.619294639921518, 
         8.636148443779636, 
         10.902335161887136, 
         12.967164107380286, 
         14.25483569210475, 
         11.883896140448755, 
         8.57836836241246, 
         12.108538628694356, 
         14.596406988325231, 
         16.899836008838136, 
         17.45983304147545, 
         15.288894988595548, 
         12.842933731423177, 
         12.634439993055487, 
         14.976273116223247, 
         13.973095461169224, 
         11.579687928823448, 
         8.722264604761211, 
         11.359478939207996, 
         9.094435881716263, 
         11.319130889988253,
         13.07468525475738, 
         14.59833565714002, 
         15.976504426364908, 
         17.14446723947686, 
         16.830891848520427, 
         15.315505333789117, 
         13.367467748270116, 
         15.839070166912261, 
         12.87843254673837, 
         10.164929113608478, 
         13.993634072088671, 
         16.962374356872935, 
         14.003942598957051, 
         10.54007037152099, 
         11.601262721432617, 
         8.402779200397974, 
         11.245974050925943, 
         14.425951003507905, 
         9.015008590849966, 
         14.564309571736706, 
         19.357772023609872, 
         21.408400176024916, 
         22.076762320615927, 
         22.213906212130908, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 
         21.680466816535866, 
         18.318708364389796, 
         18.82655443714922, 
         16.18113049469274, 
         14.009098398149911, 
         13.213893745617984, 
         13.856656820136125, 
         15.973969977325865, 
         17.548245920586588, 
         15.00980495468197, 
         11.651176735054943, 
         15.496991912175908, 
         18.56784286234689, 
         20.045827566165357, 
         20.902495718565774, 
         21.422995504684298, 
         21.740090625695657, 
         21.93237026461276, 
         22.04805031676879, 
         21.449032911053326, 
         18.37014076272792, 
         18.10276243907106, 
         14.462741105615441, 
         18.145297088036465, 
         19.812699548898404, 
         16.928681703703457, 
         18.986722142403824, 
         20.064260874597977, 
         20.715076863264, 
         21.144567800441305, 
         19.533333832340475, 
         17.022514427484026, 
         16.871544536354225, 
         12.681650425429755, 
         16.981169083729924, 
         12.070223899355673, 
         17.621894431409363, 
         20.74551066843905, 
         21.766653002813847, 
         22.113654504274805, 
         22.202216824991513, 
         22.218171082253487, 
         22.219909229400514, 
         22.219999168515447, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 22.22, 
         22.145472598606286, 
         18.58006739428217, 
         18.43113882893838, 
         15.689637672444322, 
         13.134317286682132,
         10.676381693816948, 
         8.085057874631122,
         9.174378875736755, 
         9.98018475168804, 
         12.003748839363118,
         11.960310936809222,
         9.534625930894963, 
         12.08240325969192, 
         14.153625919597872, 
         15.948497044931075, 
         17.433301895736015, 
         16.92685313686744, 
         15.14497309249129, 
         12.877644554851104, 
         15.794977074237984, 
         14.90780471652406, 
         14.860065516989591, 
         16.22542855802993, 
         19.500278937223065, 
         21.004379912565845, 
         21.73612423532988, 
         22.061354658740623, 
         22.179755727142, 
         22.21236484638815, 
         22.218935036138152, 
         22.2198938265428, 
         22.219991093964527, 
         22.219999261463496, 
         22.219999928253465, 
         22.219999990683604, 
         22.219999998377364, 
         22.219999999607847, 
         22.219999999873313, 
         18.51202873901757, 
         14.08752777298269, 
         15.37196635249142, 
         12.412120177991845, 
         9.972290580380443, 
         9.615736297801982, 
         8.35409211285864, 
         9.62441329734337, 
         10.874509878385052, 
         12.277374735275751, 
         14.046392109288826, 
         16.30065198274046, 
         16.307813604175237, 
         19.57959950050558, 
         21.332436954753742]

simulacao = SimulacaoBateria(velocidade)
simulacao.simular()
simulacao.plot1()
simulacao.plot3D()        
simulacao.plot_planos_centrais() 
'''