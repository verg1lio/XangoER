class Dynamics:
    """A Vehicle Dynamics object.

    This class will create a Vehicle Dynamics digital twin with the springs, suspension bars, damper cylinder, brakes, tires, and transmission elements provided.

    Parameters AINDA A ALTERAR, ISSO É UM EXEMPLO
    ----------
    shaft_elements : list
        List with the shaft elements.
    disk_elements : list
        List with the disk elements.
    bearing_elements : list
        List with the bearing elements.
    automeshing : boolean
        Set it True to use the automeshing method. Default is False.
        If automeshing is True, the previous shaft_elements parameter is now defined by:
            shaft_elements : list
                List with the length, inner and outter diameter and material of each element, as follows:
                    [[length, inner, outter, material],[length, inner, outter, material],...]
        For the other parameters please check the respective bearing and disk classes for more information.
    **kwargs : dict, optional
        If automeshing is True, these parameters needs to be informed.
        The possible arguments are:
            alpha : float
                Proportional damping coefficient, associated to the element Mass matrix
            beta : float
                Proportional damping coefficient, associated to the element Stiffness matrix

    Returns
    -------
    A rotor object.

    Attributes
    ----------
    MM : array
        Global mass matrix.
    KK : array
        Global stiffness matrix.
    CCgyros: array
        Global gyroscopic matrix.
    CCtotal: array
        Global damping matrix

    Examples
    --------
    >>> import lmest_rotor as lm
    >>> rotor = lm.rotor_example()
    >>> rotor.MM
    array(30x30)
    """



    def __init__(self, spring_type, spring_k, spring_F):
        self.spring_type = spring_type # Hooke, 
        self.spring_k = spring_k # rigidez da mola [N/m]
        self.spring_F = spring_F # força que a mola recebe [N]
        self.spring_non_lin_coef = spring_non_lin_coef # coeficiente de ganho não-linear
        

    def Spring(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        if self.spring_type == 'Hooke':
            spring_x = self.spring_F/self.spring_k
        if self.spring_type == 'Softening'
            spring_x = self.spring_F/(self.spring_non_lin_coef*(self.spring_k)**2)

        return spring_x



    def Suspension(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return Sus



    def Damper(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return Dam

    

    def Brake(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return Bra



    def Tire(self):
        """Description.

        Detailed description.

        Returns
        -------
        Bat : variable type
            Description.

        Examples
        --------
        >>> example
        """
        
        return Tir



    def CalculateOutputs(self, cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, cp):
        peso = massa * 9.81
        rnet = (peso * (cgx * 0.001)) / (etex * 0.001)
        rned = massa * 9.81 - rnet
        ftr = (rnet * cfat) / (1 - ((cgy * 0.001) * cfat) / (etex * 0.001))
        tdcl = (ftr * cgy * 0.001) / (etex * 0.001)
        cnet = rnet + tdcl
        ptet = cnet * rpneu * 0.001 * cfat
        cpneu = cnet / 2
        tpneu = ftr * rpneu * 0.001
        redf = redp * red1
        tpwt = ptet / (red1 * redp * cp)
        acpr = (ftr / massa) / 9.81
        acfi = acpi * 9.81
        acfr = acpr * 9.81
        fti = massa * acfi
        tpi = fti * rpneu * 0.001
        tpwti = tpi / (red1 * redp * cp)
        tci = (fti * cgy) / etex
        tcr = (ftr * cgy) / etex
        cteti = rnet + tci
        return peso, rnet, rned, ftr, tdcl, cnet, ptet, cpneu, tpneu, redf, tpwt, acpr, acfi, acfr, fti, tpi, tpwti, tci, tcr, cteti

    def CurveTorquePower(self, data_matrix):
        rpm_values = []
        torque_values = []
        power_values = []
        for data in data_matrix:
            rpm = data["rpm"]
            ptc = data["ptc"]
            trq = data["trq"]
            rpm_values.append(rpm)
            torque_values.append(trq)
            power_values.append(ptc)
        return rpm_values, torque_values, power_values
        
    # Valores experimentais para a curva de torque e potência
    matriz_dados = [
        {"rpm": 0, "ptc": 0.0, "trq": 205.89},
        {"rpm": 375, "ptc": 8.0742858188337, "trq": 205.61},
        {"rpm": 750, "ptc": 16.1100871276608, "trq": 205.12},
        {"rpm": 1125, "ptc": 24.1521716217951, "trq": 205.01},
        {"rpm": 1500, "ptc": 32.1824751434784, "trq": 204.88},
        {"rpm": 1875, "ptc": 40.1318826543315, "trq": 204.39},
        {"rpm": 2250, "ptc": 48.1182038788644, "trq": 204.22},
        {"rpm": 2625, "ptc": 56.1159133767666, "trq": 204.14},
        {"rpm": 3000, "ptc": 64.0853485407864, "trq": 203.99},
        {"rpm": 3375, "ptc": 71.9475768555021, "trq": 203.57},
        {"rpm": 3750, "ptc": 79.894628171865, "trq": 203.45},
        {"rpm": 4125, "ptc": 87.7372215324957, "trq": 203.11},
        {"rpm": 4500, "ptc": 95.5531113555708, "trq": 202.77},
        {"rpm": 4875, "ptc": 103.418873962022, "trq": 202.58},
        {"rpm": 5250, "ptc": 111.242225067649, "trq": 202.34},
        {"rpm": 5625, "ptc": 119.04083613113, "trq": 202.09},
        {"rpm": 6000, "ptc": 126.926626390747, "trq": 202.01},
        {"rpm": 6375, "ptc": 134.826161118224, "trq": 201.96},
        {"rpm": 6750, "ptc": 142.31885959706, "trq": 201.34},
        {"rpm": 7125, "ptc": 148.748128962653, "trq": 199.36},
        {"rpm": 7500, "ptc": 155.304632830716, "trq": 197.74},
        {"rpm": 7875, "ptc": 161.733902196308, "trq": 196.12},
        {"rpm": 8250, "ptc": 168.856678140183, "trq": 195.45},
        {"rpm": 8625, "ptc": 174.599116811882, "trq": 193.31},
        {"rpm": 9000, "ptc": 179.702241378574, "trq": 190.67},
        {"rpm": 9375, "ptc": 182.36945354148, "trq": 185.76},
        {"rpm": 9750, "ptc": 184.007794110332, "trq": 180.22},
        {"rpm": 10125, "ptc": 184.850526339661, "trq": 174.34},
        {"rpm": 10500, "ptc": 184.230847188738, "trq": 167.55},
        {"rpm": 10875, "ptc": 183.715625993548, "trq": 161.32},
        {"rpm": 11250, "ptc": 179.600925015495, "trq": 152.45},
        {"rpm": 11625, "ptc": 171.247037450491, "trq": 140.67},
        {"rpm": 12000, "ptc": 167.949543261456, "trq": 133.65},
        {"rpm": 12375, "ptc": 157.815550759106, "trq": 121.78},
        {"rpm": 12750, "ptc": 144.799932395241, "trq": 108.45},
        {"rpm": 13125, "ptc": 134.393406730191, "trq": 97.78},
        {"rpm": 13500, "ptc": 122.823706385146, "trq": 86.88},
        {"rpm": 13875, "ptc": 112.824409667819, "trq": 77.65},
        {"rpm": 14250, "ptc": 103.622292086342, "trq": 69.44},
        {"rpm": 14625, "ptc": 91.6771714191918, "trq": 59.86},
        {"rpm": 15000, "ptc": 74.565701633196, "trq": 47.47}
    ]
    
    def Transmission(self, cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, cp):
        peso, rnet, rned, ftr, tdcl, cnet, ptet, cpneu, tpneu, redf, tpwt, acpr, acfi, acfr, fti, tpi, tpwti, tci, tcr, cteti = self.CalculateOutputs(cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, cp)

        # Print dos resultados obtidos
        print("Resultados:")
        print("Peso:", peso, "kg")
        print("Reação no eixo traseiro:", rnet, "N")
        print("Reação no eixo dianteiro:", rned, "N")
        print("Força de trativa:", ftr, "N")
        print("Transferência de carga longitudinal:", tdcl, "N")
        print("Carga no eixo traseiro:", cnet, "N")
        print("Pico de torque no eixo traseiro:", ptet, "Nm")
        print("Carga no pneu:", cpneu, "N")
        print("Torque no pneu:", tpneu, "Nm")    
        print("Redução final:", redf)
        print("Torque necessário no motor:", tpwt, "Nm")
        print("Aceleração primária real (g):", acpr)
        print("Aceleração final ideal:", acfi, "m/s²")
        print("Aceleração final real:", acfr, "m/s²")
        print("Força trativa ideal:", fti, "N")
        print("Torque no pneu ideal:", tpi, "Nm")
        print("Torque no motor ideal:", tpwti, "Nm")
        print("Transferência de carga ideal:", tci, "N")
        print("Transferência de carga real:", tcr, "N")
        print("Carga total no eixo traseiro ideal:", cteti, "N")

        print("\nMatriz de RPM, Torque e Potência:")
        print("RPM\t\tTorque (Nm)\tPotência (kW)")
        for data in matriz_dados:
            rpm = data["rpm"]
            trq = data["trq"]
            ptc = data["ptc"]
            print("{:.2f}\t\t{:.2f}\t\t{:.2f}".format(rpm, trq, ptc))

        # Plotando o gráfico
        rpm_values, torque_values, power_values = self.CurveTorquePower(matriz_dados)
        plt.figure(figsize=(10, 6))
        plt.plot(rpm_values, torque_values, label="Torque")
        plt.plot(rpm_values, power_values, label="Potência")
        plt.title("Curva de Torque e Potência")
        plt.xlabel("RPM")
        plt.ylabel("Torque / Potência")
        plt.legend()
        plt.grid(True)
        plt.show()

        return peso, rnet, rned, ftr, tdcl, cnet, ptet, cpneu, tpneu, redf, tpwt, acpr, acfi, acfr, fti, tpi, tpwti, tci, tcr, cteti, matriz_dados

    # Criando uma instância da classe, definindo valores para os parâmetros e chamando o método Transmission
    dynamics = Dynamics()
    cgx = 853  # mm
    cgy = 294  # mm
    massa = 347  # kg
    etex = 1567  # mm
    cfat = 0.9  # coeficiente de atrito
    rpneu = 259  # mm
    acpi = 1.2  # g
    redp = 2.12  # redução primária
    red1 = 2.76  # redução da marcha única
    cp = 1.5  # relação coroa-pinhão
    dynamics.Transmission(cgx, cgy, massa, etex, cfat, rpneu, acpi, redp, red1, cp)

    def CarPerformance(self):
        
        peso = massa * 9.81 # peso
        rdt = 0.9 # rendimento da transmissão
        mtrd = 0.9 # transmissão motor-roda
        cfar = 0.54 # coeficiente de arrasto
        da = 1.162 # densidade do ar
        af = 1.06 # área frontal
        bscf = 0.015 # basic f
        spdf = 0.012 # speed f

        parametros = []

        for dado in matriz_dados:
            # Cálculo da força trativa (N)
            ftf = ((dado["trq"] * redp * red1 * cp) / (rpneu * 0.001)) * rdt

            # Cálculo da velocidade angular (rad/s)
            va = (dado["rpm"] * 2 * math.pi) / (60 * redp * red1 * cp)

            # Cálculo da velocidade linear (m/s)
            vl = (va * (rpneu * 0.001)) * mtrd

            # Cálculo da força de arrasto (N)
            fa = (da * vl ** 2 * cfar * af) / 2

            # Cálculo da resistência de rolamento (N)
            rr = (bscf + (3.24 * spdf * ((vl / 100 * 0.44704) ** 2.5))) * peso

            # Cálculo da força final (N)
            ff = ftf - fa - rr

            # Armazenar os parâmetros calculados em um dicionário
            parametro = {
                "ftf": ftf,
                "va": va,
                "vl": vl,
                "fa": fa,
                "rr": rr,
                "ff": ff
            }

            parametros.append(parametro)

        # Imprimir os parâmetros calculados com títulos
        print("Força Trativa [N]\tVelocidade Angular [rad/s]\tVelocidade Linear [m/s]\tForça de Arrasto [N]\tResistência de Rolamento [N]\tForça Final [N]")
        for param in parametros:
            print(f"{param['ftf']}\t{param['va']}\t{param['vl']}\t{param['fa']}\t{param['rr']}\t{param['ff']}")

    dynamics.CarPerformance()
            
    def HalfShaftsSizing(self, fsi=1.25, tet=786, tec=471.6, dif=1):
        # Considerando o uso do aço 4340
        # Obtendo o maior torque do motor a partir dos dados experimentais 
        tmax = max(data["trq"] for data in matriz_dados)
    
        # Calculando o torque máximo nos semieixos
        tmsx = tmax * redp * red1 * cp * dif
    
        # Calculando o torque máximo de projeto
        tmp = tmsx * fsi
    
        # Calculando o diâmetro dos semieixos (mm)
        dsx = (((2 * tmp) / (math.pi * tec * 10 ** 6)) ** (1 / 3)) * 2000
    
        # Calculando o fator de segurança obtido
        fso = (math.pi * (((dsx / 1000) / 2) ** 3) * tec * (10 ** 6)) / (2 * tmsx)
    
        # Calculando o fator de segurança para 1 polegada
        fs1p = (math.pi * ((0.0254 / 2) ** 3) * tec * (10 ** 6)) / (2 * tmsx)
    
        # Print dos resultados obtidos
        print("Dimensionamento dos Semieixos:")
        print("Torque máximo do motor:", tmax, "Nm")
        print("Torque máximo nos semieixos:", tmsx, "Nm")
        print("Torque máximo de projeto:", tmp, "Nm")
        print("Diâmetro dos semieixos:", dsx, "mm")
        print("Fator de segurança ideal:", fsi)
        print("Fator de segurança obtido:", fso)
        print("Fator de segurança para 1 polegada:", fs1p)
    
        return tmax, tmsx, tmp, dsx, fso, fs1p

    dynamics.HalfShaftsSizing()


        
       

