# Import functions
import numpy as np


# Define primitive class
class Primitives:
    def __init__(self, u, v, p, T):
        # Initalized input variables
        self.u = u  # x-direction velocity
        self.v = v  # y-direction velocity
        self.p = p  # pressure
        self.T = T  # Absolute temperature
        # Variables with constant value
        self.mu0 = 0.000017894  # Pa, Sea level dynamic viscosity
        self.T0 = 288.16 # K, Sea level temperature
        self.gm = 1.4  # Ratio of specific heat
        self.Pr = 0.71  # Prandtl number
        self.R = 287    # J/kg-K, Specific gas constant
        # Calculated variables
        self.r = p / (R * T)# Density
        self.e = cv * T # Internal energy
        self.Et = r * (e + 0.5 * (u ** 2 + v ** 2))# Total energy
        self.mu = mu0 * (T / T0) ** (3/2) * (T0 + 110) / (T + 110) # Molecular viscosity coefficient
        self.lam = - 2 / 3 * mu  # Second viscosity coefficient
        self.a = np.sqrt(gm * R * T) # Speed of sound
        self.k = cp * mu / Pr # Thermal conductivity
        self.cv = R / (gm - 1)# Specific heat for an ideal gas at constant volume
        self.cp  = gm * cv# Specific heat for an ideal gas at constant pressure

        # Define methods
    def getMuLmadaK(self):
        muOut = self.mu
        lamdaOut = -2/3 * muOut
        kOut = self.cp * muOut /self.Pr
        return muOut, lamdaOut, kOut


    def calculateReynoldsNumber(self, referenceLength):
        r_mu = self.r / self.mu
        ReX = referenceLength * self.u * r_mu
        ReY = referenceLength * self.v * r_mu
        return ReX, ReY


    def calculateTimeStep(self, dx, dy, K):
        Mu = self.mu
        vp = max(4/3 * Mu, self.gm * Mu / self.Pr) / self.r
        dtCFL = 1 / (abs(self.u) / dx + abs(self.v) / dy + self.a * np.sqrt( 1 / (dx ** 2) + 1 / (dy ** 2)) + 2 * vp * (1 / (dx **2) + 1 / (dy ** 2)))
        dt = K * min(dtCFL[:])





