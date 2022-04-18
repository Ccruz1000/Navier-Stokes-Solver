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
        self.r = p / (self.R * T)# Density
        self.mu = self.mu0 * (T / self.T0) ** (3/2) * (self.T0 + 110) / (T + 110) # Molecular viscosity coefficient
        self.lam = - 2 / 3 * self.mu  # Second viscosity coefficient
        self.a = np.sqrt(self.gm * self.R * T)  # Speed of sound
        self.cv = self.R / (self.gm - 1)  # Specific heat for an ideal gas at constant volume
        self.e = self.cv * T # Internal energy
        self.Et = self.r * (self.e + 0.5 * (u ** 2 + v ** 2))# Total energy
        self.cp  = self.gm * self.cv  # Specific heat for an ideal gas at constant pressure
        self.k = self.cp * self.mu / self.Pr  # Thermal conductivity

        # Define methods

    def deal(self):
        u = self.u
        v = self.v
        p = self.p
        T = self.T
        return u, v, p, T
    def getMuLambdaK(self):
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
        vp = max(4/3 * np.any(Mu), np.any(self.gm * Mu / self.Pr)) / self.r
        dtCFL = 1 / (np.abs(self.u) / dx + np.abs(self.v) / dy + self.a * np.sqrt( 1 / (dx ** 2) + 1 / (dy ** 2)) + 2 * vp * (1 / (dx **2) + 1 / (dy ** 2)))
        dt = K * np.amin(dtCFL)
        vp1 = np.amax(4/3 * Mu)
        vp2 = np.amax((self.gm * Mu) / self.Pr)
        if vp1 > vp2:
            vp = vp1 / self.r
        else:
            vp = vp2 / self.r
        dtCFL = 1 / (np.abs(self.u) / dx + np.abs(self.v) / dy + self.a * np.sqrt( 1 / (dx ** 2) + 1 / (dy ** 2)) + 2 * vp * (1 / (dx **2) + 1 / (dy ** 2)))
        dt = K * np.amin(dtCFL)
        return dt




