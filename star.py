import numpy as np

#class Constants():

class Star():
    def __init__(self, settings, mass = 1, X = 0.7, Y = 0.24, Pc = 1,
     Tc = 1, L = 1, R = 1):
        M_sun = 1.989e33 #cgs
        L_sun = 3.847e33

        #c The input quantities are X and Y (the hydrogen and helium
        #c mass fractions), the model mass (in Solar units), and your
        #c guesses for central pressure, central temperature, total
        #c radius, and total luminosity.

        self.solar_masses = mass
        self.H_frac = X
        self.He_frac = Y
        self.pressure_center = Pc
        self.temp_center = Tc
        self.solar_luminosities = L
        self.radius = R
        self.Me_frac = 1 - X - Y

        self.mass = M_sun * self.solar_masses

        self.pressure_field = np.empty(settings.mesh_points)
        self.temperature_field = np.empty(settings.mesh_points)
        self.mass_field = np.empty(settings.mesh_points)

