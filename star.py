import numpy as np

#class Constants():
def compnn(a, b):
    if (min(a, b)/max(a,b)) >= 0.99:
        return True
    else:
        return False

class Star():
    def __init__(self, settings, mass = 1, X = 0.7, Y = 0.24, Pc = 1,
     Tc = 1, L = 1, R = 1):
        self.M_sun = 1.989e33 #cgs
        self.L_sun = 3.847e33

        #c The input quantities are X and Y (the hydrogen and helium
        #c mass fractions), the model mass (in Solar units), and your
        #c guesses for central pressure, central temperature, total
        #c radius, and total luminosity.

        self.solar_masses = mass #mass/msun
        self.H_frac = X
        self.He_frac = Y
        self.pressure_center = Pc #cgs
        self.temp_center = Tc     #cgs
        self.luminosity = L       #cgs
        self.radius = R           #cgs
        self.Me_frac = 1 - X - Y

        self.filename = ''
        self.Teff = 0
#
#                            *****FINAL MODEL*****
# Pc: 5.6642D+16, Tc: 2.5474D+07, R: 1.8672D+11, L: 2.0190D+36
# Teff:  1.6885D+04, LOG(Teff): 4.2275, LOG(L/LSUN): 2.7200
#           1-Mr/M       LOG(r)   LOG(P)  LOG(T) LOG(RHO) LOG(L)
#    2    9.98776351D-01  9.74841 16.7439  7.4026  1.2122 34.9936

        self.coordinate_field = np.empty(settings.mesh_points)
        self.radius_field = np.empty(settings.mesh_points)
        self.pressure_field = np.empty(settings.mesh_points)
        self.temperature_field = np.empty(settings.mesh_points)
        self.density_field = np.empty(settings.mesh_points)
        self.luminosity_field = np.empty(settings.mesh_points)

    def out_parm(self):
        s =  "{0:.2f}".format(self.solar_masses) + ', '
        s += "{0:.2f}".format(self.H_frac) + ', '
        s += "{0:.2f}".format(self.He_frac) + ', '
        s += "{0:.2E}".format(self.pressure_center) + ', '
        s += "{0:.2E}".format(self.temp_center) + ', '
        s += "{0:.2E}".format(self.radius) + ', '
        s += "{0:.2f}".format(self.luminosity)
        return s

    def __str__(self):
        s = 'M/Msun=' + str(self.solar_masses) + ', '
        s += 'X= ' + str(self.H_frac) +', Y= ' + str(self.He_frac) + ', '
        s += 'Pc= ' + str(self.pressure_center) + ' dina/cm2, '
        s += 'Tc= ' + str(self.temp_center) + ' K, '
        s += 'L= ' + str(self.luminosity) + ' erg/s, '
        s += 'R= ' + str(self.radius) + ' cm '
        return s


    def __lt__(self, other):
        if not compnn(self.solar_masses, other.solar_masses):
            if self.solar_masses < other.solar_masses:
                return True
            return False

        if not compnn(self.H_frac, other.H_frac):
            if self.H_frac < other.H_frac:
                return True
            return False

        if not compnn(self.He_frac, other.He_frac):
            if self.He_frac < other.He_frac:
                return True
            return False

        if not compnn(self.pressure_center, other.pressure_center):
            if self.pressure_center < other.pressure_center:
                return True
            return False

        if not compnn(self.temp_center, other.temp_center):
            if self.temp_center < other.temp_center:
                return True
            return False


      #  if not compnn(self.luminosity, other.luminosity):
      #      if self.luminosity < other.luminosity:
      #          return True
      #      return False
      #
      #  if not compnn(self.radius, other.radius):
      #      if self.radius < other.radius:
      #          return True
      #      return False

    #def __gt__(self,other):

    def __eq__(self,other):
        b = compnn(self.solar_masses, other.solar_masses)
        b *= compnn(self.H_frac, other.H_frac)
        b *= compnn(self.He_frac, other.He_frac)
        b *= compnn(self.pressure_center, other.pressure_center)
        b *= compnn(self.temp_center, other.temp_center)
#        b *= compnn(self.luminosity, other.luminosity)
#        b *= compnn(self.radius, other.radius) #only the important parts......
        return b

    def __hash__(self):
        return(hash(str(self)))

    def get_final_teff(self):
        try:
            f = open(self.filename)
        except FileNotFoundError:
            try:
                f = open('data/'+self.filename)
            except:
                return None
        for row in f:
            if "Teff" in row:
                 self.Teff = float(row.split(',')[0].replace('Teff:','').replace('D','e'))
                 return

    def get_final_luminosity(self):
        try:
            f = open(self.filename)
        except FileNotFoundError:
            try:
                f = open('data/'+self.filename)
            except:
                return None
        flag = False
        for row in f:
            if flag:
                # Pc: 1.5624D+17, Tc: 1.5451D+07, R: 7.2280D+10, L: 4.0480D+33
                self.luminosity = float(row.split(',')[-1].replace('L:','').replace('D','e'))
                return
            else:
                if "FINAL MODEL" in row:
                    flag = True

    def log_l_lsun(self):
        return np.log(self.luminosity/self.L_sun)
