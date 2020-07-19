import numpy

#creates the file and populate w/ parms sets 
#parameter_space will be iterated by autom.py
#number of points in spaces will 
#greatly increase computational time
f = open('parameter_space.txt', 'w')

#Parameters and format accordingly to ZAMS.for
#The values below were 'inspired' by Hanser & Kawaler

#Central Temperature (K)
temperatures = numpy.logspace(5, 7, num = 3)

#Central Pressure (cgs)
pressures = numpy.logspace(13, 16, num = 4)

#Mass (M_sun)
mass = [0.1, 0.5, 1, 2, 5, 10, 50]
#mass = numpy.linspace(0.7, 0.8, num = 10)

#H %, He%
Xs = numpy.linspace(0.7, 0.75, num = 5)
Ys = numpy.linspace(0.25, 0.3, num = 5)

#Luminosity (L_sun)
Ls = numpy.logspace(-3,5, num = 9)

#Radius (cm)
Rs = numpy.logspace(9,11, num = 5)

for t in temperatures:
    for p in pressures:
        for x in Xs:
            for y in Ys:
                for l in Ls:
                    for r in Rs:
                        for m in mass:
                            #print(m, x, y, p, t, r, l)
                            s  = "{0:.2f}".format(m) + ', '
                            s += "{0:.2f}".format(x) + ', '
                            s += "{0:.2f}".format(y) + ', '
                            s += "{0:.2E}".format(p) + ', '
                            s += "{0:.2E}".format(t) + ', '
                            s += "{0:.2E}".format(r) + ', '
                            s += "{0:.2f}".format(l) + '\n'
                            #print(s)
                            f.write(s)
