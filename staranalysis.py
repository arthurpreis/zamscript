from star import Star
from settings import Settings
import matplotlib.pyplot as plt

def plot_each_metallicity(metal, group, fields, i = 0, fieldname = '',
        filenametag = '', ylabel = '', save = False):
    (x,y) = metal
    z = 1 - x - y

    chosen_stars = []
    for s in group:
        if s.H_frac == x and s.He_frac == y:
            chosen_stars.append(set_fields(s))

    packed_chosen_stars = []
    for star in chosen_stars:
        packed_chosen_stars.append(pack_star(star, *fields))


    title = 'Distribuição de ' + fieldname +' na estrutura estelar \n para X= '
    title += str(x) + ', Y=' + str(y) + ', Z=' + str(round(z,2))
    filename = filenametag + str(x).replace('.','') + str(y).replace('.','')

    plot_graf(*packed_chosen_stars, i = i, xlabel = r'$x_i = 1-M(r)/M$',
        ylabel = ylabel, xlim = (1,0), mode = 'log', labeltag = 'M=',
        title = title, fn = filename, save = save)




def pack_star(star, attr1, attr2, label = 'solar_masses', marker = 'o'):
    a = getattr(star, attr1)
    b = getattr(star, attr2)
    l = str(getattr(star, label))
    m = marker
    return (a,b,l,m)

def choose_stars(group, list_indexes):
    c = []
    for i in list_indexes:
        s = set_fields(group[i])
        c.append(s)
    return c

def plot_graf(*args, i = 0, fn = 'dummy', xlabel = '', ylabel = '',
        labeltag = '', title = 'Title', mode = 'lin', ylim = (0,0), xlim = (0,0), save = False):
    fig = plt.figure(i)
    ax = plt.gca()
    ax.set_title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    if mode =='log':
        ax.set_yscale("log")
    if mode == 'dlog':
        ax.set_xscale('log')
        ax.set_yscale('log')
    if mode == 'dloginv':
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.gca().invert_xaxis()
    if ylim != (0,0):
        ax.set_ylim(ylim)
    if xlim != (0,0):
        ax.set_xlim(xlim)
    for arg in args:
        (x, y, label, marker) = arg
        #plt.xticks(np.arange(min(x), max(x)+1, 0.1))
        #ax.set_xlim((0.0,1.0))
        #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        #print(type(x), type(y))
        ax.scatter(x, y, label = labeltag + label)

    plt.legend()
    if save:
        plt.savefig(fn)
    else:
        plt.show()
    plt.close(fig)

def set_fields(star):
    sout = Star(settings)
    sout = star
    f = open('./data/' + star.filename)
    offset = 0
    for row in f:
        if offset == 0 and "FINAL MODEL" in row:
            offset = 1
            continue

        if offset == 1:
            sout.luminosity = float((row.split(','))[-1].replace('L:', '').replace('D','e'))
            offset = 2
            continue

        if offset == 2:
            sout.Teff = float((row.split(','))[0].replace('Teff:', '').replace('D','e'))
            offset = 3
            continue

        if offset == 3:
            offset = 4
            i=1
            continue

        if offset == 4:
            s = row.replace('D','e').split()
            sout.coordinate_field[i]  = float(s[1])
            sout.radius_field[i]      = 10**float(s[2])
            sout.pressure_field[i]    = 10**float(s[3])
            sout.temperature_field[i] = 10**float(s[4])
            sout.density_field[i]     = 10**float(s[5])
            sout.luminosity_field[i]  = 10**float(s[6])
            i += 1
            if s[0] == '200':
                break

    sout.coordinate_field[0]  = sout.coordinate_field[1]
    sout.radius_field[0]      = sout.radius_field[1]
    sout.pressure_field[0]    = sout.pressure_field[1]
    sout.temperature_field[0] = sout.temperature_field[1]
    sout.density_field[0]     = sout.density_field[1]
    sout.luminosity_field[0]  = sout.luminosity_field[1]

    sout.coordinate_field[-1]  = sout.coordinate_field[-2]
    sout.radius_field[-1]      = sout.radius_field[-2]
    sout.pressure_field[-1]    = sout.pressure_field[-2]
    sout.temperature_field[-1] = sout.temperature_field[-2]
    sout.density_field[-1]     = sout.density_field[-2]
    sout.luminosity_field[-1]  = sout.luminosity_field[-2]

    return sout

settings = Settings()
star_group = []

infile = open('unique_stars_final_model.txt','r')
for row in infile:
    s = row.split(',')
    star = Star(settings, float(s[0]), float(s[1]), float(s[2]),
        float(s[3]), float(s[4]), float(s[5]), float(s[6]))
    star.filename = s[7].strip()
    star_group.append(star)

star_group.sort()
print(len(star_group))
print(star_group)

metalicities = [(0.70,0.25),
                (0.70,0.26),
                (0.70,0.28),
                (0.70,0.29),
                (0.70,0.30),
                (0.71,0.25),
                (0.71,0.26),
                (0.71,0.28),
                (0.71,0.29),
                (0.72,0.25),
                (0.72,0.26),
                (0.72,0.28),
                (0.74,0.25),
                (0.74,0.26),
                (0.75,0.25)
                ]
i = 0
for m in metalicities:
    plot_each_metallicity(m, star_group, ('coordinate_field',
        'temperature_field'), fieldname = 'temperaturas', filenametag = 'temp', ylabel = r'$T(x_i) [K]$',
        i = i + 0,  save = True)
    plot_each_metallicity(m, star_group, ('coordinate_field',
        'pressure_field'), fieldname = 'pressões', filenametag = 'press', ylabel = r'$P(x_i) [dina/cm^2]$',
        i = i + 1,  save = True)
    plot_each_metallicity(m, star_group, ('coordinate_field',
        'density_field'), fieldname = 'densidades', filenametag = 'dens', ylabel = r'$\rho(x_i) [g/cm^3]$',
        i = i + 2,  save = True)
    plot_each_metallicity(m, star_group, ('coordinate_field',
        'luminosity_field'), fieldname = 'luminosity', filenametag = 'lumn', ylabel = r'$L(x_i) [erg/s]$',
        i = i + 3,  save = True)
    i+=4
