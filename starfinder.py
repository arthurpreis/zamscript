import os
from star import Star
from settings import Settings

settings = Settings()
star_group = []


for f in os.listdir("./"):
    if f.endswith(".dat"):
        print(f)
        data = open(f,'r')
        s = ''
        nextrow = False
        for row in data:
            if "MASS/MSUN" in row:
                s+= row
            if "FINAL MODEL" in row:
                nextrow = True
                continue
            if nextrow:
                s += row
                nextrow = False

#        print(s)
#s: M...N=  2.000, X= 0.670, Y= 0.130, Pc: 6.2312D+16, Tc: 1.8647D+07, R: 1.2960D+11, L: 2.7493D+34

        mass = float(s[12:18])
        X = float(s[23:28])
        Y = float(s[33:38])
        Pc = float(s[43:54].replace('D','e'))
        Tc = float(s[59:70].replace('D','e'))
        R = float(s[74:85].replace('D','e'))
        L = float(s[89:99].replace('D','e'))

        star = Star(settings, mass, X, Y, Pc, Tc, L, R)
        star.filename = f
        star_group.append(star)

star_group.sort()

uniq = []
for i in range(len(star_group)):
    not_seen = True
    for j in range(i):
        if j != i and star_group[i] == star_group[j]:
            not_seen = False
    if not_seen:
        uniq.append(star_group[i])


print(len(uniq))

out_file = open('unique_stars_final_model.txt','w+')
for s in uniq:
    out_file.write(s.out_parm() + ', ' + s.filename + '\n')


