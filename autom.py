import os
import numpy

# open list of parms
#write to PARM.txt
#call zams
#rinse and repeat

def deleteContent(fName):
    with open(fName, "w"):
        pass

deleteContent('err_pass.txt')
input_parms = open("script_parameters.txt", "r")

i=0
for p in input_parms:
    if p[0] == '#' or not p.strip():
        pass
    else:
        if os.path.exists("PARM.txt"):
            os.remove("PARM.txt")
        out_parm = open("PARM.txt", 'w+')
        out_parm.write(p + ', teste_' + str(i) + '.dat, n')
        out_parm.close()
        os.system(r'./\zams')
    i+=1
