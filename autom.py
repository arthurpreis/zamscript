import os
import numpy
import file_filter

####### ZAMS AUTOMATION ##############
#Arthur Reis - 2020
#
#This script serves as a "wrapper" for the ZAMS executable.
#It will iterate over the parameter space and call the ZAMS.exe
#for each set of parameters.
#Be sure to execute "parameter_space.py" and have the 
#modified version of zams.for compiled before running the autom.py.
#This setup was lastly succesfuly tested on Win10 64. Some workaround w/ the os library
#may be needed on other OSs.

#Parameter space is generate by parameter_space.py
input_parms = open("parameter_space.txt", "r")
output_candidates = open('candidates.txt', 'w')

i = 0
for p in input_parms:
    if p[0] == '#' or not p.strip():
        pass
    else:
        if os.path.exists("PARM.txt"):
            #Cleaning parm.txt. Make sure to always create a fresh one-liner file
            file = open("PARM.txt",'w+')
            file.close()
            os.remove("PARM.txt")

        #PARM.txt is the file that the ZAMS.exe will read
        out_parm = open("PARM.txt", 'w+')

        outfile = str(i) + '.dat' #output file name
        out_parm.write(p + outfile + ', N')
        #parm.txt = parameters, name of the file, NO pulsation calculations

        out_parm.close()
        os.system(r'ZAMS.exe') #works on win10. May need rework on linux.
        
        #Quick and dirty filter:
        #Runs that fail to reach final model, will often be small
        if os.path.exists(outfile) and os.stat(outfile).st_size < 12000:
            os.remove(outfile)    #Test for size, remove smaller files
        elif os.path.exists(outfile) and os.stat(outfile).st_size > 12000:
            if not file_filter.check(outfile): 
                os.remove(outfile) #removes big files w/ lots of trial runs,
                #but w/ no final model

    i+=1

