# zamscript

## About

This script serves as a "wrapper" for the ZAMS executable. It will iterate over the parameter space and call the ZAMS.exe for each set of parameters. Be sure to execute "parameter_space.py" and have the modified version of zams.for compiled before running the autom.py.
The existing version of the code was built on Win10 64, gfortran on Mingw. Due to some details of the *os* library, I can't guarantee if it works on other OSs from the get-go.
The fortran 77 code comes from Hansen & Kawaler 94, modified to receive the parameters from a text file instead of prompting the user. The automation doesn't do the file handling for pulsation output.

## Instructions:

1. Compile the modified ZAMS.for to create ZAMS.exe;
2. Edit and run parameter_space.py to generate parameter_space.txt;
3. Run autom.py;
4. Wait. A list of succesfull models will be saved at unique_stars_final_model;
5. The \*.dat files remaining contain the interior model of each valid star;
6. Star_analisys provides some plots.
