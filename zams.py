import numpy as np
from settings import Settings
from star import Star

settings = Settings()
star = Star(settings)

# ========================== CONTROL ===================================
#c The independent variable is xi=ln(1-Mr/M) and the dependent
#c variables are ln P, ln T ln r, and ln L (in cgs units).
#c There are 201 points in the model and the mesh in xi is
#c determined for you.

#c The fitting point is also chosen
#c for you because convergence may depend on where it is set.
#c The fitting point is specified by ``QFIT'' which is the
#c fractional mass interior to the fitting point.

#Qfit determination:
t1 = np.log10(1.5)
t2 = 0.6/(1.0-t1)
if star.solar_masses >= 10:
    Qfit = 0.8
elif star.solar_masses <= 1.5:
    Qfit = 0.2
else:
    Qfit = 0.2 + (log10(star.solar_masses) - t1)*t2

#c Set up XI and DXI arrays. The DXI are set up as
#c two geometric series. XI(N) is undefined. XI(N-1) is
#c at the mass 1-Mr/M=1.d-10. There is no simple way
#c to set up the zoning so as to satisfy all kinds of
#c models and pulsation codes. A logarithmic mesh
#c of this sort is discussed in Kippenhahn, Weigert,
#c and Hofmeister (1967) (see chapter 7).

nf_out = 0
n = settings.mesh_points
nc = settings.mesh_centerfold

x = np.empty(n)
dx = np.empty(n)

eps1 = 0.035
epsp1 = eps1 + 1
eps2 = -0.15
epsp2 = eps2 + 1
x[n-1] = np.log(1e-10)
t1 = (epsp1**(nc -1) - 1)/eps1
t2 = (epsp1**(nc-1))*epsp2*(epsp2**(n - 1 - nc) -1)/eps2
print(t1,t2)
x[0] = 1e-10
dx[0] = x[n-1]/(t1+t2)

print(1, x[0],dx[0])
for i in range(1,n): #ta funcionando a menos de um epsilon ahahha
    if i<nc:
        dx[i] = dx[0]*(epsp1**(i-1))
        x[i] = dx[1]*(epsp1**(i-1))/eps1
    elif i == nc:
        x[nc] = x[nc-1] + dx[nc-1]
        dx[nc] = dx[nc-1]*epsp2
    elif i > nc:
        J = i-nc+1
        dx[i]=dx[nc-1]*epsp2**J
        TEMP=epsp2*(epsp2**J-1)/eps2
        x[i] = x[nc] + dx[nc-1]*TEMP
    if x[i] < np.log(1-Qfit) and nf_out == 0:
        nf_out = i

    print(i+1, x[i], dx[i])



#c The independent variable is xi=ln(1-Mr/M) and the dependent
#c variables are ln P, ln T ln r, and ln L (in cgs units).
#c There are 201 points in the model and the mesh in xi is
#c determined for you. The fitting point is also chosen
#c for you because convergence may depend on where it is set.
#c The fitting point is specified by ``QFIT'' which is the
#c fractional mass interior to the fitting point.


# These are the interior masses.
for i in range(0,n):
    if i == 0:
        mass = 0
    elif i < n:
        mass = star.mass*(1 - np.exp(x[i]))
    elif i == n:
        mass = 1

    star.mass_field[i] = mass
nf_in = n - nf_out + 1

#SUBROUTINE GOOUT (NFOUT, NGO)
#c ***********************************************
#c This is the driver routine for the outwards integration.
#c It does a first step using the central expansions.
#c There are three different integrations performed:
#c (1) vary the central pressure from its guessed value;
#c (2) vary the central temperature; (3) use the guessed
#c values for both central pressure and temperature. All
#c three are usually done (when NGO is 1) except when the
#c problem has converged (NGO is 2). Only (3) is done in
#c that case (to get a final model).
#c ***********************************************

# ==============================GOOUT===================================

# ==============================GOIN====================================

# ==============================CORRECT=================================
