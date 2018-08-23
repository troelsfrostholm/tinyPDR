# Program that tabulates the ISRF from Black 1987

from numpy import *

h = 4.13566766e-15    # Planck constant (eV s)
c = 299792458         # Speed of light (m s-1)
k = 8.61733034e-5     # Boltzmann constant (eV K-1)
eV2erg = 1.60217662e-12 # eV in ergs

mathis83_ev = loadtxt("../dat/mathis83_ev.dat")

# Planck's law in terms of temperature (K) and photon energy (eV)
def planck(T, E):
  return 2*E**3/h**2/c**2/(exp(E/(k*T))-1.0)# * eV2erg

# Soft X-ray flux in 54-500 eV range (ergs cm-2 s-1 Hz-1)
# as function of photon energy (eV)
# from Bregman and Harrington 1985 eq. 8
def F_soft_xray(E):
  F = 2.9e-25*(E/1e3)**(-2) / eV2erg
  ioutside = logical_or(E < 54.0, E > 200.0)
  F[ioutside] = 0.0
  return F

# Interpolates in tabulated intensity from Mathis 1983 (ergs cm-2 s-1 str-1 Hz-1)
# E is in eV
def J_mathis83(E):
  logJ = log(mathis83_ev)
  print logJ[:,0]
  print logJ[:,1]
  print log(E)
  logJinterp = interp(log(E), logJ[:,0], logJ[:,1], left=-1e100, right=-1e100)
  print logJinterp
  return exp(logJinterp)

# Mean intensity of ISRF from Black 1987 (ergs cm-2 s-1 str-1 Hz-1)
# E is in eV
def Black87(E):
  return planck(2.7, E) + J_mathis83(E) + F_soft_xray(E)/(4*pi)
  #return J_mathis83(E) + F_soft_xray(E)/(4*pi)

# Limits and number of points
Elow = 5e-6 # eV
Ehigh = 500 # eV
N = 50

# Logarithmic spaced energy axis
E = exp(linspace(log(Elow),log(Ehigh),N))
J = Black87(E)
data = empty((N,2))
data[:,0] = E
data[:,1] = J
savetxt("../dat/black87_eV_new.dat", data)