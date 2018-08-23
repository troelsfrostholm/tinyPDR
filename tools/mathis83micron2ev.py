# Tool for computing the interstellar radiation field in the solar vicinity

from numpy import *

h = 4.13566766e-15    # Planck constant (eV s)
c = 299792458         # Speed of light (m s-1)

data = loadtxt("../dat/mathis83_micron.dat")
l = data[:,0]  # Wavelength (lambda)

# Convert integrated intensity (erg cm-2 s-1 micrometers-1) to mean intensity (erg cm-2 s-1 str-1 Hz-1)
# Recall lambda = c/nu => d(lambda) = lambda**2/c*d(nu)
data[:,1] = data[:,1]*1e6*l**2/c/(4*pi)

# Convert wavelength (micrometers) to photon energy (eV)
data[:,0] = c/(data[:,0]/1e6)*h

i = argsort(data[:,0]) # Sort according to photon energy (increasing)
data[:,0] = data[i,0]
data[:,1] = data[i,1]

header = '# ISRF in the solar vicinity from\n\
# Mathis, Mezger and Panagia 1983 table A3 and B1\n\
# First column: photon energy (eV)\n\
# Second column: J (erg cm-2 s-1 str-1 Hz-1)\n\n'

with open("../dat/mathis83_eV.dat", 'w') as f:
  f.write(header)
  savetxt(f,data)
