import sys
sys.path.append('../utils')
import numpy as np
from matplotlib import pyplot as plt
import reader

species = "H2"
N = 5
if len(sys.argv) > 1:
  species = sys.argv[1]

if len(sys.argv) > 2:
  N = int(sys.argv[2])

print "Plotting top", N, "reactions with", species

isnapshot = 100

fluxfile = reader.FluxFile("fluxes.dat")
reactions = fluxfile.getReactions()
with_H2 = reactions.selectByMolecule(species).selectByTimeStep(99,100).top(N)
for (name, flux) in with_H2.iter():
  plt.loglog(fluxfile.av, flux[0,:], label=name)
plt.legend()
plt.show()
