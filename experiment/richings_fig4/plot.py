import sys
sys.path.append('../../utils')
import numpy as np
from matplotlib import pyplot as plt
from cycler import cycler
import richings as r
import KromeInfo as I
import reader

# Conversion factor from total H column density -> Av
fav = 1e0/4e-22

info = I.KromeInfo()
run = reader.Run("output.dat")
abundances = run.getMolecules(itime=run.ntime-1)
axis = run.getAv(itime=run.ntime-1)

datadir = "output_black"
if(len(sys.argv)>1):
  datadir = sys.argv[1]

# Plot number density of species
f, axarr = plt.subplots(3, sharex=True, figsize=(8,9))

def totalNumberDensity(atom):
  ids, counts = info.countAtomInAllSpecies(atom)
  ids = map(lambda x: x - 1, ids)
  return np.sum(abundances[:,ids]*counts, axis=1)

atoms = ["H", "He", "C", "O"]
ntot_atom = { a : totalNumberDensity(a) for a in atoms}

# Correct for ices
ice = False
if ice:
  iH2O       = info.symbol2i("H2O")
  iH2O_total = info.symbol2i("H2O_total")
  iCO        = info.symbol2i("CO")
  iCO_total  = info.symbol2i("CO_total")

  ntot_atom["H"] += 2*(abundances[:,iH2O_total] - abundances[:,iH2O])
  ntot_atom["O"] += abundances[:,iH2O_total] - abundances[:,iH2O] + abundances[:,iCO_total] - abundances[:,iCO]
  ntot_atom["C"] += abundances[:,iCO_total] - abundances[:,iCO]

#colors = ['black', 'brown', 'red', 'green', 'purple', 'orange', 'blue', 'magenta', 'lightblue']
#species =       ["E", "H", "H+", "He", "He+", "HCO+", "C+", "C"]
#normalized_to = ["H", "H", "H" , "He", "He" , "C"   , "C" , "C"]
colors = ['black', 'brown', 'blue']
species =       ["E", "H", "C+"]
normalized_to = ["H", "H", "C" ]
axarr[1].set_prop_cycle(cycler('color', colors))
for s, n in zip(species, normalized_to):
  axarr[1].loglog(axis, abundances[:,info.symbol2i(s)]/ntot_atom[n], label=s)

# Plot series from R14II fig. 4 for comparison
R14II_dir = "R14IIfig4/"
R_HI = np.loadtxt(R14II_dir+"HI.csv")
axarr[1].semilogy(R_HI[:,0]/fav, R_HI[:,1], color="brown", linestyle="--", label="Richings")
R_CII = np.loadtxt(R14II_dir+"CII.csv")
axarr[1].semilogy(R_CII[:,0]/fav, R_CII[:,1], color="blue", linestyle="--")
R_E = np.loadtxt(R14II_dir+"E.csv")
axarr[1].semilogy(R_E[:,0]/fav, R_E[:,1], color="black", label="Richings", linestyle="--")

plt.xlim(min(axis),max(axis))
axarr[1].set_ylim(3e-5,2)
plt.xlabel(r"$A_v$")
axarr[1].set_ylabel("Ionization Fraction")
axarr[1].legend(loc="best")

# Plot molecular fractions
species = ["H2", "CO", "OH"]
normalized_to = ["H", "C", "O"]
colors = ['black', 'blue', 'red', 'brown', 'purple', 'orange']
axarr[0].set_prop_cycle(cycler('color', colors))
abundances[:,info.symbol2i("H2")] *= 2.0         # Like Richings, we plot 2*n_H2/n_Htot
for s, n in zip(species, normalized_to):
  axarr[0].loglog(axis, abundances[:,info.symbol2i(s)]/ntot_atom[n], label=s)

if ice:
  nH2OIce = abundances[:,info.symbol2i("H2O_total")] - abundances[:,info.symbol2i("H2O")]
  axarr[0].loglog(axis, nH2OIce/ntot_atom["O"], label="H2O ice")
  nCOIce = abundances[:,info.symbol2i("CO_total")] - abundances[:,info.symbol2i("CO")]
  axarr[0].loglog(axis, nCOIce/ntot_atom["C"], label="CO ice")
  xH2O = abundances[:,info.symbol2i("H2O")] / abundances[:,info.symbol2i("H2O_total")]
  axarr[0].loglog(axis, xH2O, label="H2O")

# Plot series from R14II fig. 4 for comparison
R_H2 = np.loadtxt(R14II_dir+"H2.csv")
axarr[0].semilogy(R_H2[:,0], R_H2[:,1], label="Richings", color="black", linestyle="--")
R_CO = np.loadtxt(R14II_dir+"CO.csv")
axarr[0].semilogy(R_CO[:,0], R_CO[:,1], color="blue", linestyle="--")
R_OH = np.loadtxt(R14II_dir+"OH.csv")
axarr[0].semilogy(R_OH[:,0], R_OH[:,1], color="red", linestyle="--")

axarr[0].set_ylim(1e-9,1)
axarr[0].set_ylabel("Molecular Fraction")
axarr[0].legend(loc='best')

# Plot density and temperature
T = np.loadtxt(R14II_dir+"/T.csv")
# axarr[2].loglog(snapshot_left[:,2], snapshot_left[:,3])
# axarr[2].loglog(snapshot_right[:,2], snapshot_right[:,3])
axarr[2].loglog(T[:,0]/fav, T[:,1], color='black')
axarr[2].set_ylabel("T (K)")

# And density
n = np.loadtxt(R14II_dir+"/n.csv")
ax_n = axarr[2].twinx()
ax_n.loglog(n[:,0]/fav, n[:,1], label=r"$n_{H_{tot}}$", color='lightseagreen')
ax_n.loglog(0,0,color='black',label="T")   # Hack: Add empty series to get a label for the temperature
ax_n.set_ylabel(r"$cm^{-3}$")
ax_n.legend(loc=3)

f.subplots_adjust(hspace=0)

plt.savefig("plots/R14II-4.pdf")

plt.figure()
tau = run.getTau(itime=run.ntime-1)
for i in xrange(run.nbins):
  plt.loglog(axis, tau[:,i], label=i)
plt.xlabel(r"$A_v$")
plt.ylabel(r"$\tau$")
plt.legend()
plt.savefig("plots/tau.pdf")

iHejMax = np.argmax(abundances[:,info.symbol2i("He+")])

plt.figure()
plt.loglog(axis, ntot_atom["He"], label="nHetot")
plt.loglog(axis, abundances[:,info.symbol2i("He")], label="He")
plt.loglog(axis, abundances[:,info.symbol2i("He+")], label="He+")
plt.ylim(1e-1,20)
plt.legend()
plt.savefig("plots/He.pdf")

plt.figure()
plt.loglog(axis, run.getTgas(itime=run.ntime-1))
plt.loglog(axis, run.getTdust(itime=run.ntime-1))
plt.savefig("plots/Tgas.pdf")

