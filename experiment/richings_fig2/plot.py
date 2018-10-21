import sys
sys.path.append('../../utils')
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from cycler import cycler
import richings as r
import KromeInfo as I
import reader

# Parameters
ice = False
show_dust_T = False

info = I.KromeInfo()
run = reader.Run("output.dat")
abundances = run.getMolecules(itime=run.ntime-1)
axis = run.getAv(itime=run.ntime-1)

datadir = "output_black"
if(len(sys.argv)>1):
  datadir = sys.argv[1]

# Plot number density of species
Nplots = 2
if show_dust_T: Nplots += 1

fig = plt.figure(figsize=(8,10))
f = fig
hratios = [2,2,1]
subplots = gridspec.GridSpec(Nplots, 1, height_ratios=hratios[0:Nplots])
axarr = range(3)
axarr[0] = fig.add_subplot(subplots[0])
for i in xrange(Nplots-1):
  axarr[i+1] = fig.add_subplot(subplots[i+1], sharex=axarr[0])
#top_ax = fig.add_subplot(top_pos)
#bot_ax = fig.add_subplot(bot_pos, sharex=top_ax)

#f, axarr = plt.subplots(Nplots, sharex=True, figsize=(8,9))

def totalNumberDensity(atom):
  ids, counts = info.countAtomInAllSpecies(atom)
  ids = map(lambda x: x - 1, ids)
  return np.sum(abundances[:,ids]*counts, axis=1)

atoms = ["H", "He", "C", "O"]
ntot_atom = { a : totalNumberDensity(a) for a in atoms}

# Correct for ices
if ice:
  iH2O       = info.symbol2i("H2O")
  iH2O_total = info.symbol2i("H2O_total")
  iCO        = info.symbol2i("CO")
  iCO_total  = info.symbol2i("CO_total")

  ntot_atom["H"] += 2*(abundances[:,iH2O_total] - abundances[:,iH2O])
  ntot_atom["O"] += abundances[:,iH2O_total] - abundances[:,iH2O] + abundances[:,iCO_total] - abundances[:,iCO]
  ntot_atom["C"] += abundances[:,iCO_total] - abundances[:,iCO]

colors = ['black', 'brown', 'red', 'purple', 'orange', 'blue', 'magenta', 'lightblue', 'grey']
species =       ["E", "H", "H+", "He+",  "C", "C+", "O", "O+", "HCO+"]
normalized_to = ["H", "H", "H" , "He" ,  "C", "C" , "O", "O", "C"   ]
axarr[1].set_prop_cycle(cycler('color', colors))
for s, n in zip(species, normalized_to):
  axarr[1].loglog(axis, abundances[:,info.symbol2i(s)]/ntot_atom[n], label=s)

r.plot_abundance('electrons', color="black", linestyle="--", plotting_function=axarr[1].loglog)
r.plot_abundance('HI', color="brown", linestyle="--", plotting_function=axarr[1].loglog)
r.plot_abundance('HII', color="red", linestyle="--", plotting_function=axarr[1].loglog)
#r.plot_abundance('HeI', relative_to="He", color="green", linestyle="--", plotting_function=axarr[1].loglog)
r.plot_abundance('HeII', relative_to="He", color="purple", linestyle="--", plotting_function=axarr[1].loglog)
r.plot_abundance('CI', relative_to="C", color="orange", linestyle="--", plotting_function=axarr[1].loglog)
r.plot_abundance('CII', relative_to="C", color="blue", linestyle="--", plotting_function=axarr[1].loglog)
r.plot_abundance('OI', relative_to="O", color="magenta", linestyle="--", plotting_function=axarr[1].loglog)
r.plot_abundance('OII', relative_to="O", color="lightblue", linestyle="--", plotting_function=axarr[1].loglog)
r.plot_abundance('HCOp', relative_to="C", color="grey", linestyle="--", plotting_function=axarr[1].loglog)

axarr[1].set_ylim(3e-8,2)
#plt.xlim(1e-6,1e2)
plt.xlabel(r"$A_v$")
axarr[1].set_ylabel("Ionization Fraction")
axarr[1].legend(loc="best")

# Plot molecular fractions
species = ["H2", "CO", "OH", "H2O"]
normalized_to = ["H", "C", "O", "O"]
colors = ['black', 'blue', 'red', 'brown', 'purple', 'orange']
axarr[0].set_prop_cycle(cycler('color', colors))
abundances[:,info.symbol2i("H2")] *= 2.0         # Like Richings, we plot 2*n_H2/n_Htot
for s, n in zip(species, normalized_to):
  axarr[0].loglog(axis, abundances[:,info.symbol2i(s)]/ntot_atom[n], label=s)

if ice:
  nCOIce = abundances[:,info.symbol2i("CO_total")] - abundances[:,info.symbol2i("CO")]
  axarr[0].loglog(axis, nCOIce/ntot_atom["C"], label="CO ice")
  nH2OIce = abundances[:,info.symbol2i("H2O_total")] - abundances[:,info.symbol2i("H2O")]
  axarr[0].loglog(axis, nH2OIce/ntot_atom["O"], label="H2O ice")

r.plot_abundance('H2', factor=2.0, color="black", label="Richings", linestyle="--", plotting_function=axarr[0].loglog)
r.plot_abundance('CO', relative_to="C", color="blue", linestyle="--", plotting_function=axarr[0].loglog)
r.plot_abundance('OH', relative_to="O", color="red", linestyle="--", plotting_function=axarr[0].loglog)
r.plot_abundance('H2O', relative_to="O", color="brown", linestyle="--", plotting_function=axarr[0].loglog)

#Plot series from R14II fig. 2 for comparison
cloudystyles = {
  "big_H2" : ("black", ":"),
  "small_H2" : ("black", "-."),
  "big_CO" : ("blue", ":"),
  "small_CO" : ("blue", "-."),
  "big_OH" : ("red", ":"),
  "small_OH" : ("red", "-.")
}
for name, data in r.cloudy.iteritems():
  s = cloudystyles[name]
  label = ""
  if name=="big_H2": label="Cloudy big"
  if name=="small_H2": label="Cloudy small"
  axarr[0].loglog(data[:,0], data[:,1], color=s[0], linestyle=s[1], label=label)

axarr[0].set_ylim(1e-6,2)
axarr[0].set_ylabel("Molecular Fraction")
axarr[0].legend(loc='best')

if show_dust_T:
  dustT = run.getTdust(itime=run.ntime-1)
  axarr[2].semilogx(axis, dustT, label="dust temperature")
  axarr[2].set_ylabel(r"$T_d$ (K)")
  axarr[2].legend(loc='best')

f.subplots_adjust(hspace=0)

#plt.show()
plt.savefig("plots/R14II-2.pdf")

plt.figure()
tau = run.getTau(itime=run.ntime-1)
for i in xrange(run.nbins):
  plt.loglog(axis, tau[:,i], label=i)
plt.xlabel(r"$A_v$")
plt.ylabel(r"$\tau$")
plt.legend()
plt.savefig("plots/tau.pdf")

iHejMax = np.argmax(abundances[:,info.symbol2i("He+")])
print "iHejMax", iHejMax

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
plt.tight_layout()
plt.savefig("plots/Tgas.pdf")

