#!/usr/local/bin/python
import sys
sys.path.append('../../utils')
import numpy as np
from matplotlib import pyplot as plt
from cycler import cycler
import richings as r
import KromeInfo as I
import reader

run = reader.Run("output.dat")
itime = run.ntime-1

if len(sys.argv) > 1:
  itime = int(sys.argv[1])

species = ["H"]

if len(sys.argv) > 2:
  species = sys.argv[2:]

info = I.KromeInfo()
abundances = run.getMolecules(itime=itime)
axis = run.getR(itime=itime)

ntot = np.sum(abundances, axis=1)
T = run.getTgas(itime=itime)

plt.figure()
plt.semilogy(axis, ntot, label=r"$N_{tot}$")
plt.semilogy(axis, T, label=r"$T_{gas}$")
plt.legend()

plt.savefig("plots/N_T.pdf")