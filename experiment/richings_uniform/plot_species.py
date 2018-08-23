import sys
sys.path.append('../../utils')
import numpy as np
from matplotlib import pyplot as plt
from cycler import cycler
import richings as r
import KromeInfo as I
import reader

info = I.KromeInfo()
run = reader.Run("output.dat")
abundances = run.getMolecules(itime=run.ntime-1)
axis = run.getAv(itime=run.ntime-1)

s = "H"
nHtot = 1e2

if len(sys.argv) > 1:
  s = sys.argv[1]

plt.loglog(axis, abundances[:,info.symbol2i(s)]/nHtot, label="krome")

def isNumber(x):
  try:
    int(x)
    return True
  except:
    return False

if s == "E":
  richsym = "electrons"
elif "-" in s:
  richsym = s
else:
  s = s.replace("HE", "He")
  base = s.replace("+", "")
  if(len(base)==1 or (len(base)==2 and (base[-1]==base[-1].lower() and not isNumber(base[-1])))):
    richsym = base + "I"*(s.count("+")+1)
  else:
    richsym = s.replace("+", "p")

try:
  r.plot_abundance(richsym, linestyle="--", plotting_function=plt.loglog, label="richings")
except:
  print "Warning. Species missing in Richings data: " + richsym + " ("+s+")"

plt.legend()
plt.show()