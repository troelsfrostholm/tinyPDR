import sys
sys.path.append('../../utils')
import numpy as np
from matplotlib import pyplot as plt
from cycler import cycler
import richings as r
import KromeInfo as I

def ionized(name):
  return "II" in name or \
         "V"  in name or \
         "X"  in name or \
         "p"  in name or \
         ("I" in name and "V" in name) or \
         ("I" in name and "X" in name)

threshold = 1e-7
names = r.getNames()

for name in names:
  if not ionized(name): continue
  n = r.get_abundance(r.ilast, name)
  if np.all(n < threshold): continue
  r.plot_abundance(name, linestyle="--", plotting_function=plt.loglog, label=name)
plt.legend()
plt.show()
exit()


info = I.KromeInfo()

def isNumber(x):
  try:
    int(x)
    return True
  except:
    return False

for s in map(lambda x: x.symbol, info.species):
  if not "+" in s: continue
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
    r.plot_abundance(richsym, linestyle="--", plotting_function=plt.loglog, label=s)
  except:
    print "Warning. Species missing in Richings data: " + richsym + " ("+s+")"

plt.legend()
plt.show()