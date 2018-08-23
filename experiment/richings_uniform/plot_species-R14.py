import sys
sys.path.append('../../utils')
import numpy as np
from matplotlib import pyplot as plt
from cycler import cycler
import richings as r
import KromeInfo as I


threshold = 1e-5
names = r.getNames()

name = "MgI"

if len(sys.argv) > 1:
  name = sys.argv[1]

r.plot_abundance(name, linestyle="--", plotting_function=plt.loglog, label=name)
plt.legend()
plt.show()
