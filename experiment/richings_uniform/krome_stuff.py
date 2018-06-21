import numpy as np
from matplotlib import pyplot as plt

idx_E=1
idx_Hk=2
idx_Ck=3
idx_Ok=4
idx_H=5
idx_HE=6
idx_H2=7
idx_C=8
idx_O=9
idx_OH=10
idx_CO=11
idx_CH=12
idx_CH2=13
idx_C2=14
idx_HCO=15
idx_H2O=16
idx_O2=17
idx_Hj=18
idx_HEj=19
idx_H2j=20
idx_Cj=21
idx_Oj=22
idx_HOCj=23
idx_HCOj=24
idx_H3j=25
idx_CHj=26
idx_CH2j=27
idx_COj=28
idx_CH3j=29
idx_OHj=30
idx_H2Oj=31
idx_H3Oj=32
idx_O2j=33
idx_HEjj=34
idx_CR=35
idx_g=36
idx_Tgas=37
idx_dummy=38

# idx_E=1
# idx_Hk=2
# idx_Ck=3
# idx_Ok=4
# idx_H=5
# idx_HE=6
# idx_H2=7
# idx_C=8
# idx_O=9
# idx_OH=10
# idx_CO=11
# idx_CH=12
# idx_CH2=13
# idx_C2=14
# idx_HCO=15
# idx_H2O=16
# idx_O2=17
# idx_N=18
# idx_Hj=19
# idx_HEj=20
# idx_H2j=21
# idx_Cj=22
# idx_Oj=23
# idx_HOCj=24
# idx_HCOj=25
# idx_H3j=26
# idx_CHj=27
# idx_CH2j=28
# idx_COj=29
# idx_CH3j=30
# idx_OHj=31
# idx_H2Oj=32
# idx_H3Oj=33
# idx_O2j=34
# idx_Nj=35
# idx_HEjj=36
# idx_CR=37
# idx_g=38
# idx_Tgas=39
# idx_dummy=40

pc = 3.08567758e18                 # parsec in cm
units_density=3.4667880e-21        # 1e3 atoms / cm^3
units_length=1.5900779138177769e19 # corresponding to an opacity of 10 and to ~ 5.15 parsec
m_H = 1.6726219e-24/units_density
m_Hj = m_H
m_H2 = 2*m_H
m_C = 12*m_H
m_Cj = 12*m_H
m_O = 16*m_H
m_CO = m_C + m_O
m_He = 4*m_H

ntot=1e2
n_H = ntot/3.
n_H2 = ntot/3.*2.
n_C = 1e-4*ntot
n_O = 3e-4*ntot
n_He = 1e-1*ntot

o = 11

m_tot = (n_H*m_H + n_H2*m_H2 + n_C*m_C + n_O*m_O + n_He*m_He)

x_H = n_H*m_H/m_tot
x_H2 = n_H2*m_H2/m_tot
x_C = n_C*m_C/m_tot
x_O = n_O*m_O/m_tot
x_He = n_He*m_He/m_tot

o = 3

def plot_series(filename, labelled, linestyle, marker, axis=1):
  A = np.loadtxt(filename)
  xaxis = A[:,axis]
  plt.semilogy(xaxis, A[:,idx_H+o], label="H" if labelled else "", linestyle=linestyle, marker=marker)
  plt.semilogy(xaxis, A[:,idx_Hj+o], label="H+" if labelled else "", linestyle=linestyle, marker=marker)
  plt.semilogy(xaxis, A[:,idx_H2+o], label="H2" if labelled else "", linestyle=linestyle, marker=marker)
  plt.semilogy(xaxis, A[:,idx_C+o], label="C" if labelled else "", linestyle=linestyle, marker=marker)
  plt.semilogy(xaxis, A[:,idx_Cj+o], label="C+" if labelled else "", linestyle=linestyle, marker=marker)
  plt.semilogy(xaxis, A[:,idx_CO+o], label="CO" if labelled else "", linestyle=linestyle, marker=marker)
  plt.semilogy(xaxis, A[:,idx_O+o], label="O" if labelled else "", linestyle=linestyle, marker=marker)

def plot_series_log(filename, labelled, linestyle, marker, axis=1):
  A = np.loadtxt(filename)
  xaxis = A[:,axis]
  n_Htot = A[:,idx_H+o] + 2*A[:,idx_H2+o] + A[:,idx_Hj+o] + A[:,idx_H2j+o] + A[:,idx_H3j+o] + A[:,idx_Hk+o]
  n_C_tot = A[:,idx_C+o] + A[:,idx_Cj+o] + A[:,idx_CO+o] + A[:,idx_CH3j+o] + A[:,idx_COj+o] + A[:,idx_CH2j+o] + A[:,idx_CH3j+o] + A[:,idx_CHj+o]
  n_O_tot = A[:,idx_O+o] + A[:,idx_Ok+o] + A[:,idx_OH+o] + A[:,idx_H2O+o] + A[:,idx_CO+o] + A[:,idx_COj+o] + 2.*A[:,idx_O2+o] + A[:,idx_HCO+o] + A[:,idx_Oj+o] + A[:,idx_HOCj+o] + A[:,idx_HCOj+o] + A[:,idx_OHj+o] + A[:,idx_H2Oj+o]
  plt.loglog(xaxis, A[:,idx_H+o]/n_Htot, label="H" if labelled else "", linestyle=linestyle, marker=marker)
  plt.loglog(xaxis, A[:,idx_Hj+o]/n_Htot, label="H+" if labelled else "", linestyle=linestyle, marker=marker)
  plt.loglog(xaxis, 2.*A[:,idx_H2+o]/n_Htot, label="H2" if labelled else "", linestyle=linestyle, marker=marker)
  plt.loglog(xaxis, A[:,idx_E+o]/n_Htot, label="E" if labelled else "", linestyle=linestyle, marker=marker)
  plt.loglog(xaxis, A[:,idx_C+o]/n_C_tot, label="C" if labelled else "", linestyle=linestyle, marker=marker)
  plt.loglog(xaxis, A[:,idx_Cj+o]/n_C_tot, label="C+" if labelled else "", linestyle=linestyle, marker=marker)
  plt.loglog(xaxis, A[:,idx_CO+o]/n_C_tot, label="CO" if labelled else "", linestyle=linestyle, marker=marker)
  plt.loglog(xaxis, A[:,idx_O+o]/n_O_tot, label="O" if labelled else "", linestyle=linestyle, marker=marker)

def plot_molecular_fraction(filename, axis=1, pltfun=plt.loglog, **kwdargs):
  A = np.loadtxt(filename)
  xaxis = A[:,axis]
  n_Htot = A[:,idx_H+o] + 2*A[:,idx_H2+o] + A[:,idx_Hj+o] + A[:,idx_H2j+o] + A[:,idx_H3j+o] + A[:,idx_Hk+o]
  H2_fraction = 2*A[:,idx_H2+o]/n_Htot
  n_C_tot = A[:,idx_C+o] + A[:,idx_Cj+o] + A[:,idx_CO+o] + A[:,idx_CH3j+o] + A[:,idx_COj+o] + A[:,idx_CH2j+o] + A[:,idx_CH3j+o] + A[:,idx_CHj+o]
  CO_fraction = A[:,idx_CO+o]/n_C_tot
  n_O_tot = A[:,idx_O+o] + A[:,idx_Ok+o] + A[:,idx_OH+o] + A[:,idx_H2O+o] + A[:,idx_CO+o] + A[:,idx_COj+o] + 2.*A[:,idx_O2+o] + A[:,idx_HCO+o] + A[:,idx_Oj+o] + A[:,idx_HOCj+o] + A[:,idx_HCOj+o] + A[:,idx_OHj+o] + A[:,idx_H2Oj+o]
  OH_fraction = A[:,idx_OH+o]/n_O_tot
  pltfun(xaxis, H2_fraction, label="H2", color="black", **kwdargs)
  pltfun(xaxis, CO_fraction, label="CO", color="blue", **kwdargs)
  pltfun(xaxis, OH_fraction, label="OH", color="red", **kwdargs)

def plot_ionization_degree(filename, axis=1, pltfun=plt.loglog, **kwdargs):
  A = np.loadtxt(filename)
  xaxis = A[:,axis]
  n_Htot = A[:,idx_H+o] + 2*A[:,idx_H2+o] + A[:,idx_Hj+o] + A[:,idx_H2j+o] + A[:,idx_H3j+o] + A[:,idx_Hk+o]
  n_C_tot = A[:,idx_C+o] + A[:,idx_Cj+o] + A[:,idx_CO+o] + A[:,idx_CH3j+o] + A[:,idx_COj+o] + A[:,idx_CH2j+o] + A[:,idx_CH3j+o] + A[:,idx_CHj+o]
  n_O_tot = A[:,idx_O+o] + A[:,idx_Ok+o] + A[:,idx_OH+o] + A[:,idx_H2O+o] + A[:,idx_CO+o] + A[:,idx_COj+o] + 2.*A[:,idx_O2+o] + A[:,idx_HCO+o] + A[:,idx_Oj+o] + A[:,idx_HOCj+o] + A[:,idx_HCOj+o] + A[:,idx_OHj+o] + A[:,idx_H2Oj+o]
  n_HE_tot = A[:,idx_HE+o] + A[:,idx_HEj+o] + A[:,idx_HEjj+o]
  H_fraction = A[:,idx_H+o]/n_Htot
  Hj_fraction = A[:,idx_Hj+o]/n_Htot
  Cj_fraction = A[:,idx_Cj+o]/n_C_tot
  HE_fraction = A[:,idx_HE+o]/n_HE_tot
  HEj_fraction = A[:,idx_HEj+o]/n_HE_tot
  pltfun(xaxis, A[:,idx_E+o]/n_Htot, label="E", color="black", **kwdargs)
  pltfun(xaxis, H_fraction, label="H", color="brown", **kwdargs)
  pltfun(xaxis, Hj_fraction, label="H+", color="red", **kwdargs)
  pltfun(xaxis, HE_fraction, label="HE", color="green", **kwdargs)
  pltfun(xaxis, HEj_fraction, label="HE+", color="purple", **kwdargs)
  pltfun(xaxis, Cj_fraction, label="C+", color="blue", **kwdargs)

def plot_series_convergence(filename, label, style):
  A = np.loadtxt(filename)
  xaxis = A[:,1]
  plt.semilogy(xaxis, A[:,idx_H+o],  style, label=label)
  plt.semilogy(xaxis, A[:,idx_Hj+o], style)
  plt.semilogy(xaxis, A[:,idx_H2+o], style)
  plt.semilogy(xaxis, A[:,idx_C+o],  style)
  plt.semilogy(xaxis, A[:,idx_Cj+o], style)
  plt.semilogy(xaxis, A[:,idx_CO+o], style)
  plt.semilogy(xaxis, A[:,idx_O+o],  style)

def plot_0D_time_series(filename, **kwdargs):
  data = np.loadtxt(filename)
  xaxis = data[:,0]
  oo=1
  n_Htot = data[:,idx_H+oo] + 2*data[:,idx_H2+oo] + data[:,idx_Hj+oo] + data[:,idx_H2j+oo] + data[:,idx_H3j+oo] + data[:,idx_Hk+oo]
  H2_fraction = 2*data[:,idx_H2+oo]/n_Htot
  n_C_tot = data[:,idx_C+oo] + data[:,idx_Cj+oo] + data[:,idx_CO+oo] + data[:,idx_CH3j+oo] + data[:,idx_COj+oo] + data[:,idx_CH2j+oo] + data[:,idx_CH3j+oo] + data[:,idx_CHj+oo]
  CO_fraction = data[:,idx_CO+oo]/n_C_tot
  plt.semilogy(xaxis, H2_fraction, label="H2", color="black", **kwdargs)
  plt.semilogy(xaxis, CO_fraction, label="CO", color="blue", **kwdargs)

def plot_ionization_degree_0D(filename, **kwdargs):
  A = np.loadtxt(filename)
  xaxis = A[:,0]
  oo=1
  n_Htot = A[:,idx_H+oo] + 2*A[:,idx_H2+oo] + A[:,idx_Hj+oo] + A[:,idx_H2j+oo] + A[:,idx_H3j+oo] + A[:,idx_Hk+oo]
  n_C_tot = A[:,idx_C+oo] + A[:,idx_Cj+oo] + A[:,idx_CO+oo] + A[:,idx_CH3j+oo] + A[:,idx_COj+oo] + A[:,idx_CH2j+oo] + A[:,idx_CH3j+oo] + A[:,idx_CHj+oo]
  n_oo_tot = A[:,idx_O+oo] + A[:,idx_Ok+oo] + A[:,idx_OH+oo] + A[:,idx_H2O+oo] + A[:,idx_CO+oo] + A[:,idx_COj+oo] + 2.*A[:,idx_O2+oo] + A[:,idx_HCO+oo] + A[:,idx_Oj+oo] + A[:,idx_HOCj+oo] + A[:,idx_HCOj+oo] + A[:,idx_OHj+oo] + A[:,idx_H2Oj+oo]
  n_HE_tot = A[:,idx_HE+oo] + A[:,idx_HEj+oo] + A[:,idx_HEjj+oo]
  Hj_fraction = A[:,idx_Hj+oo]/n_Htot
  Cj_fraction = A[:,idx_Cj+oo]/n_Htot
  HE_fraction = A[:,idx_HE+oo]/n_Htot
  HEj_fraction = A[:,idx_HEj+oo]/n_Htot
  HEjj_fraction = A[:,idx_HEjj+oo]/n_Htot
  plt.loglog(xaxis, A[:,idx_E+oo]/n_Htot, label="E", color="black", **kwdargs)
  plt.loglog(xaxis, Hj_fraction, label="H+", color="brown", **kwdargs)
  plt.loglog(xaxis, HE_fraction, label="HE", color="green", **kwdargs)
  plt.loglog(xaxis, HEj_fraction, label="HE+", color="purple", **kwdargs)
  plt.loglog(xaxis, HEjj_fraction, label="HE++", color="orange", **kwdargs)
  plt.loglog(xaxis, Cj_fraction, label="C+", color="blue", **kwdargs)
