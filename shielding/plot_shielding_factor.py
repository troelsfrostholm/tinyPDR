# Plot of shielding factor of H2 dissociation rate from
# Draine and Bertoldi 1996

from numpy import *
from matplotlib.pyplot import *

omega_H = 0.035
alpha = 2.0
N_ref = 5e14 # Reference number density of H2 (cm^-3)
b = 3e5 # Doppler broadening parameter (cm s^-1)
b5 = b/1e5
N_H2_max = 1e24

N_H2 = 10**(arange(10,log10(N_H2_max)))
x = N_H2/N_ref

S = (1 - omega_H)/(1+x/b5)**alpha + omega_H/(1+x)**0.5*exp(-8.5e-4*(1+x)**0.5)
loglog(N_H2, S, label="Draine and Bertoldi")
loglog(N_H2[6:12], N_H2[6:12]**(-0.5)*1e6, linestyle="--", label=r"$N_{H2}^{-0.5}$")

# Also plot text file dump of S_H2 for Richtings et. al. 2014 paper II
data_R14II = loadtxt("richtings_S_H2.txt")
loglog(data_R14II[8:,0], data_R14II[8:,2], label="Richtings et. al. ")

title("Shielding fractor for H2 dissociation rate")
xlabel(r"$N_{H_2} (cm^{-2})$")
ylabel(r"$S_{H_2}$")
xlim([1e13,3e21])
ylim([1e-6,1])
legend()
show()