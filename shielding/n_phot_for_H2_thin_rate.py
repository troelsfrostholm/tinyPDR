# This script calculates the number density of photons in the
# Black ISRF in the energy range relevant for H2 dissociation
# as used by Richings 2014 to parameterize the H2 dissociation
# rate
from numpy import *

h_planck = 4.13566766e-15 # Planck constant (eV s)
cs = 29979245800.0 # Speed of light (cm s^-1)

table = loadtxt("black87_eV.dat")
low = 12.24 # eV
high = 13.51 # eV

i_l = 39
i_c = 40
i_h = 41

l = table[i_l,0]
c = table[i_c,0]
h = table[i_h,0]

f_l = table[i_l,1]
f_c = table[i_c,1]
f_h = table[i_h,1]

print "e_low, e_high", low, high
print "ditto, from file", l, h

f_low = f_l + (f_c - f_l)/(c - l)*(low - l)
A_l = 0.5*(f_low + f_c)*(c - low)
f_high = f_c + (f_h - f_c)/(h - c)*(c - high)
A_h = 0.5*(f_c + f_high)*(high - c)

A_tot = A_l + A_h
print "Square approximation of A", f_c*(high-low)
print "Triangle approximation of A using 3 points", A_tot

n = 4*pi*A_tot / (cs * h_planck * c)
print "n", n
