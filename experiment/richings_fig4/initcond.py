# This script computes the distance from column- and number density in 
# R14II fig 4, and interpolates n and T to a grid logarithmic in space. 
from numpy import *
from matplotlib import pyplot as plt

do_plot = True

datadir = "R14IIfig4/"
fav = 1e0/4e-22
pc = 3.08567758e18  # parsec in cm
n = loadtxt(datadir+"/n.csv")
T = loadtxt(datadir+"/T.csv")

# Find distance to first point with known density, 
# assuming constant density below this point
n0 = n[0,1]   # First density point
N0 = 3e20     # Column depth to this point, read off from R14II fig 4
x0 = N0/n0    # Estimated distance

# Find column depth and distance axis for the number density series
N = n[:,0]                          # Column depth
n = n[:,1]                          # total number density of H
dN = N[1:]-N[:-1]                   # Column depth elements
dx = empty(N.shape)                 
dx[0] = 0                           
dx[1:] = dN/(0.5*(n[0:-1]+n[1:]))   # Distance element is column depth element over average number density
x = cumsum(dx) + x0                 # Distance axis is the cumulative sum + distance to the first point

# Find distance for the temperature series by interpolating the axis above
xx = empty(len(x)+1)                # Create new distance axis that includes zero
xx[0] = 0e0
xx[1:] = x
NN = empty(len(N)+1)                # same for column depth axis
NN[0] = 0e0
NN[1:] = N
N_T = T[:,0]                        # column depth for temperature series
x_T = interp(N_T, NN, xx)           # interpolate to find distance axis for temperature series

# Interpolate to logarithmic grid
Ngrid = 256
x_u = logspace(log10(2e2),log10(6e3),num=Ngrid)*pc
n_u = interp(x_u, x, n, left=n0, right=n[-1])
T_u = interp(x_u, x_T, T[:,1], left=T[0,1], right=T[-1,1])

# Optionally plot the series
if do_plot:
  fig, ax1 = plt.subplots()
  ax2 = ax1.twinx()
  ax1.loglog(x_T/pc, T[:,1], linestyle="", marker=".", color="black")
  ax1.loglog(x_u/pc, T_u, linestyle="--", color="black")
  ax2.loglog(x/pc, n, linestyle="", marker=".", color="blue")
  ax2.loglog(x_u/pc, n_u, linestyle="--", color="blue")
  ax1.set_ylim(1e1,1e4)
  ax2.set_ylim(1e-1,1e2)
  plt.show()

# Write series to file
M_u = empty((3,Ngrid))
M_u[0,:] = x_u
M_u[1,:] = n_u
M_u[2,:] = T_u
savetxt("initcond.dat",M_u)