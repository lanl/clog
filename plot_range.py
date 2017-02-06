import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#######################################
# dedx
data = np.genfromtxt('dedx.001.out')
E = data[:,1] / 1000. # [MeV]
dedx = data[:,2]      # [MeV/mu-m]
dedxi = data[:,3]     # [MeV/mu-m]
dedxe = data[:,4]     # [MeV/mu-m]
#
Emax = 3.54 # [MeV]
ymax = 0.26
plt.plot(E,dedx, label=r'$dE/dx$')
plt.plot(E,dedxi, label=r'$dE_I/dx$')
plt.plot(E,dedxe, label=r'$dE_e/dx$')
plt.xlim(-0.1,Emax)
plt.ylim(0,ymax)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Stopping Power: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_dedx_001.pdf')
plt.show()

#######################################
# range
data = np.genfromtxt('range.001.out')
x = data[:,2] * 1.e4 # [mu-m] cm * 10^4 = micron
dedx = data[:,5]     # [MeV/mu-m]
dedxi = data[:,6]    # [MeV/mu-m]
dedxe = data[:,7]    # [MeV/mu-m]
#
rmax =  32.0
ymax = 0.26
plt.plot(x,dedx, label=r'$dE/dx$')
plt.plot(x,dedxi, label=r'$dE_I/dx$')
plt.plot(x,dedxe, label=r'$dE_e/dx$')
plt.xlim(-0.1,rmax)
plt.ylim(0,ymax)
plt.xlabel(r'$x \,\, {\rm \mu m}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Range: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_range_001.pdf')
plt.show()
