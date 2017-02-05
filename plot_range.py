import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

# Read in data from an ASCII data table
data = np.genfromtxt('dedx.001.out')
E = data[:,1] / 1000. # [keV]
dedx = data[:,2]      # [MeV/mu-m]
dedxi = data[:,3]      # [MeV/mu-m]
dedxe = data[:,4]      # [MeV/mu-m]

print "***", dedx

# Create a loglog plot of data
plt.plot(E,dedx)
plt.plot(E,dedxi)
plt.plot(E,dedxe)
plt.xlim(-0.1,3.54)
plt.ylim(0,0.3)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Stopping Power: BPS')
plt.legend(loc=2)
plt.grid(True)
#plt.savefig('plot_range.pdf')
plt.show()


