import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams


#######################################
# output from main.range.f90
# DT equimolar
data1 = np.genfromtxt('range_DT.out')
x1 = data1[:,2] * 1.e4 # [mu-m] cm * 10^4 = micron
dedx1 = data1[:,5]     # [MeV/mu-m]
dedxi1 = data1[:,6]    # [MeV/mu-m]
dedxe1 = data1[:,7]    # [MeV/mu-m]
# DTH p=1
data2 = np.genfromtxt('range_Hp1.out')
x2 = data2[:,2] * 1.e4 # [mu-m] cm * 10^4 = micron
dedx2 = data2[:,5]     # [MeV/mu-m]
dedxi2 = data2[:,6]    # [MeV/mu-m]
dedxe2 = data2[:,7]    # [MeV/mu-m]
# DTH p=05
data3 = np.genfromtxt('range_Hp05.out')
x3 = data3[:,2] * 1.e4 # [mu-m] cm * 10^4 = micron
dedx3 = data3[:,5]     # [MeV/mu-m]
dedxi3 = data3[:,6]    # [MeV/mu-m]
dedxe3 = data3[:,7]    # [MeV/mu-m]
# DTH p=05
data4 = np.genfromtxt('range_Hp01.out')
x4 = data4[:,2] * 1.e4 # [mu-m] cm * 10^4 = micron
dedx4 = data4[:,5]     # [MeV/mu-m]
dedxi4 = data4[:,6]    # [MeV/mu-m]
dedxe4 = data4[:,7]    # [MeV/mu-m]

#
rmax =  35.0
ymax = 0.40
plt.plot(x1,dedx1, label=r'$dE^{DT}/dx$  $p=0$')
plt.plot(x1,dedxe1, label=r'$dE^{DT}_e/dx$')
plt.plot(x2,dedx2, label=r'$dE^{H}/dx$  $p=1$')
plt.plot(x3,dedx3, label=r'$dE^{H}/dx$  $p=0.5$')
plt.plot(x4,dedx4, label=r'$dE^{H}/dx$  $p=0.1$')
plt.xlim(-0.1,rmax)
plt.ylim(0,ymax)
plt.xlabel(r'$x \,\, {\rm \mu m}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Range: BPS')
plt.legend(loc=2)
plt.grid(True)
plt.savefig('plot_dedxR_DTH.pdf')
plt.show()
