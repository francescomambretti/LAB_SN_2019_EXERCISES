import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#01.1.1
#plt.suptitle('Pseudo-Random Numbers generator test: <r>')
#
#x, f, error = np.loadtxt("data.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
#plt.errorbar(x,f,yerr=error,color='C0')
#plt.xlabel('$n_{throws}$')
#plt.ylabel('<r>')
#
#plt.savefig('01.1.1.png')
#
#plt.clf()
##01.1.2
#
#plt.suptitle('Pseudo-Random Numbers generator test: $<\sigma^2>$')
#x1, f1, error = np.loadtxt("data.dat", usecols=(0,3,4), delimiter=' ', unpack='true')
#plt.errorbar(x1,f1,yerr=error,color='C1')
#plt.xlabel('$n_{throws}$')
#plt.ylabel('$<\sigma^2>$')
#
#plt.savefig('01.1.2.png')
#
#plt.clf()

#01.1.3
plt.suptitle('Pearson $\chi^2$ test')
x,f = np.loadtxt("chisquare.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,f,color='C2')
plt.xlabel('j')
plt.ylabel('$\chi_j^2$')

#plt.savefig('01.1.3.png')
plt.show()
