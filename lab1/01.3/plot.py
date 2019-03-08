import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#01.3
plt.suptitle('Buffon\'s needle: PI estimation')
x,f,error = np.loadtxt("PI.dat", usecols=(0,1,2), delimiter=' ', unpack='true')

plt.xscale('log')
plt.errorbar(x,f,error,fmt='o',markersize='4',capsize=4,color='C2')
plt.xlabel('$N_throws$')
plt.ylabel('$PI \: estimation$')

plt.hlines(3.1415,0,5000000)

plt.savefig('01.3.png')
plt.show()
