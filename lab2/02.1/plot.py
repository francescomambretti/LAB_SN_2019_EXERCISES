import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#02.1
plt.suptitle('Integral evalution via Monte Carlo method')
x,f,error = np.loadtxt("I_uniform.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(x,f,error,color='C3',fmt='o',markersize='2',capsize=3,label="Uniform sampling")
x,f,error = np.loadtxt("I_importance.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(x,f,error,color='C2',fmt='o',markersize='2',capsize=3,label="Importance sampling")
plt.xlabel('Block number')
plt.ylabel('Integral value')
plt.legend() #show legend
plt.savefig('02.1.png')
plt.show()
