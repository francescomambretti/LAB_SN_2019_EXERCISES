import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def f_fit(x,k):  # defining the fitting function
    return k * np.sqrt(x)

k=1 #trial value

#02.2
plt.suptitle('Random walk in 3D')
x,f,error = np.loadtxt("RW3D.dat", usecols=(0,1,2), delimiter=' ', unpack='true')

p_opt, p_cov = curve_fit(f_fit, x, f)
y_fit = f_fit(x,p_opt[0])
plt.errorbar(x,y_fit,error,color='C1',fmt='o',markersize='2',capsize=3,label="Discrete lattice")

x,f,error = np.loadtxt("RW3D_continuum.dat", usecols=(0,1,2), delimiter=' ', unpack='true')

p_opt, p_cov = curve_fit(f_fit, x, f)
y_fit = f_fit(x,p_opt[0])
print (y_fit/x)
plt.errorbar(x,f,error,color='C2',fmt='o',markersize='2',capsize=3,label="Continuous")

plt.xlabel('Step')
plt.ylabel('$\sqrt{<r^2>}$')
plt.legend() #show legend
plt.savefig('02.2.png')
plt.show()
