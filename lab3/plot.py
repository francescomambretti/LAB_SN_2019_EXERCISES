import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#03.1.1
x=np.arange(100)
plt.suptitle('Call prices: direct and stepwise sampling')
f,error = np.loadtxt("call_final.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.errorbar(x,f,error,color='C3',fmt='o',markersize='2',capsize=3,label="Call; direct final sampling")
f,error = np.loadtxt("call_discrete.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.errorbar(x,f,error,color='C4',fmt='o',markersize='2',capsize=3,label="Call; stepwise discrete sampling")
plt.xlabel('Block number')
plt.ylabel('Price')
plt.legend() #show legend
plt.savefig('call.png')
plt.show()

plt.clf()

plt.suptitle('Put prices: direct and stepwise sampling')
f,error = np.loadtxt("put_final.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.errorbar(x,f,error,color='C8',fmt='o',markersize='2',capsize=3,label="Put; direct final sampling")
f,error = np.loadtxt("put_discrete.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.errorbar(x,f,error,color='C1',fmt='o',markersize='2',capsize=3,label="Put; stepwise discrete sampling")
plt.xlabel('Block number')
plt.ylabel('Price')
plt.legend() #show legend
plt.savefig('call.png')
plt.show()


