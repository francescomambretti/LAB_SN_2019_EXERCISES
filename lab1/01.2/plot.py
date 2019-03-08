import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#here I need subplots
#capire come assegnare i vari blocchi di dati ai vari subplots

#01.2.1
fig, ax = plt.subplots(2, 2)
plt.suptitle('Uniform RNG')

plt.subplots_adjust(wspace=0.3,hspace=0.3)

#metterci anche titoletti o legende
N, x, f = np.loadtxt("uniform.dat", usecols=(0,1,2), delimiter=' ', unpack='true')

ax[0, 0].plot(x[:10000], f[:10000], 'r') #row=0, col=0
ax[0, 0].set_title('N=1')
ax[0, 0].set_ylim(0,1)
ax[1, 0].plot(x[10000:20000], f[10000:20000], 'b') #row=1, col=0
ax[1, 0].set_title('N=2')
ax[1, 0].set_ylim(0,1)
ax[0, 1].plot(x[20000:30000], f[20000:30000], 'g') #row=0, col=1
ax[0, 1].set_title('N=10')
ax[0, 1].set_ylim(0,1)
ax[1, 1].plot(x[30000:40000], f[30000:40000], 'k') #row=1, col=1
ax[1, 1].set_title('N=100')
ax[1, 1].set_ylim(0,1)

plt.savefig('01.2.1.png')

#01.2.2
fig, ax = plt.subplots(2, 2)
plt.suptitle('Exponential RNG')

plt.subplots_adjust(wspace=0.3,hspace=0.3)

#metterci anche titoletti o legende
N, x, f = np.loadtxt("expo.dat", usecols=(0,1,2), delimiter=' ', unpack='true')

ax[0, 0].plot(x[:10000], f[:10000], 'r') #row=0, col=0
ax[0, 0].set_title('N=1')
ax[0, 0].set_ylim(0,8)
ax[1, 0].plot(x[10000:20000], f[10000:20000], 'b') #row=1, col=0
ax[1, 0].set_title('N=2')
ax[1, 0].set_ylim(0,8)
ax[0, 1].plot(x[20000:30000], f[20000:30000], 'g') #row=0, col=1
ax[0, 1].set_title('N=10')
ax[0, 1].set_ylim(0,8)
ax[1, 1].plot(x[30000:40000], f[30000:40000], 'k') #row=1, col=1
ax[1, 1].set_title('N=100')
ax[1, 1].set_ylim(0,8)

plt.savefig('01.2.2.png')

#01.2.3
fig, ax = plt.subplots(2, 2)
plt.suptitle('Cauchy-Lorentz RNG')

plt.subplots_adjust(wspace=0.3,hspace=0.3)

#metterci anche titoletti o legende
N, x, f = np.loadtxt("lorentz.dat", usecols=(0,1,2), delimiter=' ', unpack='true')

ax[0, 0].plot(x[:10000], f[:10000], 'r') #row=0, col=0
ax[0, 0].set_title('N=1')
ax[0, 0].set_yscale('log')
ax[1, 0].plot(x[10000:20000], f[10000:20000], 'b') #row=1, col=0
ax[1, 0].set_title('N=2')
ax[1, 0].set_yscale('log')
ax[0, 1].plot(x[20000:30000], f[20000:30000], 'g') #row=0, col=1
ax[0, 1].set_title('N=10')
ax[0, 1].set_yscale('log')
ax[1, 1].plot(x[30000:40000], f[30000:40000], 'k') #row=1, col=1
ax[1, 1].set_title('N=100')
ax[1, 1].set_yscale('log')
plt.savefig('01.2.3.png')
