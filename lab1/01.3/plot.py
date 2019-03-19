import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#01.3
plt.suptitle('Buffon\'s needle: PI estimation')
plt.title('L=1, d=1.2, blk_size=500, tot_steps=1000')
x,f,error = np.loadtxt("PI.dat", usecols=(0,1,2), delimiter=' ', unpack='true')

plt.errorbar(x,f,error,fmt='o',markersize='2',capsize=3,color='C2',ecolor='C2',elinewidth=1,errorevery=2)
plt.xlabel('$N_{blocks}$')
plt.ylabel('PI estimation')

plt.hlines(3.14159,0,200)

plt.savefig('01.3.png')
plt.show()
