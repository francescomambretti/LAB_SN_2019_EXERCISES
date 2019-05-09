import os
import string
import subprocess
import numpy as np

tot_steps=100000
nblk = [10,20,25,40,50,100,200,250,400,500,1000,2000,2500,4000,5000]

ene_arr=np.loadtxt('gas/output.epot.0', usecols=(1))

count=0
ene_glob=0
ene2_glob=0

for i in nblk:
	ene_glob=0
	ene2_glob=0
	nstep_per_blk=tot_steps/i
	for count in np.arange(1,i+1):
		ene_blk=0
		#data blocking
		for j in np.arange(0,nstep_per_blk):
			index=int(j+(nstep_per_blk*(count-1)))
			ene_blk+=ene_arr[index]
		ene_glob+=ene_blk/nstep_per_blk
		ene2_glob+=(ene_blk/nstep_per_blk)*(ene_blk/nstep_per_blk)
		if (count==i):
			print(nstep_per_blk,np.sqrt(ene2_glob/count-(ene_glob/count)*(ene_glob/count))/np.sqrt(count-1))
