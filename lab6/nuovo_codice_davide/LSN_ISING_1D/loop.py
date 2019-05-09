import os
import string
import subprocess
import numpy as np

T = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]

os.system('rm -r T_*')
os.system('rm -r *.0')
os.system('rm -r *.temp')


for i in T:
	path='T_{}'.format(i)
	if not os.path.exists(path):
        	os.mkdir(path)
        	os.system('cp input.dat '+path)
        	os.system('cp Monte_Carlo_ISING_1D.exe '+path)
	
	s = open(path+"/input.dat").read()
	s = s.replace('$temp', format(i))
	f = open(path+"/input.dat", 'w')
	f.write(s)
	f.close()

	os.chdir(path)
	os.system('./Monte_Carlo_ISING_1D.exe')
	os.system('cat output.ene.final >> ../energy.temp')
	os.system('cat output.heat.final >> ../heat.temp')
	os.system('cat output.mag.final >> ../magn.temp')
	os.system('cat output.chi.final >> ../chi.temp')
	os.chdir('..')

