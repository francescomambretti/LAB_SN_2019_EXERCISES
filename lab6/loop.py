import os
import string
import subprocess

T = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

os.system('rm -r T_*')
os.system('rm -r *dat')

for i in T:
	path='T_{}'.format(i)
	if not os.path.exists(path):
        	os.mkdir(path)
        	os.system('cp param.in '+path)
        	os.system('cp ising1D.x '+path)
	
	s = open(path+"/param.in").read()
	s = s.replace('$temp', format(i))
	f = open(path+"/param.in", 'w')
	f.write(s)
	f.close()

	os.chdir(path)
	os.system('./ising1D.x')
	os.chdir('..')

os.system('list=(T*); for i in ${list[*]}; do tail -n 1 $i/ene.output >> energy.dat; done')
os.system('list=(T*); for i in ${list[*]}; do tail -n 1 $i/heat.output >> heat.dat; done')
os.system('list=(T*); for i in ${list[*]}; do tail -n 1 $i/chi.output >> chi.dat; done')
os.system('list=(T*); for i in ${list[*]}; do tail -n 1 $i/magn.output >> magn.dat; done')
os.system('list=(T*); for i in ${list[*]}; do tail -n 50 $i/correl.output >> correl.dat; done')
