/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Modified by Francesco Mambretti, 19/03/2019

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

//3D random walk

int tot_times = 10000; //each random walk must be repeated tot_times times
int tot_steps = 100; //each RW lasts for tot_steps steps
int index=0,sign;
double lattice=1.0; //lattice step for next move
vector<double > r2_ave;
vector<double > r4_ave;
r2_ave.resize(tot_steps); //this is the observable we want to monitor (we will plot its sqrt in the end)
r4_ave.resize(tot_steps); //necessary for error evaluation
double r2=0;

ofstream out;
out.open ("RW3D.dat");

double pos [3] ={0,0,0}; //keep track of the actual global position

for (int i=0; i < tot_times; i++){
	pos[0]=0; //reset actual position for the next experiment
	pos[1]=0;
	pos[2]=0;
	r2_ave[0]=0;
	r4_ave[0]=0;
	for (int t=0; t < tot_steps; ++t){
		index = int(rnd.Rannyu (0,3)); //choose randomly the axis
		if (rnd.Rannyu() < 0.5) 
			sign = -1;
		else
			sign=1;	
		pos[index]+=lattice*sign; //move by a single lattice step in a random (discrete) direction
		r2=pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]; //square modulus of the actual position array
		r2_ave[t+1]+=r2; //accumulate
		r4_ave[t+1]+=r2*r2;
	}//end t-loop
}//end i-loop

//for the error, you have to propagate, because you are applying a sqrt to the measured observable values

double error=0;

for (int t=0; t < tot_steps; t++)
{
	error=sqrt (r4_ave[t]/double(tot_times)-(r2_ave[t]*r2_ave[t])/(double(tot_times)*double(tot_times)))/sqrt(double(tot_times-1)) ; //this is the error on r2
	out << t << ' ' << sqrt(r2_ave[t]/double(tot_times)) << ' ' << 0.5*error/sqrt(r2_ave[t]/double(tot_times)) << '\n';

}
	out.close();

// Continuum random walk
out.open ("RW3D_continuum.dat");

double theta, phi;
r2_ave.clear();
r4_ave.clear();
r2_ave.resize(tot_steps);
r4_ave.resize(tot_steps);

for (int i=0; i < tot_times; i++){
	pos[0]=0;
	pos[1]=0;
	pos[2]=0;
	r2_ave[0]=0;
	r4_ave[0]=0;
	for (int t=0; t < tot_steps; ++t){
		phi=rnd.Rannyu(0,2.*M_PI); //Uniformly between 0 and 2*PI
		theta=acos(1.-2.*rnd.Rannyu()); //Mind the Jacobian!
		pos[0]+=lattice*sin(theta)*cos(phi);
		pos[1]+=lattice*sin(theta)*sin(phi);
		pos[2]+=lattice*cos(theta);
		r2=pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]; //square modulus of the actual position array
		r2_ave[t+1]+=r2; //accumulate
		r4_ave[t+1]+=r2*r2;

	}//end t-loop
}//end i-loop

for (int t=0; t < tot_steps; t++)
{
	error=sqrt (r4_ave[t]/double(tot_times)-(r2_ave[t]*r2_ave[t])/(double(tot_times)*double(tot_times)))/sqrt(double(tot_times-1)) ; //this is the error on r2
	out << t << ' ' << sqrt(r2_ave[t]/double(tot_times)) << ' ' << 0.5*error/sqrt(r2_ave[t]/double(tot_times)) << '\n';

}out.close();

 	
rnd.SaveSeed();

return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
