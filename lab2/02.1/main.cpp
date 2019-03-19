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

double eval (double);
 
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

//compute integral of a given function

int blk_size=100; //arbitrary, can be changed. Less than 100 would probably be too low
int tot_steps = 10000; // M
int nblk=tot_steps/blk_size; //number of blocks
double x; 
double I_tot=0., I2_tot=0., I_blk=0.;

ofstream out;
out.open ("I_uniform.dat");

for (int i=0; i < nblk; ++i){
	I_blk=0.;
	for (int t=0; t < blk_size; ++t){
		x=rnd.Rannyu(); //uniformly generate a random number between 0 and 1
		I_blk+=eval(x); //accumulate evaluation of the integrand function
	}//end t-loop

	I_tot+=I_blk/double(blk_size); //all the measures of I within a block contribute to a single block average
	I2_tot+=(I_blk/double(blk_size))*(I_blk/double(blk_size));

	if (i>0)
		out << i << ' ' << I_tot/double(i+1) << ' ' << sqrt( I2_tot/double(i+1) -I_tot*I_tot/(double((i+1)*(i+1)) ))/sqrt(double(i)) <<endl;
	else
		out << i << ' ' << I_tot/double(i+1) << ' ' << 0 << endl;
}//end i-loop

	out.close();

I_tot=0;
I2_tot=0;

out.open ("I_importance.dat");

for (int i=0; i < nblk; ++i){
	I_blk=0.;
	for (int t=0; t < blk_size; ++t){
		//extract x according to weight function
		//inverse of the cumulative function is used
		//this should improve the quality of the integral evaluation by quicker reducing the errors
		x=1-sqrt(1-rnd.Rannyu());	
		I_blk+=eval(x)/(2*(1-x)); // 2*(1-x) is the normalized probability density
	}//end t-loop

        I_tot+=I_blk/double(blk_size);
        I2_tot+=(I_blk/double(blk_size))*(I_blk/double(blk_size));

        if (i>0)
                out << i << ' ' << I_tot/double(i+1) << ' ' << sqrt( I2_tot/double(i+1) -I_tot*I_tot/(double((i+1)*(i+1))) )/sqrt(double(i)) <<endl;
        else
                out << i << ' ' << I_tot/double(i+1) << ' ' << 0 << endl;

}//end i-loop

	out.close();

   	rnd.SaveSeed();

return 0;
}

double eval (double x){
	return 0.5*M_PI*cos(M_PI*x/2.);
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
