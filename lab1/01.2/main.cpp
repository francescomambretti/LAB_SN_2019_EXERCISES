/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Modified by Francesco Mambretti, 03/03/2019

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

ofstream out;
out.open("uniform.dat");

int N [4]={1,2,10,100};
double sum=0;

for (int i = 0; i < 4; i++){
	for (int j=1; j <= 10000; ++j){
		for (int k=0; k < N[i]; ++k){
			sum+=rnd.Rannyu();
		}//end for over k
		out << N[i] << ' ' << j << ' ' << sum/double(N[i]) << endl;
		sum=0;
	}//end for over j
}//end for over i

out.close();

out.open("expo.dat");

sum=0;

for (int i = 0; i < 4; i++){
	for (int j=1; j <= 10000; ++j){
		for (int k=0; k < N[i]; ++k){
			sum+=rnd.Expo(1.);
		}//end for over k
		out << N[i] << ' ' << j << ' ' << sum/double(N[i]) << endl;
		sum=0;
	}//end for over j
}//end for over i

out.close();

out.open("lorentz.dat");

sum=0;

for (int i = 0; i < 4; i++){
	for (int j=1; j <= 10000; ++j){
		for (int k=0; k < N[i]; ++k){
			sum+=rnd.Cauchy_Lorentz(0.,1.);
		}//end for over k
		out << N[i] << ' ' << j << ' ' << sum/double(N[i]) << endl;
		sum=0;
	}//end for over j
}//end for over i

out.close();

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
