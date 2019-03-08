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

float ave_r=0, ave2_r=0, ave_v=0, ave2_v=0;
float tmp, counter=0;
ofstream out;
out.open("data.dat");

for (int j=1; j <= 10000; j++){
	ave_r=0;
	ave2_r=0;
	ave_v=0;	
	ave2_v=0;
	counter=0;
	for(int i=0; i<10*j; i++){
		tmp=rnd.Rannyu();
		ave_r+=tmp;
		ave2_r+=tmp*tmp;
		ave_v+=(tmp-0.5)*(tmp-0.5);
		ave2_v+=(tmp-0.5)*(tmp-0.5)*(tmp-0.5)*(tmp-0.5);
		counter++;
 //     cout << rnd.Rannyu() << endl;
   	}
	out << counter << ' ' << ave_r/counter << ' ' << sqrt (ave2_r/counter-(ave_r/counter)*(ave_r/counter)) /sqrt(counter)<< ' ' << ave_v/counter << ' ' << sqrt (ave2_v/counter-(ave_v/counter)*(ave_v/counter)) /sqrt(counter) << endl; //error=std dev of the average
}
	out.close();
//End points 1) and 2)

	// Chi-2 test
	float M = 100;
	float n = 10000;
	int times = 100;
	float chisquare = 0.;
	vector<float> histo(M);	

	float expected = n/M; 

	out.open("chisquare.dat");
	
	for (int i=0; i < times; i++){ //for a certain number of times
		histo.resize(M);
		for (int j=0; j < n; j++){ //generate n random numbers
			tmp=rnd.Rannyu();
			histo[floor(tmp*100)]++;
		}
		for (int k = 0; k < M; k++) chisquare+=(histo[k]-expected)*(histo[k]-expected)/(expected);
		out << i << ' ' << chisquare << endl;
		chisquare = 0.;
		histo.clear();
	}

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
