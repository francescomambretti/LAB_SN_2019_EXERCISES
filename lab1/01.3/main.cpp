/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Modified by Francesco Mambretti, 08/03/2019

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

//Buffon's needle
double L = 1.0;
double d = 2.1;
int N_hit=0;
double x_c, x_t, y_t;
double dist2=0,dist=0;
double radius=L/2.;
ofstream out;
out.open("PI.dat");
double PI, PI2;

vector<int> N_throws = {10,100,1000,2500,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000};

for (int j=0; j<N_throws.size(); ++j){
	PI=0.;
	PI2=0.;
	for (int t=0; t < 10; t++){
		N_hit=0;
		for (int i=0; i < N_throws[j]; i++){
			x_c=rnd.Rannyu(0,d);
			//now I have to extract needle's direction
			//I generate random pairs of (x_t, y_t); these coordinates are accepted (and indicate a direction) as long as they are inside a circle centred in (x_c,0) whose radius is L/2
			do {
				x_t=rnd.Rannyu(-1,1);
				y_t=rnd.Rannyu(-1,1);	
				dist2=(x_t)*(x_t)+(y_t)*(y_t);
				dist=sqrt(dist2);
			} while( dist > 1); //stop when the extracted point falls inside the circle
	
			if (x_t>0){
				if (x_c+x_t*radius/dist>d || x_c-x_t*radius/dist<0)
					N_hit++;
			}
	
			else{
				if (x_c+x_t*radius/dist<0 || x_c-x_t*radius/dist>d)
					N_hit++;
			}	

		}//end i-for
		PI+=2*L/d*(double(N_throws[j])/double(N_hit));
		PI2+=(2*L/d*(double(N_throws[j])/double(N_hit)))*(2*L/d*(double(N_throws[j])/double(N_hit)));
		//out << N_throws[j] << '\t' << 2*L/d*(double(N_throws[j])/double(N_hit)) << endl;
		//cout << "With " <<  N_throws[j] << " throws, PI is equal to = " << 2*L/d*(double(N_throws[j])/double(N_hit)) << endl;
	}//end t-for
		out << N_throws[j] << ' ' << PI/10.0 << ' ' << sqrt(PI2/10.0-PI/10.0*PI/10.0) << endl;
}//end j-for
   	rnd.SaveSeed();

out;

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
