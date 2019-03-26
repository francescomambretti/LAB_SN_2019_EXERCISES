/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Modified by Francesco Mambretti, 25/03/2019

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
//#include <algorithm>
#include "random.h"

using namespace std;

double M=100000.;
int nblk=100;
double S0=100.;
double T=1.;
double K=100.;
double r=0.1;
double sigma=0.25;
int nstep=100; 

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

//03.1.1
	ofstream out,out1;
	double final_C=0, w, final_C2=0,appo,S, final_P=0, final_P2=0, appo1, blk_C, blk_P;

	int nstep_per_blk=M/double(nblk);

	out.open("call_final.dat");
	out1.open("put_final.dat");
	for (int i=0; i < nblk; ++i){
		blk_C=0;
		blk_P=0;
		for (int j=0; j < nstep_per_blk; ++j){
			w= rnd.Gauss(0,sigma);
			S=S0*exp((r-sigma*sigma/2.0)*T+w);//call
			appo=exp(-r*T)*max(0.,S-K);
			w= rnd.Gauss(0,sigma);
			S=S0*exp((r-sigma*sigma/2.0)*T+w);//put
			appo1=exp(-r*T)*max(0.,-S+K);
			blk_C+=appo;
			blk_P+=appo1;
		}//end j-for
		blk_C/=double(nstep_per_blk); //questa è la mia osservabile di blocco
		final_C+=blk_C;
		final_C2+=blk_C*blk_C;
		blk_P/=double(nstep_per_blk); //questa è la mia osservabile di blocco
		final_P+=blk_P;
		final_P2+=blk_P*blk_P;
		if (i==0){
			out << final_C/double(i+1) << ' ' << 0 << endl;
			out1 << final_P/double(i+1) << ' ' << 0 << endl;
		}
		else{
			out << final_C/double(i+1) << ' '  << (final_C2/double(i+1)-(final_C/double(i+1))*(final_C/double(i+1)))/sqrt(double(i)) << endl;
			out1 << final_P/double(i+1) << ' '  <<(final_P2/double(i+1)-(final_P/double(i+1))*(final_P/double(i+1)))/sqrt(double(i)) << endl;
		}
	}//end i-for

	cout << "The direct final call price is " << final_C/double(nblk) << " " << sqrt(final_C2/double(nblk)-final_C*final_C/double(nblk*nblk))/sqrt(double(nblk-1)) << endl;
	cout << "The direct final put price is " << final_P/double(nblk) << " " << sqrt(final_P2/double(nblk)-final_P*final_P/double(nblk*nblk))/sqrt(double(nblk-1)) << endl;

	out.close();
	out1.close();
//03.1.2

	cout << "Using discrete step approximation..." << endl;

	double dt=T/double(nstep);
	cout << "...with dt=" <<dt << endl;
	final_C=0;
	final_C2=0;
	final_P=0;
	final_P2=0;
	double Sc,Sp;

	out.open("call_discrete.dat");
	out1.open("put_discrete.dat");
	for (int i=0; i < nblk; ++i){
		blk_C=0;
		blk_P=0;
		for (int j=0; j < nstep_per_blk; ++j){
			Sc=S0;
			Sp=S0;
			//generate 100 discrete steps
			for (int k=0; k < nstep; ++k){ //from 0 to T
				w= rnd.Gauss(0,sigma);
				Sc=Sc*exp((r-sigma*sigma/2.0)*dt+w*sqrt(dt));//call
				if (k==nstep-1) appo=exp(-r*T)*max(0.,Sc-K);
				w= rnd.Gauss(0,sigma);
				Sp=Sp*exp((r-sigma*sigma/2.0)*dt+w*sqrt(dt));//put
				if (k==nstep-1) appo1=exp(-r*T)*max(0.,-Sp+K);
			}//end k-for
			blk_C+=appo;
			blk_P+=appo1;		
		}//end j-for
		blk_C/=double(nstep_per_blk); //questa è la mia osservabile di blocco
		final_C+=blk_C;
		final_C2+=blk_C*blk_C;
		blk_P/=double(nstep_per_blk); //questa è la mia osservabile di blocco
		final_P+=blk_P;
		final_P2+=blk_P*blk_P;
		if (i==0){
			out << final_C/double(i+1) << ' ' << 0 << endl;
			out1 << final_P/double(i+1) << ' ' << 0 << endl;
		}
		else{
			out << final_C/double(i+1) << ' '  << (final_C2/double(i+1)-(final_C/double(i+1))*(final_C/double(i+1)))/sqrt(double(i)) << endl;
			out1 << final_P/double(i+1) << ' '  <<(final_P2/double(i+1)-(final_P/double(i+1))*(final_P/double(i+1)))/sqrt(double(i)) << endl;
		}
	}//end i-for

	cout << "The discrete step final call price is " << final_C/double(nblk) << " " << sqrt(final_C2/double(nblk)-final_C*final_C/double(nblk*nblk))/sqrt(double(nblk-1)) << endl;
	cout << "The discrete step final put price is " << final_P/double(nblk) << " " << sqrt(final_P2/double(nblk)-final_P*final_P/double(nblk*nblk))/sqrt(double(nblk-1)) << endl;

	out.close();
	out1.close();

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
