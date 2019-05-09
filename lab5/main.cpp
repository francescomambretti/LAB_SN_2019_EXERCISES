#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "random.h"

#define BOHR 0.0529 //nm

using namespace std;

int M=1E6;
double delta1, delta2;

//functions
void ReadInput();
double prob(double,double,double,bool);

int main(){

	ofstream out;

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

	int nblk = 1000;
	int nstep_per_blk=M/nblk;

	ReadInput();

	double r[3]={0.,0.,0.};
	double ave_r=0., blk_r=0., blk_r2=0.;
	double dx,dy,dz;
	int accept;

	out.open("gs.dat");

	for (int i=0; i < nblk; ++i){
		ave_r=0.;
		accept = 0;
		for (int j=0; j < nstep_per_blk; ++j){
			//misura
			ave_r+=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			//proponi update
			dx=rnd.Rannyu(-0.5,0.5)*delta1; 
			dy=rnd.Rannyu(-0.5,0.5)*delta1;
			dz=rnd.Rannyu(-0.5,0.5)*delta1;
			/*dx=rnd.Gauss(0,0.5*delta1);
			dy=rnd.Gauss(0,0.5*delta1);
			dz=rnd.Gauss(0,0.5*delta1);*/
			//accettalo eventualmente
			if (rnd.Rannyu() < min ( prob(r[0]+dx,r[1]+dy,r[2]+dz, true)/(prob(r[0],r[1],r[2], true)) , 1.) ){
				r[0]+=dx;
				r[1]+=dy;
				r[2]+=dz;
				accept++;
			}	
		}
		blk_r+=ave_r/double(nstep_per_blk);
		blk_r2+=(ave_r/double(nstep_per_blk))*(ave_r/double(nstep_per_blk));
		if (i==0)
			out << BOHR*blk_r/double(i+1) << '\t' << 0 << endl;
		else
			out << BOHR*blk_r/double(i+1) << '\t' << sqrt (BOHR*BOHR*blk_r2/double(i+1)-BOHR*BOHR*(blk_r/double(i+1))*(blk_r/double(i+1)))/sqrt(double(i)) << endl;
		cout << "Acceptance ratio: " << double(accept)/double(nstep_per_blk) << endl;
	}	
	cout << "The average value of r on the ground state is: " << BOHR*blk_r/double(nblk) << " +- " << sqrt (BOHR*BOHR*blk_r2/double(nblk)-BOHR*BOHR*(blk_r/double(nblk))*(blk_r/double(nblk)))/sqrt(double(nblk-1)) << endl;	
	cout << "In units of Bohr radius: " << blk_r/double(nblk) << " +- " << sqrt (blk_r2/double(nblk)-(blk_r/double(nblk))*(blk_r/double(nblk)))/sqrt(double(nblk-1)) << endl;

	out.close();
	out.open("210.dat");

	for (int i=0; i < 3; ++i) r[i]=0.5; //in 0,0,0 it is not defined!!!
	blk_r=0;
	blk_r2=0;

	for (int i=0; i < nblk; ++i){
		ave_r=0.;
		accept = 0;
		for (int j=0; j < nstep_per_blk; ++j){
			//misura
			ave_r+=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			//proponi update
			/*dx=rnd.Rannyu(-0.5,0.5)*delta2;
			dy=rnd.Rannyu(-0.5,0.5)*delta2;
			dz=rnd.Rannyu(-0.5,0.5)*delta2;*/
			dx=rnd.Gauss(0,0.5*delta2);
			dy=rnd.Gauss(0,0.5*delta2);
			dz=rnd.Gauss(0,0.5*delta2);

			//accettalo eventualmente
			if (rnd.Rannyu() < min ( prob(r[0]+dx,r[1]+dy,r[2]+dz, false)/(prob(r[0],r[1],r[2], false)) , 1.) ){
				r[0]+=dx;
				r[1]+=dy;
				r[2]+=dz;
				accept++;
			}	
		}
		blk_r+=ave_r/double(nstep_per_blk);
		blk_r2+=(ave_r/double(nstep_per_blk))*(ave_r/double(nstep_per_blk));
		cout << "Acceptance ratio: " << double(accept)/double(nstep_per_blk) << endl;
		if (i==0)
        		out << BOHR*blk_r/double(i) << '\t' << 0 << endl;
		else	
        		out << BOHR*blk_r/double(i) << '\t' << sqrt (BOHR*BOHR*blk_r2/double(i+1)-BOHR*BOHR*(blk_r/double(i+1))*(blk_r/double(i+1)))/sqrt(double(i)) << endl;
	}	
	cout << "The average value of r on the ground state is: " << BOHR*blk_r/double(nblk) << " +- " << sqrt (BOHR*BOHR*blk_r2/double(nblk)-BOHR*BOHR*(blk_r/double(nblk))*(blk_r/double(nblk)))/sqrt(double(nblk-1)) << endl;	
	 cout << "In units of Bohr radius: " << blk_r/double(nblk) << " +- " << sqrt (blk_r2/double(nblk)-(blk_r/double(nblk))*(blk_r/double(nblk)))/sqrt(double(nblk-1)) << endl;

	out.close();

	return 0;
}

//functions implementation

void ReadInput(){
	fstream in ("input.dat");
	in >> delta1;
	in >> delta2;

	in.close();
}

double prob (double a, double b, double c, bool gs){
	if (gs==true)
		return (exp(-2*sqrt(a*a+b*b+c*c)));
		//return (exp(-2*sqrt(a*a+b*b+c*c)))/M_PI;
	else
	{
		double r = sqrt(a*a+b*b+c*c);
		double theta = acos (c/(sqrt(a*a+b*b+c*c)));
		return cos (theta)*cos(theta)*exp(-r)*r*r;
		//return cos (theta)*cos(theta)*exp(-r)*r*r/(32.*M_PI);
	}
}
