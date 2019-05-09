#include "functions.h"

double J = 1.0;

//global variables
int N, steps, blocks, steps_equil, nsteps_per_block;
double T, h;
vector<int> lattice;
double ene_blk, ene2_blk, ene_glob=0, ene2_glob=0, cv_glob=0, cv2_glob=0, magn_blk, magn_glob=0, magn2_glob=0, chi_blk, chi_glob=0, chi2_glob=0;
vector<double> gofx_blk, gofx_glob, gofx2_glob;

//functions

//===============
// flip
//===============

void flip(){

	//select a random index
	int index, left, right;
	double deltaE;
	double gibbs;

	for (int t=0; t<N; ++t){ //try to flip all the spins
		index=rand() % N; //choose index between 0 and N-1
		//locate left and right neighbors of the selected lattice site, using PBC
		
		PBC(left,right,index);
		#ifndef GIBBS
		//final - initial energy is returned as deltaE
		deltaE = compute_deltaE (lattice[left],lattice[right],lattice[index]); //the initial state is passed as function argument

		if (rand()/double(RAND_MAX) <= exp(-deltaE/T) )	//accept the flip
		{
			lattice[index]=-lattice[index];
		}
		#else
		gibbs = 1./(1.+exp(2.*J/T*(lattice[left]+lattice[right])));
		if (rand()/double(RAND_MAX) <= gibbs)
			lattice[index]=+1;
		else
			lattice[index]=-1;
		#endif
	}//end t-for
	return;
}

//============
// PBC
//============

void PBC (int& left, int& right, int index){

	if (index==0) 
		left = N-1;
	else
		left = index-1;

	if (index==N-1)
		right = 0;
	else
		right = index+1;

	return;
}

//=================
//compute_energy
//=================

double compute_deltaE (int a, int b, int c){

	if(h==0) //no external field
		return 2.0*J*double(c)*(double(a)+double(b));
	else 
		return 2.0*double(c)*(h+J*(double(a)+double(b)));
}

//===================
// measure_energy
//===================

double measure_energy (){ //measure energy after a single overall flip
	double ene;

	ene = -J*lattice[0]*lattice[N-1]; //PBC
	if(h==1)
		ene += -h/2.*(lattice[N-1]+lattice[0]); //PBC 
	for (int i=0; i < N-1; ++i){
		ene+= -J*lattice[i]*lattice[i+1];
		if(h==1) 
			ene+=-h/2.*(lattice[i]+lattice[i+1]);
	}

	return ene;
}

//============
// reset
//============

void reset(){

	gofx_blk.clear();
	gofx_blk.resize(N/2);
	ene_blk=0;
	ene2_blk=0;
	chi_blk=0;
	if(h==1){
		magn_blk=0;
	}
	return;
}

//============
// accumulate
//============

void accumulate(){

	double appo=measure_energy()/double(N); //energy per spin
	ene_blk+=appo;
	ene2_blk+=appo*appo;

	appo=0;

	//susceptivity
	for (auto& n : lattice)
        	appo+=n;

	chi_blk+=appo*appo/(T*double(N)*double(N)); //if the susc.is not per spin, you must delete double(N)

	//measure spin/spin correlation function
	for (int i=0; i < N/2; ++i){
        	for (int x=0; x < N/2; ++x){
                	gofx_blk[x]+=(lattice[i]*lattice[i+x])*2./double(N);
        	}//end x-for
	}//end i-for


	if (h==1)//then measure M 
	{
		for (auto& n : lattice)
			appo+=n;
		magn_blk+=appo;
		magn_blk/=double(N);
	}

	return;
}

//=================
// averages
//=================

void averages(int iblk){

	double ene_ave_blk=ene_blk/double(nsteps_per_block);
	double ene2_ave_blk=ene2_blk/double(nsteps_per_block);
	double cv_ave_blk = N*(ene2_ave_blk-ene_ave_blk*ene_ave_blk)/T; //you must multiply by N because ene_blk is divided by N
	double chi_ave_blk = chi_blk/double(nsteps_per_block);

	ene_glob += ene_ave_blk;
	ene2_glob += ene_ave_blk*ene_ave_blk;

	cv_glob += cv_ave_blk;
	cv2_glob += cv_ave_blk*cv_ave_blk;
	//perform block averages
	ofstream write("ene.output",ios::app);
	write << T << '\t' << iblk << '\t' << ene_glob/double(iblk+1) << '\t' << error (ene_glob, ene2_glob, iblk) << endl;
	write.close();

	write.open("heat.output",ios::app);
	write << T<< '\t' <<iblk << '\t' << cv_glob/double(iblk+1) << '\t' << error (cv_glob, cv2_glob, iblk) << endl;
	write.close();	

	chi_glob += chi_ave_blk;
	chi2_glob += chi_ave_blk*chi_ave_blk;

	write.open("chi.output",ios::app);
	write << T<< '\t' <<iblk << '\t' << chi_glob/double(iblk+1) << '\t' << error (chi_glob, chi2_glob, iblk) << endl;
	write.close();	

	write.open("correl.output",ios::app);
	vector<double> gofx_ave_blk;
	gofx_ave_blk.resize(N/2);

	for (int x=0; x < N/2; ++x){// spin/spin correlation function
		gofx_ave_blk[x]=gofx_blk[x]/double(nsteps_per_block);
		gofx_glob[x]+=gofx_ave_blk[x];
		gofx2_glob[x]+=gofx_ave_blk[x]*gofx_ave_blk[x];
		write << T<< '\t' <<iblk << '\t' << x << '\t' << gofx_glob[x]/double(iblk+1) << '\t' << error (gofx_glob[x], gofx2_glob[x], iblk) << endl;
	}//end x-for
	write.close();
	
	if (h==1)
	{
		double magn_ave_blk = magn_blk/double(nsteps_per_block);
		magn_glob+= magn_ave_blk;
		magn2_glob+=magn_ave_blk*magn_ave_blk;
		
		write.open("magn.output",ios::app);
		write << T<< '\t' <<iblk << '\t' << magn_glob/double(iblk+1) << '\t' << error (magn_glob, magn2_glob, iblk) << endl;
		write.close();
	}

	return;	
}

//=============
// error
//=============

double error(double sum, double sum2, int iblk)
{
	if (iblk!=0)
        	return sqrt(sum2/((double)iblk+1.) - (sum/(double)(iblk+1.))*(sum/(double)(iblk+1.)))/((double)iblk);
	else
		return 0;
}//end Error

//===================
// print_config
//===================

void print_config(){
	ofstream write("spins.config");
	for (int i=0; i < N; ++i)
		write << lattice[i] << '\n';
	write.close();
	return;
}
