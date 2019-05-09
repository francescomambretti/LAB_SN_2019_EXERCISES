/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"
//#define RESTART //if it is defined, read velocities from config.0 file, otherwise ignore them and generate random velocities

using namespace std;

int main(){ 
	system("./clean.sh");
	Input();             //Inizialization
	int nconf = 1;
  	nstep_per_block=nstep/nblocks;
  	for(int j=0; j < nblocks; ++j){
		cout << j << endl;
    		Reset();
    		for(int istep=1; istep <= nstep_per_block; ++istep){
       			Move();           //Move particles with Verlet algorithm
       			if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
       			if(istep%10 == 0){
          			Measure();     //Properties measurement
		          	ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          			nconf += 1;
       			}
			Accumulate();
     		}
		PrintAverages(j+1);
  	}

	ConfFinal();         //Write final configuration to restart
	delete[] gofr;
	delete[] gofr_blk;
	delete[] gofr_tot;
	delete[] gofr_tot2;
	return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//measurement of g(r)
  nbins = 100;
  bin_size = (box/2.0)/(double)nbins;
  gofr = new double [nbins];
  gofr_blk = new double [nbins];
  gofr_tot = new double [nbins];
  gofr_tot2 = new double [nbins];
//Read initial configuration

  double sumv[3] = {0.0, 0.0, 0.0};

  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
#ifdef RESTART
    ReadConf >> x[i] >> y[i] >> z[i]  >> vx[i] >> vy[i] >> vz[i];
#else 
    ReadConf >> x[i] >> y[i] >> z[i];
#endif
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
#ifndef RESTART //If I am equilibrating
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl; 
   for (int i=0; i<npart; ++i){
     vx[i] = ((double) rand() / (RAND_MAX)) - 0.5;
     vy[i] = ((double) rand() / (RAND_MAX)) - 0.5;
     vz[i] = ((double) rand() / (RAND_MAX)) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
#else
   cout << "Read also velocities from file config.0 " << endl << endl;
   for (int i=0; i<npart; ++i){
	sumv[0] += vx[i];
	sumv[1] += vy[i];
	sumv[2] += vz[i];
   }
#endif

   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = x[i] - vx[i] * delta;
     yold[i] = y[i] - vy[i] * delta;
     zold[i] = z[i] - vz[i] * delta;
   }
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, w, t, vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press, Gofr;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Press.open("output_press.dat",ios::app);
  Gofr.open("output_gofr.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.0;

//reset the hystogram of g(r)
  for (int k=0; k<nbins; ++k) gofr[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //update of the histogram of g(r)
     bin = (int)(dr/bin_size);
     if(bin < nbins) gofr[bin] += 2.0;

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
//Potential energy
       v += vij;
//Virial
      w+=wij;
     }
    }          
  }

  stima_press = rho * (2.0 / 3.0) * t/(double)npart + 48.0/3.0*w/vol;

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    if (fabs(stima_temp-temp)>0.5){
    //rescale velocities
	for (int i=0; i<npart; ++i){
		vx[i]=vx[i]*sqrt(temp/stima_temp);
		vy[i]=vy[i]*sqrt(temp/stima_temp);
		vz[i]=vz[i]*sqrt(temp/stima_temp);
	}
    }

    stima_etot = (t+v)/(double)npart; //Total enery

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press << endl;

    for (int k=0; k < nbins; ++k)
	Gofr << k*bin_size << '\t' << gofr[k] << endl; 

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();
    Gofr.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration and velocities to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << "  " << vx[i] << "   " << vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void Reset () {
	epot_blk = 0;
	ekin_blk = 0;
	etot_blk = 0;
	temp_blk = 0;
	press_blk = 0;
	for (int k=0; k< nbins;++k)
		gofr_blk[k]=0;
}

void Accumulate() //Update block averages
{
	epot_blk+=stima_pot;
	ekin_blk+=stima_kin;
	etot_blk+=stima_etot;
	temp_blk+=stima_temp;
	press_blk+=stima_press;
	for (int k=0; k< nbins;++k)
		gofr_blk[k]+=gofr[k];
} //end Accumulate

void PrintAverages(int iblk){
	ofstream Epot, Ekin, Etot, Temp, Press, Gave;

	//update
	epot_blk=epot_blk/double(nstep_per_block);
	ekin_blk=ekin_blk/double(nstep_per_block);
	etot_blk=etot_blk/double(nstep_per_block);
	temp_blk=temp_blk/double(nstep_per_block);
	press_blk=press_blk/double(nstep_per_block);
	epot_tot+=epot_blk;
	ekin_tot+=ekin_blk;
	etot_tot+=etot_blk;
	temp_tot+=temp_blk;
	press_tot+=press_blk;
	epot2_tot+=epot_blk*epot_blk;
	ekin2_tot+=ekin_blk*ekin_blk;
	etot2_tot+=etot_blk*etot_blk;
	temp2_tot+=temp_blk*temp_blk;
	press2_tot+=press_blk*press_blk;

	Gave.open("g_ave.out",ios::app);

	double r, gdir;

    	for (int k=0; k< nbins; ++k)
    	{
        	r = k * bin_size;
        	gdir = gofr_blk[k]/double(nstep_per_block);
        	gdir *= 3.0/((4.0)*M_PI * (pow(r + bin_size,3) - pow(r,3)) * rho * (double)npart);
        	gofr_tot[k] += gdir;
       		gofr_tot2[k] += gdir*gdir;
        	if(iblk == nblocks)
            		Gave << setw(12) << scientific << r << "  " << gofr_tot[k]/(double)iblk << "  " << sqrt (gofr_tot2[k]/double(iblk)-(gofr_tot[k]/double(iblk))*(gofr_tot[k]/double(iblk)))/sqrt(double(iblk-1)) << endl;
    	}//end g(r)

	Epot.open("epot_ave.out",ios::app);
	Ekin.open("ekin_ave.out",ios::app);
	Etot.open("etot_ave.out",ios::app);
	Temp.open("temp_ave.out",ios::app);
	Press.open("press_ave.out",ios::app);	

	if (iblk>1){
		Epot << epot_tot/double(iblk) << '\t' << sqrt (epot2_tot/double(iblk)-(epot_tot/double(iblk))*(epot_tot/double(iblk)))/sqrt(double(iblk-1)) << endl;
		Ekin << ekin_tot/double(iblk) << '\t' << sqrt (ekin2_tot/double(iblk)-(ekin_tot/double(iblk))*(ekin_tot/double(iblk)))/sqrt(double(iblk-1)) << endl;
		Etot << etot_tot/double(iblk) << '\t' << sqrt (etot2_tot/double(iblk)-(etot_tot/double(iblk))*(etot_tot/double(iblk)))/sqrt(double(iblk-1)) << endl;
		Temp << temp_tot/double(iblk) << '\t' << sqrt (temp2_tot/double(iblk)-(temp_tot/double(iblk))*(temp_tot/double(iblk)))/sqrt(double(iblk-1)) << endl;
		Press << press_tot/double(iblk) << '\t' << sqrt (press2_tot/double(iblk)-(press_tot/double(iblk))*(press_tot/double(iblk)))/sqrt(double(iblk-1)) << endl;
	}
	else {
                Epot << epot_tot/double(iblk) << '\t' << 0.0 << endl;
                Ekin << ekin_tot/double(iblk) << '\t' << 0.0 << endl;
                Etot << etot_tot/double(iblk) << '\t' << 0.0 << endl;
                Temp << temp_tot/double(iblk) << '\t' << 0.0 << endl;
                Press << press_tot/double(iblk) << '\t' << 0.0 << endl;
	}

	Epot.close();
	Ekin.close();
	Etot.close();
	Temp.close();
	Press.close();
	Gave.close();
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
