#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

//#define GIBBS

using namespace std;

//global variables

extern double J;
extern int N, steps, blocks, steps_equil, nsteps_per_block;
extern double T, h;
extern vector<int> lattice;
extern double ene_blk, ene2_blk, ene_glob, ene2_glob, cv_glob, cv2_glob, magn_blk, magn_glob, magn2_glob, chi_blk, chi_glob, chi2_glob;
extern vector<double> gofx_blk, gofx_glob, gofx2_glob;

//functions
void flip();
void PBC (int&, int&, int);
double compute_deltaE (int, int, int);
double measure_energy();
void reset();
void accumulate();
void averages(int);
double error(double,double,int);
void print_config();
