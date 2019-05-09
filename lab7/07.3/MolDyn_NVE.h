/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
double* gofr, *gofr_blk, *gofr_tot, *gofr_tot2;
int nbins;
double bin_size;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_press;
int nblocks=10, nstep_per_block;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//averages
double epot_blk, ekin_blk, etot_blk, temp_blk, press_blk;
double epot_tot=0, ekin_tot=0, etot_tot=0, temp_tot=0, press_tot=0, epot2_tot=0, ekin2_tot=0, etot2_tot=0, temp2_tot=0, press2_tot=0;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Reset();
void Accumulate();
void PrintAverages(int);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
