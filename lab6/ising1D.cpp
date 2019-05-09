#include "functions.h"

int main (int argc, char** argv){

	system("rm *.output");

	//read parameters
	fstream Read("param.in");
	Read >> N;
	Read >> steps;
	Read >> blocks;
	Read >> steps_equil;
	Read >> T;
	Read >> h;

	if (h!=0 and h!=1)
	{
		cerr << "Error! h value not acceptable, please change it." << endl;
		return -1;
	}

	Read.close();

	gofx_blk.resize(N/2);
	gofx_glob.resize(N/2);
	gofx2_glob.resize(N/2);

	//setup initial configuration
	lattice.assign(N,1); //all spins are up, at the beginning. This could be improved

	nsteps_per_block = steps/blocks;

	cout << "There will be " << nsteps_per_block << " steps in each block" << endl;

	srand(0); //initialize RNG

	//equilibration
	for (int j=0; j < steps_equil; ++j)
		flip();
	//simulation
	for (int i=0; i < blocks; ++i){
		reset();
		for (int j=0; j < nsteps_per_block; ++j){
			flip();	
			accumulate();
		}//end j-for
		averages(i);
	}//end i-for

	print_config();

	return 0;
}
