#include "mpi.h"
#include <iostream>

using namespace std;
int main(int argc, char* argv[])  
{
	MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
        int rank = MPI::COMM_WORLD.Get_rank();

	double tstart = MPI_Wtime();

	int n = 50000;
	int* imesg = new int[n];
	int sum=0;

	for(int i=0;i<n;i++){
		imesg[i]=rank;
		if(rank==1) MPI::COMM_WORLD.Send(&imesg[0],n,MPI::INTEGER,0,i);
		else if(rank==0) MPI::COMM_WORLD.Recv(&imesg[0],n,MPI::INTEGER,1,i);
		sum += imesg[i];
	}

	double tend = MPI_Wtime();
	double dt = tend - tstart;
	cout<<"io sono "<<rank<<"; somma = "<<sum<<"; tempo = "<<dt<<endl;

	MPI::Finalize();

	return 0;
}
