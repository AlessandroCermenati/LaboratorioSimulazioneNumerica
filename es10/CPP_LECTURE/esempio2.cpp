#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])  
{
	MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
        int rank = MPI::COMM_WORLD.Get_rank();
	
	int my_values[3];

	for(int i=0;i<3;i++)
		if(rank==0)
			my_values[i]=i+1;
		else
			my_values[i]=0;
	
	cout<<"Prima: "<<my_values[0]<<" "<<my_values[1]<<" "<<my_values[2]<<" per il processo "<<rank<<endl;
	
	MPI_Bcast ( my_values, 3, MPI_INTEGER, 0, MPI::COMM_WORLD);
	
	cout<<"Dopo: "<<my_values[0]<<" "<<my_values[1]<<" "<<my_values[2]<<" per il processo "<<rank<<endl;
	
	
	MPI::Finalize();

        return 0;

}
