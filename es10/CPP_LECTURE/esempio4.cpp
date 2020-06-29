#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])  
{
	MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
        int rank = MPI::COMM_WORLD.Get_rank();
	
	int isend[2],irecv[2];
	for(int i=0;i<2;i++)
		isend[i]=rank+i+1;
	
	MPI_Reduce(&isend[0],&irecv[0],1,MPI_INTEGER,MPI_SUM,0,MPI::COMM_WORLD);
	MPI_Reduce(&isend[1],&irecv[1],1,MPI_INTEGER,MPI_PROD,0,MPI::COMM_WORLD);
		
	
	if(rank==0)
		cout<<"irecv: "<<irecv[0]<<" "<<irecv[1]<<endl;			
	
	MPI::Finalize();

        return 0;

}
