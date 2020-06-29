#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])  
{
	MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
        int rank = MPI::COMM_WORLD.Get_rank();
	
	int itag=1;
	int imesg = rank;
	
	if(rank==1)
		MPI::COMM_WORLD.Send(&imesg,1,MPI::INTEGER,0,itag);
	else if(rank==0)
		MPI::COMM_WORLD.Recv(&imesg,1,MPI::INTEGER,1,itag);
				
	cout<<"messaggio = "<<imesg<<endl;
	
	MPI::Finalize();

        return 0;
}
