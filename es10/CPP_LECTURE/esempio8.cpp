#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])  
{
	MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
        int rank = MPI::COMM_WORLD.Get_rank();
	
	
	
	int itag=1;
	int itag2=1;
	
	int imesg = rank;
	int imesg2 = rank+1;
	
	
	if(rank==1){
		MPI::COMM_WORLD.Recv(&imesg2,1,MPI::INTEGER,0,itag2);
		MPI::COMM_WORLD.Send(&imesg,1,MPI::INTEGER,0,itag);
		cout<<"messaggio = "<<imesg2<<endl;
	}
	else if(rank==0){
		MPI::COMM_WORLD.Send(&imesg2,1,MPI::INTEGER,1,itag2);
		MPI::COMM_WORLD.Recv(&imesg,1,MPI::INTEGER,1,itag);
		// scrivere sempre i send-recv in modo speculare...altrimenti...
//		MPI::COMM_WORLD.Recv(&imesg,1,MPI::INTEGER,1,itag);
//		MPI::COMM_WORLD.Send(&imesg2,1,MPI::INTEGER,1,itag2);
		cout<<"messaggio = "<<imesg<<endl;
	}
	

	
	MPI::Finalize();

        return 0;

}
