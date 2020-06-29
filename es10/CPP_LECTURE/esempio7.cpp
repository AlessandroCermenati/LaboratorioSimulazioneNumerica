#include "mpi.h"
#include <iostream>

using namespace std;
const int n = 10000;
//const int n = 20000;
int main(int argc, char* argv[])  
{
	MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
        int rank = MPI::COMM_WORLD.Get_rank();
	
	int* imesg = new int[n];
	int* imesg2 = new int[n];
	int itag=1;
	int itag2=2;
	for(int i=0;i<n;i++)
	{
		imesg[i]=rank;
		imesg2[i]=rank+1;
	}
	
	if(rank==1){
		MPI::COMM_WORLD.Send(&imesg[0],n,MPI::INTEGER,0,itag);
		MPI::COMM_WORLD.Recv(&imesg2[0],n,MPI::INTEGER,0,itag2);
		cout<<"messaggio = "<<imesg2[0]<<endl;

	}
	else if(rank==0){
		MPI::COMM_WORLD.Send(&imesg2[0],n,MPI::INTEGER,1,itag2);
		MPI::COMM_WORLD.Recv(&imesg[0],n,MPI::INTEGER,1,itag);
//		MPI::COMM_WORLD.Send(&imesg2[0],n,MPI::INTEGER,1,itag2);
		cout<<"messaggio = "<<imesg[0]<<endl;
	}			
	
	MPI::Finalize();

        return 0;

}
