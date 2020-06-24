#include "random.h"
#include "random_walk.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(){

    ofstream mean;                                                  //creo i file di output
    ofstream sigma;
    ofstream sigmaL;
    ofstream meanL;
    meanL.open("Lattice_RW_mean.txt");
    sigmaL.open("Lattice_RW_sigma.txt");
    mean.open("Space_RW_mean.txt");
    sigma.open("Space_RW_sigma.txt");

	Random *rnd = new Random();										//creo l'oggetto random per i random walk
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd->SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    int stepRW = 1;                                                      //parto con RW da uno step
    int M = 100;                                                         //massimo di step RW
    int N = 10000;                                                       //per ogni numero di step, eseguo N RW su cui medio il valore rms finale

    Random_Walk *RW = new Random_Walk(stepRW);                           //creo un oggetto per eseguire random walk da stepRW step

    double y[N];                                                         //arrays per calcolare la media del valore rms della posizione finale dei RW
    double y2[N];
    double sum_prog[M];
    double sum2_prog[M];
    double devstd[M];

    for(int l = 0; l<M; l++){
        sum_prog[l]=0;
        sum2_prog[l]=0;
        devstd[l]=0;
	}
    //parte 1: svolgo l'esercizio in un lattice 3D di spaziatura 1
    for(int l = 0; l<M; l++){                                            //ciclo sul numero di step del RW da 1 a 100
    
        for(int i = 0; i<N; i++){                                        //ciclo per eseguire 10000 RW per ogni numero di step
            RW->RW_lattice(rnd, 1.);                                     //eseguo il RW e calcolo il valore rms
    		y[i] = sqrt(pow(RW->GetX(),2) + pow(RW->GetY(),2) + pow(RW->GetZ(),2));
	    	y2[i] = pow(y[i],2);
            RW->nullorigin();                                            //reimposto a (0,0,0) le coordinate dell'origine prima di ripetere il RW
            sum_prog[l] = sum_prog[l] + y[i];                            //calcolo il valor medio del valore rms
            sum2_prog[l] = sum2_prog[l] + y2[i];
        }

        sum_prog[l] = sum_prog[l]/N;
        sum2_prog[l] = sum2_prog[l]/N;
        devstd[l] = (sqrt(sum2_prog[l] - pow(sum_prog[l],2)));           //calcolo l'incertezza sulla media del valore rms della posizione finale
        meanL<<sum_prog[l]<<endl;                                        //salvo in output i dati
		sigmaL<<devstd[l]<<endl;
        stepRW++;                                                        //aumento di uno il numero di step
        RW->setsteps(stepRW);
        }

    for(int l = 0; l<M; l++){
        sum_prog[l]=0;
        sum2_prog[l]=0;
        devstd[l]=0;
	}
    //ripeto l'esercizio nello spazio invece che in un lattice; l'unica differenza è il metodo di RW utilizzato
    stepRW = 1;
    RW->setsteps(stepRW);

    for(int l = 0; l<M; l++){
    
        for(int i = 0; i<N; i++){
            RW->RW_continuum(rnd, 1.);
    		y[i] = sqrt(pow(RW->GetX(),2) + pow(RW->GetY(),2) + pow(RW->GetZ(),2));
	    	y2[i] = pow(y[i],2);
            sum_prog[l] = sum_prog[l] + y[i];
            sum2_prog[l] = sum2_prog[l] + y2[i];
            RW->nullorigin();
    	}
        sum_prog[l] = sum_prog[l]/N;
        sum2_prog[l] = sum2_prog[l]/N;
        devstd[l] = (sqrt(sum2_prog[l] - pow(sum_prog[l],2)));
        mean<<sum_prog[l]<<endl;
		sigma<<devstd[l]<<endl;
        stepRW++;
        RW->setsteps(stepRW);
    }

    meanL.close();
    sigmaL.close();
    mean.close();
    mean.close();

    delete rnd;
    delete RW;

return 0;
}