#include "funzionebase.h"
#include "coseno_es2.h"
#include "random.h"
#include "integral.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(){
	ofstream mean;													//file di output
	mean.open("integral_mean_value.txt");
	ofstream sigma;
	sigma.open("integral_sigma_mean.txt");
	ofstream NUmean;
	NUmean.open("non_uniform_integral_mean.txt");
	ofstream NUsigma;
	NUsigma.open("non_uniform_integral_sigma_mean.txt");

	coseno_es2 *F = new coseno_es2(M_PI/2);							//costruttori per la funzione da valutare e per il calcolo del suo integrale
	integral *I = new integral(0., 1.);

	Random *rnd = new Random();										//creo l'oggetto random per produrre i campionamenti della variabile
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


	int hitpoint = 10000;													//numero di campionamenti della variabile per valutare l'integrale ad ogni blocco
	int N = 100;															//numero di blocchi -> il totale dei campionamenti è N*hitpoint
	double y[N];
	double y2[N];
	double sum_prog[N];
	double sum2_prog[N];
	double devstd[N];

	for(int i=0; i<N; i++){
		sum_prog[i]=0;
		sum2_prog[i]=0;
		y[i] = I->integration_MeanValue(F, rnd, hitpoint);					//l'array y è riempito con N stime dell'integrale, utilizzando ad ogni stima hitpoint punti
		y2[i] = pow(y[i],2);
	}

	for(int i=0; i<N; i++){													//calcola valor medio e incertezza statistica ad ogni blocco
		for(int j=0; j<i+1; j++){
			sum_prog[i]= sum_prog[i] + y[j];
			sum2_prog[i]= sum2_prog[i] + y2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		if(i==0){
			devstd[i] = 0;
		}
		else{
		devstd[i] = sqrt((sum2_prog[i] - pow((sum_prog[i]),2))/i);
		}
	}

	for(int i=0; i<N; i++){
		mean<<sum_prog[i]<<endl;
		sigma<<devstd[i]<<endl;
	}
																			//RIPETO L'ESERCIZIO PER UN CAMPIONAMENTO NON UNIFORME
	for(int i=0; i<N; i++){
		sum_prog[i]=0;
		sum2_prog[i]=0;
		y[i] = I->integration_importance(F, rnd, hitpoint);					//l'array y è riempito con N stime dell'integrale, utilizzando ad ogni stima hitpoint punti campionati con importance sampling
		y2[i]= pow(y[i],2);
	}

	for(int i=0; i<N; i++){													//calcola valor medio e incertezza statistica ad ogni blocco
		for(int j=0; j<i+1; j++){
			sum_prog[i]= sum_prog[i] + y[j];
			sum2_prog[i]= sum2_prog[i] + y2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		if(i==0){
			devstd[i] = 0;
		}
		else{
		devstd[i] = sqrt((sum2_prog[i] - pow((sum_prog[i]),2))/i);
		}
	}

	for(int i=0; i<N; i++){
		NUmean<<sum_prog[i]<<endl;
		NUsigma<<devstd[i]<<endl;
	}

	delete F;																//dealloco la memoria
	delete I;
	delete rnd;

	return 0;
}