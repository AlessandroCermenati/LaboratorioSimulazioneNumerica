#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;																//genero l'oggetto random per svolgere il moto Browniano geometrico
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
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    ofstream mean_call;                                                  //creo i file di output
    ofstream sigma_call;
	ofstream mean_put;
    ofstream sigma_put;
    ofstream mean_call_discrete;
	ofstream sigma_call_discrete;
	ofstream mean_put_discrete;
	ofstream sigma_put_discrete;
    mean_call.open("mean_call.txt");
    sigma_call.open("sigma_call.txt");
    mean_put.open("mean_put.txt");
    sigma_put.open("sigma_put.txt");
	mean_call_discrete.open("mean_call_discrete.txt");
    sigma_call_discrete.open("sigma_call_discrete.txt");
    mean_put_discrete.open("mean_put_discrete.txt");
    sigma_put_discrete.open("sigma_put_discrete.txt");

    int M = 100000;                                     //totale GBM
    int N = 100;                                        //blocchi
    int L = M/N;                                        //GBM per blocco
    double C[M];										//array per contenere le M stime del prezzo dell'opzione call
	double P[M];										//array per contenere le M stime del prezzo dell'opzione put
	double S;											//variabile per contenere la stima dell'asset prize simulata
	double AVE[N];										//arrays per l'analisi dati
	double AV2[N];
	double sum_prog[N];
	double sum2_prog[N];
	double devstd[N];

	double drift = 0.1;									//parametri dell'esercizio
	double volatility = 0.25;
	double T = 1.;
	int iter = 100;
	double S0 = 100.;
	double r = 0.1;
	double K = 100;

    for(int i=0; i<N; i++){
		AVE[i]=0;
		AV2[i]=0;
		sum_prog[i]=0;
		sum2_prog[i]=0;
	}

	for(int i=0; i<M; i++){								//eseguo le stime di GBM e calcolo call e put
		S = rnd.GBM(S0, drift, volatility, T);			//stimo l'evoluzione dell'asset price
		if(S >= K ){									//calcolo il prezzo delle opzioni in funzione del profitto
			C[i] = exp(-1.*r*T)*(S - K);
			P[i] = 0;
		}
		else{
			C[i] = 0;
			P[i] = exp(-1.*r*T)*(K - S);
		}
	}

	double sum = 0;

	for(int i=0; i<N; i++){								//calcolo del valor medio e della sua incertezza statistica della call option
		sum = 0;
		for(int j=0; j<(L); j++){
			int k = j + i*(L);
			sum = sum + C[k];
		}
		AVE[i]=sum/L;
		AV2[i]=pow(AVE[i],2);
	}

	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i]= sum_prog[i] + AVE[j];
			sum2_prog[i]= sum2_prog[i] + AV2[j];
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

	for(int i=0; i<N; i++){									//scrittura dei risultati per valor medio e incertezza associata ad ogni blocco per la call option
		mean_call<<sum_prog[i]<<endl;
		sigma_call<<devstd[i]<<endl;
	}

	for(int i=0; i<N; i++){
		AVE[i]=0;
		AV2[i]=0;
		sum_prog[i]=0;
		sum2_prog[i]=0;
	}

	for(int i=0; i<N; i++){									//calcolo del valor medio e della sua incertezza statistica per la put option
		sum = 0;
		for(int j=0; j<(L); j++){
			int k = j + i*(L);
			sum = sum + P[k];
		}
		AVE[i]=sum/L;
		AV2[i]=pow(AVE[i],2);
	}

	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i]= sum_prog[i] + AVE[j];
			sum2_prog[i]= sum2_prog[i] + AV2[j];
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

	for(int i=0; i<N; i++){									//scrittura dei risultati per valor medio e incertezza associata ad ogni blocco per la put option
		mean_put<<sum_prog[i]<<endl;
		sigma_put<<devstd[i]<<endl;
	}

	mean_call.close();
	sigma_call.close();
	mean_put.close();
	sigma_put.close();

	for(int i=0; i<N; i++){										//ripeto l'esercizio per path discreto del GBM che simula l'evoluzione dell'asset price.
		AVE[i]=0;
		AV2[i]=0;
		sum_prog[i]=0;
		sum2_prog[i]=0;
	}

	for(int i=0; i<M; i++){
		S = rnd.GBM_step(S0, drift, volatility, T, iter);		//l'unica differenza rispetto alla prima parte è il metodo random invocato per la stima dell'asset price al tempo T
		if(S >= K ){
			C[i] = exp(-1.*r*T)*(S - K);
			P[i] = 0;
		}
		else{
			C[i] = 0;
			P[i] = exp(-1.*r*T)*(K - S);
		}
	}

	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<(L); j++){
			int k = j + i*(L);
			sum = sum + C[k];
		}
		AVE[i]=sum/L;
		AV2[i]=pow(AVE[i],2);
	}

	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i]= sum_prog[i] + AVE[j];
			sum2_prog[i]= sum2_prog[i] + AV2[j];
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
		mean_call_discrete<<sum_prog[i]<<endl;
		sigma_call_discrete<<devstd[i]<<endl;
	}

	for(int i=0; i<N; i++){
		AVE[i]=0;
		AV2[i]=0;
		sum_prog[i]=0;
		sum2_prog[i]=0;
	}

	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<(L); j++){
			int k = j + i*(L);
			sum = sum + P[k];
		}
		AVE[i]=sum/L;
		AV2[i]=pow(AVE[i],2);
	}

	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i]= sum_prog[i] + AVE[j];
			sum2_prog[i]= sum2_prog[i] + AV2[j];
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
		mean_put_discrete<<sum_prog[i]<<endl;
		sigma_put_discrete<<devstd[i]<<endl;
	}

	mean_call_discrete.close();
	sigma_call_discrete.close();
	mean_put_discrete.close();
	sigma_put_discrete.close();

return 0;
}