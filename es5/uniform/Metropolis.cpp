#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "Autofunzione.h"
#include <cmath>

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;																//genero l'oggetto random per svolgere l'esercizio
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

    ofstream r100;
    ofstream sr100;
    ofstream r210;
    ofstream sr210;
    ofstream X100;
    ofstream Y100;
    ofstream Z100;
    ofstream X210;
    ofstream Y210;
    ofstream Z210;
    r100.open("mean_100.dat");
    sr100.open("sigma_100.dat");
    r210.open("mean_210.dat");
    sr210.open("sigma_210.dat");
    X100.open("X_100.dat");
    Y100.open("Y_100.dat");
    Z100.open("Z_100.dat");
    X210.open("X_210.dat");
    Y210.open("Y_210.dat");
    Z210.open("Z_210.dat");

    Autofunzione* Psi = new Autofunzione();

    double origin[3] = {1., -1., 0.5};                          //punto di partenza
    double move[3];                                             //array per le coordinate della mossa metropolis
    int step = 1000000;                                         //mosse Metropolis
    int accepted = 0;                                           //contatore per le mosse accettate
    int blocks = 500;                                           //numero di blocchi per data blocking
    int L = step/blocks;                                        //step per blocco
    double acceptance;                                          //salva il valore di accettanza
    double fstep;                                               //variabili per il calcolo del valore di accettanza
    double fmove;
    double AVE[blocks];                                         //array per data blocking
	double AV2[blocks];
	double sum_prog[blocks];
	double sum2_prog[blocks];
	double devstd[blocks];
    double sum=0;

    for(int i=0; i< blocks; i++){
        sum=0;
        for(int j=0; j<L; j++){
            move[0] = origin[0] + (rnd.Rannyu(-1.225 , 1.225));        //calcola la mossa x' con prob. di trans. uniforme su un intervallo di più o meno un'unità di raggio di Bohr
            move[1] = origin[1] + (rnd.Rannyu(-1.225 , 1.225));
            move[2] = origin[2] + (rnd.Rannyu(-1.225 , 1.225));
            fstep = Psi->n1l0m0(origin[0], origin[1], origin[2]);       //calcola il valore delle pdf nel punto di partenza e di arrivo dello step
            fmove = Psi->n1l0m0(move[0], move[1], move[2]);
            if(fmove/fstep <1){
                acceptance = fmove/fstep;                       //calcola l'accettanza supponendo T(x'|x) = T(x|x')
            }
            else{
                acceptance = 1.;
            }
            double reject = rnd.Rannyu();                       //check per accettazione/reiezione della mossa
            if(reject <= acceptance){
                origin[0] = move[0];
                origin[1] = move[1];
                origin[2] = move[2];
                accepted++;                                     //conto il numero di mosse accettate per portare a 50% l'accettanza
            }
            sum += sqrt(pow(origin[0],2) + pow(origin[1],2) + pow(origin[2],2));        //sommatoria per la media di r in ogni blocco
	    if(j%10 == 0){
            X100<<origin[0]<<endl;                              //salvo le coordinate
            Y100<<origin[1]<<endl;
            Z100<<origin[2]<<endl;
	    }
        }
        AVE[i] = sum/L;                                         //calcolo il valor medio di r in ogni blocco
        AV2[i]=pow(AVE[i],2);
    }
    
    X100.close();
    Y100.close();
    Z100.close();

    cout<<"la percentuale di accettazione dell'algoritmo di Metropolis, per l'autofunzione 1,0,0, è di: "<<(double)accepted/step<<endl<<endl;
    accepted=0;

    for(int i=0; i<blocks; i++){                                //data-blocking
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

    for(int i=0; i<blocks; i++){
        r100<<sum_prog[i]<<endl;
        sr100<<devstd[i]<<endl;
    }

    r100.close();
    sr100.close();

    for(int i=0; i<blocks; i++){
        AVE[i]=0;
        AV2[i]=0;
        sum_prog[i]=0;
        sum2_prog[i]=0;
        devstd[i]=0;
    }
                                        //l'esercizio si ripete, cambia solo la pdf valutata e la prob. di trans.

    origin[0] = 3.;
    origin[1] = -3.;
    origin[2] = 2.;

    for(int i=0; i< blocks; i++){
        sum=0;
        for(int j=0; j<L; j++){
            move[0] = origin[0] + (rnd.Rannyu(-3. , 3.));
            move[1] = origin[1] + (rnd.Rannyu(-3. , 3.));
            move[2] = origin[2] + (rnd.Rannyu(-3. , 3.));
            fstep = Psi->n2l1m0(origin[0], origin[1], origin[2]);                   //uso la pdf per l'orbitale 2p
            fmove = Psi->n2l1m0(move[0], move[1], move[2]);
            if(fmove/fstep <1){
                acceptance = fmove/fstep;
            }
            else{
                acceptance = 1.;
            }
            double reject = rnd.Rannyu();
            if(reject <= acceptance){
                origin[0] = move[0];
                origin[1] = move[1];
                origin[2] = move[2];
                accepted++;
            }
            sum += sqrt(pow(origin[0],2) + pow(origin[1],2) + pow(origin[2],2));
	    if(j%10 == 0){
            X210<<origin[0]<<endl;
            Y210<<origin[1]<<endl;
            Z210<<origin[2]<<endl;
	    }
        }
        AVE[i] = sum/L;
        AV2[i]=pow(AVE[i],2);
    }
    
    cout<<"la percentuale di accettazione dell'algoritmo di Metropolis, per l'autofunzione 2,1,0, è di: "<<(double)accepted/step<<endl<<endl;

    X210.close();
    Y210.close();
    Z210.close();

    for(int i=0; i<blocks; i++){
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

    for(int i=0; i<blocks; i++){
        r210<<sum_prog[i]<<endl;
        sr210<<devstd[i]<<endl;
    }

    r210.close();
    sr210.close();

    delete Psi;

    return 0;
}