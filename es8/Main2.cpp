#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Main2.h"

using namespace std;

int main()
{

  Input(); //Inizialization
  Optimize();   //ottimizzazione dei parametri
  for(int iblk=1; iblk <= nblocks; ++iblk) //Misura dell'energia di ground state con l'autofunzione trial a parametri ottimizzati
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nsteps; ++istep)
    {
      Move();       //esegue la mossa Metropolis
      Measure();    //misura l'energia e fa l'update del data-blocking e dell'istogramma dell'autofunzione
    }
    Averages(iblk);   //Print results for current block
  }

  return 0;
}

void Optimize(void){

  ifstream parameter;
  parameter.open("parameter_optimized.dat");                //leggo da file i parametri ottimizzati se il file esiste e se l'ottimizzazione non è richiesta
  if(parameter.is_open() && optimize != 1){                  //per richiedere l'ottimizzazione mettere =1 l'ultimo dato in input.dat
      cout<<"leggo i parametri ottimizzati da file"<<endl;
      parameter >> mu;
      parameter >> sigma;
      parameter.close();
      cout<<"mu = "<<mu<<endl;
      cout<<"sigma = "<<sigma<<endl;
	}

  else{

    double Tstart = 1.0 + cool;
    int Nopt = Tstart/cool;

    cout << "OTTIMIZZAZIONE DEI PARAMETRI" << endl << endl;
    cout << "ricerca dei parametri nell'intervallo con l'algoritmo Simulated Annealing"<<endl;
    cout << "numero di raffreddamenti = "<< Nopt << endl;
    cout << "L'algoritmo parte da temperatura = "<< Tstart << endl;
    cout << "raffreddamento ad ogni step = "<< cool << endl;
    cout << "calcolo dell'energia effettuato con "<<nsteps<<" punti di sampling del modulo quadro della funzione d'onda"<<endl<<endl;

    double beta;
    double stima_old;
    double stima_new;
    double mu_old, sigma_old;
    double mu_new, sigma_new;
    double p;

    for(int i=0; i<Nopt; i++){
        beta = 1/Tstart;
        for(int j=0; j<nsteps; j++){
            Move();
            Measure();
		}
        stima_old = blk_av[0]/blk_norm;
        Reset(0);
        mu_old = mu;
        sigma_old = sigma;
        mu_new = mu + rnd.Rannyu(-deltamu, deltamu);
        sigma_new = sigma + rnd.Rannyu(-deltasig, deltasig);
        mu = mu_new;
        sigma = sigma_new;
        for(int j=0; j<nsteps; j++){
            Move();
            Measure();
		}
        stima_new = blk_av[0]/blk_norm;
        p = rnd.Rannyu();
        if(p>exp(-beta*(stima_new - stima_old))){
              mu = mu_old;
              sigma = sigma_old;
		}
        Tstart = Tstart - cool;
        if((i+1)%10 == 0){
            cout<<"ciclo di ottimizzazione numero "<<i+1<<"/"<<Nopt<<endl;
            cout<<"-------------------------------------"<<endl;
		}
	}

    cout<<"mu = "<<mu<<endl;
    cout<<"sigma = "<<sigma<<endl<<endl;

    ofstream par;
    par.open("parameter_optimized.dat");
    par << mu << endl << sigma << endl;
    par.close();
  }

return;
}

void Input(void){
    ifstream ReadInput, ParInput;
    ReadInput.open("input.dat");

    cout<<"Calcolo dell'energia di ground state con metodi variazionali"<<endl;
    cout<<"potenziale: V = x^4 - 5/2 x^2"<<endl;
    cout<<"Autofunzione di test psi(x) = exp(-(x-mu)^2/2s^2) + exp(-(x+mu)^2/2s^2)"<<endl;
    cout<<"massa e costante di Plank poste uguali a 1"<<endl;
    
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();

    //Informazioni di input
    ReadInput >> delta;
    cout<<"Mosse Metropolis scelte con probabilità uniforme in un intervallo [x - "<<delta<<" ; x + "<<delta<<"]"<<endl<<endl;

    ReadInput >> nblocks;
    cout<<"numero di blocchi = "<<nblocks<<endl;

    ReadInput >> nsteps;
    cout<<"numero di step in ogni blocco = "<<nsteps<<endl<<endl;

    ReadInput >> nbins;
    binl = 5.0/(double)nbins;
    cout<<"Il programma produce un istogramma per il campionamento del modulo quadro dell'autofunzione di test nell'intervallo [-2.5;2.5]"<<endl;
    cout<<"numero di bin dell'istogramma = "<<nbins<<endl;
    cout<<"dimensione di ogni bin = "<<binl<<endl<<endl;

    ReadInput >> cool;
    ReadInput >> deltamu;
    ReadInput >> deltasig;
    ReadInput >> optimize;

    ReadInput.close();
    //azzera la variabile per le misure
    hamiltonian = 0;
    //inizializza i parametri
    mu = 0.9;
    sigma = 0.6;

return;
}

void Measure(void){
    hamiltonian = Hamiltonian(X);                           //misuratore per il valor medio dell'energia
    blk_av[0] = blk_av[0] + hamiltonian;                    // accumulatore per il data blocking
    for(int i=1; i<=nbins; ++i){                            //accumulatore dell'istogramma
        if(X > hist_min + (double)(i-1)*binl && X <= hist_min + (double)i*binl){
            blk_av[i] = blk_av[i] + 1.0;
        }
    }
    blk_norm = blk_norm + 1.0;

return;
}

void Averages(int iblk){

    const int wd=12;

    ofstream EGS, histogram;
    EGS.open("ave_EGS.dat",ios::app);
    histogram.open("wavefunction.dat",ios::app);

    cout<<"block number: "<<iblk<<endl;
    cout<<"accptance rate = "<<accepted/attempted<<endl<<endl;

    stima_H = blk_av[0]/blk_norm;                                 //calcola la stima dell'hamiltoniana ed esegue il data blocking
    glob_av[0] += stima_H;
    glob_av2[0] += stima_H*stima_H;
    err_H = Error(glob_av[0],glob_av2[0],iblk);

    EGS << setw(wd) << iblk << setw(wd) << glob_av[0]/(double)iblk << setw(wd) << err_H << endl;
    EGS.close();

    for(int k = 1; k <=nbins; k++){
        stima_PSI = blk_av[k]/blk_norm;
        glob_av[k] += stima_PSI;
        glob_av2[k] += stima_PSI*stima_PSI;
        if(iblk == nblocks){                                          //calcolo il valor medio dell'autofunzione nei bin solo all'ultimo blocco
            err_PSI = Error(glob_av[k],glob_av2[k],iblk);
            histogram << setw(wd) << k << setw(wd) << glob_av[k]/(double)iblk << setw(wd) << err_PSI << endl;
        }
		}

    histogram.close();

return;
}

void Move(void){
    Xold = X;
    Xnew = Xold + rnd.Rannyu(-delta, delta);
    PSIold = PSI_mod2(Xold);
    PSInew = PSI_mod2(Xnew);
    double acc = rnd.Rannyu();
    if(PSInew/PSIold >= acc){
        accepted += 1.0;
        X = Xnew;
		}
    attempted += 1.0;

return;
}

double Hamiltonian(double x){
    return (-0.5*PSI_derivative2(x) + Potential(x)*PSI(x))/PSI(x);
}

double Potential(double x){
    return pow(x,4) - (5.0/2.0)*x*x;
}

double PSI_mod2(double x){   //valuta il modulo quadro dell'autofunzione
    double f = PSI(x);
    return f*f;
}

double PSI(double x){
    double A = exp(-((pow(x-mu,2))/(2*sigma*sigma)));
    double B = exp(-((pow(x+mu,2))/(2*sigma*sigma)));
    return A + B;
}

double PSI_derivative2(double x){
    double A = exp(-((pow(x-mu,2))/(2*sigma*sigma)));
    double B = exp(-((pow(x+mu,2))/(2*sigma*sigma)));
    return pow((x-mu)/(sigma*sigma),2)*A + pow((x+mu)/(sigma*sigma),2)*B - (A+B)/(sigma*sigma);
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<=nbins; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<=nbins; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;

return;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1){
     return 0;
	}
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/((double)iblk -1.));
}