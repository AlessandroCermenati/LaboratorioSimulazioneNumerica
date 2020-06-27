#ifndef __es8code_
#define __es8code_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//osservabili
const int m_bins = 1000;														//numero massimo di bin per l'istogramma dell'autofunzione
int nbins;																					//numero di bin per l'istogramma dell'autofunzione
double hamiltonian;																	//misuratore per l'energia
const double hist_min = -2.5;												//valori minimo e massimo dell'istogramma
const double hist_max = 2.5;
double binl;																				//lunghezza dei bin

//data-blocking
double blk_av[m_bins+1],blk_norm,accepted,attempted;
double glob_av[m_bins+1],glob_av2[m_bins+1];
double stima_H, stima_PSI, err_H, err_PSI;

//simulazione
int nsteps, nblocks;																//numero di blocchi e step per blocco

//Metropolis
double Xold, Xnew, X = 0.1;													//punto di partenza e finale delle mosse Metropolis
double PSIold, PSInew;															//per l'algoritmo Metropolis
double delta;																				//step per la mossa Metropolis

//trial wavefunction
double mu, sigma;							//parametri dell'autofunzione (impostati in input)
double cool;								//raffreddamento del simulated annealing (da input)
double deltamu, deltasig;					//estremi dell'intervallo di ricerca dei parametri ottimizzati (da input)
int optimize;								//se =1, fa l'ottimizzazione, altrimenti, se esiste il file contenente i parametri già ottimizzati, non la esegue e carica i valori da file

//funzioni
void Input(void);
void Optimize(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
double Error(double,double,int);
double PSI_mod2(double);
double PSI(double);
double PSI_derivative2(double);
double Hamiltonian(double);
double Potential(double);

#endif
