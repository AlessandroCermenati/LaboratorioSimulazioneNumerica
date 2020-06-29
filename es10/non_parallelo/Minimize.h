#ifndef __es10code_
#define __es10code_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//coordinate delle città
const int MNoC = 50;								//numero massimo di città
int NoC;											//numero di città (letto da input)
double Xcd[MNoC], Ycd[MNoC];						//coordinate delle città

//Metropolis
double acceptance;									//variabili per la misura dell'accettanza del Metropolis
int accepted;

//Simulated Annealing
double Tinit, Tstep;								//temperatura iniziale (da input) del sistema e raffreddamento ad ogni step di raffreddamento
int freezingstep;									//numero di step di raffreddamento (passato da input)
double Teff,beta;									//temperatura "attuale"

// array per l'utilizzo dei percorsi
int individual[MNoC];								//percorso salvato
int newindividual[MNoC];							//nuovo percorso generato per la mossa Metropolis
double lenght;										//lunghezza del cammino della configurazione salvata
double newlenght;									//lunghezza del cammino della nuova configurazione
int elite[MNoC];									//per salvare il miglior percorso trovato
double elitelenght;

//check dei cromosomi
int somma, prodotto;								//inizializzati in input (noto il numero di città), sono somma e prodotto degli interi fino a NoC per il check

//algoritmo
int iterations;										//numero di iterazioni del Metropolis (passato da input) per ogni step di raffreddamento
int generate;										//passata da input, se uguale a zero genera le città sulla circonferenza, altrimenti in un quadrato

//funzioni
void Input(void);									//legge l'input
void GenerateCirc(void);							//genera le città su una circonferenza
void GenerateSquare(void);							//genera le città in un quadrato
void Measure(void);									//calcola la lunghezza del percorso creato
void Copy(void);									//salva le informazioni del percorso creato
void Exchange(int);									//scambia l'indice passato ad argomento con uno scelto casualmente
void Savepath(void);								//salva le coordinate delle città nell'ordine indicato dal percorso salvato
void Savedata(void);								//salva il percorso e la sua lunghezza in funzione della temperatura del sistema
int Pbc(int);										//implementa le periodic boundary condition per il calcolo del percorso
void Check(void);									//controlla che il percorso rispetti le boundary conditions
void NewPath(void);									//crea un nuovo percorso
void Elitair(void);									//salva sempre il miglior percorso visitato
#endif
