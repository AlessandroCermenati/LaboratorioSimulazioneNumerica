#ifndef __es9code_
#define __es9code_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//coordinate delle città
const int MNoC = 50;								//numero massimo di città
int NoC;											//numero di città (letto da input)
double Xcd[MNoC], Ycd[MNoC];						//coordinate delle città

//popolazione
const int Mpopulation = 500;					    //numero massimo di individui della popolazione
int population;										//numero di individui della popolazione (passato da input)
int individual[Mpopulation][MNoC];					//individui nella popolazione; il primo indice corre sugli individui, il secondo sull'ordine delle città da testare
int exchangecity;									//serve per riordinare la popolazione
double exchangelenght;								//serve per riordinare la popolazione

//Selezione e riproduzione
double lenght[Mpopulation];							//contiene la lunghezza del percorso per ogni individuo della popolazione
int selectedindex1, selectedindex2;					//contengono gli indici dei genitori selezionati per la riproduzione
int new1[MNoC], new2[MNoC];							//contengono i figli da inserire poi nella popolazione
double Crossprob;									//variabile per la probabilità di fare crossover
int Crossing1[MNoC], Crossing2[MNoC];				//contengono le città che fanno crossover
int index1[MNoC], index2[MNoC];						//contengono la posizione nel cammino dell'altro genitore delle città che fanno crossover
int newindividual[Mpopulation][MNoC];				//contiene la nuova popolazione generata

//parallelo
int migrations;										//ogni quante iterazione faccio migrazione dei cromosomi
int send[MNoC];										//cromosoma da scambiare
double lenghtsend;									//lunghezza del cromosoma da scambiare

//check dei cromosomi
int somma, prodotto;								//inizializzati in input (noto il numero di città), sono somma e prodotto degli interi fino a NoC per il check

//algoritmo
int iterations;										//numero di iterazioni (generazioni di nuova popolazione con crossover e mutazione) (passato da input)
int generate;										//passata da input, se uguale a zero genera le città sulla circonferenza, altrimenti in un quadrato
int print;											//ogni quante iterazione avviene la stampa a terminale (da input)

//funzioni
void Input(int, int);								//legge l'input
void GenerateCirc(void);							//genera le città su una circonferenza
void GenerateSquare(void);							//genera le città in un quadrato
void StartGenetic(void);							//crea la popolazione iniziale
void Test(int);										//calcola la lunghezza del percorso dell'individuo passato come indice all'argomento
void Ordering();									//riordina la popolazione in funzione crescente della lunghezza del percorso
void Selection();									//selezione due genitori per la riproduzione
void Crossover();									//genera due figli riordinando parte dei cromosomi dei genitori nell'ordine in cui i geni si trovano nell'altro genitore
void NotCrossover();								//genera due figli identici ai genitori
void Elitarism();									//genera due figli identici ai due migliori cromosomi della popolazione
void Place(int,int);								//mette i nuovi cromosomi per creare la successiva popolazione in newindividual
void Replace();										//sostituisce la vecchia popolazione con la nuova
void Mutation(int);									//esegue, con una certa probabilità, le funzione di mutazione genetica per un individuo passato ad argomento
void Permutation(int, int);							//mutazione genetica: scambia un gene con uno scelto a caso
void Flip(int);										//mutazione genetica: scambia due geni successivi
void FallBack(int);									//mutazione genetica: sceglie casualmente un gene e lo "porta indietro" di un numero casuale di posizioni
void FallForw(int);									//mutazione genetica: sceglie casualmente un gene e lo "porta avanti" di un numero casuale di posizioni
void Reverse(int);									//mutazione genetica: scaglie due geni distanti massimo 5 posizioni e inverte il loro ordine
void Savepath(int);									//salva le coordinate delle città nell'ordine indicato dal percorso migliore della popolazione
void Savedata(int);									//salva il percorso migliore e la sua lunghezza ad ogni iterazione
int Pbc(int);										//implementa le periodic boundary condition per il calcolo del percorso
void Check(int);									//controlla che nel cromosoma il primo indice sia 0 (città di partenza fissata) e che le città siano contenute tutte e una sola volta
void NewGeneration(int);							//se il check trova un errore, rimpiazza il cromosoma sbagliato con uno completamente nuovo

//funzioni per il calcolo in parallelo
void Migrazione(int, int);							//chiama la migrazione
void nodes2(int);									//eseguono la migrazione (da 2 a 4 nodi)
void nodes3(int);
void nodes4(int);
void CopySend();									//copia le variabili da inviare nelle variabili send per la migrazione
void CopyRecv();									//copia le variabili ricevute nelle variabili utilizzate dall'algoritmo
void SendCoordinates(int, int);						//chiama l'invio delle coordinate tra i nodi
#endif
