#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main(){

	ifstream data;									//carico i file contenenti i numeri generati secondo le tre distribuzioni (unif_dice.txt, exp_dice.txt, Lor_dice.txt)
	data.open("Lor_dice.txt");						//modifico il nome del file di input per le tre distribuzioni

	ofstream mean1;									//creo i file di output necessari
	ofstream mean2;
	ofstream mean10;
	ofstream mean100;

	mean1.open("Lorentzian_1.txt");						//modifico il nome dei file di output per le 3 distribuzioni
	mean2.open("Lorentzian_2.txt");
	mean10.open("Lorentzian_10.txt");
	mean100.open("Lorentzian_100.txt");

	int M = 1000000;										//totale dei numeri casuali generati
	double sum2 = 0;
	double sum10 = 0;
	double sum100 = 0;
	double y[M];

	for(int i=0; i<M; i++){
		data>>y[i];
	}

	for(int i=0; i<10000; i++){								//10000 medie su
		mean1<<y[i]<<endl;									//N=1 numeri casuali; ovviamente la media su un solo numero coincide col numero stesso.
	}

	for(int i=10000; i<30000; i = i+2){						//N=2 numeri casuali; non utilizzo gli stessi numeri per la media con N=1;
		sum2=0;
		for(int j=0; j<2; j++){
			sum2 = sum2 + y[i+j];
		}
		mean2<<sum2/2<<endl;
	}

	for(int i=100000; i<200000; i = i+10){					//N=10 numeri casuali; non utilizzo gli stessi numeri per la media con N=1,2;
		sum10=0;
		for(int j=0; j<10; j++){
			sum10 = sum10 + y[i+j];
		}
		mean10<<sum10/10<<endl;
	}

	for(int i=0; i<1000000; i=i+100){						//N=100 numeri casuali; utilizzo tutti i numeri casuali generati;
		sum100=0;
		for(int j=0; j<100; j++){
			sum100 = sum100 + y[i+j];
		}
		mean100<<sum100/100<<endl;
	}

	mean1.close();
	mean2.close();
	mean10.close();
	mean100.close();

	return 0;
}
