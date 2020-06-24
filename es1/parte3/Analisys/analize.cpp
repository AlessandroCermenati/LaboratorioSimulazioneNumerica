#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main(){

	ifstream data;									//carico il file contenente 2 valori pseudocasuali: xc uniforme in (0-1), theta uniforme in (-pi/2 ; pi/2).
	data.open("Buffon.txt");	

	ofstream mean;									//creo i file di output necessari
	ofstream sigma;

	mean.open("mean_pi.txt");
	sigma.open("sigma_pi.txt");

	int M = 100000;										//totale dei lanci casuali generati
	int N = 100;										//numero di blocchi
	int L = M/N;										//numero di lanci in ogni blocco
	double AVE[N];										//vettori per la stima del valore di pi greco e della sua incertezza ad ogni blocco
	double AV2[N];
	double sum_prog[N];
	double sum2_prog[N];
	double devstd[N];
	int counter=0;										//ci servirà un contatore per determinare la probabilita dell'intersezione dell'ago con la griglia: p=counter/L
														//vettori contenenti i dati sul lancio
	double xc[M];
	double theta[M];

	for(int i=0; i<M; i++){
		data>>xc[i];
		data>>theta[i];
	}

	double l = 0.8;										//lunghezza dell'asta lanciata
	double d = 1.;										//distanza tra le maglie della griglia; in questo modo i valori di xc coprono tutta la griglia.
	double Xdx;											//coordinate x degli estremi dell'asta; se la griglia è verticale, l'ago la interseca se i suoi estremi escono dall'intervallo (0,1)
	double Xsx;											//per costruzione, l'estremo destro può intersecare la griglia a destra, viceversa l'estremo sinistro

	for(int i=0; i<N; i++){								//stimo pi greco in ogni blocco
		counter=0;
		for(int j=0; j<L; j++){
			int k= j + i*L;
			Xdx = xc[k] + (l/2)*cos(theta[k]);			//calcolo la posizione x delle estremità dell'ago
			Xsx = xc[k] - (l/2)*cos(theta[k]);
			if(Xdx>1.){									//conta un hit se l'ago tocca la griglia, cioè se uno dei due estremi esce dall'intervallo (0,1)
				counter++;
			}
			if(Xsx<0.){
				counter++;
			}
		}
		AVE[i] = (2*l*L)/(counter*d);
		AV2[i] = pow(AVE[i],2);
	}

	for(int i=0; i<N; i++){								//calcolo del valor medio e della sua incertezza statistica in funzione del numero di blocchi
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
		mean<<sum_prog[i]<<endl;
		sigma<<devstd[i]<<endl;
	}

	mean.close();
	sigma.close();

	return 0;
}
