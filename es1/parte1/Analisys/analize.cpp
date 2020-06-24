#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main(){

	ifstream file;
	file.open("random_values.txt");							//leggo il file contenente i numeri pseudocasuali da studiare
	ofstream mean;											//creo i file di output, che contengono i dati da esporre nel notebook
	ofstream stddev;
	ofstream variance;
	ofstream stddevvariance;
	ofstream chiquadro;
	ofstream results;
	mean.open("mean.txt");
	stddev.open("stddev.txt");
	variance.open("variance.txt");
	stddevvariance.open("stddev_variance.txt");
	chiquadro.open("chi_quadro.txt");
	results.open("results.txt");

	int M = 100000;										//imposta il numero di step montecarlo
	int N=100;											//imposta il numero di blocchi per l'analisi e di intervalli per il calcolo del chi quadro
	int L = M/N;
	double y[M];
	double AVE[N];
	double AV2[N];
	double sum_prog[N];
	double sum2_prog[N];
	double devstd[N];
	double chi2[N];
	int counter = 0;
	double meanvalue = 0.5;								//definisco il valore atteso per il calcolo della varianza
	int expected[N];									//vettore che deve essere riempito coi valori di aspettazione in ogni intervallo per il calcolo del chi quadro
	double a=0.;
	double b=1.;									//estremi dell'intervallo della distribuzione di numeri casuali, per il calcolo del chi quadro
	double chi2tot = 0.;

	for(int i=0; i<N; i++){
		AVE[i]=0;
		AV2[i]=0;
		sum_prog[i]=0;
		sum2_prog[i]=0;
	}

	for(int i=0; i<M; i++){								//leggo i valori da analizzare dal file
		file >> y[i];
	}

	file.close();

	double sum = 0;

	for(int i=0; i<N; i++){								//calcolo del valor medio e della sua incertezza statistica
		sum = 0;
		for(int j=0; j<(L); j++){
			int k = j + i*(L);
			sum = sum + y[k];
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

	for(int i=0; i<N; i++){									//scrittura dei risultati per valor medio e incertezza associata ad ogni blocco, o stima
		mean<<sum_prog[i]<<endl;
		stddev<<devstd[i]<<endl;
	}

	for(int i=0; i<N; i++){
		AVE[i]=0;
		AV2[i]=0;
		sum_prog[i]=0;
		sum2_prog[i]=0;
		devstd[i]=0;
	}

	for(int i=0; i<N; i++){											//calcolo della varianza e della sua incertezza statistica
		sum = 0;
		for(int j=0; j<(L); j++){
			int k = j + i*(L);
			sum = sum + pow((y[k]-meanvalue),2);
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
	for(int i=0; i<N; i++){										//scrittura dei dati sulla varianza e incertezza associata ad ogni blocco, o stima
		variance<<sum_prog[i]<<endl;
		stddevvariance<<devstd[i]<<endl;
	}

	for(int i=0; i<N; i++){										//fisso il valore di aspettazione in ogni intervallo per il calcolo del chi quadro
		expected[i]=L/N;
	}


	for(int i=0; i<N; i++){										//conto, per ogni sottointervallo, il numero di valori di un set degli step che cadono nel sottointervallo
		counter=0;
		for(int j=0; j<L; j++){
			int k = j + i*L;
			if(y[k] >= a + ((i + 1. - 1.)/N)*(b-a) && y[k] < a + ((i+1.)/N)*(b-a)){
				counter++;
			}
		}
		chi2[i] = (pow((counter-expected[i]),2))/expected[i];			//calcolo il chi quadro in ogni intervallo
		chi2tot = chi2tot + chi2[i];
	}

	for(int i=0; i<N; i++){										//scrivo il chi quadro in ogni intervallo
		chiquadro<<chi2[i]<<endl;
	}

	results<<chi2tot<<endl;										//salvo il valore del chi quadro nel file results

	mean.close();
	stddev.close();
	variance.close();
	stddevvariance.close();
	chiquadro.close();
	results.close();

return 0;
}
	
