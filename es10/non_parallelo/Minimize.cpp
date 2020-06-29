#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Minimize.h"

using namespace std;

int main()
{
    Input();                                    //leggo le informazioni di input
    for(int i = 0; i < freezingstep; i++){
        cout << "freezing step numero " << i+1 <<"/" << freezingstep << endl << "Temperatura del sistema = " << Teff << endl;
        accepted = 0;
        for(int j = 0; j < iterations; j++){
            NewPath();                          //creo un nuovo percorso
            Check();                            //controllo che rispetti le boundary conditions
            Measure();                          //calcolo la lunghezza del percorso creato
            if(newlenght < elitelenght){
                Elitair();     
			}
            double prob = rnd.Rannyu();
            if(prob < exp((-1.)*beta*(newlenght - lenght))){
                Copy();                         //con probabilità Metropolis accetto il nuovo percorso
                accepted++;
			}
		}
        acceptance = (double)accepted/(double)iterations;
        Savedata();
        Teff -= Tstep;                          //raffreddo il sistema
        beta = 1/Teff;
        cout << "Lunghezza del cammino = " << lenght << endl << "Accettanza del Metropolis = " << acceptance << endl << "-------------------------------------" << endl << endl;
	}
    Savepath();
  return 0;
}

void Input(){
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();

    ifstream ReadInput;
    ReadInput.open("input.dat");

    cout<<"Algoritmo Simulated Annealing per la minimizzazione del cammino del commesso viaggiatore"<<endl<<endl;

    ReadInput >> NoC;           //numero di città
    somma = 0;
    prodotto = 1;

    for(int i = 1; i<NoC; i++){
        somma += i;
        prodotto *= i;
	}

    cout<<"Valori di check:"<<endl<<"somma = "<< somma <<endl<<"prodotto = "<< prodotto <<endl<<endl;
    cout<<"Numero di città da visitare: "<<NoC<<endl;

    ReadInput >> Tinit;         //temperatura iniziale del sistema
    Teff = Tinit;
    beta = 1/Teff;

    cout<<"Temperatura iniziale del sistema: "<<Tinit<<endl;

    ReadInput >> freezingstep;    //numero di raffreddamenti del sistema
    Tstep = Tinit/(double)freezingstep;

    cout<<"Numero di step di raffreddamento dell'algoritmo : "<<freezingstep<<endl;
    cout<<"Raffreddamento ad ogni step = "<<Tstep<<endl;

    ReadInput >> iterations;

    cout<<"Numero di mosse Metropolis ad ogni step di temperatura = "<<iterations<<endl;

    ReadInput >> generate;      //opzione per generare le città su una circonferenza o in un quadrato

    if(generate == 0){
        GenerateCirc();
        cout<<"Città generate in maniera casuale su una circonferenza di raggio 10"<<endl;
	}
    else{
        GenerateSquare();
        cout<<"Città generate in maniera casuale in un quadrato di lato 10"<<endl;
	}
    cout << "------------------------------------------------------------"<<endl<<endl;
    cout << "Genero il percorso iniziale"<<endl;

    for(int i = 0; i < NoC; i++){
        individual[i] = i;
	}
    NewPath();             //crea il percorso iniziale
    Check();
    Measure();             //misura la lunghezza del percorso iniziale
    Elitair();             //salva il percorso iniziale come miglior percorso
    Copy();                //salva il percorso creato come percorso iniziale

    return;
}

void Elitair(){
    for(int i=0; i<NoC; i++){
        elite[i] = newindividual[i];
	}
    elitelenght = newlenght;
    return;
}

void Check(){                               //esegue il check delle boundary conditions
    int sumcheck = 0;
    int prodcheck = 1;
    for(int i = 1; i<NoC; i++){
        sumcheck += newindividual[i];
        prodcheck *= newindividual[i];
	}
    if(sumcheck != somma || prodcheck != prodotto || newindividual[0] != 0){
        cout<<"errore nel percorso!"<<endl;
        cout<<"indice di partenza: "<<newindividual[0]<<endl;
        cout<<"somma degli indici = "<<sumcheck<<"    somma attesa = "<<somma<<endl;
        cout<<"prodotto degli indici = "<<prodcheck<<"    prodotto atteso = "<<prodotto<<endl;
        NewPath();
	}
    return;
}

void Measure(){                               //calcola la lunghezza del nuovo percorso proposto
    double sum = 0.;
    double l;
    for(int i=0; i<NoC; i++){
        l = sqrt(pow(Xcd[newindividual[i]] - Xcd[newindividual[Pbc(i+1)]],2) + pow(Ycd[newindividual[i]] - Ycd[newindividual[Pbc(i+1)]],2));
        sum += l;
	}
    newlenght = sum;
    return;
}

void NewPath(){                             //genera un nuovo percorso
    for(int j=0; j<NoC; j++){
        newindividual[j] = individual[j];               //parto in maniera ordinata
    }
    int j = 1 + rnd.Rannyu()*(NoC - 1);
    Exchange(j);
    return;
}

void Copy(){                                //copia le informazioni del nuovo percorso al posto di quello vecchio
    for(int i=0; i<NoC; i++){
        individual[i] = newindividual[i];
	}
    lenght = newlenght;
    return;
}

void Exchange(int j){                                   //scambia la città all'indice j con un indice casuale
    int ind = 1 + rnd.Rannyu()*(NoC - 1);
    double l1 = sqrt(pow(Xcd[individual[ind - 1]]-Xcd[individual[ind]],2) + pow(Ycd[individual[ind - 1]]-Ycd[individual[ind]],2)) + sqrt(pow(Xcd[individual[ind]]-Xcd[individual[Pbc(ind+1)]],2) + pow(Ycd[individual[ind]]-Ycd[individual[Pbc(ind+1)]],2));
    double l2 = sqrt(pow(Xcd[individual[j - 1]]-Xcd[individual[j]],2) + pow(Ycd[individual[j - 1]]-Ycd[individual[j]],2)) + sqrt(pow(Xcd[individual[j]]-Xcd[individual[Pbc(j+1)]],2) + pow(Ycd[individual[j]]-Ycd[individual[Pbc(j+1)]],2));
    double prob = rnd.Rannyu();
    if(prob < l2/l1){
        int a = newindividual[j];
        newindividual[j] = newindividual[ind];
        newindividual[ind] = a;
	}
    return;
}

void GenerateCirc(){                            //genera le coordinate delle città da visitare su una circonferenza di raggio 10
    double R = 10.;
    for(int i=0; i<NoC; i++){
        double s = rnd.Rannyu(0., 2*M_PI);
        Xcd[i] = R*cos(s);
        Ycd[i] = R*sin(s);
	}
    return;
}

void GenerateSquare(){                          //genera le coordinate delle città da visitare in un quadrato di lato 10
    for(int i=0; i<NoC; i++){
        Xcd[i] = rnd.Rannyu(0., 10.);
        Ycd[i] = rnd.Rannyu(0., 10.);
	}
    return;
}

void Savedata(){                                //salva il percorso scelto e la sua lunghezza in funzione della temperatura del sistema
    ofstream Cities;
    ofstream Lenght;
    Cities.open("cooling_cities.dat",ios::app);
    Lenght.open("cooling_lenght.dat",ios::app);
    const int wd = 3;
    const int dwd = 12;
    Cities << Teff << setw(dwd);
    for(int j = 0; j < NoC; j++){
        Cities << individual[j] << setw(wd);
	}
    Cities << endl;
    Lenght << Teff << setw(dwd) << lenght << endl;
    Cities.close();
    Lenght.close();
    return;
}

void Savepath(){                                //salva le coordinate delle citàà nell'ordine del percorso salvato
    ofstream Output;
    ofstream Loutput;
    const int wd = 12;
    if(generate==0){
        Output.open("CircumferencePath.dat",ios::app);
	}
    else{
        Output.open("SquarePath.dat",ios::app);
	}
    Loutput.open("Minimum_path_lenght.dat",ios::app);
    Loutput << elitelenght;
    Loutput.close();
    for(int j = 0; j<=NoC; j++){
        Output << Xcd[elite[Pbc(j)]] << setw(wd) << Ycd[elite[Pbc(j)]] << endl;
	}
    Output.close();
    return;
}

int Pbc(int i){                                     //Periodic boundary condition per il calcolo del percorso
    if(i==NoC){
        return 0;
	}
    else{
        return i;
	}
}