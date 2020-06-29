#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Minimize.h"

using namespace std;

int main()
{
    Input();
    for(int a=1; a<=iterations; a++){
        if(a%print == 0){
            cout <<"Iterazione numero "<<a<<"/"<<iterations<<endl;  
		}
        for(int j=0; j<population; j++){
            Check(j);
            Test(j);
	    }
        Ordering();
        Savedata();
        if(a%print == 0){
            cout <<"Lunghezza del percorso migliore = "<<lenght[0]<<endl<<"inizio il crossover"<<endl;
		}
        Elitarism();            //rende l'algoritmo elitario: i primi due individui della nuova popolazione sono i due migliori individui della vecchia popolazione
        Place(0,1);
        for(int b = 2; b<population; b=b+2){
            Selection();
            Crossprob = rnd.Rannyu();
            if(Crossprob<=0.7){
                Crossover();
            }
            else{
                NotCrossover();
			}
            Place(b,b+1);
		}
        Replace();
        if(a%print == 0){
            cout<<"Mutazione dei geni"<<endl;  
		}
        for(int c = 1; c<population; c++){
            Mutation(c);  
		}
        if(a%print == 0){
            cout<<"Fine dell'iterazione"<<endl<<"-----------------------------------------------------------"<<endl<<endl;  
		}
    }
    Savepath();
  return 0;
}

void Replace(){
    for(int k = 0; k<NoC; k++){
        for(int j = 0; j<population; j++){
            individual[j][k] = newindividual[j][k];  
		}
	}
    return;
}

void Mutation(int x){
    double p = rnd.Rannyu();
    if(p<0.4){                                                      //poichè Permutation la uso per costruire la popolazione iniziale provando a scambiare tutti i geni con probabilità del 15%, qui
        int y = 1 + rnd.Rannyu()*(NoC-1);                           //la chiamo con probabilità del 40% (complessivamente <10%) e per un solo gene scelto anch'esso a caso
        Permutation(x,y);
	}
    p = rnd.Rannyu();
    if(p<0.1){
        Flip(x);
	}
    p = rnd.Rannyu();
    if(p<0.1){
        Reverse(x);
	}
    p = rnd.Rannyu();
    if(p<0.8){
        FallBack(x);
	}
    p = rnd.Rannyu();
    if(p<0.8){
        FallForw(x);
	}
    return;
}

void Elitarism(){
    for(int k=0; k<NoC; k++){
        new1[k] = individual[0][k];
        new2[k] = individual[1][k];
	}
    return;
}

void NotCrossover(){
    for(int k=0; k<NoC; k++){
        new1[k] = individual[selectedindex1][k];
        new2[k] = individual[selectedindex2][k];
	}
    return;
}

void Crossover(){
    int Crossind = 1 + (NoC/2)*rnd.Rannyu();                        //stabilisce da dove parte il crossover
    int Crosslenght = rnd.Rannyu()*(NoC/2 -1);                      //stabilisce la lunghezza del crossover
    int i=0;    //indice
    int a;      //variabile per il riordinamento
    for(int k = 0; k< Crossind; k++){                               //copio nei figli le parti che non fanno crossover
        new1[k] = individual[selectedindex1][k];
        new2[k] = individual[selectedindex2][k];
	}
    for(int k = Crosslenght; k<NoC; k++){
        new1[k] = individual[selectedindex1][k];
        new2[k] = individual[selectedindex2][k];
	}
    for(int k = Crossind; k < Crossind+Crosslenght; k++){           //copio in vettori esterni le parti che fanno crossover
        Crossing1[i] = individual[selectedindex1][k];
        Crossing2[i] = individual[selectedindex2][k];
        i++;
	}
    for(i=0; i<Crosslenght; i++){                                   //inserisco in index la posizione a cui si trovano gli elementi che fanno crossover nell'altro genitore
        for(int k=0; k<NoC; k++){
            if(individual[selectedindex1][k] == Crossing2[i]){
                index2[i] = k;
			}
            if(individual[selectedindex2][k] == Crossing1[i]){
                index1[i] = k;     
			}
	    }
    }
    for(i=0; i<Crosslenght-1; i++){                                 //riordino gli indici, quindi i valori che fanno crossover, secondo l'ordine in cui si trovano nell'altro genitore
        if(index1[i]>index1[i+1]){
              a = index1[i];
              index1[i] = index1[i+1];
              index1[i+1] = a;
              a = Crossing1[i];
              Crossing1[i] = Crossing1[i+1];
              Crossing1[i+1] = a;
		}
        if(index2[i]>index2[i+1]){
              a = index2[i];
              index2[i] = index2[i+1];
              index2[i+1] = a;
              a = Crossing2[i];
              Crossing2[i] = Crossing2[i+1];
              Crossing2[i+1] = a;
		}
        for(int j = i; j>0; j--){
             if(index1[j]<index1[j-1]){
              a = index1[j];
              index1[j] = index1[j-1];
              index1[j-1] = a;
              a = Crossing1[j];
              Crossing1[j] = Crossing1[j-1];
              Crossing1[j-1] = a;
		     }
             if(index2[j]<index2[j-1]){
                 a = index2[j];
                index2[j] = index2[j-1];
                 index2[j-1] = a;
                a = Crossing2[j];
                Crossing2[j] = Crossing2[j-1];
                Crossing2[j-1] = a;
		     } 
		}
	}
    i=0;
    for(int k= Crossind; k<Crossind+Crosslenght; k++){              //inserisco nei figli le parti di crossover riordinate secondo l'ordine di apparizione nell'altro genitore
        new1[k] = Crossing1[i];
        new2[k] = Crossing2[i];
        i++;
	}
    return;
}

void Place(int x, int y){                                         //sostituisce gli individui x e y della popolazione con i due nuovi individui generati da crossover
    for(int k=0; k<NoC; k++){
        newindividual[x][k] = new1[k];
        newindividual[y][k] = new2[k];
	}
    return;
}

void Selection(){                                                   //seleziona 2 individui da riprodurre per crossover e mutazione
    int i = (double)population*pow(rnd.Rannyu(),2.5);
    int j = (double)population*pow(rnd.Rannyu(),1.2);
    selectedindex1 = i;
    selectedindex2 = j;
    return;
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

    cout<<"Algoritmo Genetico per la minimizzazione del cammino del commesso viaggiatore"<<endl<<endl;

    ReadInput >> NoC;           //numero di città

    cout<<"Numero di città da visitare: "<<NoC<<endl;

    ReadInput >> population;    //numero di individui nella popolazione

    cout<<"Numero di cromosomi nella popolazione iniziale: "<<population<<endl;

    ReadInput >> iterations;    //numero di iterazioni di crossover + mutazione

    cout<<"Numero di iterazioni dell'algoritmo genetico (Crossover + Mutazioni): "<<iterations<<endl;

    ReadInput >> generate;      //opzione per generare le città su una circonferenza o in un quadrato

    ReadInput >> print;

    if(generate == 0){
        GenerateCirc();
        cout<<"Città generate in maniera casuale su una circonferenza di raggio 10"<<endl;
	}
    else{
        GenerateSquare();
        cout<<"Città generate in maniera casuale in un quadrato di lato 10"<<endl;
	}
    cout << "------------------------------------------------------------"<<endl<<endl;
    cout << "Genero la popolazione iniziale"<<endl<<endl;

    StartGenetic();             //crea la popolazione iniziale

    somma = 0;
    prodotto = 1;

    for(int i = 1; i<NoC; i++){
        somma += i;
        prodotto *= i;
	}

    cout<<"valori di check:"<<endl<<"somma = "<< somma <<endl<<"prodotto = "<< prodotto <<endl<<endl;
    return;
}

void NewGeneration(int j){
    for(int i = 0; i<NoC; i++){
        individual[j][i] = i;
	}
    for(int i = 1; i<NoC; i++){
        Permutation(j,i);
	}
    return;
}

void Check(int j){
    int sumcheck = 0;
    int prodcheck = 1;
    for(int i = 1; i<NoC; i++){
        sumcheck += individual[j][i];
        prodcheck *= individual[j][i];
	}
    if(sumcheck != somma || prodcheck != prodotto || individual[j][0] != 0){
        cout<<"errore nel cromosoma "<<j<<":"<<endl;
        cout<<"indice di partenza: "<<individual[j][0]<<endl;
        cout<<"somma degli indici = "<<sumcheck<<"    somma attesa = "<<somma<<endl;
        cout<<"prodotto degli indici = "<<prodcheck<<"    prodotto atteso = "<<prodotto<<endl;
        NewGeneration(j);
	}
    return;
}

void Ordering(){                                               //riordina la popolazione in funzione crescente della lunghezza del percorso
    for(int j = 0; j<population-1; j++){
        if(lenght[j] > lenght[j+1]){
            for(int k=0; k<NoC; k++){
                exchangecity = individual[j][k];
                individual[j][k] = individual[j+1][k];
                individual[j+1][k] = exchangecity;
			}
            exchangelenght = lenght[j];
            lenght[j] = lenght[j+1];
            lenght[j+1] = exchangelenght;
		}
        for(int i = j; i>0; i--){
            if(lenght[i]<lenght[i-1]){
                for(int k=0; k<NoC; k++){
                exchangecity = individual[i][k];
                individual[i][k] = individual[i-1][k];
                individual[i-1][k] = exchangecity;
			}
            exchangelenght = lenght[i];
            lenght[i] = lenght[i-1];
            lenght[i-1] = exchangelenght;     
			}  
		}
	}
    return;
}

void Test(int j){                               //calcola la funzione costo per l'individuo j
    double sum = 0.;
    double l;
    for(int i=0; i<NoC; i++){
        l = sqrt(pow(Xcd[individual[j][i]] - Xcd[individual[j][Pbc(i+1)]],2) + pow(Ycd[individual[j][i]] - Ycd[individual[j][Pbc(i+1)]],2));
        sum += l;
	}
    lenght[j] = sum;
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

void StartGenetic(){                            //genera la popolazione iniziale
    for(int j=0; j<population; j++){
        for(int i=0; i<NoC; i++){
            individual[j][i] = i;               //genero la configurazione di partenza in maniera ordinata
	    }
    }
    for(int j=0; j<population; j++){            //creo configurazioni casuali (con lo stesso punto di partenza)
        for(int i=1; i<NoC; i++){
              Permutation(j,i);
		}
	}
    return;
}

void Permutation(int j, int i){                //scambia il gene i con un gene scelto a caso tra gli altri (tranne il primo) per l'individuo j, con probabilità del 15%
    double prob = rnd.Rannyu();
    if(prob <= 0.3){
        int ind = 1 + rnd.Rannyu()*(double)(NoC-1);
        int a = individual[j][i];
        individual[j][i] = individual[j][ind];
        individual[j][ind] = a;
	}
    return;
}

void Flip(int x){
    int y = 1 + rnd.Rannyu()*(NoC-2);
    int a = individual[x][y];
    individual[x][y] = individual[x][y+1];
    individual[x][y+1] = a;
    return;
}

void Reverse(int x){
    int rev = 1 + 4*rnd.Rannyu();
    int y = 1 + rnd.Rannyu()*(NoC-rev);
    for(int i = 0; i<(double)rev/2.; i++){
        int a = individual[x][y+i];
        individual[x][y+i] = individual[x][y+rev-1-i];
        individual[x][y+rev-1-i] = a;
	}
    return;
}

void FallForw(int x){
    int forw = 1 + (NoC - 8)*rnd.Rannyu();
    int y = 1+ rnd.Rannyu()*(NoC-forw);
    for(int i = y; i<y+forw-1; i++){
        int a = individual[x][i];
        individual[x][i] = individual[x][i+1];
        individual[x][i+1] = a;
	}
    return;
}

void FallBack(int x){
    int back = 1 + (NoC-8)*rnd.Rannyu();
    int y = 1 + back +(NoC-back-1)*rnd.Rannyu();
    for(int i = y; i>y-back; i--){
        int a = individual[x][i];
        individual[x][i] = individual[x][i-1];
        individual[x][i-1] = a;
	}
    return;
}

void Savedata(){
    ofstream Cities;
    ofstream Lenght;
    ofstream MeanLenght;
    Cities.open("Best_Path_Iterator.dat",ios::app);
    Lenght.open("Lenght_of_the_best_path.dat",ios::app);
    MeanLenght.open("Mean_Path_half_population.dat",ios::app);
    const int wd = 3;
    for(int j = 0; j < NoC; j++){
        Cities << individual[0][j] << setw(wd);
	}
    Cities << endl;
    Lenght << lenght[0] << endl;
    Cities.close();
    Lenght.close();
    double sum = 0;
    for(int k=0; k<population/2; k++){
        sum += lenght[k];
    }
    sum = sum*2/population;
    MeanLenght << sum << endl;
    MeanLenght.close();
    return;
}

void Savepath(){
    ofstream Output;
    const int wd = 12;
    if(generate==0){
        Output.open("CircumferencePath.dat",ios::app);
	}
    else{
        Output.open("SquarePath.dat",ios::app);
	}
    Ordering();
    for(int j = 0; j<=NoC; j++){
        Output << Xcd[individual[0][Pbc(j)]] << setw(wd) << Ycd[individual[0][Pbc(j)]] << endl;
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