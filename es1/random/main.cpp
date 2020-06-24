/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
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

//inizio modifica//
   ofstream output;                                 //parte 1: scrivo un file con 100000 numeri casuali con distribuzione uniforme in (0-1)
   output.open("random_values.txt");

   for(int i=0; i<100000; i++){
      output << rnd.Rannyu() << endl;
   }
   output.close();

   ofstream output2;                                //parte 2: genero 3 file con 1000000 valori casuali secondo le distribuzioni date, implementando i metodi della trasformata
   ofstream output3;
   ofstream output4;
   output2.open("unif_dice.txt");
   output3.open("exp_dice.txt");
   output4.open("Lor_dice.txt");

   for(int i=0; i<1000000; i++){
       output2<<rnd.Rannyu()<<endl;
       output3<<rnd.Exp(1.)<<endl;
       output4<<rnd.Lorentz(0. , 1.)<<endl;
   }
   output2.close();
   output3.close();
   output4.close();

   ofstream data;                               //parte 3: genero un file contenente coppie di numeri pseudocasuali, per la posizione del centro dell'ago e per la sua inclinazione
   data.open("Buffon.txt");
   double xc;
   double x1;
   double x2;
   double y1;
   double y2;
   double theta;

   for(int i=0; i<100000; i++){
        xc = rnd.Rannyu();                      //posizione del centro dell'ago, distribuita uniformemente fra 0 e 1 (utilizzo 1 come distanza fra le maglie della griglia)
        x1 = rnd.Rannyu();
        x2 = rnd.Rannyu();
        y1 = rnd.Rannyu();
        y2 = rnd.Rannyu();
        theta = atan((y2-y1)/(x2-x1));          //genera l'inclinazione in modo uniforme fra meno pi/2 e pi/2
        data<<xc<<" "<<theta<<endl;
   }
   data.close();
//fine modifica//

   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
