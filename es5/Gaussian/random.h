/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

#include "funzionebase.h"

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
  double Ranpolar();                                                        //campiona uniformemente un angolo fra 0 e pi greco (angolo polare)
  double Ranazim();                                                         //campiona uniformemente un angolo fra 0 e 2pi greco (angolo azimutale)
  double Exp(double A);                                                     //campionamento esponenziale col metodo della trasformata
  double Lorentz(double mu, double gamma);                                  //campionamento Lorentziano col metodo della trasformata
  double hitmiss(funzionebase *F, double a, double b);                      //campionamento, col metodo di reiezione, secondo una funzione passata come classe funzionebase, nell'intervallo (a,b)
  double importance_sampling(funzionebase *F, double a, double b);          //importance sampling di una funzione; la pdf da usare è contenuta in funzionebase, nell'intervallo (a,b)
  double GBM(double start_value, double mean, double sigma, double T);      //esegue GBM(mean, sigma) per un tempo T senza discretizzazione
  double GMB_step(double start_value, double mean, double sigma, double T, int N);     //esegue GBM(mean, sigma) per un tempo T, dividendo l'intervallo (0,T) in N sottointervalli
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
