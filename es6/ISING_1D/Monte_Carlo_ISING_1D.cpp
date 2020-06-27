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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  ofstream Energy, Magnetization, Susceptibility, Capacity;
  Energy.open("Energy.dat");
  Magnetization.open("Magnetization.dat");
  Susceptibility.open("Susceptibility.dat");
  Capacity.open("Heat_Capacity.dat");
  Input(); //Inizialization
  while(temp<=2.0){
      cout<< endl << "-------------------------" << endl <<"Simulazione a temperatura  "<<temp<< endl << "-------------------------" << endl << endl;
      for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
         Reset(iblk);   //Reset block averages
         for(int istep=1; istep <= nstep; ++istep){
            Move(metro);
            Measure();
            Accumulate(); //Update block averages
         }
         Averages(iblk);   //Print results for current block
      }
      ConfFinal(); //Write final configuration
      const int wd = 12;
      if(h == 0){
        Energy << setw(wd) << temp << setw(wd) << mean_u_final << setw(wd) << err_u_final << endl;                    //salva i dati finali
        Capacity << setw(wd) << temp << setw(wd) << mean_c_final << setw(wd) << err_c_final << endl;
        Susceptibility << setw(wd) << temp << setw(wd) << mean_x_final << setw(wd) << err_x_final << endl;
      }
      else{
        Magnetization << setw(wd) << temp << setw(wd) << mean_m_final << setw(wd) << err_m_final << endl;
      }
      temp += 0.05;
      beta = 1.0/temp;
      cout << "re - equilibrating the system to the new temperature"<<endl<<endl;
      for(int i=0; i< 50000; i++){
        move(metro);
	  }
  }
  Energy.close();
  Capacity.close();
  Susceptibility.close();
  Magnetization.close();
  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration

  ifstream OldConf;
  OldConf.open("config.old");
  
  if(OldConf.is_open()){        //inizializzazione degli spin da una configurazione precedente (richiede che il file config.final sia copiato in un file config.old)
    cout<< "read the old configuration saved in config.old" <<endl<<endl;
    for(int i=0; i<nspin; i++){
        OldConf >> s[i];
	}
  }

  else{     //inizializzazione random degli spin
    for (int i=0; i<nspin; ++i)
     {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
      }
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for all the observables
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
  cout << "Initial heat capacity = " << walker[ic]/(double)nspin << endl;
  cout << "Initial magnetization = " << walker[im]/(double)nspin << endl;
  cout << "Initial magnetic susceptibility = " << walker[ix]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
    // INCLUDE YOUR CODE HERE
        attempted++;
        if(s[o] == +1){             //cambia lo spin alla posizione scelta randomicamente
            sm = -1;  
		}
        if(s[o] == -1){
            sm = +1;
		}
        energy_old = Boltzmann(s[o], o);            //energia dello stato prima della mossa
        energy_new = Boltzmann(sm, o);              //energia della stato dopo la mossa
        p = rnd.Rannyu();                           //check per l'accettazione
        if(p < exp(-1.*beta*(energy_new - energy_old))){
            s[o] = sm;
            accepted++;
		}
    }
    else //Gibbs sampling
    {
    // INCLUDE YOUR CODE HERE
        energy_up = Boltzmann(+1, o);
        energy_down = Boltzmann(-1, o);
        p = rnd.Rannyu();
        if(p <= 1/(1 + exp(-1.*beta*(energy_down - energy_up)))){               //campiono la probabilità che lo spin sia +1
            s[o] = +1;
	    }
        else{
            s[o] = -1;
		}
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
// INCLUDE YOUR CODE HERE
  }
  walker[iu] = u;
  walker[ic] = pow(beta,2)*pow(u,2);
  walker[im] = m;
  walker[ix] = beta*m*m;
// INCLUDE YOUR CODE HERE
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    Heat.open("output.heat.0",ios::app);
    Mag.open("output.mag.2",ios::app);
    Chi.open("output.chi.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    stima_c = (blk_av[ic]/blk_norm - pow(beta,2)*pow(stima_u*(double)nspin,2))/nspin; //Heat capacity
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    stima_x = blk_av[ix]/blk_norm/(double)nspin; //Susceptibility
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Ene.close();
    Heat.close();
    Mag.close();
    Chi.close();
    if(iblk == nblk -1){                                //salva i dati finali
        mean_u_final = glob_av[iu]/(double)iblk;
        err_u_final = err_u;
        mean_c_final = glob_av[ic]/(double)iblk;
        err_c_final = err_c;
        mean_m_final = glob_av[im]/(double)iblk;
        err_m_final = err_m;
        mean_x_final = glob_av[ix]/(double)iblk;
        err_x_final = err_x;
	}
// INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout<<"Print final configuration to file config.final"<< endl <<"If you want to restart your simulation from this configuration, make a copy of config.final in a file named config.old"<<endl<<endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk == 1){
        return 0.;
	}
    else{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/((double)iblk-1));
    }
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
