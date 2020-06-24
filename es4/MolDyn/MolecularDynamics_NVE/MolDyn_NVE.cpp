/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){

    ofstream meanT;
    ofstream meanK;
    ofstream meanU;
    ofstream meanE;
    ofstream sigmaT;
    ofstream sigmaK;
    ofstream sigmaU;
    ofstream sigmaE;
    ofstream meanP;
    ofstream sigmaP;
    meanT.open("ave_temp.dat");
    meanK.open("ave_ekin.dat");
    meanU.open("ave_epot.dat");
    meanE.open("ave_etot.dat");
    meanP.open("ave_pres.dat");
    sigmaT.open("err_temp.dat");
    sigmaK.open("err_ekin.dat");
    sigmaU.open("err_epot.dat");
    sigmaE.open("err_etot.dat");
    sigmaP.open("err_pres.dat");

  Input();             //Inizialization
  int nconf = 1;
  int j=0;
  int L = nstep/(10*nblocks);                   //poichè si ha una misura ogni 10 step, L è il numero di misure in ogni blocco
  double sumT=0, sumU=0, sumK=0, sumE=0, sumP=0;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
        sumT += stima_temp;                     //sommatorie, misura per misura, per effettuare il data-blocking
        sumK += stima_kin;
        sumU += stima_pot;
        sumE += stima_etot;
        sumP += stima_pres;
     }
     if(istep%(10*L) == 0){                       //L è il numero di misure per blocco; poichè si ha una misura ogni 10 step, il blocco si conclude ogni 10*L step
        AVT[j]=sumT/L;
		AV2T[j]=pow(AVT[j],2);
        AVK[j]=sumK/L;
		AV2K[j]=pow(AVK[j],2);
        AVU[j]=sumU/L;
		AV2U[j]=pow(AVU[j],2);
        AVE[j]=sumE/L;
		AV2E[j]=pow(AVE[j],2);
        AVP[j]=sumP/L;
		AV2P[j]=pow(AVP[j],2);
        j += 1;                                 //l'indice j scorre sui blocchi;
        sumT=0;                                 //le sommatorie sono azzerate per il calcolo delle medie del blocco seguente
        sumU=0;
        sumK=0;
        sumE=0;
        sumP=0;
	 }
  }
  ConfOld();
  ConfFinal();         //Write final configuration to restart

  for(int i=0; i<nblocks; i++){                               //avendo salvato le medie a blocchi nei vettori AV, il calcolo del data blocking si può fare in seguito;
		for(int j=0; j<i+1; j++){
			progT[i]= progT[i] + AVT[j];
			prog2T[i]= prog2T[i] + AV2T[j];
            progK[i]= progK[i] + AVK[j];
			prog2K[i]= prog2K[i] + AV2K[j];
            progU[i]= progU[i] + AVU[j];
			prog2U[i]= prog2U[i] + AV2U[j];
            progE[i]= progE[i] + AVE[j];
			prog2E[i]= prog2E[i] + AV2E[j];
            progP[i]= progP[i] + AVP[j];
			prog2P[i]= prog2P[i] + AV2P[j];
		}
		progT[i] = progT[i]/(i+1);
		prog2T[i] = prog2T[i]/(i+1);
        progK[i] = progK[i]/(i+1);
		prog2K[i] = prog2K[i]/(i+1);
        progU[i] = progU[i]/(i+1);
		prog2U[i] = prog2U[i]/(i+1);
        progE[i] = progE[i]/(i+1);
		prog2E[i] = prog2E[i]/(i+1);
        progP[i] = progP[i]/(i+1);
		prog2P[i] = prog2P[i]/(i+1);
		if(i==0){
			errT[i] = 0;
            errK[i] = 0;
            errU[i] = 0;
            errE[i] = 0;
            errP[i] = 0;
		}
		else{
		errT[i] = sqrt((prog2T[i] - pow((progT[i]),2))/i);
        errK[i] = sqrt((prog2K[i] - pow((progK[i]),2))/i);
        errU[i] = sqrt((prog2U[i] - pow((progU[i]),2))/i);
        errE[i] = sqrt((prog2E[i] - pow((progE[i]),2))/i);
        errP[i] = sqrt((prog2P[i] - pow((progP[i]),2))/i);
		}
	}

    for(int i=0; i<nblocks; i++){
     meanT<<progT[i]<<endl;
     sigmaT<<errT[i]<<endl;
     meanK<<progK[i]<<endl;
     sigmaK<<errK[i]<<endl;
     meanU<<progU[i]<<endl;
     sigmaU<<errU[i]<<endl;
     meanE<<progE[i]<<endl;
     sigmaE<<errE[i]<<endl;
     meanP<<progP[i]<<endl;
     sigmaP<<errP[i]<<endl;
	}

    meanT.close();
    sigmaT.close();
    meanK.close();
    sigmaK.close();
    meanU.close();
    sigmaU.close();
    meanE.close();
    sigmaE.close();
    meanP.close();
    sigmaP.close();

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf,ReadOld;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  ReadOld.open("old.0");
  if(ReadOld.is_open()){
     //Read old configuration                                          Legge la configurazione al tempo -dt
    cout << "Read old configuration from file old.0 " << endl << endl;

    for (int i=0; i<npart; ++i){
        ReadOld >> xold[i] >> yold[i] >> zold[i];
        xold[i] = xold[i] * box;
        yold[i] = yold[i] * box;
        zold[i] = zold[i] * box;
    }
    ReadOld.close();

    Move();                                   //step di algoritmo di Verlet per calcolare r(t+dt) (inserito in x,y,z) e pone in x,y,z old la configurazione iniziale (di config.0);
                                            //calcola anche le velocità al tempo t, che utilizzo per calcolare le velocità al tempo t + dt/2: v(t+dt/2)=v(t)+(dv/dt)*(dt/2) al primo ordine di Taylor
    double vxnew[npart], vynew[npart], vznew[npart];
    double sumv[3] = {0.0, 0.0, 0.0};

    for (int i=0; i<npart; ++i){                                  //calcolo v(t + dt/2) al primo ordine di Taylor
        ReadOld >> xold[i] >> yold[i] >> zold[i];
        vxnew[i] = vx[i] + ((x[i] - xold[i])/delta);
        vynew[i] = vy[i] + ((y[i] - yold[i])/delta);
        vznew[i] = vz[i] + ((z[i] - zold[i])/delta);
    
        sumv[0] += vxnew[i];
        sumv[1] += vynew[i];
        sumv[2] += vznew[i];
    }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;                //calcola il fattore di scala
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
        vxnew[i] = vxnew[i] - sumv[0];
        vynew[i] = vynew[i] - sumv[1];
        vznew[i] = vznew[i] - sumv[2];

        sumv2 += vxnew[i]*vxnew[i] + vynew[i]*vynew[i] + vznew[i]*vznew[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vxnew[i] *= fs;                                   //imposta le velocità al tempo t + dt/2 secondo il fattore di scala
     vynew[i] *= fs;
     vznew[i] *= fs;

     xold[i] = Pbc(x[i] - vxnew[i] * delta);           //calcola x(t) = x(t+dt) - v(t+dt/2)*dt 
     yold[i] = Pbc(y[i] - vynew[i] * delta);
     zold[i] = Pbc(z[i] - vznew[i] * delta);
   }
 }

 else{
//Prepare initial velocities
   cout << "File Old.0 does not exist " << endl << endl;
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
 }
 return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pres.open("output_pres.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;                      //accumulatore per il calcolo della pressione

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       pij = 48*((1/pow(dr,12)) - (1/2)*(1/pow(dr,6)));

//Potential energy
       v += vij;
//Pressure
       p += pij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_pres = rho*stima_temp + p/(3*vol);      //pressure of the system

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

//void ConfXYZ(int nconf){ //Write configuration in .xyz format
  //ofstream WriteXYZ;

  //WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  //WriteXYZ << npart << endl;
  //WriteXYZ << "This is only a comment!" << endl;
  //for (int i=0; i<npart; ++i){
    //WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  //}
  //WriteXYZ.close();
//}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

//funzione per salvare le penultime coordinate, in modo da ripartire con la simulazione senza doverle ricreare. Il file deve essere letto dalla funzione Input

void ConfOld(void){ //Write old configuration to restart the system without creating the coordinates at time -dt
  ofstream WriteConfOld;

  cout << "Print old configuration to restart the simulation without creating the old coordinates for the Verlet algorithm to file old.final " << endl << endl;
  WriteConfOld.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConfOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConfOld.close();
  return;
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
