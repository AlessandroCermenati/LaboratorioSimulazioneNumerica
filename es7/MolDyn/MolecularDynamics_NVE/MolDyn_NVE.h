/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
int nblocks = 100;                         //numero di blocchi per il data-blocking
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
double AVT[1000], AV2T[1000], progT[1000], prog2T[1000], errT[1000];       //vettori per il data blocking
double AVK[1000], AV2K[1000], progK[1000], prog2K[1000], errK[1000];       //vettori per il data blocking
double AVU[1000], AV2U[1000], progU[1000], prog2U[1000], errU[1000];       //vettori per il data blocking
double AVE[1000], AV2E[1000], progE[1000], prog2E[1000], errE[1000];       //vettori per il data blocking
double AVP[1000], AV2P[1000], progP[1000], prog2P[1000], errP[1000];       //vettori per il data blocking
double radialf[100];
double binl;
double AVG[100][1000], AV2G[100][1000], progG[100][1000], prog2G[100][1000], errG[100][1000];


// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfOld(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
