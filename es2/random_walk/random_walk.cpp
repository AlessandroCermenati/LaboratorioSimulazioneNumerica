#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random_walk.h"
#include "random.h"

using namespace std;

Random_Walk :: Random_Walk(){										//costruttore
	_N = 1;
	origin[0] = 0.;
	origin[1] = 0.;
	origin[2] = 0.;
}

Random_Walk :: ~Random_Walk(){}										//distruttore

Random_Walk :: Random_Walk(int N){									//costruttore con numero di step
	_N = N;
	origin[0] = 0.;
	origin[1] = 0.;
	origin[2] = 0.;
}

void Random_Walk :: setsteps(int N){								//imposta il numero di step del random walk
	_N = N;
}

void Random_Walk :: setorigin(double x, double y, double z){		//imposta le coordinate della particella
	origin[0] = x;
	origin[1] = y;
	origin[2] = z;
}

void Random_Walk :: nullorigin(){									//mette la particella nell'origine
	origin[0] = 0.;
	origin[1] = 0.;
	origin[2] = 0.;
}
double Random_Walk :: GetX(){										//accede alle coordinate della particella
	return origin[0];
}

double Random_Walk :: GetY(){
	return origin[1];
}

double Random_Walk :: GetZ(){
	return origin[2];
}

void Random_Walk :: RW_lattice(Random *rnd, double a){				//esegue _N step di RW in lattice 3D di passo a; le coordinate finali sono salvate nella classe
	for(int i = 0; i < _N; i++){
		double r = rnd->Rannyu();									//per decidere la direzione dello step suddivido in 3 l'intervallo (0,1) e campiono r uniformemente in (0,1)
		double s = rnd->Rannyu();									//per decidere se lo spostamento è nel senso positivo o negativo suddivido (0,1) in due sottointervalli di uguale probabilità
		if(r<(1./3.) && s<0.5){
			origin[0] = origin[0] + a;
		}
		if(r<(1./3.) && s>=0.5){
			origin[0] = origin[0] - a;
		}
		if(r>= (1./3.) && r < (2./3.) && s<0.5){
			origin[1] = origin[1] + a;
		}
		if(r>= (1./3.) && r < (2./3.) && s>= 0.5){
			origin[1] = origin[1] - a;
		}
		if(r>= (2./3.) && s<0.5){
			origin[2] = origin[2] + a;
		}
		if(r>= (2./3.) && s>=0.5){
			origin[2] = origin[2] - a;
		}
	}
}

void Random_Walk :: RW_continuum(Random *rnd, double r){			//esegue _N step di RW nello spazio con spostamenti di modulo r in direzione qualunque
	for(int i = 0; i < _N; i++){
		double theta = rnd->Ranpolar();								//campiono uniformemente l'angolo solido per produrre la direzione dello spostamento
		double phi = rnd->Ranazim();
		origin[2] = origin[2] + r*cos(theta);
		origin[1] = origin[1] + r*sin(theta)*sin(phi);
		origin[0] = origin[0] + r*sin(theta)*cos(phi);
	}
}