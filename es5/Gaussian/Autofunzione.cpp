#include "Autofunzione.h"
#include <math.h>
#include <iostream>

Autofunzione :: Autofunzione(){}

Autofunzione :: ~Autofunzione(){}

double Autofunzione :: n1l0m0(double x, double y, double z) {			//valuta l'autofunzione dell'orbitale 1s

	return exp(-2.*sqrt(x*x+y*y+z*z));
}

double Autofunzione :: n2l1m0(double x, double y, double z) {			//valuta l'autofunzione dell'orbitale 2p

	return z*z*exp(-sqrt(x*x+y*y+z*z));
}