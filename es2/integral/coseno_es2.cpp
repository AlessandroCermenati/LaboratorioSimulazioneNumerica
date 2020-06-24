#include "funzionebase.h"
#include "coseno_es2.h"
#include <math.h>
#include <iostream>

coseno_es2 :: coseno_es2(){
		_A = 1;
}

coseno_es2 :: ~coseno_es2(){}

coseno_es2 :: coseno_es2(double x){
		_A = x;
}

void coseno_es2 :: setparameter(double x){
		_A = x;
}

double coseno_es2 :: getparameter(){
		return _A;
}

double coseno_es2 :: Eval(double x) const{								//funzione richiesta dall'esercizio
		double y;
		y = _A*(cos(_A*x));
		return y;
}

double coseno_es2 :: Max() const{										//il massimo della funzione si ha per x=0 (cos(x)=1) e vale quanto il parametro
	return _A;
}

double coseno_es2 :: importance_pdf(double x) const{					//funzione p(x) normalizzata, usata come pdf per l'importance sampling
	return (3./2.)*(1 - pow(x,2));
}
double coseno_es2 :: importance_func(double x) const{					//funzione riscritta con la pdf per l'importance sampling
	return Eval(x)/importance_pdf(x);
}
double coseno_es2 :: Max_pdf() const{									//massimo della pdf scelta per l'importance sampling
	return 1.5;
}