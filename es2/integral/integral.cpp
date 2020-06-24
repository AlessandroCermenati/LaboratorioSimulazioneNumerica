#include "funzionebase.h"
#include "random.h"
#include "integral.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

integral :: integral(){
    _a=0;
    _b=1;
}

integral :: ~integral(){}

integral :: integral(double inf, double sup){
    _a=inf;
    _b=sup;
}

void integral :: setdomain(double inf, double sup){
    _a=inf;
    _b=sup;
}

double integral :: integration_MeanValue(funzionebase *F, Random *rnd, int Nhit){          //calcolo dell'integrale per campionamento uniforme della variabile in (a,b)
    double sum=0;
    double x;
    for(int i=0; i<Nhit; i++){
        x = _a + (_b - _a)*rnd->Rannyu();                                                  //campiona x uniformemente in (a,b)
        sum = sum + F->Eval(x);                                                            //calcola il valor medio della funzione
    }
    double I = ((_b - _a)/Nhit)*sum;
    return I;
}

double integral :: integration_importance(funzionebase *F, Random *rnd, int Nhit){         //calcolo dell'integrale con importance sampling della variabile in (a,b)
    double sum=0;
    double x;
    for(int i=0; i<Nhit; i++){
        x = rnd->importance_sampling(F, _a, _b);                                           //col metodo importance_sampling, x è campionata secondo la pdf scelta
        sum = sum + F->importance_func(x);
    }
    double I = ((_b - _a)/Nhit)*sum;
    return I;
}