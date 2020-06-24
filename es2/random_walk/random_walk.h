#ifndef _Random_Walk_h
#define _Random_Walk_h

#include "random.h"

class Random_Walk {
    private:
        int _N;                                                 //fisso il numero di step di ogni random walk
        double origin[3];                                       //coordinate della particella che compie random walk
    public:
        Random_Walk();                                          //costruttore
        ~Random_Walk();                                         //distruttore
        Random_Walk(int N);                                     //costruttore con numero di passi del RW
        void setsteps(int N);                                   //imposta il numero di step del RW
        void setorigin(double x, double y, double z);           //imposta l'origine
        void nullorigin();                                      //imposta l'origine al punto (0,0,0)
        double GetX();                                          //metodi per accedere alle coordinate
        double GetY();
        double GetZ();
        void RW_lattice(Random *rnd, double a);                 //esegue un RW in un lattice 3D, con _N passi, partendo dalle coord. in origin, con step a in una delle direzioni
        void RW_continuum(Random *rnd, double r);               //esegue un RW nello spazio, con _N passi, partendo dalle coord. in origin, con step di modulo r in una direzione casuale nell'angolo solido
};

#endif 