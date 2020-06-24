#ifndef integral_h_
#define integral_h_
#include "funzionebase.h"
#include "random.h"

class integral {
	private:
		double _a, _b;																//intervallo di integrazione (a,b)
	public:
		integral();																	//costruttore di default
		~integral();																//distruttore
		integral(double inf, double sup);											//costruttore a partire dall'intervallo di integrazione
		void setdomain(double inf, double sup);										//metodo per impostare un dominio di integrazione
		double integration_MeanValue(funzionebase *F, Random *rnd, int Nhit);		//calcola l'integrale col metodo del valor medio della funzione producendo Nhit campionamenti della variabile
		double integration_importance(funzionebase *F, Random *rnd, int Nhit);		//calcola l'integrale col metodo del valor medio della funzione producendo Nhit importance sampling della variabile
};
#endif