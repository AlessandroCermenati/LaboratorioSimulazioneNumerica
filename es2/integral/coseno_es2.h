#ifndef coseno_es2_h_
#define coseno_es2_h_
#include "funzionebase.h"

class coseno_es2 : public funzionebase {
	private:
		double _A;												//parametro della funzione
	public:
		coseno_es2();											//costruttore
		~coseno_es2();											//distruttore
		coseno_es2(double x);									//costruttore che imposta il parametro
		void setparameter(double x);							//metodo per impostare il parametro
		double getparameter();									//metodo per accedere al parametro
		virtual double Eval(double x) const;					//metodi virtuali ereditati da funzionebase, chiamati nei metodi di integral e Random
		virtual double Max() const;
		virtual double importance_pdf(double x) const;
        virtual double importance_func(double x) const;
		virtual double Max_pdf() const;
};
#endif