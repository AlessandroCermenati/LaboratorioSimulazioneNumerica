#ifndef funzione_base_h_
#define funzione_base_h_

class funzionebase {
    public:
        virtual double Eval(double x) const =0;                 //valuta la funzione
        virtual double Max() const =0;                          //valuta il massimo della funzione per il metodo di reiezione
        virtual double importance_pdf(double x) const=0;        //valuta la pdf per l'importance sampling
        virtual double importance_func(double x) const=0;       //valuta la funzione riscritta per l'importance sampling
        virtual double Max_pdf() const =0;                      //valuta il massimo della pdf per l'importance sampling da campionare col metodo di reiezione
};
#endif