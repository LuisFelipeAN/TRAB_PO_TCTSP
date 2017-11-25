#ifndef VERTICE_H
#define VERTICE_H
#include <stdio.h>
#include <math.h>
#include <cstdlib>
class Vertice
{
    private:
        int id;
        int idTabu;
        int idCluster;
        Vertice * prox;
        double x;
        double y;

    public:
        Vertice();
        Vertice(int indice, double cordX, double cordY,int indiceCluster,int indiceTabu);
        void setProximo(Vertice *v){prox=v;};
        Vertice* getProximo(){return prox;};

        int getIndiceTabu(){return idTabu;}
        int getIndiceCluster(){return idCluster;}
        double calculaCusto(Vertice *v);
        int getIDVertice(){return id;}
        double getCordX(){return x;}
        double getCordY(){return y;}
        virtual ~Vertice();

    protected:

};

#endif // VERTICE_H
