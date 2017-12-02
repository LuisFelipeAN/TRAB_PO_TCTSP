#include "Vertice.h"

Vertice::Vertice(int indice, double cordX, double cordY,int indiceCluster,int indiceTabu)
{
    id=indice;
    x=cordX;
    y=cordY;
    idTabu=indiceTabu;
    idCluster = indiceCluster;
    prox=NULL;
}
double Vertice::calculaCusto(Vertice *v){
       return  sqrt( (v->getCordX()-x)*(v->getCordX()-x)+(v->getCordY()-y)*(v->getCordY()-y));
};
Vertice::~Vertice()
{
     prox=NULL;
}





