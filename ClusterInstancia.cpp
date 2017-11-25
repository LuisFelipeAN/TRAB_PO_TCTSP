#include "ClusterInstancia.h"
#include <cstdlib>
ClusterInstancia::ClusterInstancia()
{
    //ctor
    primeiroTabu=NULL;
    numTabus = 0;
    visitado=false;
    numPresenteSolucao=0;
    numVertices=0;
}
void ClusterInstancia::adicionaTabu(Tabu *t){
    NoCluster *novo = new NoCluster();
    novo->tabu=t;
    novo->proximo=NULL;

    if(primeiroTabu==NULL){
        numTabus++;
        primeiroTabu=novo;
    }
    else{
        bool jaAdicionado=false;
        for(NoCluster *p=primeiroTabu;p!=NULL;p=p->proximo){
            if(p->tabu==t){
                delete novo;
                jaAdicionado=true;
                break;
            }
        }
        if(!jaAdicionado){
            numTabus++;
            novo->proximo=primeiroTabu;
            primeiroTabu = novo;
        }

    }

}

int ClusterInstancia::getDifVisitados_Total(){
    int numVisitados=0;
    for(NoCluster *p=primeiroTabu;p!=NULL;p=p->proximo){
            if(p->tabu->visitado){
                numVisitados++;
            }
    }
    return numTabus-numVisitados;
}
Tabu*ClusterInstancia::getTabuRandom(){
     int cont=0;
     for(NoCluster *p=primeiroTabu;p!=NULL;p=p->proximo){
            if(!p->tabu->visitado){
                cont++;
            }
    }
    int tRand;
    if(cont!=0)
        tRand = rand() % cont;
    else
        tRand=0;

    cont = 0;
    for(NoCluster *p=primeiroTabu;p!=NULL;p=p->proximo){
            if(!p->tabu->visitado){
                if(cont==tRand){
                    return p->tabu;
                }
                cont++;
            }
    }
    return NULL;
}
ClusterInstancia::~ClusterInstancia()
{
    if(primeiroTabu!=NULL){
        NoCluster *p=primeiroTabu->proximo;
         delete primeiroTabu;
         while(p=NULL){
                primeiroTabu=p;
                p=p->proximo;
                delete primeiroTabu;
        }
        primeiroTabu=NULL;
    }

}
