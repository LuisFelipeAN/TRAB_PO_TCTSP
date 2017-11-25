#ifndef CLUSTERINSTANCIA_H
#define CLUSTERINSTANCIA_H
#include "Tabu.h"
typedef struct NoCluster{
    Tabu* tabu;
    NoCluster* proximo;
}NoCluster;

class ClusterInstancia
{
    public:
       int numTabus;
       int numVertices;
       bool visitado;
       int numPresenteSolucao;
       void adicionaTabu(Tabu* t);
       ClusterInstancia();
       virtual ~ClusterInstancia();
       Tabu* getTabuRandom();
       int getDifVisitados_Total();
       void imprimeSeusTabus(){

       }

    protected:

    private:
        bool tabuAdicionado(int id);///verifica se o id do tabu ja foi adicionado ao cluster
        NoCluster *primeiroTabu;

};

#endif // CLUSTERINSTANCIA_H
