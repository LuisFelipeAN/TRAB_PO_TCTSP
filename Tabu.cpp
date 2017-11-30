#include "Tabu.h"
using namespace std;
Tabu::Tabu()
{
    primeiroNo=NULL;
    ultimoNo  = NULL;
    visitado=false;
    NumTrocasVertices=0;
    numVertices=0;
}
///insere os Nos ordenados na lista encadeada
void Tabu::insereVertice(Vertice* v){
    NoTabu *novo = new NoTabu();
    novo->vertice=v;
    novo->proximo=NULL;
    numVertices++;
    if(primeiroNo==NULL){
        primeiroNo=novo;
        ultimoNo=novo;
    }
    else{
        ultimoNo->proximo=novo;
        ultimoNo = novo;
    }
}
///retorna o primeiro vertice diferente do vertice passado
Vertice* Tabu::outroVertice(Vertice *v){
    Vertice* retorno=primeiroNo->vertice;
    NoTabu* percorre;
     NumTrocasVertices++;
    if(primeiroNo->vertice==v&&primeiroNo->proximo!=NULL){///retorna o segundo vertice da lista encadeada

        retorno=primeiroNo->vertice;

        ultimoNo->proximo=primeiroNo;

        ultimoNo=primeiroNo;

        primeiroNo=primeiroNo->proximo;

        ultimoNo->proximo=NULL;




    }else if(primeiroNo->proximo!=NULL&&ultimoNo->vertice!=v){
        retorno=primeiroNo->vertice;
        percorre = primeiroNo;

        while(percorre->proximo!=NULL&&percorre->proximo->vertice!=v){///percorre ate que o proximo tenha o vertice procurado
            percorre=percorre->proximo;
        }
        if(percorre->proximo!=NULL&&percorre->proximo->vertice==v){
            ///insere o no do vertice procurado no final da lista encadeada

            ultimoNo->proximo=percorre->proximo;
            ultimoNo=percorre->proximo;
            percorre->proximo=percorre->proximo->proximo;
            ultimoNo->proximo=NULL;

        }

     }
    return retorno;
}

void Tabu::imprimeVertices(FILE* arquivo){
    for(NoTabu *aux = primeiroNo;aux!=NULL;aux=aux->proximo){
            fprintf(arquivo," %d",aux->vertice->getIDVertice());
    }

}
///retorna um vertice randomico do tabu
Vertice* Tabu::getRandom(){

    int vRand = rand() % numVertices;
    NoTabu* aux=primeiroNo;
    for(int i=0;i<vRand;i++){
        aux=aux->proximo;
    }
    return aux->vertice;
}

Vertice * Tabu::getRamdomCluster(int IdCluster){
    int cont =0;
     for(NoTabu *aux = primeiroNo;aux!=NULL;aux=aux->proximo){
           if(aux->vertice->getIndiceCluster()==IdCluster){
                cont++;
           }
    }
   int vRand;
    if(cont!=0)
        vRand = rand() % cont;
    else
        vRand=0;
    cont=0;
    for(NoTabu *aux = primeiroNo;aux!=NULL;aux=aux->proximo){
           if(aux->vertice->getIndiceCluster()==IdCluster){
                if(vRand==cont){
                    return aux->vertice;
                }
                cont++;
           }
    }
    return NULL;

}
bool Tabu::efetuouTodasTrocas(){
    if(NumTrocasVertices<numVertices){
        return false;
    }else{
        NumTrocasVertices=0;
        return true;
    }
}
Tabu::~Tabu()
{
    NoTabu* aux;
    aux=primeiroNo->proximo;
    delete primeiroNo;

    while(aux!=NULL){
        primeiroNo = aux;
        aux = aux->proximo;
        delete primeiroNo;
    }
    primeiroNo=NULL;
    ultimoNo =NULL;
}
int Tabu::getNumVertices(){return numVertices;}
