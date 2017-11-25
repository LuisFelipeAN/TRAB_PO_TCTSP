#include "PesquisaOperacional.h"
///estrutura para encadear os clusters presentes na solucao
typedef struct Cluster{ ///estrutura para fragmentar a solucao em clusters com vertices de entrada e saida do cluster
    No*inicio;
    No*fim;
    Cluster* proximo;
    Cluster* anterior;///duplamente encadeada para efetuar a troca tanto na solucao como tambem na lista
}Cluster;///utulizada nas funcoes BuscaLocal2 e BuscaLocal4

static FILE *arqClusters;
static FILE *arqYr;

void inicializaArquivosEscrita(char* nomeArquivoClusters,char * nomeArquivoYr){
    arqClusters = fopen(nomeArquivoClusters, "w");
    arqYr = fopen(nomeArquivoYr, "w");
}
void finalizarArquivosEscrita(){
    if(arqClusters)
        fclose(arqClusters);
    if(arqYr)
        fclose(arqYr);
}
void inicializaArquivosLeitura(char* nomeArquivoClusters,char * nomeArquivoYr){
    arqClusters = fopen(nomeArquivoClusters, "r");
    arqYr = fopen(nomeArquivoYr, "r");
    if(!arqClusters){
        fprintf(stderr,"ERRO: arqivo não encontrado %s\n",nomeArquivoClusters);
    }
     if(!arqYr){
        fprintf(stderr,"ERRO: arqivo não encontrado %s\n",nomeArquivoYr);
    }
    if(!arqClusters|!arqYr){
        exit(1);
    }
}
static double calculaCustoIntraCluster(Cluster *c){
    No* no = c->inicio;
    if(c->inicio!=c->fim){
        double custo=0;
        while(no!=c->fim){
            custo+= no->vertice->calculaCusto(no->proximo->vertice);
            no=no->proximo;
        }
        return custo;
    }else return 0;
}

void salvarSolucaoArquivosPO(No* solucao){
    Cluster* clusterFinal = NULL;///ultimo cluster da solucao
    int controle=1;
    int idClusterAtual=-1;

    No* ultimo=solucao;
    while(ultimo->proximo!=NULL){///encontra o ultimo No da solucao;
        ultimo=ultimo->proximo;
    }
    No* aux = solucao;
    No* anterior=solucao;

    ///Cria uma lista de clusters e encadeia pelos anteriores
    while(aux!=NULL){
        if(aux->vertice->getIndiceCluster()!=idClusterAtual){
            if(controle==1){
                controle=0;
                Cluster* no= new Cluster();
                no->inicio=anterior;
                no->anterior=clusterFinal;
                clusterFinal=no;
                idClusterAtual=aux->vertice->getIndiceCluster();
                anterior=aux;
                aux=aux->proximo;

            }else{
                clusterFinal->fim=anterior;
                controle=1;
                anterior=aux;
                idClusterAtual=-1;
            }

        }else{
            anterior=aux;
            aux=aux->proximo;
        }

    }
    clusterFinal->fim=ultimo;
    Cluster *clusterInicial=clusterFinal;
    Cluster *auxcc=NULL;

    ///percorre a lista encadeada por anteriores ligando os proximos
    while(clusterInicial->anterior!=NULL){
            auxcc=clusterInicial;
            clusterInicial=clusterInicial->anterior;
            clusterInicial->proximo=auxcc;

    }
    clusterInicial->proximo=auxcc;

    Cluster *c = clusterInicial;
    while(c!=NULL){
        Vertice *vEntrada = c->inicio->vertice;
        Vertice *vSaida = c->fim->vertice;
        fprintf(arqClusters,"%d\t %d\t %d\t %lf\n",vEntrada->getIndiceCluster()+1,vEntrada->getIDVertice(),vSaida->getIDVertice(),calculaCustoIntraCluster(c));
        c=c->proximo;
    }
    c=clusterInicial->proximo;

    ///Faz a lista de clusters ficar circular
    clusterFinal->proximo=clusterInicial;
    clusterInicial->anterior=clusterFinal;

    while(c!=clusterInicial){
        Vertice *vEntrada =  c->inicio->vertice;
        Vertice *vSaida = c->anterior->fim->vertice;
        fprintf(arqYr,"%d\t %d\t %lf\n",vEntrada->getIDVertice(),vSaida->getIDVertice(),vEntrada->calculaCusto(vSaida));
        c=c->proximo;
    }
     fprintf(arqYr,"#\n");

    c = clusterInicial; ///desaloca a lista encadeada de clusters
    while(c!=clusterFinal){
        c=c->proximo;
        delete clusterInicial;
        clusterInicial=NULL;
        clusterInicial=c;
    }
    delete clusterFinal;
    clusterInicial=NULL;
    clusterFinal=NULL;

}
