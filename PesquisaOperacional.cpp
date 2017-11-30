#include "PesquisaOperacional.h"
///estrutura para encadear os clusters presentes na solucao
typedef struct Cluster{ ///estrutura para fragmentar a solucao em clusters com vertices de entrada e saida do cluster
    No*inicio;
    No*fim;
    Cluster* proximo;
    Cluster* anterior;///duplamente encadeada para efetuar a troca tanto na solucao como tambem na lista
}Cluster;///utulizada nas funcoes BuscaLocal2 e BuscaLocal4

typedef struct X{
    int idX;
    int idCluster;
    int i;
    int j;
    double custo;
    int emUso;
    X* proximo;
}X;

typedef struct ArestaInterC{
    int vSaida;
    int vEntrada;
    double custo;
    ArestaInterC *proxima;
}ArestaInterC;

typedef struct Yr{
    ArestaInterC *primeira;
    Yr* proximo;
}Yr;

static Yr* primeiroYr;
static X * primeiroX;

static int contaX=0;
static FILE *arqClusters;
static FILE *arqYr;
static int contYr=0;
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

void inicializaLeitura(char* nomeArquivoClusters,char * nomeArquivoYr){
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
    int leitura=0;
    primeiroX = new X();
    primeiroX->proximo=NULL;

    int idCluster;
    int vEntrada;
    int vSaida;
    double custo;
    while(leitura!=-1){
        X* aux= new X();
        leitura=fscanf(arqClusters,"%d\t%d\t%d\t%lf\n",&idCluster,&vEntrada,&vSaida,&custo);
        aux->idCluster= idCluster;
        aux->i=vEntrada;
        aux->j=vSaida;
        aux->custo=custo;
        aux->emUso=0;
        aux->proximo=primeiroX->proximo;
        aux->idX=contaX;
        contaX++;
        primeiroX->proximo=aux;
    };
    X *pri = primeiroX;
    primeiroX=primeiroX->proximo;
    delete pri;
    /*X * p =primeiroX;
    while(p!=NULL){
        fprintf(stdout,"%d\t%d\t%d\t%lf\n",p->idCluster,p->i,p->j,p->custo);
        p=p->proximo;
    }*/


    leitura=0;
    int cont;
    int numArestas;

    primeiroYr = new Yr();
    primeiroYr->proximo=NULL;
    while(leitura!=-1){
       leitura=fscanf(arqYr,"%d\t%d\n",&cont,&numArestas);
       //fprintf(stdout,"%d\t%d\n",cont,numArestas);
       if(leitura == -1) break;
       Yr * yi = new Yr();
        yi->proximo=primeiroYr->proximo;
        primeiroYr->proximo=yi;
        ArestaInterC* primeira = new ArestaInterC();
        fscanf(arqYr,"%d\t%d\t%lf\n",&vSaida,&vEntrada,&custo);
       // fprintf(stdout,"%d\t%d\t%lf\n",vSaida,vEntrada,custo);
        primeira->vSaida=vSaida;
        primeira->vEntrada=vEntrada;
        primeira->custo=custo;
        primeira->proxima=NULL;
        for(int i=0;i<(numArestas-1);i++){
           ArestaInterC * aux =  new ArestaInterC();
            leitura=fscanf(arqYr,"%d\t%d\t%lf\n",&vSaida,&vEntrada,&custo);
           // fprintf(stdout,"%d\t%d\t%lf\n",vSaida,vEntrada,custo);
            aux->vSaida=vSaida;
            aux->vEntrada=vEntrada;
            aux->custo=custo;
            aux->proxima=primeira;
            primeira=aux;
        }
        yi->primeira=primeira;
    }

    Yr * prj = primeiroYr;
    primeiroYr = primeiroYr->proximo;
    delete prj;

   /* Yr *perc=primeiroYr;
    while(perc!=NULL){
        ArestaInterC *a = perc->primeira;
        while(a!=NULL){
            fprintf(stdout,"%d\t%d\t%lf\n",a->vSaida,a->vEntrada,a->custo);
            a=a->proxima;
        }
        fprintf(stdout,"\n");
        perc=perc->proximo;
    }*/
    fclose(arqYr);
    fclose(arqClusters);
}
void emitirSistemaLinear(char* nomeArquivoSl){
    FILE * arq;
    arq = fopen(nomeArquivoSl,"w");
    X* p= primeiroX;
    ///Funcao objetivo
    fprintf(arq," 0.000 ");
    Yr * pr = primeiroYr;
    while(pr!=NULL){
          fprintf(arq," %lf ",p->custo);
          pr=pr-> proximo;
    }

    int cAtual=1;
    X * x = primeiroX;
    while(x!=NULL){
        fprintf(arq," %lf ",x->custo);
        x=x->proximo;

    }
    fprintf(arq,"\n\n\n");


    int numClusters = getNumTotalClusters();
    ///Primeira restriçao
    while(cAtual<=numClusters){
        fprintf(arq,"1 ");
        pr = primeiroYr;
        while(pr!=NULL){
              fprintf(arq," 0 ");
              pr=pr-> proximo;
        }
        p = primeiroX;
        while(p!=NULL){
            if(p->idCluster==cAtual){
                fprintf(arq," -1 ");
            }else{
                fprintf(arq," 0 ");
            }
            p=p->proximo;
        }
        fprintf(arq,"\n\n");
        cAtual++;
    }
    pr = primeiroYr;
    fprintf(arq," 1 ");
    while(pr!=NULL){
      fprintf(arq," -1 ");
      pr=pr-> proximo;
    }
    p = primeiroX;
    while(p!=NULL){
        fprintf(arq," 0 ");
        p=p->proximo;
    }

    fclose(arq);
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
    fprintf(arqYr,"%d\t%d\n",contYr,getNumTotalClusters());///numero de arestas inter cluster
    while(c!=clusterInicial){///salva as areatas inter cluster
        Vertice *vEntrada =  c->inicio->vertice;
        Vertice *vSaida = c->anterior->fim->vertice;
        fprintf(arqYr,"%d\t %d\t %lf\n",vEntrada->getIDVertice(),vSaida->getIDVertice(),vEntrada->calculaCusto(vSaida));
        c=c->proximo;
    }

    ///salva a ultima aresta inter cluster ligando o final ao inicial
    Vertice *vEntrada =  c->inicio->vertice;
    Vertice *vSaida = c->anterior->fim->vertice;
    if(vEntrada->getIndiceCluster()!=vSaida->getIndiceCluster()){
        fprintf(arqYr,"%d\t %d\t %lf\n",vEntrada->getIDVertice(),vSaida->getIDVertice(),vEntrada->calculaCusto(vSaida));
    }

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
    contYr++;
}
