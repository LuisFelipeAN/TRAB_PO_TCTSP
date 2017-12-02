#include "PesquisaOperacional.h"
#include "SimplexSolver.h"
#include "Eigen/Dense"
#include <iostream>

///estrutura para encadear os clusters presentes na solucao
typedef struct Cluster{ ///estrutura para fragmentar a solucao em clusters com vertices de entrada e saida do cluster
    No*inicio;
    No*fim;
    Cluster* proximo;
    Cluster* anterior;///duplamente encadeada para efetuar a troca tanto na solucao como tambem na lista
}Cluster;///utulizada nas funcoes BuscaLocal2 e BuscaLocal4

typedef struct ClusterLista{
    Cluster *cluster;
    ClusterLista *proximo;
}ClusterLista;

typedef struct X{
    int idX;
    int idCluster;
    Cluster *cluster;
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

static Yr* primeiroYr=NULL;
static X * primeiroX=NULL;
static double** matriz;

static ClusterLista *primeiroCluster = NULL;


static int contaX=0;
static FILE *arqClusters;
static FILE *arqYr;
static int contaYr=0;
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
static bool buscaTabu(Cluster *atual,int idTabu){
    Cluster *c =atual;
    No *p=c->inicio;
     /*fprintf(stdout,"buscando: %d\n",idTabu);///imprime para confirmar a busca
     while(p!=c->fim){
        fprintf(stdout,"%d ",p->vertice->getIndiceTabu());
        p=p->proximo;
     }
     fprintf(stdout,"%d ",p->vertice->getIndiceTabu());
     fprintf(stdout,"\n");*/

     p=atual->inicio;
     while(p!=atual->fim){
        if(p->vertice->getIndiceTabu()==idTabu){
            return true;
        }
        p=p->proximo;
     }
     if(p->vertice->getIndiceTabu()==idTabu) return true;
     return false;
}

void inicializaLeitura(char* nomeArquivoClusters,char * nomeArquivoYr){
    arqYr = fopen(nomeArquivoYr, "r");
    if(!arqYr){
        fprintf(stderr,"ERRO: arqivo não encontrado %s\n",nomeArquivoYr);
    }
    if(!arqYr){
        exit(1);
    }
    int leitura=0;
    int idCluster;
    int vEntrada;
    int vSaida;
    double custo;
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
        contaYr++;
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
}

static bool buscaEntrada(X* xAtual, Yr* yAtual){
    ArestaInterC* arestaAtual = yAtual->primeira;
    while(arestaAtual){
        if(arestaAtual->vEntrada == xAtual->j){
            return true;
        }
        arestaAtual = arestaAtual->proxima;
    }
    return false;
}

static bool buscaSaida(X* xAtual, Yr* yAtual){
    ArestaInterC* arestaAtual = yAtual->primeira;
    while(arestaAtual){
        if(arestaAtual->vSaida == xAtual->i){
            return true;
        }
        arestaAtual = arestaAtual->proxima;
    }
    return false;
}

void emitirSistemaLinear(char* nomeArquivoSl){
    FILE * arq;
    arq = fopen(nomeArquivoSl,"w");
    X* p= primeiroX;

    ///Dimensao
    int num_restricoes= 2*contaX + getNumTotalClusters() + 2+ getNumTotalTabus();
    int num_variaveis=  contaX + contaYr + 1;

    fprintf(arq, "%d %d\n", num_restricoes, num_variaveis);

    matriz = new double *[num_restricoes];
    for(int i=0;i<num_restricoes;i++){
        matriz[i]= new double[num_variaveis];
    }
    int i=0,j=0;

    ///Funcao objetivo
    Yr * pr = primeiroYr;
    while(pr!=NULL){
        matriz[i][j] = pr->primeira->custo;
        j++;
        fprintf(arq,"%.0lf ",pr->primeira->custo);
        pr= pr->proximo;
    }

    int cAtual=1;
    X * x = primeiroX;
    while(x){
        matriz[i][j] = x->custo;
        j++;
        fprintf(arq,"%.0lf ",x->custo);
        x=x->proximo;

    }
    j=0;
    i++;
    fprintf(arq,"\n\n");

    VectorXd funcaoObjetivo(num_variaveis-1);
    for(int k=0; k<num_variaveis-1; k++){
        funcaoObjetivo(k) = matriz[0][k];
    }

    //std::cout << funcaoObjetivo << std::endl;


    int numClusters = getNumTotalClusters();
    /// Primeira Restriçao
    while(cAtual<=numClusters){
        pr = primeiroYr;
        while(pr!=NULL){
            matriz[i][j] = 0;
            j++;
            fprintf(arq," 0\t");
            pr=pr-> proximo;
        }
        p = primeiroX;
        while(p!=NULL){
            if(p->idCluster==cAtual){
                matriz[i][j] = 1;
                j++;
                fprintf(arq,"1 ");
            }else{
                matriz[i][j] = 0;
                j++;
                fprintf(arq,"0 ");
            }
            p=p->proximo;
        }
        matriz[i][j] = 1;
        j++;
        fprintf(arq,"1\t");
        j=0;
        i++;
        fprintf(arq,"\n");
        cAtual++;
    }

    /// Segunda Restricao
    fprintf(arq,"\n");
    pr = primeiroYr;
    while(pr!=NULL){
        matriz[i][j] = 1;
        j++;
        fprintf(arq," 1\t");
        pr=pr-> proximo;
    }
    p = primeiroX;
    while(p!=NULL){
        matriz[i][j] = 0;
        j++;
        fprintf(arq,"0 ");
        p=p->proximo;
    }
    matriz[i][j] = 1;
    j++;
    fprintf(arq,"1\t");
    j=0;
    i++;

    /// Terceira Restricao
    fprintf(arq,"\n\n");
    x = primeiroX;
    while(x){
        Yr* y = primeiroYr;
        while(y){
            if(buscaEntrada(x, y)){
                matriz[i][j] = -1;
                j++;
                fprintf(arq, "-1\t");
            }
            else{
                matriz[i][j] = 0;
                j++;
                fprintf(arq, " 0\t");
            }
            y = y->proximo;
        }
        X* xAux = primeiroX;
        while(xAux){
            if(xAux == x){
                matriz[i][j] = 1;
                j++;
                fprintf(arq, "1 ");
            }
            else{
                matriz[i][j] = 0;
                j++;
                fprintf(arq, "0 ");
            }
            xAux = xAux->proximo;
        }
        matriz[i][j] = 0;
        j++;
        fprintf(arq, "0\t");
        fprintf(arq,"\n");
        j=0;
        i++;
        x = x->proximo;
    }
    fprintf(arq,"\n\n");
    x = primeiroX;
    while(x){
        Yr* y = primeiroYr;
        while(y){
            if(buscaSaida(x, y)){
                matriz[i][j] = -1;
                j++;
                fprintf(arq, "-1\t");
            }
            else{
                matriz[i][j] = 0;
                j++;
                fprintf(arq, " 0\t");
            }
            y = y->proximo;
        }
        X* xAux = primeiroX;
        while(xAux){
            if(xAux == x){
                matriz[i][j] = 1;
                j++;
                fprintf(arq, "1 ");
            }
            else{
                matriz[i][j] = 0;
                j++;
                fprintf(arq, "0 ");
            }
            xAux = xAux->proximo;
        }
        matriz[i][j] = 0;
        j++;
        fprintf(arq, "0\t");
        j=0;
        i++;
        fprintf(arq,"\n");
        x = x->proximo;
    }

    /// Quarta Restricao
    fprintf(arq,"\n");
    for(int idTabu=0;idTabu<getNumTotalTabus();idTabu++){
        pr = primeiroYr;
        while(pr!=NULL){
            fprintf(arq,"0\t");
            matriz[i][j] = 0;
            j++;
            pr=pr-> proximo;
        }
        x=primeiroX;
        while(x){
            if(buscaTabu(x->cluster,idTabu)){
                fprintf(arq,"1 ");
                 matriz[i][j] = 1;
                j++;
            }else{
                fprintf(arq,"0 ");
                 matriz[i][j] = 0;
                j++;
            }
            x=x->proximo;
        }
        fprintf(arq,"1\n");
         matriz[i][j] = 1;
        i++;
        j=0;
    }
    MatrixXd restricoes(num_restricoes-1, num_variaveis);
    for(int i=1; i<num_restricoes; i++){
        for(int j=0; j<num_variaveis; j++){
            restricoes(i-1, j) = matriz[i][j];
        }
    }

    FILE * arqMat = fopen("matriz.txt","w");
    fprintf(arqMat, "%d %d\n", num_variaveis, num_restricoes);
    for(int i=0; i<num_restricoes; i++){
        for(int j=0; j<num_variaveis; j++){
            fprintf(arqMat, "%.0f ", matriz[i][j]);
        }
        fprintf(arqMat, "\n");
    }



/*

    for(int i=0;i<num_variaveis;i++){
       delete matriz[i];
    }
    delete [] matriz;*/
    fclose(arq);

    //std::cout << funcaoObjetivo << std::endl;
    SimplexSolver *simplex = new SimplexSolver(SIMPLEX_MINIMIZE, funcaoObjetivo, restricoes);

    std::cout << funcaoObjetivo.transpose() << std::endl << simplex->getSolution().transpose() << std::endl << "CUSTO: "<<simplex->getOptimum() << std::endl;
}
static void verificaTabus(){
    for(int i =1;i<=getNumTotalClusters();i++){
        X* x =primeiroX;
        while(x){
            if(x->idCluster==i){
                No*inicio = x->cluster->inicio;
                fprintf(stdout,"%d ,",inicio->vertice->getIndiceTabu());
                while(inicio!= x->cluster->fim){
                    fprintf(stdout,"%d ,",inicio->vertice->getIndiceTabu());
                    inicio=inicio->proximo;
                }
            }
        }
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

void salvarSolucaoArquivosPO(No* s){
    ///faz uma copia da solução
    No *anterior = NULL;
    No *solucao=NULL;
    No *copia = s;
    while(copia->proximo!=NULL){
        copia->anterior=anterior;
        anterior=copia;
        copia = copia->proximo;
    }
    copia->anterior=anterior;
    fprintf(stdout,"ultimo atualizado: %d\n",copia->vertice->getIDVertice());
    No*proximo=NULL;
    while(copia){
        if(solucao==NULL){
            solucao= new No();
            solucao->tabu=copia->tabu;
            solucao->vertice=copia->vertice;
            solucao->cluster=copia->cluster;
            solucao->proximo=proximo;
            proximo=solucao;
            fprintf(stdout,"copiado: %d\n",solucao->vertice->getIDVertice());
        }else{
            No *novo;
            novo= new No();
            novo->tabu=copia->tabu;
            novo->vertice=copia->vertice;
            novo->cluster=copia->cluster;
            novo->proximo=proximo;
            solucao->anterior=novo;
            solucao=novo;
            proximo=novo;
        }
        copia=copia->anterior;
    }

    Cluster* clusterFinal = NULL;///ultimo cluster da solucao
    int controle=1;
    int idClusterAtual=-1;

    No* ultimo=solucao;
    while(ultimo->proximo!=NULL){///encontra o ultimo No da solucao;
        ultimo=ultimo->proximo;
    }
    No* aux = solucao;
    anterior=solucao;

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
         if(primeiroX==NULL){
            contaX=1;
            primeiroX= new X();
            primeiroX->i = vEntrada->getIDVertice();
            primeiroX->j = vSaida->getIDVertice();
            primeiroX->idCluster = vEntrada->getIndiceCluster()+1;
            primeiroX->cluster = c;
            primeiroX->proximo=NULL;
            primeiroX->custo= calculaCustoIntraCluster(c);
        }else{
            X* aux;
            contaX++;
            aux= new X();
            aux->i = vEntrada->getIDVertice();
            aux->j = vSaida->getIDVertice();
            aux->idCluster = vEntrada->getIndiceCluster()+1;
            aux->cluster = c;
            aux->custo=calculaCustoIntraCluster(c);
            aux->proximo=primeiroX;
            primeiroX=aux;
        }
        c=c->proximo;

    }
    fprintf(stdout,"\n");
   /* X *x=primeiroX;
    while(x!=NULL){
        fprintf(stdout,"x :%d\t %d\t %d\t %lf\n",x->idCluster,x->i,x->j,x->custo);
        x=x->proximo;
    }*/

    c=clusterInicial->proximo;
    ///Faz a lista de clusters ficar circular
    clusterFinal->proximo=clusterInicial;
    clusterInicial->anterior=clusterFinal;
    fprintf(arqYr,"%d\t%d\n",contaYr,getNumTotalClusters());///numero de arestas inter cluster
     while(c!=clusterInicial){
        Vertice *vEntrada =  c->inicio->vertice;
        Vertice *vSaida = c->anterior->fim->vertice;
        fprintf(arqYr,"%d\t %d\t %lf\n",vEntrada->getIDVertice(),vSaida->getIDVertice(),vEntrada->calculaCusto(vSaida));
        c=c->proximo;
    }

    ///salva a ultima aresta inter cluster ligando o final ao inicial
    Vertice *vEntrada =  c->inicio->vertice;
    Vertice *vSaida = c->anterior->fim->vertice;
    if(vEntrada->getIndiceCluster()!=vSaida->getIndiceCluster()){
        fprintf(arqYr,"%d\t %d\t %lf\n",vEntrada->getIDVertice(),vSaida->getIDVertice(),vEntrada->calculaCusto(vSaida)*10*getPenalizacao());
    }

   c = clusterInicial; ///corta a lista circular
   salvarSolucao(solucao);///inprime para verificar os tabus
   while(c!=clusterFinal){
        No *p=c->inicio;
         while(p!=c->fim){
            fprintf(stdout,"%d ",p->vertice->getIndiceTabu()+1);
            p=p->proximo;
         }
         fprintf(stdout,"%d ",p->vertice->getIndiceTabu()+1);
         fprintf(stdout,"\n");
         c=c->proximo;
    }
    No *p=clusterFinal->inicio;
     while(p!=clusterFinal->fim){
        fprintf(stdout,"%d ",p->vertice->getIndiceTabu()+1);
        p=p->proximo;
     }
     fprintf(stdout,"%d ",p->vertice->getIndiceTabu()+1);
     fprintf(stdout,"\n");

    if(primeiroCluster==NULL){
        primeiroCluster = new ClusterLista();
        primeiroCluster->proximo=NULL;
        primeiroCluster->cluster = clusterInicial;
    }else{
        ClusterLista *aux =  new ClusterLista();
        aux->proximo=primeiroCluster;
        aux->cluster=clusterInicial;
        primeiroCluster=aux;
    }
}

void imprimeX(){
    /*fprintf(stdout,"ImprimindoX");
    X *x=primeiroX;
    while(x!=NULL){
        fprintf(stdout,"x :%d\t %d\t %d\t %lf\n",x->idCluster,x->i,x->j,x->custo);
        x=x->proximo;
    }*/
}
