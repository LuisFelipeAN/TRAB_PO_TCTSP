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
    int id;
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
    No *p=atual->inicio;
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
       if(leitura == -1) break;
       Yr * yi = new Yr();

        yi->proximo=primeiroYr->proximo;
        primeiroYr->proximo=yi;
        ArestaInterC* primeira = new ArestaInterC();
        fscanf(arqYr,"%d\t%d\t%lf\n",&vSaida,&vEntrada,&custo);
        primeira->vSaida=vSaida;
        primeira->vEntrada=vEntrada;
        primeira->custo=custo;
        primeira->proxima=NULL;
        for(int i=0;i<(numArestas-1);i++){
           ArestaInterC * aux =  new ArestaInterC();
            leitura=fscanf(arqYr,"%d\t%d\t%lf\n",&vSaida,&vEntrada,&custo);
            aux->vSaida=vSaida;
            aux->vEntrada=vEntrada;
            aux->custo=custo;
            aux->proxima=primeira;
            primeira=aux;
        }
        yi->primeira=primeira;
        contaYr++;
        yi->id=contaYr;
    }

    Yr * prj = primeiroYr;
    primeiroYr = primeiroYr->proximo;
    delete prj;
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

    fprintf(arq, "Minimize\n", num_restricoes, num_variaveis);

    matriz = new double *[num_restricoes];
    for(int i=0;i<num_restricoes;i++){
        matriz[i]= new double[num_variaveis];
    }
    int i=0,j=0;

    ///Funcao objetivo
    Yr * pr = primeiroYr;
    while(pr!=NULL){
        ArestaInterC *a=pr->primeira;
        double soma=0;
        while(a){
            soma+=a->custo;
            a=a->proxima;
        }
        matriz[i][j] = soma;
        j++;
        fprintf(arq,"+ %.0lf Y%d ",soma,pr->id);
        pr= pr->proximo;
    }

    int cAtual=1;
    X * x = primeiroX;
    while(x){
        matriz[i][j] = x->custo;
        j++;
        fprintf(arq," + %.0lf X%d ",x->custo,x->idX);
        x=x->proximo;

    }
    j=0;
    i++;
    fprintf(arq,"\nSubject To\n");

    VectorXd funcaoObjetivo(num_variaveis-1);
    for(int k=0; k<num_variaveis-1; k++){
        funcaoObjetivo(k) = matriz[0][k];
    }
    int numClusters = getNumTotalClusters();
    /// Primeira Restriçao
    while(cAtual<=numClusters){
        pr = primeiroYr;
        fprintf(arq," R%d: ",i);
        while(pr!=NULL){
            matriz[i][j] = 0;
            j++;
            pr=pr-> proximo;
        }
        p = primeiroX;
        while(p!=NULL){
            if(p->idCluster==cAtual){
                matriz[i][j] = 1;
                j++;
                fprintf(arq,"+ X%d ",p->idX);
            }else{
                matriz[i][j] = 0;
                j++;
            }
            p=p->proximo;
        }
        matriz[i][j] = 1;
        j++;
        fprintf(arq," = 1");
        j=0;
        i++;
        fprintf(arq,"\n");
        cAtual++;
    }

    /// Segunda Restricao
    pr = primeiroYr;
    fprintf(arq," R%d:",i);
    while(pr!=NULL){
        matriz[i][j] = 1;
        j++;
        fprintf(arq,"+ Y%d ",pr->id);
        pr=pr-> proximo;
    }
    p = primeiroX;
    while(p!=NULL){
        matriz[i][j] = 0;
        j++;
        p=p->proximo;
    }
    matriz[i][j] = 1;
    j++;
    fprintf(arq," = 1\n");
    j=0;
    i++;

    /// Terceira Restricao

    x = primeiroX;
    while(x){

        fprintf(arq," R%d: ",i);
        Yr* y = primeiroYr;
        while(y){
            if(buscaEntrada(x, y)){
                matriz[i][j] = -1;
                j++;
                fprintf(arq, "- Y%d ",y->id);
            }
            else{
                matriz[i][j] = 0;
                j++;
            }
            y = y->proximo;
        }
        X* xAux = primeiroX;
        while(xAux){
            if(xAux == x){
                matriz[i][j] = 1;
                j++;
                fprintf(arq, "+ X%d ",xAux->idX);
            }
            else{
                matriz[i][j] = 0;
                j++;
            }
            xAux = xAux->proximo;
        }
        matriz[i][j] = 0;
        j++;

        j=0;
        i++;
        x = x->proximo;
        fprintf(arq," = 0 \n");
    }
    x = primeiroX;
    while(x){

        fprintf(arq," R%d: ",i);
        Yr* y = primeiroYr;
        while(y){
            if(buscaSaida(x, y)){
                matriz[i][j] = -1;
                j++;
                fprintf(arq, "- Y%d ",y->id);
            }
            else{
                matriz[i][j] = 0;
                j++;
            }
            y = y->proximo;
        }
        X* xAux = primeiroX;
        while(xAux){
            if(xAux == x){
                matriz[i][j] = 1;
                j++;
                 fprintf(arq, "+ X%d ",xAux->idX);
            }
            else{
                matriz[i][j] = 0;
                j++;
            }
            xAux = xAux->proximo;
        }
        matriz[i][j] = 0;
        j++;
        j=0;
        i++;
        fprintf(arq," = 0 \n");
        x = x->proximo;
    }

    /// Quarta Restricao
    for(int idTabu=0;idTabu<getNumTotalTabus();idTabu++){
        pr = primeiroYr;
        fprintf(arq," R%d: ",i);
        while(pr!=NULL){
            matriz[i][j] = 0;
            j++;
            pr=pr-> proximo;
        }
        x=primeiroX;
        while(x){
            if(buscaTabu(x->cluster,idTabu)){
                fprintf(arq,"+ X%d ",x->idX);
                 matriz[i][j] = 1;
                j++;
            }else{
                 matriz[i][j] = 0;
                j++;
            }
            x=x->proximo;
        }
        fprintf(arq," = 1\n");
         matriz[i][j] = 1;
        i++;
        j=0;
    }    fprintf(arq,"Bounds\nBinary \n");
    pr = primeiroYr;
    while(pr!=NULL){
        fprintf(arq,"Y%d ",pr->id);
        pr= pr->proximo;
    }
    x = primeiroX;
    while(x){
        fprintf(arq,"X%d ",x->custo,x->idX);
        x=x->proximo;

    }

    fprintf(arq,"\nEnd\n");
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
    fclose(arqMat);
    fclose(arq);
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
    if(c->inicio->proximo==c->fim) return c->inicio->vertice->calculaCusto(c->fim->vertice);
    if(c->inicio!=c->fim){
        double custo=0;
        while(no->proximo!=c->fim){
            custo+= no->vertice->calculaCusto(no->proximo->vertice);
            no=no->proximo;
        }
        custo+= no->vertice->calculaCusto(c->fim->vertice);
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
    No*proximo=NULL;
    while(copia){
        if(solucao==NULL){
            solucao= new No();
            solucao->tabu=copia->tabu;
            solucao->vertice=copia->vertice;
            solucao->cluster=copia->cluster;
            solucao->proximo=proximo;
            proximo=solucao;
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
    if(clusterInicial->inicio->vertice->getIndiceCluster()==clusterFinal->inicio->vertice->getIndiceCluster()){
        clusterFinal->fim->proximo=clusterInicial->inicio;
        clusterInicial->inicio=clusterFinal->inicio;

        Cluster * cexcluir=clusterFinal;
        clusterFinal=clusterFinal->anterior;
        clusterFinal->proximo=NULL;
        clusterFinal->fim->proximo=NULL;
        delete cexcluir;
    }

    while(c!=NULL){
        Vertice *vEntrada = c->inicio->vertice;
        Vertice *vSaida = c->fim->vertice;
        fprintf(arqClusters,"%d\t %d\t %d\t %lf\n",vEntrada->getIndiceCluster()+1,vEntrada->getIDVertice(),vSaida->getIDVertice(),calculaCustoIntraCluster(c));
         if(primeiroX==NULL){
            contaX=1;
            primeiroX= new X();
            primeiroX->idX=contaX;
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
            aux->idX=contaX;
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

    c=clusterInicial;
    ///Faz a lista de clusters ficar circular
    clusterFinal->proximo=clusterInicial;
    clusterInicial->anterior=clusterFinal;

    fprintf(arqYr,"%d\t%d\n",contaYr,getNumTotalClusters());///numero de arestas inter cluster
     while(c!=clusterFinal){
        Vertice *vEntrada =  c->inicio->vertice;
        Vertice *vSaida = c->anterior->fim->vertice;
        fprintf(arqYr,"%d\t %d\t %lf\n",vEntrada->getIDVertice(),vSaida->getIDVertice(),vEntrada->calculaCusto(vSaida)+getPenalizacao());
        c=c->proximo;
    }

    ///salva a ultima aresta inter cluster ligando o final ao inicial
    Vertice *vEntrada =  c->inicio->vertice;
    Vertice *vSaida = c->anterior->fim->vertice;
    fprintf(arqYr,"%d\t %d\t %lf\n",vEntrada->getIDVertice(),vSaida->getIDVertice(),vEntrada->calculaCusto(vSaida)+getPenalizacao());


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
void desalocaClusters(){
    ClusterLista * pri = primeiroCluster;
    ClusterLista *aux;
    fprintf(stdout,"Excluindo CLUSTERS\n");
    while(pri!=NULL){
        aux=pri->proximo;
        Cluster *c= pri->cluster->proximo;
        while(c!=pri->cluster){
                Cluster *c2 = c->proximo;
                No* inicio = c->inicio;
                while(inicio!=c->fim){
                    No * n=inicio->proximo;
                    delete inicio;
                    inicio=n;
                }
                delete inicio;
                delete c;
                c=c2;
        }
        No* inicio = c->inicio;
        while(inicio!=c->fim){
            No * n=inicio->proximo;
            delete inicio;
            inicio=n;
        }
        delete inicio;
        delete pri->cluster;
        delete pri;
        pri = aux;
    }
    primeiroCluster = NULL;
     fprintf(stdout,"Excluindo X\n");
    X *x=primeiroX->proximo;
    while(primeiroX!=NULL){
         x=primeiroX->proximo;
        delete primeiroX;
        primeiroX=x;

    }
    fprintf(stdout,"Excluindo Yr\n");
    primeiroX=NULL;

    Yr *y;
    while(primeiroYr!=NULL){
        y=primeiroYr->proximo;

        delete primeiroYr;
        primeiroYr=y;
    }
    primeiroYr=NULL;
}
