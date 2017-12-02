#include "Programa.h"
#include "Vertice.h"

///vertor de tabus do problema
static Tabu* tabus;
static ClusterInstancia* clusters;

///variaveis para a quantidade de tabus,vertices e clusters do problema
static int numTotalTabus;
static int numTotalClusters;
static int numTotalVertices;

///Ponteiro para a lista encadeada de vertices
static Vertice *primeiro;
static Vertice *ponteiro;

static FILE* arqSaida;
static double penalisacao = 0;
///calcula o custo da solucao


static int z=0;

static void atualizaPenalizacao(){
    double maiorCusto=0;
    for(Vertice * v = primeiro;v!=NULL;v=v->getProximo()){
                for(Vertice * v2 = primeiro;v2!=NULL;v2=v2->getProximo()){
                if(v!=v2){
                    double d = v->calculaCusto(v2);
                    if(d>maiorCusto){
                        maiorCusto = d;
                    }
                }
        }
    }
    penalisacao = 10*maiorCusto;
}
double getPenalizacao(){
    return penalisacao;
}
int calculaCustoSolucao(No* solucao){

    double somaCusto=0;
    No * aux ;
    int k=0;
    No* anterior;
    if(solucao!=NULL){
        for(aux = solucao; aux->proximo!=NULL;aux=aux->proximo){
                somaCusto+=aux->vertice->calculaCusto(aux->proximo->vertice);
                if(aux->vertice->getIndiceCluster()!= aux->proximo->vertice->getIndiceCluster()){
                    somaCusto+=penalisacao;
                }
                anterior=aux;

        }
        ///calcula a distancia do penultimo ao ultimo
        somaCusto+=anterior->vertice->calculaCusto(aux->vertice);

        if(anterior->vertice->getIndiceCluster()!= aux->vertice->getIndiceCluster()){
            somaCusto+=penalisacao;
        }

        ///calcula a distancia do primeiro ao ultimo

            somaCusto+=solucao->vertice->calculaCusto(aux->vertice);
        if(solucao->vertice->getIndiceCluster()!= aux->vertice->getIndiceCluster()){
            somaCusto+=penalisacao;
        }


    }

    return floor(somaCusto);
}

void imprimeSolucao(No* solucao){
     fprintf(arqSaida,"Rota -> ");
     for(No* aux = solucao; aux!=NULL;aux=aux->proximo){
            fprintf(arqSaida,"%d ", aux->vertice->getIDVertice());
    }
    fprintf(arqSaida,"\n");
}
No* solucao=NULL;
///Funcao calcula aondo o novo Vertice sera inserido na lista de solucao
void calculaInsercao(Vertice *v){
     fprintf(stdout,"inserindo: %d do tabu: %d\n",v->getIDVertice(),v->getIndiceTabu());
     No *anterior;
     No *novo = new No();

     novo->vertice=v;
     novo->proximo=NULL;
     novo->tabu=&tabus[v->getIndiceTabu()];
     novo->cluster=&clusters[v->getIndiceCluster()];
     if(solucao==NULL){
        novo->proximo=novo;
        solucao=novo;
     }else{

    if(v->getIDVertice()==35){
        int k;
    }

     double custoMinimo = 9999999;
     bool controleInsecaoIntraClusters=false;
     double custo=0;
     No *noAtual;
     for(noAtual=solucao;noAtual->proximo!=solucao;noAtual=noAtual->proximo){
        ///verifica se o cluster do vertice a ser inserido é igual a pelo menos um dos extremos (vertices de um mesmo cluster tem que estar consecutivos)
        if(noAtual->vertice->getIndiceCluster()==v->getIndiceCluster()||noAtual->proximo->vertice->getIndiceCluster()==v->getIndiceCluster()){
            custo=0;
            custo+=noAtual->vertice->calculaCusto(v);
            custo+=v->calculaCusto(noAtual->proximo->vertice);
            custo-=noAtual->vertice->calculaCusto(noAtual->proximo->vertice);

            if(custo<custoMinimo){
                custoMinimo=custo;
                anterior=noAtual;
            }
            controleInsecaoIntraClusters=true;
        }
     }
    ///verifica se o cluster do vertice a ser inserido é igual a pelo menos um dos extremos (vertices de um mesmo cluster tem que estar consecutivos)
    if(noAtual->vertice->getIndiceCluster()==v->getIndiceCluster()){
        custo=0;
        custo+=solucao->vertice->calculaCusto(v);
        custo+=v->calculaCusto(noAtual->vertice);
        custo-=noAtual->vertice->calculaCusto(solucao->vertice);
        if(custo<custoMinimo){
            anterior=noAtual;
        }
        controleInsecaoIntraClusters=true;
    }
    ///caso se tenha encontrado pelo menos um vertice de mesmo cluster insere na posicao calculada
    if(controleInsecaoIntraClusters){
        novo->proximo = anterior->proximo;
        anterior->proximo=novo;
        fprintf(stdout,"inserindo %d depois do %d\n",v->getIDVertice(),anterior->vertice->getIDVertice());
    }else{ /// senao insira o vertice no final da solucao
        noAtual->proximo= novo;
        novo->proximo=solucao;
     }
    }
}
int clusterMenorDiferenca(){
    ClusterInstancia* menor;
    int id=-1;
    int dif=9999999;
    for(int i=0;i<numTotalClusters;i++){
            if(clusters[i].getDifVisitados_Total()!=0){
                if(clusters[i].getDifVisitados_Total()<dif){
                    dif=clusters[i].getDifVisitados_Total();
                    menor = &clusters[i];
                    id=i;
                }
            }
    }
    return id;
}
No* construtivo(){
    int numTabusVisitados = 0;
    solucao;
    ///inicia a solucao com vertices aleatorios de 3 tabus aleatorios
    solucao = NULL;
    for(int i=0;i<numTotalClusters;i++){
        clusters[i].visitado=false;
        clusters[i].numPresenteSolucao=0;
    }
    numTabusVisitados=0;
    ///Perrcorre a lista de vertices inserindo na solução
    while(numTabusVisitados<numTotalTabus){
        Tabu *sorteado=NULL;
        int idC=0;
        while(sorteado==NULL){
            idC=rand() % numTotalClusters;
            ClusterInstancia *c = &clusters[idC];
            sorteado=c->getTabuRandom();
            if(sorteado&&sorteado->visitado==true){
                    sorteado=NULL;
            }
        }
        clusters[idC].numPresenteSolucao++;
        clusters[idC].visitado=true;
        sorteado->visitado=true;
        numTabusVisitados++;
        Vertice* v = sorteado->getRamdomCluster(idC);
        calculaInsercao(v);
        fprintf(stdout,"tabu: %d\n",v->getIndiceTabu());


    }
    for(int i=0;i<numTotalTabus;i++){
        tabus[i].visitado=false;
    }
    No * corta=solucao->proximo;
    while(corta->proximo!=solucao){
        corta=corta->proximo;
    }
    corta->proximo=NULL;

    bool ok=true;

    for(int i=0;i<numTotalClusters;i++){
       if(clusters[i].visitado==false){
            ok=false;
            break;

       }
    }
    /*if(ok){
        for(int i=0;i<numTotalTabus;i++){
           if(tabus[i].visitado==false){
                ok=false;
                fprintf(stdout,"Nem teodo mudo foi visitado\n");
                break;

           }
        }
    }*/
    if(!ok){
         desalocaSolucao(solucao);
         return NULL;
    }else{
         return solucao;
         }
}
ClusterInstancia* getCluster(int id){
    return &clusters[id];
}
///percorre o vetor de tabus imprimindo os seus respectivos vertices
void imprimeTabus(No* solucao){
    fprintf(arqSaida,"Tabus e seus Vertices:\n");
    for(int i=0;i<numTotalTabus;i++){
      fprintf(arqSaida,"T:%d ->",i+1);
      tabus[i].imprimeVertices(arqSaida);
      fprintf(arqSaida,"\n");
    }
}
 ///imprime os ids dos vertices seus tabus e clusters
void imprimeVerificacaoViabilidade(No* solucao){
    fprintf(arqSaida,"Solucao com tabus e clusters:\n");
    for(No * aux = solucao; aux!=NULL;aux=aux->proximo){
        fprintf(arqSaida,"Vertice: %d \tTabu: %d\t Cluster:%d\n",aux->vertice->getIDVertice(),aux->vertice->getIndiceTabu()+1,aux->vertice->getIndiceCluster());
    }
}
///efetua a leitura do fluxo passado como parametro
void inicializa(FILE* arquivoEntrada, FILE* arquivoSaida){
    arqSaida=arquivoSaida;

    fscanf(arquivoEntrada,"%d %d %d",&numTotalVertices,&numTotalClusters,&numTotalTabus);
    ///Aloca o vetor de Tabus
    tabus= new Tabu[numTotalTabus];
    clusters = new ClusterInstancia[numTotalClusters];
    double x,y;
    int id,cont=0,idT,idC;

    ///cria o primeiro vertice
    fscanf(arquivoEntrada, "%d\t%lf\t%lf\t%d\t%d\n",&id,&x,&y,&idC,&idT);
    ///Converte o tabus para a forma de represetnacao
    idT=idT-1;
    idC=idC-1;
    Vertice* inicial = new Vertice(id,x,y,idC,idT);

    tabus[idT].insereVertice(inicial);
    clusters[idC].adicionaTabu(&tabus[idT]);
    clusters[idC].numVertices++;
    ponteiro = inicial;
    primeiro = inicial;

    ///encadeia os demais vertices
    while(cont<numTotalVertices-1){
        fscanf(arquivoEntrada, "%d\t%lf\t%lf\t%d\t%d\n",&id,&x,&y,&idC,&idT);
        idT=idT-1;
        idC=idC-1;
        Vertice* aux = new Vertice(id,x,y,idC,idT);
        tabus[idT].insereVertice(aux);
        clusters[idC].adicionaTabu(&tabus[idT]);
        clusters[idC].numVertices++;
        ponteiro->setProximo(aux);
        ponteiro = aux;
        cont++;
    }
    fclose(arquivoEntrada);
    atualizaPenalizacao();

};

///salva a solucao no arquivo de saida
void salvarSolucao(No* solucao){
    int i = 0;
    Vertice *v;
    fprintf(arqSaida,"#Resultado\n");
    fprintf(arqSaida,"#Custo: %d\n",calculaCustoSolucao(solucao));
    fprintf(arqSaida,"#Penalizacao: %lf\n",penalisacao);
    fprintf(arqSaida,"#Rota: (idV|X|Y|idC|idT)\n");

    for(No *aux=solucao;aux!=NULL;aux=aux->proximo){
        v=aux->vertice;
        fprintf(arqSaida,"%d\t%lf\t%lf\t%d\t%d\n",v->getIDVertice(),v->getCordX(),v->getCordY(),v->getIndiceCluster()+1,v->getIndiceTabu()+1);
        i++;
    }
    v=solucao->vertice;
    fprintf(arqSaida,"%d\t%lf\t%lf\t%d\t%d\n",v->getIDVertice(),v->getCordX(),v->getCordY(),v->getIndiceCluster()+1,v->getIndiceTabu()+1);
    if(i==numTotalTabus){
        fprintf(arqSaida,"#Todos os Tabus foram visitados\n");
    }
}
int getNumTotalClusters(){
    return numTotalClusters;
}
int getNumTotalTabus(){
    return numTotalTabus;
}
///desaloca uma solucao com o ponteiro inicial dado pelo parametro solucao
void desalocaSolucao(No* solucao){
   No* aux;
    if(solucao!=NULL){
        aux=solucao->proximo;
        delete solucao;
        while(aux!=NULL){
            solucao = aux;
            aux = aux->proximo;
            delete solucao;
        }
        solucao=NULL;
    }

}

void desalocaMemoria(No* solucao){
    ///Desaloca a lista encadeada de vertices
    ponteiro=primeiro;
    while(ponteiro!=NULL){
        ponteiro=primeiro->getProximo();
        primeiro->~Vertice();
        primeiro=ponteiro;
    }

    ///Desaloca o vetor de tabus
    delete ponteiro;
    for(int i=0;i<numTotalTabus;i++){
        tabus[i].~Tabu();
    }

    ///Desaloca os elementos da solucao
    No* aux;
    if(solucao!=NULL){
        aux=solucao->proximo;
        delete solucao;
        while(aux!=NULL){
            solucao = aux;
            aux = aux->proximo;
            delete solucao;
        }
        solucao=NULL;
    }
}
