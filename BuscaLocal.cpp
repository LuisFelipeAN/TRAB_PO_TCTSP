#include "BuscaLocal.h"
///estrutura para encadear os clusters presentes na solucao
typedef struct Cluster{ ///estrutura para fragmentar a solucao em clusters com vertices de entrada e saida do cluster
    No*inicio;
    No*fim;
    Cluster* proximo;
    Cluster* anterior;///duplamente encadeada para efetuar a troca tanto na solucao como tambem na lista
}Cluster;///utulizada nas funcoes BuscaLocal2 e BuscaLocal4

typedef struct ClusterOrdenado{///estrutura para ordenar os clusters de acordo com o peso de suas arestas de entrada e saida
    Cluster* clusterAtual;
    ClusterOrdenado *proximo;
    double custo;
}ClusterOrdenado;///utulizada nas funcoes BuscaLocal2 e BuscaLocal4

typedef struct NoOrdenado{///estruturas para ordenar os tabus da solução  de acordo com o peso das arestas de entrada e saida
    No* noAnterior;
    No* noAtual;        ///sao necessarios os vertices anterior e proximo tanto para garantir a viabilidade da alteracaom
    No* noProximo;          ///como para possibilitar o calculo dos pesos se vai melhorar ou nao com a alteracao
    NoOrdenado *proximo;
    double custo;///utilizada na funcao de BuscaLocal
}NoOrdenado;

typedef struct TabuOrdenado{///estrura para ordenar os tabus da solução de forma a efetuar troca 2 a 2
    No *noAtual;
    TabuOrdenado *proximo;
    double custo;
}TabuOrdenado;///utilizada na funcao de BuscaLocal3

///Busca o vertice de menor custo dentro do tabu
bool melhorVertice(NoOrdenado *noOrdenado){
    No* anterior=noOrdenado->noAnterior;
    No* atual =noOrdenado->noAtual;
    No* proximo=noOrdenado->noProximo;
    ClusterInstancia * cBkp=atual->cluster;

    Vertice *melhor = atual->vertice;
    bool atualizou=false;
    double custoAtual=anterior->vertice->calculaCusto(atual->vertice)+atual->vertice->calculaCusto(proximo->vertice);///calcula o peso atual das aresatas
    Vertice *candidato=atual->vertice;
    while(!atual->tabu->efetuouTodasTrocas()){///verifica se nao testou todas as possibilidades de troca
        candidato = atual->tabu->outroVertice(candidato);
        if(candidato->getIndiceCluster()==atual->vertice->getIndiceCluster()|candidato->getIndiceCluster()==proximo->vertice->getIndiceCluster()|candidato->getIndiceCluster()==anterior->vertice->getIndiceCluster()){
            double novoCusto=anterior->vertice->calculaCusto(candidato)+candidato->calculaCusto(proximo->vertice);

             if (candidato->getIndiceCluster()!=atual->vertice->getIndiceCluster()){

            }///verificar se o cluster so esta sendo visitado por este vertice
            if(novoCusto<custoAtual){

                atualizou=true;
                melhor=candidato;
                custoAtual=novoCusto;
            }

        }
    }
    if(atualizou){
        atual->vertice=melhor;
        atual->cluster=getCluster(melhor->getIndiceCluster());
        if(atual->cluster!=cBkp){///se trocou de cluster atualiza os valores de vertices presentes na solucao
            cBkp->numPresenteSolucao--;
            atual->cluster->numPresenteSolucao++;
        }
    }
    return atualizou;

}
///Escolhe o melhor vertice de cada tabu
void buscaLocal(No* solucao){
    Vertice* candidato; ///vertice candidado a ser o substitudo na solução
    Vertice* bkp; ///armazena o vertice que foi alterado na solução
    int custoAtual = calculaCustoSolucao(solucao);///armazena o custo atual da solução
    int novoCusto;///armazena o custo da nova soluçã
    No * ultimo=solucao;
    No * anterior;
    while(ultimo->proximo!=NULL){///encontra o ultimo No da solucao;
        ultimo=ultimo->proximo;
    }
    ///primeiro no ordenado
    NoOrdenado *noOrdenadoInicial = new NoOrdenado();
    noOrdenadoInicial->noAnterior=ultimo;
    noOrdenadoInicial->noAtual=solucao;
    noOrdenadoInicial->noProximo=solucao->proximo;
    noOrdenadoInicial->custo=solucao->tabu->getNumVertices();
    noOrdenadoInicial->proximo=NULL;

    ///cria uma lista ordenadas de nos da solucao
    anterior=solucao;
    No* atual = solucao->proximo;
    while(atual!=NULL){///insere os tabus em uma lista ordenando pelo maior numero de vertices
        NoOrdenado *aux = new NoOrdenado();
        aux->noAnterior=anterior;
        aux->noAtual=atual;
        if(atual->proximo==NULL){
             aux->noProximo=solucao;
        }else{
             aux->noProximo=atual->proximo;
        }
        aux->custo=atual->tabu->getNumVertices();
        aux->proximo=NULL;

        if(aux->custo>noOrdenadoInicial->custo){
            aux->proximo=noOrdenadoInicial;
            noOrdenadoInicial=aux;
        }else{
            NoOrdenado *noOrdAnterior=noOrdenadoInicial;
            NoOrdenado *noOrdProximo=noOrdenadoInicial->proximo;
            if(noOrdProximo==NULL){
                noOrdAnterior->proximo=aux;
            }
            else{
                bool atualizou=false;
                while(noOrdProximo!=NULL){
                     if(aux->custo>noOrdProximo->custo){
                        noOrdAnterior->proximo=aux;
                        aux->proximo=noOrdProximo;
                        atualizou=true;
                        break;
                    }
                    noOrdAnterior=noOrdProximo;
                    noOrdProximo=noOrdProximo->proximo;

                }
                if(!atualizou){
                    noOrdAnterior->proximo=aux;
                }

            }
        }
        anterior=atual;
        atual=atual->proximo;
    }
    NoOrdenado *p;
    p=noOrdenadoInicial;
    while(p!=NULL){
        if(p->noAtual->cluster->numPresenteSolucao>1){///se o cluster atual tem mais de um vertice presente na solucao efetua a troca se possivel
                if(melhorVertice(p)){
                    p=noOrdenadoInicial;
                }else{
                    p=p->proximo;
                }
        }else{
            p=p->proximo;
        }

    }
    p=noOrdenadoInicial->proximo;
    delete noOrdenadoInicial;
    while(p!=NULL){
        noOrdenadoInicial=p;
        p=p->proximo;
        delete noOrdenadoInicial;
    }
}

///Funcao que efetua a troca de dois clustes adjacentes tanto na solucao quanto na lista duplamente encadeada atual
bool trocaClusters(ClusterOrdenado *coAtual){

    Cluster* cAtual=coAtual->clusterAtual;
    Cluster* cAnt  =coAtual->clusterAtual->anterior;
    Cluster* cProx =coAtual->clusterAtual->proximo;
    ///primeiro verifica se a troca valerá a pena
    bool custoAtual = cAnt->fim->vertice->calculaCusto(cAtual->inicio->vertice)+cAtual->fim->vertice->calculaCusto(cProx->inicio->vertice)+cProx->fim->vertice->calculaCusto(cProx->proximo->inicio->vertice);
    bool novoCusto  = cAnt->fim->vertice->calculaCusto(cProx->inicio->vertice)+cProx->fim->vertice->calculaCusto(cAtual->inicio->vertice)+cAtual->fim->vertice->calculaCusto(cProx->proximo->inicio->vertice);
    if(novoCusto<custoAtual){

        No* bkpProxAnt=cAnt->fim->proximo;

        if(cAtual->fim->proximo!=NULL){
            cAnt->fim->proximo=cProx->inicio;
        }
        else
            cAnt->fim->proximo=NULL;

        ///o fim do atual aponta para o fim do proximo
        cAtual->fim->proximo=cProx->fim->proximo;


        ///o fim do proximo aponta para o inicio do atual
        if(bkpProxAnt==NULL){
            cProx->fim->proximo=NULL;
        }else{
             cProx->fim->proximo=cAtual->inicio;
        }

        cProx->anterior=cAnt;
        cAnt->proximo=cProx;

        cAtual->proximo=cProx->proximo;
        cAtual->anterior=cProx;

        cProx->proximo=cAtual;
        return true;
    }
    else return false;

}

///Efetua uma busca local trocando clusters adjacentes baseado em uma lista duplamente encadeada de clusters
void buscaLocal2(No* solucao){

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

    Cluster *ultimoAdicionado=NULL;

    ///emenda a lista de clusters crianda uma lista circular
    clusterInicial->anterior=clusterFinal;
    clusterFinal->proximo=clusterInicial;

    ClusterOrdenado *clusterOrdenadoInicial = new ClusterOrdenado();
    clusterOrdenadoInicial->clusterAtual=clusterInicial;
    clusterOrdenadoInicial->custo=solucao->vertice->calculaCusto(ultimo->vertice)+clusterFinal->anterior->fim->vertice->calculaCusto(clusterFinal->inicio->vertice);
    clusterOrdenadoInicial->proximo=NULL;

    Cluster* atual=clusterInicial->proximo;

    while(atual!=clusterInicial){///insere os clusters em uma lista ordenando pelo maior peso das arestas de entrada e saida
        ClusterOrdenado *aux = new ClusterOrdenado();
        aux->clusterAtual=atual;
        aux->custo=atual->anterior->fim->vertice->calculaCusto(atual->inicio->vertice)+atual->fim->vertice->calculaCusto(atual->proximo->inicio->vertice);
        aux->proximo=NULL;

        ///insere ordenado na lista
        if(aux->custo>clusterOrdenadoInicial->custo){
            aux->proximo=clusterOrdenadoInicial;
            clusterOrdenadoInicial=aux;
        }else{
            ClusterOrdenado *clusterOrdAnterior=clusterOrdenadoInicial;
            ClusterOrdenado *clusterOrdProximo=clusterOrdenadoInicial->proximo;
            if(clusterOrdProximo==NULL){
                clusterOrdAnterior->proximo=aux;
            }
            else{
                bool atualizou=false;
                while(clusterOrdProximo!=NULL){
                     if(aux->custo>clusterOrdProximo->custo){
                        clusterOrdAnterior->proximo=aux;
                        aux->proximo=clusterOrdProximo;
                        atualizou=true;
                        break;
                    }
                    clusterOrdAnterior=clusterOrdProximo;
                    clusterOrdProximo=clusterOrdProximo->proximo;

                }
                if(!atualizou){
                    clusterOrdAnterior->proximo=aux;
                }

            }
        }
        atual=atual->proximo;
    }
    ClusterOrdenado * ponteiro = clusterOrdenadoInicial;
    ClusterOrdenado *p;
    while(ponteiro!=NULL){
        if(trocaClusters(ponteiro)){
            ponteiro=clusterOrdenadoInicial;
        }else{
          ponteiro=ponteiro->proximo;
        }

    }
    ///desaloca a lista ordenada de clusters
    p=clusterOrdenadoInicial->proximo;
    delete clusterOrdenadoInicial;
    while(p!=NULL){
        clusterOrdenadoInicial=p;
        p=p->proximo;
        delete clusterOrdenadoInicial;

    }
    ///desaloca a lista de clusters criada
    Cluster *c;
    c = clusterInicial;
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

bool trocaVertices(TabuOrdenado *toAtual){///efetua a troca de dois vertices na solucao

    No* nAtual=toAtual->noAtual;
    No* nAnt  =toAtual->noAtual->anterior;
    No* nProx =toAtual->noAtual->proximo;
    if(nProx->vertice->getIndiceCluster()==nAtual->vertice->getIndiceCluster())	{
	    double custoAtual=nAnt->vertice->calculaCusto(nAtual->vertice)+nAtual->vertice->calculaCusto(nProx->vertice)+nProx->vertice->calculaCusto(nProx->proximo->vertice);
	    double novoCusto =nAnt->vertice->calculaCusto(nProx->vertice)+nProx->vertice->calculaCusto(nAtual->vertice)+nAtual->vertice->calculaCusto(nProx->proximo->vertice);
	    if (novoCusto< custoAtual){

	        nProx->proximo->anterior=nAtual;
	        nAtual->proximo=nProx->proximo;
	        nProx->proximo=nAtual;

	        nProx->anterior=nAnt;
	        nAtual->anterior=nProx;
	        nAnt->proximo=nProx;

	        return true;
	    }
	    return false;
	}
    return false;
}

///Efetua uma busca local trocando vertices adjacentes que possuam um mesmo cluster
void buscaLocal3(No* solucao){

    No * ultimo=solucao;
    No * auxAnt=NULL;
    ///encontra o ultimo vertice atualizando os ponteiros da lista duplamente encadeada
    while(ultimo->proximo!=NULL){///encontra o ultimo No da solucao;
        ultimo->anterior=auxAnt;
        auxAnt=ultimo;
        ultimo=ultimo->proximo;

    }

    ///emenda a lista solucao criando uma lista circular
    ultimo->anterior=auxAnt;
    ultimo->proximo=solucao;
    solucao->anterior=ultimo;


    TabuOrdenado *tabuOrdenadoInicial = new TabuOrdenado();
    tabuOrdenadoInicial->custo=solucao->anterior->vertice->calculaCusto(solucao->vertice)+solucao->vertice->calculaCusto(solucao->proximo->vertice);
    tabuOrdenadoInicial->noAtual=solucao;
    tabuOrdenadoInicial->proximo=NULL;


    No* atual = solucao->proximo;
    while(atual!=solucao){///insere os tabus ordenando pelo peso das arestas de entrada e saida
        TabuOrdenado *aux = new TabuOrdenado();
        aux->noAtual=atual;
        aux->custo=atual->anterior->vertice->calculaCusto(atual->vertice)+atual->vertice->calculaCusto(atual->proximo->vertice);
        aux->proximo=NULL;

        if(aux->custo>tabuOrdenadoInicial->custo){
            aux->proximo=tabuOrdenadoInicial;
            tabuOrdenadoInicial=aux;
        }else{
            TabuOrdenado *noOrdAnterior=tabuOrdenadoInicial;
            TabuOrdenado *noOrdProximo=tabuOrdenadoInicial->proximo;
            if(noOrdProximo==NULL){
                noOrdAnterior->proximo=aux;
            }
            else{
                bool atualizou=false;
                while(noOrdProximo!=NULL){
                     if(aux->custo>noOrdProximo->custo){
                        noOrdAnterior->proximo=aux;
                        aux->proximo=noOrdProximo;
                        atualizou=true;
                        break;
                    }
                    noOrdAnterior=noOrdProximo;
                    noOrdProximo=noOrdProximo->proximo;

                }
                if(!atualizou){
                    noOrdAnterior->proximo=aux;
                }
            }
        }
        atual=atual->proximo;
    }
    TabuOrdenado *p;
    TabuOrdenado *ponteiro= tabuOrdenadoInicial;
    while(ponteiro!=NULL){///efetua a busca local
        if(trocaVertices(ponteiro)){
            ponteiro=tabuOrdenadoInicial;
        }else{
            ponteiro=ponteiro->proximo;
        }
    }

    ///desemenda a lista circular
    ultimo=solucao;
    while(ultimo->proximo!=solucao){///encontra o ultimo No da solucao;
        ultimo=ultimo->proximo;

    }
    ultimo->proximo=NULL;
    ///desaloca a lista ordenada de tabus
    p=tabuOrdenadoInicial->proximo;
    delete tabuOrdenadoInicial;
    while(p!=NULL){
        tabuOrdenadoInicial=p;
        p=p->proximo;
        delete tabuOrdenadoInicial;
    }

}
///efetivamente inverte a ordem de visita de todos os vertices do clueter
bool inverteOrdemVisita(Cluster *cAtual){
    No* bkpProxFim=cAtual->fim->proximo;

    No* ant= cAtual->inicio;
    No* prox= cAtual->inicio->proximo;
    if(prox->vertice->getIndiceCluster()!=ant->vertice->getIndiceCluster()){
        return false;
    }

    while(prox!=cAtual->fim){
        No* aux=prox->proximo;
        prox->proximo=ant;
        ant=prox;
        prox=aux;
    }
    prox->proximo=ant;

    No* aux=cAtual->inicio;
    cAtual->inicio=cAtual->fim;
    cAtual->fim=aux;

    cAtual->fim->proximo= bkpProxFim;
    cAtual->anterior->fim->proximo=cAtual->inicio;

    return true;

}
///verifica se inverter a ordem de visita melhora a solucao caso sim chama a funcao que inverte
bool inverteOrdemVisita(ClusterOrdenado *clusterOrdenado){
    Cluster *cAtual = clusterOrdenado->clusterAtual;
    Cluster *cAnterior = clusterOrdenado->clusterAtual->anterior;
    Cluster *cProximo= clusterOrdenado->clusterAtual->proximo;
    ///verifica se compensa inverter a ordem de visita do cluster
    double custoLigacaAtual=cAnterior->fim->vertice->calculaCusto(cAtual->inicio->vertice)+cAtual->fim->vertice->calculaCusto(cProximo->inicio->vertice);
    double custoNovaLigacao=cAnterior->fim->vertice->calculaCusto(cAtual->fim->vertice)+cAtual->inicio->vertice->calculaCusto(cProximo->inicio->vertice);
    if(custoNovaLigacao<custoLigacaAtual){
             if(inverteOrdemVisita(cAtual)){
                 return true;
             }
        }else return false;
}

///Efetua uma busca local trocando clusters adjacentes baseado em uma lista duplamente encadeada de clusters
void buscaLocal4(No* solucao){

    Cluster* clusterFinal = NULL;///ultimo cluster da solucao
    int controle=1;
    int idClusterAtual=-1;

    No* ultimo=solucao;
    No* auxAnt=NULL;
     while(ultimo->proximo!=NULL){///encontra o ultimo No da solucao;
        ultimo->anterior=auxAnt;
        auxAnt=ultimo;
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
    Cluster *ultimoAdicionado=NULL;


    ClusterOrdenado *clusterOrdenadoInicial = new ClusterOrdenado();
    clusterOrdenadoInicial->clusterAtual=clusterInicial;
    clusterOrdenadoInicial->custo=solucao->tabu->getNumVertices();
    clusterOrdenadoInicial->proximo=NULL;

    Cluster* atual=clusterInicial->proximo;
    while(atual!=NULL){
        ClusterOrdenado *aux = new ClusterOrdenado();
        aux->clusterAtual=atual;
        if(atual->proximo==NULL){
            aux->custo=atual->anterior->fim->vertice->calculaCusto(atual->inicio->vertice)+atual->fim->vertice->calculaCusto(clusterInicial->inicio->vertice);
        }else{
            aux->custo=atual->anterior->fim->vertice->calculaCusto(atual->inicio->vertice)+atual->fim->vertice->calculaCusto(atual->proximo->inicio->vertice);
        }
        aux->proximo=NULL;

        ///insere ordenado na lista
        if(aux->custo>clusterOrdenadoInicial->custo){
            aux->proximo=clusterOrdenadoInicial;
            clusterOrdenadoInicial=aux;
        }else{
            ClusterOrdenado *clusterOrdAnterior=clusterOrdenadoInicial;
            ClusterOrdenado *clusterOrdProximo=clusterOrdenadoInicial->proximo;
            if(clusterOrdProximo==NULL){
                clusterOrdAnterior->proximo=aux;
            }
            else{
                bool atualizou=false;
                while(clusterOrdProximo!=NULL){
                     if(aux->custo>clusterOrdProximo->custo){
                        clusterOrdAnterior->proximo=aux;
                        aux->proximo=clusterOrdProximo;
                        atualizou=true;
                        break;
                    }
                    clusterOrdAnterior=clusterOrdProximo;
                    clusterOrdProximo=clusterOrdProximo->proximo;

                }
                if(!atualizou){
                    clusterOrdAnterior->proximo=aux;
                }

            }
        }
        atual=atual->proximo;


    }

    clusterInicial->anterior=clusterFinal;
    clusterFinal->proximo=clusterInicial;

    ClusterOrdenado * ponteiro = clusterOrdenadoInicial;
    double custoAtual=calculaCustoSolucao(solucao);
    ///imprime o resultado
    ClusterOrdenado *p;
    p=clusterOrdenadoInicial;

    ///emenda a lista encadeada trasformando em uma lista circular
    ultimo->anterior=auxAnt;
    ultimo->proximo=solucao;
    solucao->anterior=ultimo;

    while(ponteiro!=NULL){
        if(inverteOrdemVisita(ponteiro)){
            ponteiro=clusterOrdenadoInicial;
        }else{
            ponteiro=ponteiro->proximo;
        }
    }
    ultimo=solucao;
    ///desemenda a lista circular
    while(ultimo->proximo!=solucao){///encontra o ultimo No da solucao;
        ultimo=ultimo->proximo;

    }
    ultimo->proximo=NULL;

    ///desaloca a lista ordenada de clusters
    p=clusterOrdenadoInicial->proximo;
    delete clusterOrdenadoInicial;
    while(p!=NULL){
        clusterOrdenadoInicial=p;
        p=p->proximo;
        delete clusterOrdenadoInicial;

    }


    ///desaloca a lista de clusters criada
    Cluster *c;
    c = clusterInicial;
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
