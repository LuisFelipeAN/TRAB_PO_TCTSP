#ifndef PROGRAMA_H_INCLUDED
#define PROGRAMA_H_INCLUDED
#include <fstream>
#include <string.h>
#include <stdio.h>
#include "Tabu.h"
#include "ClusterInstancia.h"
using namespace std;

///estrutura para encadear os vertices presentes na solucao
typedef struct No{
    No *anterior;///ponteiro do anterior e atualizado apenas para efetuar uma busca local
    No *proximo;
    Tabu* tabu;
    Vertice* vertice;
    ClusterInstancia* cluster;
} No;
double getPenalizacao();

void inicializa(FILE* arquivoEntrada,FILE* arquivoSaida);

No* construtivo();

void imprimeSolucao(No* solucao);

void imprimeTabus(No* solucao);

void imprimeVerificacaoViabilidade(No* solucao);

ClusterInstancia* getCluster(int id);


int calculaCustoSolucao(No* solucao);

void desalocaMemoria(No* solucao);

void desalocaSolucao(No* solucao);

void salvarSolucao(No* solucao);
#endif // PROGRAMA_H_INCLUDED
