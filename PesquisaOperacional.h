#ifndef PESQUISAOPERACIONAL_H_INCLUDED
#define PESQUISAOPERACIONAL_H_INCLUDED
#include "Programa.h"
void inicializaArquivosEscrita(char* nomeArquivoClusters,char * nomeArquivoYr);
void inicializaLeitura(char* nomeArquivoClusters,char * nomeArquivoYr);
void salvarSolucaoArquivosPO(No* solucao);
void emitirSistemaLinear(char *nomeArquivoSl);
void finalizarArquivosEscrita();
void imprimeX();
void desalocaClusters();
#endif // PESQUISAOPERACIONAL_H_INCLUDED
