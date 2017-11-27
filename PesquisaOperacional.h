#ifndef PESQUISAOPERACIONAL_H_INCLUDED
#define PESQUISAOPERACIONAL_H_INCLUDED
#include "Programa.h"
void inicializaArquivosEscrita(char* nomeArquivoClusters,char * nomeArquivoYr);
void inicializaLeitura(char* nomeArquivoClusters,char * nomeArquivoYr);
void salvarSolucaoArquivosPO(No* solucao);
void finalizarArquivosEscrita();
#endif // PESQUISAOPERACIONAL_H_INCLUDED
