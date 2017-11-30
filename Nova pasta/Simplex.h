#ifndef SIMPLEX_H_INCLUDED
#define SIMPLEX_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define M 20
#define N 20

int equal(double a, double b);

typedef struct {
    int m, n; // m=rows, n=columns, mat[m x n]
    double mat[M][N];
} Tableau;

void nl(int k);

void print_tableau(Tableau *tab, const char* mes);

/* Example input file for read_tableau:
     4 5
      0   -0.5  -3 -1  -4
     40    1     1  1   1
     10   -2    -1  1   1
     10    0     1  0  -1
*/
void read_tableau(Tableau *tab, const char * filename);

void pivot_on(Tableau *tab, int row, int col);

// Find pivot_col = most negative column in mat[0][1..n]
int find_pivot_column(Tableau *tab);

// Find the pivot_row, with smallest positive ratio = col[0] / col[pivot]
int find_pivot_row(Tableau *tab, int pivot_col);

void add_slack_variables(Tableau *tab);

void check_b_positive(Tableau *tab);

// Given a column of identity matrix, find the row containing 1.
// return -1, if the column as not from an identity matrix.
int find_basis_variable(Tableau *tab, int col);

void print_optimal_vector(Tableau *tab, char *message);

void simplex(Tableau *tab);


#endif // SIMPLEX_H_INCLUDED
