#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

double signum(double k);
double frobeniusNorm(double **mat, int rows, int cols);
double** transpose(int n, int m, double ** matrix);
void jacobi(double **A, int max_sweeps, double toler, int m, int n, double **V, double **S, double **U);
void sort(int k, double**U, double**V, double**S, int m, int n);
void free_double_matrix(double **matrix, int rows);
void free_int_matrix(int **matrix, int rows);

#endif


