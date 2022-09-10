#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>
#include <string.h>


double** create_empty_mat(int n, int d);

double** create_W_mat(double **points, int n, int d);

double** create_DD_mat(double **W_mat, int n, int d);

double** mat_mult(double **mat1, double **mat2, int n);

double** lnorm_calc(double **points, int n, int d);

double*** apply_Jacobi(double **A, int n, int d);

int compare_vectors(const void *vector1, const void *vector2);

int get_k(double **A, double **V, int n);