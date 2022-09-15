#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>
#include <string.h>

typedef struct
{
    double eigenvalue;
    double *eigenvector;
    
}eigenvector_object;

double** create_I_mat(int n);

double** convert_arr_to_mat(double *points_arr, int n, int d);

double** create_empty_mat(int n, int d);

double** create_W_mat(double **points, int n, int d);

double** create_DD_mat(double **W_mat, int n);

double** mat_mult(double **mat1, double **mat2, int n);

double** lnorm_calc(double **points, int n, int d);


double*** apply_Jacobi(double **A, int n);

int compare_vectors(const void *vector1, const void *vector2);

int calc_k(eigenvector_object *eigenvectors, int n);


double** calc_U_mat(eigenvector_object *eigenvectors_lst, int n, int k);

double get_row_norm(double *row_vec,int vecDim);

void normalize(double **uMat, int cnt, int clust_num);

double** calc_T_mat(eigenvector_object *eigenvector_lst, int n, int k);


void free_mat(double** mat, int n);

eigenvector_object* create_eigenlist(double **A, double **V, int n);

void print_mat(double **mat, int n, int d);

void print_eigenvalues(double **eigenvalues_mat, int n, int d);