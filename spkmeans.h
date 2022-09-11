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

double** convert_1D_arr_to_mat(double *points_arr, int n, int d);

double** create_empty_mat(int n, int d);

double** create_W_mat(double **points, int n, int d);

double** create_DD_mat(double **W_mat, int n, int d);

double** mat_mult(double **mat1, double **mat2, int n);

double** lnorm_calc(double **points, int n, int d);

double*** apply_Jacobi(double **A, int n, int d);

int compare_vectors(const void *vector1, const void *vector2);

int get_k(double **A, double **V, int n);

double* get_eigen_val(double **eigenValuesMat, int cnt);

double** get_U_mat(eigenvector_object *eigen_pairs, int cnt, int clust_num);

double get_row_norm(double *row_vec,int vecDim);

void normalize(double **uMat, int cnt, int clust_num);

double** get_T_mat(eigenvector_object *eigen_pairs, int cnt, int clust_num);

void free_mat(double** mat, int n);




