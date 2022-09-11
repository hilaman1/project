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

int get_dimention(FILE *ifp){
    /* a function that returns the dimention of the points */
    char c;
    int size;
    size = 0;

    while((c = getc(ifp)) != '\n'){
        if (c == ','){
            size ++;
        }
    }
    size ++;
    fseek(ifp, 0L, SEEK_SET);
    rewind(ifp);

    return (size);
}

int get_points_num(FILE *ifp){
    /* a function that returns the number of points. each line in input file represents a points.
    thus, we'll return the number of lines on input file. */
    char c;
    int line_counter;

    line_counter = 0;

    while( (c = getc(ifp)) !=  EOF )
    {
        if(c == '\n'){
            line_counter++;
        }
    }
    fseek(ifp, 0L, SEEK_SET);
    rewind(ifp);
    return (line_counter);
}

void get_points_from_input(FILE *ifp, double* points_arr){
    /* reads the points from the file and save them in a given array */
    double point = (double) 0;
    int i = 0;

    while (fscanf(ifp, "%lf,", &point) != EOF)
    {
        points_arr[i] = point;
        i++;
    }
}

double** create_empty_mat(int n, int d){
    /* creates an empty matrix from order nxd */
    double **mat; 
    double *p; 
    int i;
    p = calloc(n*d, sizeof(double));
    if (p == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    mat = calloc(n, sizeof(double*));
    if (mat == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    for ( i=0; i<n; i++){
        mat[i] = p + i*d;
    }
    return mat;
}

double** convert_1D_arr_to_mat(double *points_arr, int n, int d){
    int i,j;
    double **points_mat;
    points_mat = create_empty_mat(n,d);
    if (points_mat == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    for( i = 0; i < n ;i++){
        for( j = 0; j < d; j++){
            points_mat[i][j] = points_arr[i * d + j];
        }
    }
    return points_mat;
}



double L2_norm(double *x1, double *x2, int d){
    /* return L2 norm of x1 and x2 */
    double sum = 0;
    double delta;
    int i;
    for( i = 0; i < d ;i++)
    {
        delta = pow((x1[i] - x2[i]), 2);
        sum += delta;
    }
    return sum;
}


double** create_W_mat(double **points, int n, int d){
    /* return a mat W s.t W_mat[i][j] = Radial basis function kernel */
    double **W;
    double delta;
    double temp_val;
    int i,j;

    W = create_empty_mat(n,n);
    if (W == NULL){
        return NULL;
    }

    for ( i = 0; i < n; ++i){
        for( j = 0; j < n; ++j){
            if (i == j){
                W[i][j] = 0;
            }
            else{
                delta = L2_norm(points[i], points[j], d);
                temp_val = -1 * (sqrt(delta) / 2);
                W[i][j] = exp(temp_val);
            }
        }
    }
    return W;
}

void print_mat(double **mat, int n, int d){
    /* prints the matrix values */
    int i,j;
    for ( i = 0; i < n; i++){
        for ( j = 0; j < d; j++){
            printf("%.4f", mat[i][j]);
            if (j < (d - 1)){
                printf(",");
            }
            else{
                printf("\n");
            }
        }
    }
}

void print_eigenvalues(double **eigenvalues_mat, int n, int d){
    /* prints the eigenvalues */
    int i;
    for ( i = 0; i < n; i++){
        printf("%.4f", eigenvalues_mat[i][i]);
        if (i < (d - 1)){
            printf(",");
        }
    }
    printf("\n");
}

double** create_DD_mat(double **W_mat, int n, int d){
    /* return a mat D s.t D_mat[i][j] = Diagonal Degree Matrix */
    double **D_mat;
    double d_i;
    int i,j;

    D_mat = create_empty_mat(n, n);
    if (D_mat == NULL){
        return NULL;
    }

    for ( i = 0; i < n; i++){
        d_i = 0;
        for ( j = 0; j < d; j++){
            d_i += W_mat[i][j];
        }
        D_mat[i][i] = 1/sqrt(d_i);
    }
    return D_mat;
}

double** mat_mult(double **mat1, double **mat2, int n){
    double **mult;
    int j;
    int k;
    int i;
    double tmp;
    mult = calloc(n, sizeof(double*));
    for (i = 0; i < n; i++){
        mult[i] = calloc(n, sizeof(double));
    }
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            tmp=0;
            for (k=0; k<n; k++){
                tmp+= mat1[i][k]*mat2[k][j];
            }
            mult[i][j]= tmp;
        }
    }
    return mult;

}

void free_mat(double** mat, int n){
    int i;
    for( i = 0; i < n; i++){
        free(mat[i]);
    }
    free(mat);
}


double** lnorm_calc(double **points, int n, int d){
    int i;
    double **D_mat;
    double **mult;
    double **W_mat;
    
    W_mat = create_W_mat(points, n, d);
    if (W_mat == NULL){
        return NULL;
    }
    
    D_mat = create_DD_mat(W_mat, n, d);
    D_mat = calloc(n, sizeof(double*));
    if (D_mat == NULL){
        return NULL;
    }

    mult = calloc(n, sizeof(double*));
    for (i = 0; i <n; i++){
        mult[i] = calloc(n, sizeof(double));
    }
    
    mult = mat_mult(mat_mult(D_mat, W_mat, n), D_mat, n);

    for (i = 0; i < n; i++){
        mult[i][i] = 1 - mult[i][i];
    }

    free_mat(W_mat, n);
    free_mat(D_mat, n);

    return mult;
}

double** create_I_mat(int n){
    /* return an identity matrix of size n */
    double ** I_mat;
    int i,j;
    I_mat = create_empty_mat(n,n);
    if (I_mat == NULL)
    {
        return NULL;
    }

    for ( i = 0; i < n; i++){
        for ( j = 0; j < n; j++){
            if (i == j){
                I_mat[i][j] = 1;
            }else{
                I_mat[i][j] = 0;
            }
        }
    }
    return I_mat;
}

int* get_max_not_in_diag(double** sym_mat, int n){
    /* calculates the Aij entery for pivoting */
    int* indexes_to_return;
    int coord_i, coord_j;
    double max_abs = 0;
    int i,j;

    indexes_to_return = (int*) calloc(2, sizeof(int));
    if (indexes_to_return == NULL){
        return NULL;
    }

    for( i = 0; i < n; i++){
        for( j = i; j < n; j++){
            if ((i != j) && (fabs(sym_mat[i][j]) > fabs(max_abs))){
                max_abs = sym_mat[i][j];
                coord_i = i;
                coord_j = j;
            }
        }
    }

    indexes_to_return[0] = coord_i;
    indexes_to_return[1] = coord_j;

    return indexes_to_return;
}

double sign(double num){
    /* sign function */
    if (num >= 0){
        return 1;
    }else{
        return -1;
    }
}

double** transpose(double** mat, int n){
    /* recieves a matrix dimantion and calculates the transpose of it */
    double** transpose;
    int i,j;
    transpose = create_empty_mat(n, n);
    if (transpose == NULL){
        return NULL;
    }
    for( i = 0; i < n; i++){
        for( j = 0; j < n; j++){
            transpose[i][j] = mat[j][i];
        }
    }
    return transpose;
}

double** calc_P_mat(double **A, int n){
    double **P;
    int i,j; /* indices for Aij */
    double theta,t,c,s;
    int *indices_max_abs;
    
    indices_max_abs = get_max_not_in_diag(A, n); /* returns an array with 2 indices */
    i = indices_max_abs[0];
    j = indices_max_abs[1];

    theta = (A[j][j] - A[i][i])/ (2*A[i][j]);
    t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    c = 1 / (sqrt(pow(t, 2) + 1));
    s = t * c;

    P = create_I_mat(n);
    if (P == NULL){
        return NULL;
    }
    P[i][j] = s;
    P[i][i] = c;
    P[j][j] = c;
    P[j][i] = (-1) * s;

    free(indices_max_abs);
    return P;
}

double calc_off(double **A, int n){
    int i,j;
    double off = 0;

    for ( i = 0; i < n; i++){
        for ( j = 0; j < n; j++){
            if (i != j){
                off = off + A[i][j] * A[i][j];
            }
        }
    }
    return off;
}

int does_conv(double **A, double **PTAP, int n){
    double off_A, off_PTAP;
    double epsilon = 0.00001;

    off_A = calc_off(A, n);
    off_PTAP = calc_off(PTAP, n); 

    if ((off_A - off_PTAP) <= epsilon){
        return 1;
    }else{
        return 0;
    }
}

double*** apply_Jacobi(double **A, int n){
    double **P, **P_T, **PTAP, **last_V, **V; 
    int conv_flag = 0;
    int iteration_counter = 0;
    int max_iter = 100;
    double ***return_arr;

    return_arr = (double ***) calloc(2, sizeof(double **)); /* will contain A and V to return */
    if (return_arr == NULL){
        return NULL;
    }

    last_V = create_I_mat(n);

    while (iteration_counter < max_iter && conv_flag != 1){
        /* calculate P mat */
        P = calc_P_mat(A,n);

        /* calculate P transpose */
        P_T = transpose(P ,n);
        if (P_T == NULL){
            return NULL;
        }
        PTAP = mat_mult(mat_mult(P_T, A, n), P, n);
        if (PTAP == NULL){
            return NULL;
        }
        last_V = V;
        V = mat_mult(last_V, P, n);
        if (V == NULL){
            return NULL;
        }

        conv_flag = does_conv(A, PTAP, n);
        A = PTAP; /* update A = A' */
        iteration_counter++;
    }

    /* free all matrixes */
    free_mat(last_V, n);
    free_mat(P, n);
    free_mat(P_T, n);
    free_mat(PTAP, n);
    return_arr[0] = A; /* The diagonal of the final A' is the eigenvalues of the original A */
    return_arr[1] = V; /* The eigenvectors of A */
    return return_arr;
}

int compare_vectors(const void *vector1, const void *vector2){
    /* sorting two vectors in a descending order by eigenvalues */
    const eigenvector_object *vec1 = vector1;
    const eigenvector_object *vec2 = vector2;

    if (vec1->eigenvalue > vec2->eigenvalue){
        return -1;
    }
    if (vec1->eigenvalue < vec2->eigenvalue){
        return 1;
    }
    return 0;

}


int get_k(double **A, double **V, int n){
    /* calculates k from eigenvectors and eigenvalues of lnorm */
    double delta;
    double max_delta = 0;
    int k;
    int i,j;

    eigenvector_object *eigenvectors;
    eigenvectors = (eigenvector_object*) calloc(n, sizeof(eigenvector_object));
    for ( i = 0; i < n; i++){
        eigenvectors[i].eigenvector = V[i];
        eigenvectors[i].eigenvalue = A[i][i];
    }
    
    qsort(eigenvectors, n, sizeof(eigenvector_object), compare_vectors);
    for ( j = 0; j < n - 1; j++){
        delta = fabs(eigenvectors[j].eigenvalue - eigenvectors[j+1].eigenvalue);
        if (delta > max_delta){
            max_delta = delta;
        }
    }
    k = (int) max_delta;
    return k;
}

int main(int argc, char** argv){
    /* if there are command-line arguments, they are interpered as filenames, and processed in order */
    
    FILE *ifp; /* input_file */
    
    int d; /* the dimention of the given points */
    int n;  /* the number of points, equals to the number of raws in input_file */
    char* input_filename;
    double *points_arr;
    double **points;
    double **W_mat, **D_mat, **lnorm_mat;
    double **eigenvalues , **eigenvectors;
    double ***eigen;
    n = 0;

    /* checking if the input is valid */
    if (argc == 3){ /* we include the program's name and there are 4 arguments if max_iter is provided */
        input_filename = argv[2];
    }else{
        printf("Invalid Input!");
        return 1;                                                                                                                                                  
    }

    ifp = fopen(input_filename,"r");
    if (ifp == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    /* saving the points in a 1-D array */
    d = get_dimention(ifp);
    n = get_points_num(ifp);
    points_arr = calloc(d * (n), sizeof(double));
    if (points_arr == NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    get_points_from_input(ifp,points_arr);
    fclose(ifp);

    /* convert 1-D array to matrix */
    points = convert_1D_arr_to_mat(points_arr, n, d);
    if (points == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    if (strcmp(argv[1],"wam") == 0){
        W_mat = create_W_mat(points, n, d);
        if (W_mat == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
        print_mat(W_mat, n, n);
    }
    
    else if (strcmp(argv[1],"ddg") == 0){
        W_mat = create_W_mat(points, n, d);
        if (W_mat == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
        D_mat = create_DD_mat(W_mat, n, d);
        if (D_mat == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
        print_mat(D_mat, n, n);
    }

    else if (strcmp(argv[1],"lnorm") == 0){
        lnorm_mat = lnorm_calc(points, n, d);
        if (lnorm_mat == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
        print_mat(lnorm_mat, n, n);
    }
    else if (strcmp(argv[1],"jacobi") == 0){
        eigen = apply_Jacobi(points, n);
        eigenvalues = eigen[0];
        eigenvectors = eigen[1];
        print_eigenvalues(eigenvalues, n, n);
        print_mat(eigenvectors, n, n);
    }
    else{
        printf("An Error Has Occurred");
        return 1;
    }
    
    free(input_filename);
    free(points_arr);
    free_mat(points, n);
    free_mat(W_mat, n);
    free_mat(D_mat, n);
    free_mat(lnorm_mat, n);
    free_mat(eigenvalues, n);
    free_mat(eigenvectors, n);

    return 0;
}


double* get_eigen_val(double **eigenValuesMat, int cnt){
    double* res;
    int i;

    res = (double*) calloc(cnt, sizeof(double));
    if (res == NULL)
    {
        return NULL;
    }
    for (i = 0; i < cnt; i++){
        res[i] = eigenValuesMat[i][i];
    }
    return res;
}



double** get_U_mat(eigenvector_object *eigen_pairs, int cnt, int clust_num){
    /* return U matrix containing getting the first k columns of V by eigen value order*/
    int i,j;
    double **uMat = create_empty_mat(cnt, clust_num);
    if (uMat == NULL){
        return NULL;
    }
    for(i = 0; i < cnt; i++){
        for(j = 0; j < clust_num; j++){
            uMat[i][j] = eigen_pairs[j].eigenvector[i];
        }
    }
    return uMat;
}


double get_row_norm(double *row_vec,int vecDim){
    /* calc the norm of a vector*/
    int i;
    double sum = 0;
    for(i = 0; i < vecDim; i++){
        sum += pow(row_vec[i], 2);
    }
    return sqrt(sum);
}


void normalize(double **uMat, int cnt, int clust_num){
    /* normalizing U matrix using get_row_norm function*/
    int i, j;
    double rowNorm;
    for(i = 0; i < cnt; i++){
        rowNorm = get_row_norm(uMat[i], clust_num); /*changed*/
        if (rowNorm != 0) {
            for(j = 0; j < clust_num; j++){
                uMat[i][j] /= rowNorm;
            }
        }
    }
}


double** get_T_mat(eigenvector_object *eigen_pairs, int cnt, int clust_num){
    double** tMat;
    tMat = get_U_mat(eigen_pairs, cnt, clust_num);
    if (tMat == NULL)
    {
        return NULL;
    }
    normalize(tMat, cnt, clust_num);
    return tMat;
}