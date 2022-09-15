#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>
#include <string.h>

/* a class that saves each eigenvector with its eigenvalue */
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
    int i;

    i = 0;
    while (fscanf(ifp, "%lf,", &point) != EOF)
    {
        points_arr[i] = point;
        i++;
    }
}

double** create_empty_mat(int n, int d){
    /* creates an empty matrix from order nxd */
    double **mat; 
    int i;
    mat = (double**) calloc(n, sizeof(double*));
    if (mat == NULL){
        return NULL;
    }
    for (i = 0; i < n; i++){
        mat[i] = (double*) calloc(d, sizeof(double));
        if (mat[i] == NULL){
            return NULL;
        }
    }
    return mat;
}

double** convert_arr_to_mat(double *points_arr, int n, int d){
    int i,j;
    double **points_mat;
    points_mat = create_empty_mat(n,d);
    if (points_mat == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    for(i = 0; i < n ;i++){
        for(j = 0; j < d; j++){
            points_mat[i][j] = points_arr[i * d + j];
        }
    }
    return points_mat;
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

double** create_W_mat(double **points, int n, int d){
    /* return a mat W s.t W_mat[i][j] = Radial basis function kernel */
    double **W_mat;
    double delta = 0;
    double temp_val;
    int i,j, k;
    double sum = 0;

    W_mat = create_empty_mat(n,n);
    if (W_mat == NULL){
        return NULL;
    }
    for (i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if (i == j){
                W_mat[i][j] = 0;
            }else{
                sum = 0;
                for(k = 0; k < d ;k++){
                    delta = pow((points[i][k] - points[j][k]), 2);
                    sum += delta;
                }
                temp_val = ((-1) * sqrt(sum)) / 2;
                W_mat[i][j] = exp(temp_val);
            }
        }
    }
    return W_mat;
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

double** create_DD_mat(double **W_mat, int n){
    /* return a mat D s.t D_mat[i][j] = Diagonal Degree Matrix */
    double **D_mat;
    double sum;
    int i ,j;

    D_mat = create_empty_mat(n, n);

    if (D_mat == NULL){
        return NULL;
    }

    sum = 0.0;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            sum = sum + W_mat[i][j];
        }
        D_mat[i][i] = sum;
        sum = 0.0;
    }
    return D_mat;
}

double** DD_mat_for_lnorm(double **W_mat, int n){
    /* return a mat D s.t D_mat[i][j] = Diagonal Degree Matrix */
    double **D_mat;
    double sum;
    int i ,j;

    D_mat = create_empty_mat(n, n);

    if (D_mat == NULL){
        return NULL;
    }

    sum = 0.0;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            sum = sum + W_mat[i][j];
        }
        D_mat[i][i] = 1/sqrt(sum);
        sum = 0.0;
    }
    return D_mat;
}

double** mat_mult(double **mat1, double **mat2, int n){
    double **mult;
    int j;
    int k;
    int i;
    double sum;
    double temp;
    mult = create_empty_mat(n,n);
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            sum = 0;
            for (k = 0; k < n; k++){
                temp = mat1[i][k]*mat2[k][j];
                sum += temp;
            }
            mult[i][j]= sum;
        }
    }
    return mult;
}

void free_mat(double **mat, int n){
    int i;
    for(i = 0; i < n; i++){
        free(mat[i]);
    }
    free(mat);
}


double** lnorm_calc(double **points, int n, int d){
    int i, j;
    double **D_mat;
    double **first_mult, **second_mult;
    double **W_mat;
    double temp_val;
    
    W_mat = create_W_mat(points, n, d);
    if (W_mat == NULL){
        return NULL;
    }
    
    D_mat = DD_mat_for_lnorm(W_mat, n);
    if (D_mat == NULL){
        return NULL;
    }

    first_mult = mat_mult(D_mat, W_mat, n);
    if (first_mult == NULL){
        return NULL;
    }
    second_mult = mat_mult(first_mult, D_mat, n);
    if (second_mult == NULL){
        return NULL;
    }

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i == j){
                temp_val = second_mult[i][j];
                second_mult[i][j] = 1 - temp_val;
            }else{
                temp_val = second_mult[i][j];
                second_mult[i][j] = (-1) * temp_val;
            }
        }
    }

    free_mat(W_mat, n);
    free_mat(D_mat, n);
    free_mat(first_mult, n);
    return second_mult;
}

double** create_I_mat(int n){
    /* return an identity matrix of size n */
    double **I_mat;
    int i;
    I_mat = create_empty_mat(n,n);
    if (I_mat == NULL)
    {
        return NULL;
    }

    for (i = 0; i < n; i++){
        I_mat[i][i] = 1;
    }
    return I_mat;
}

int* get_max_not_in_diag(double **sym_mat, int n){
    /* calculates the Aij entery for pivoting and returns its index*/
    int *coordinates;
    int coord_i = 0;
    int coord_j = 0;

    double max_abs = 0;
    int i,j;

    coordinates = (int*) calloc(2, sizeof(int));
    if (coordinates == NULL){
        return NULL;
    }

    for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
            /* A is a symmetric matrix, check only above the diagonal */
            if ((i != j) && (fabs(sym_mat[i][j]) > fabs(max_abs))){
                max_abs = sym_mat[i][j];
                coord_i = i;
                coord_j = j;
            }
        }
    }

    coordinates[0] = coord_i;
    coordinates[1] = coord_j;

    return coordinates;
}

double sign(double num){
    /* sign function */
    if (num >= 0){
        return 1;
    }else{
        return -1;
    }
}

double** transpose(double **mat, int n){
    /* recieves a matrix dimantion and calculates the transpose of it */
    double **mat_T;
    double temp;
    int i,j;
    mat_T = create_empty_mat(n, n);
    if (mat_T == NULL){
        return NULL;
    }
    for(i = 0; i < n; i++){
        for( j = 0; j < n; j++){
            temp = mat[i][j];
            mat_T[j][i] = temp;
        }
    }
    return mat_T;
}

double** calc_P_mat(double **A, int n){
    double **P;
    int i,j; /* indices for Aij */
    double theta,t,c,s;
    int *max_abs_coordinates;
    
    max_abs_coordinates = get_max_not_in_diag(A, n); /* returns an array with 2 indices */
    if (max_abs_coordinates == NULL){
        return NULL;
    }
    i = max_abs_coordinates[0];
    j = max_abs_coordinates[1];

    theta = (A[j][j] - A[i][i])/ (2*A[i][j]);
    t = sign(theta) / ((fabs(theta) + sqrt(pow(theta, 2) + 1)));
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

    free(max_abs_coordinates);
    return P;
}

double calc_off(double **A, int n){
    int i,j;
    double off = 0;

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i != j){
                off += A[i][j] * A[i][j];
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

int is_it_diag(double** A, int n, int d){
    /* cheacks if a matrix is diagonal */
    int i;
    int j;

    for (i = 0; i < n; i++){
        for (j = 0; j < d; j++){
            if (i != j) {
                if (A[i][j] != 0){
                    return 0;
                }
            }

        }
    }
    return 1;
}

double*** apply_Jacobi(double **A, int n){
    double **P, **P_T, **last_V, **V, **A_new; 
    double **PTA, **PTAP;
    int conv_flag = 0;
    int iteration_counter = 0;
    int max_iter = 100;
    double ***return_arr;
    A_new = A;
    return_arr = (double***) calloc(2, sizeof(double**)); /* will contain A and V to return */
    if (return_arr == NULL){
        return NULL;
    }
    
    V = create_I_mat(n);


    if (is_it_diag(A, n, n) == 1){
        return_arr[0] = A_new; /* The diagonal of the final A' is the eigenvalues of the original A */
        return_arr[1] = V; /* The eigenvectors of A */
        return return_arr;
    }


    while (conv_flag != 1 && iteration_counter < max_iter){
        P = calc_P_mat(A_new, n);
        if (P == NULL){
            return NULL;
        }

        P_T = transpose(P ,n);
        if (P_T == NULL){
            return NULL;
        }

        PTA = mat_mult(P_T, A_new, n);
        if (PTA == NULL){
            return NULL;
        }

        PTAP = mat_mult(PTA, P, n);
        if (PTAP == NULL){
            return NULL;
        }

        last_V = V;
        V = mat_mult(last_V, P, n);
        if (V == NULL){
            return NULL;
        }
        conv_flag = does_conv(A_new, PTAP, n);

        /* free current matrixes */
        if (iteration_counter > 0){
            free_mat(A_new, n);
        }
        free_mat(last_V, n);
        free_mat(P, n);
        free_mat(P_T, n);
        free_mat(PTA, n);
        A_new = PTAP; /* update A = A' */
        iteration_counter++;
    }

    return_arr[0] = A_new; /* The diagonal of the final A' is the eigenvalues of the original A */
    return_arr[1] = V; /* The eigenvectors of A */
    return return_arr;
}

int compare_vectors(const void *vector1, const void *vector2){
    /* sorting two vectors in a descending order by eigenvalues */
    const eigenvector_object *vec1 = (eigenvector_object*) vector1;
    const eigenvector_object *vec2 = (eigenvector_object*)vector2;

    if (vec1->eigenvalue > vec2->eigenvalue){
        return -1;
    }
    if (vec1->eigenvalue < vec2->eigenvalue){
        return 1;
    }
    return 0;

}

eigenvector_object* create_eigenlist(double **A, double **V, int n){
    int i, j;
    double temp;
    eigenvector_object *eigenvectors;
    eigenvectors = (eigenvector_object*) calloc(n, sizeof(eigenvector_object));
    if (eigenvectors == NULL){
        return NULL;
    }
    for (i = 0; i < n; i++){
        eigenvectors[i].eigenvector = (double*) calloc(n, sizeof(double));
        if (eigenvectors[i].eigenvector == NULL){
            return NULL;
        }
        temp = A[i][i];
        eigenvectors[i].eigenvalue = temp;
        for (j = 0; j < n; j++){
            eigenvectors[i].eigenvector[j] = V[j][i];
        }
    }
    return eigenvectors;
}

int calc_k(eigenvector_object *eigenvectors, int n){
    /* calculates k from eigenvectors and eigenvalues of lnorm */
    double delta;
    double max_delta = 0;
    int i;
    int k = 1;

    for (i = 0; i < floor((n/2)); i++){
        delta = fabs(eigenvectors[i].eigenvalue - eigenvectors[i+1].eigenvalue);
        if (delta > max_delta){
            max_delta = delta;
            k = i;
        }
    }
    return k+1;
}

int main(int argc, char** argv){
    /* if there are command-line arguments, they are interpered as filenames, and processed in order */
    
    FILE *ifp; 
    
    int d; 
    int n;  
    char *input_filename;
    double *points_arr;
    double **points;
    double **W_mat, **D_mat, **lnorm_mat;
    double **eigenvalues, **eigenvectors;
    double ***eigen;

    /* checking if the input is valid */
    if (argc == 3){ 
        input_filename = argv[argc - 1];
    }else{
        printf("Invalid Input!");
        return 1;                                                                                                                                                  
    }
    if (strcmp(argv[1], "wam") != 0  && strcmp(argv[1], "ddg") != 0 && strcmp(argv[1], "lnorm") != 0 && strcmp(argv[1], "jacobi") != 0){
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
    points_arr = calloc(n*d, sizeof(double));
    if (points_arr == NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    get_points_from_input(ifp,points_arr);
    fclose(ifp);

    points = convert_arr_to_mat(points_arr,n,d);
 
    if (points == NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    /* all goals */
    if (strcmp(argv[1],"wam") == 0){
        W_mat = create_W_mat(points, n, d);
        if (W_mat == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
        print_mat(W_mat, n, n);
        free_mat(W_mat, n);
        free_mat(points,n);
    }
    if (strcmp(argv[1],"ddg") == 0){
        W_mat = create_W_mat(points, n, d);
        if (W_mat == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
        D_mat = create_DD_mat(W_mat, n);
        if (D_mat == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
        print_mat(D_mat, n, n);
        free_mat(D_mat, n);
        free_mat(W_mat, n);
        free_mat(points,n);
    }
    if (strcmp(argv[1],"lnorm") == 0){
        lnorm_mat = lnorm_calc(points, n, d);
        if (lnorm_mat == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
        print_mat(lnorm_mat, n, n);
        free_mat(lnorm_mat, n);
        free_mat(points,n);
    }
    if (strcmp(argv[1],"jacobi") == 0){
        eigen = apply_Jacobi(points, n);
        eigenvalues = eigen[0];
        eigenvectors = eigen[1];
        print_eigenvalues(eigenvalues, n, n);
        print_mat(eigenvectors, n, n);
        free(eigen);
        free_mat(eigenvalues, n);
        free_mat(eigenvectors, n);
        free_mat(points, n);
    }  
    free(points_arr);
    return 0;
}



double** calc_U_mat(eigenvector_object *eigenvectors_lst, int n, int k){
    /* returns U matrix containing the k vectors corresponding to the k biggest eigenvalues */
    int i,j;
    double **U_mat = create_empty_mat(n, k);
    if (U_mat == NULL){
        return NULL;
    }
    for(i = 0; i < n; i++){
        for(j = 0; j < k; j++){
            U_mat[i][j] = eigenvectors_lst[j].eigenvector[i];
        }
    }
    return U_mat;
}


double get_row_norm(double *eigenvector,int k){
    /* calc the norm of a vector*/
    int i;
    double sum = 0;
    double norm;
    for(i = 0; i < k; i++){
        sum += eigenvector[i] * eigenvector[i];
    }
    norm = sqrt(sum);
    return norm;
}


void normalize(double **U_mat, int n, int k){
    /* normalizing U matrix using get_row_norm function*/
    int i, j;
    double vec_norm;
    for(i = 0; i < n; i++){
        vec_norm = get_row_norm(U_mat[i], k); 
        if (vec_norm != 0) {
            for(j = 0; j < k; j++){
                U_mat[i][j] /= vec_norm;
            }
        }
    }
}


double** calc_T_mat(eigenvector_object *eigenvector_lst, int n, int k){
    /* calculates T_mat from U_mat */
    double **U_mat;

    U_mat = calc_U_mat(eigenvector_lst, n, k);
    if (U_mat == NULL){
        return NULL;
    }
    /* normalized U_mat = T_mat */
    normalize(U_mat, n, k);
    return U_mat;
}