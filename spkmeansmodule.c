#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>
#include "spkmeans.h"


/* creating a class for creating clusters */
typedef struct
{
    double *centroid; /* the centroid of this cluster */
    int size; /* number of points in this cluster */
    double *points_sum; /* the sum of the points in this cluster - it's an array because we sum over each dimension */
} cluster;



double* convert_Py_to_C_mat(PyObject *py_points, int n, int d){
    double *C_points;
    C_points = (double*) calloc(n * d, sizeof(double));
    /* get input points from python */
    for (int i = 0; i < n; i++) {
        point = PyList_GetItem(py_points, i); /* Return the object at position index in the list pointed to by list */
        if (!PyList_Check(point)){
            continue;
        }
        /* save py_points in C_points */
        for (int j = 0; j < d; j++) {
            index = PyList_GetItem(point, j);
            points[i*d + j] = PyFloat_AsDouble(index);
            if (PyErr_Occurred() && points[d*i + j]  == -1.0){
                return NULL;
            }
        }
    }
    return C_points;
}

PyObject* convert_C_to_Py_mat(double** c_mat, int mat_size){
    PyObject *py_mat, *row;
    py_mat =  PyList_New(mat_size);
    if (py_mat == NULL)
        return NULL;
    
    for (int i = 0; i < mat_size; i++){
        row = PyList_New(mat_size);
        if (row == NULL){
            return NULL;
        }
        for (int j = 0; j < mat_size; j++){
            PyList_SET_ITEM(row, j, Py_BuildValue("d", c_mat[i][j]));
        }
        PyList_SetItem(py_mat, i, Py_BuildValue("O", row));
    }
    return py_mat;
}

void free_mat(double** mat, int row){
    for(int i = 0; i < row; i++){
        free(mat[i]);
    }
    free(row);
}

static double** create_empty_mat(int mat_order){
    /* creates an empty matrix from order nxn */
    double **mat; 
    int mat_order;

    mat = calloc(mat_order, sizeof(double*));
    if (W == NULL){
        return NULL;
    }
    for (int i = 0; i< mat_order; i++){
        mat[i] = (double*) calloc(mat_order, sizeof(double));
        if (mat[i] == NULL)
        {
            return NULL;
        }
    }
    return mat;
}

static double L2_norm(double *x1, double *x2, int d){
    /* return L2 norm of x1 and x2 */
    double delta = 0;
    double current_coor;
    for(int i = 0; i < d ;i++)
    {
        current_coor = pow((x1[i] - x2[i]), 2);
        delta += delta;
    }
    return delta;
}



static double** create_W_mat(double *C_points, int n, int d){
    /* return a mat W s.t W_mat[i][j] = Radial basis function kernel */
    double **W;
    double delta;
    double RBF;

    W = create_empty_mat(n);
    if (W == NULL){
        return NULL;
    }
    for (int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if (i == j){
                W[i][j] = 0;
            }
            else{
                delta = L2_norm(C_points + (i * d), C_points + (j * d), d);
                RBF = -1 * (sqrt(delta) / 2);
                W[i][j] = exp(val);
            }
        }
    }
    return W;
}

static double** create_DD_mat(double **W_mat, int n, int d){
    /* return a mat D s.t D_mat[i][j] = Diagonal Degree Matrix */
    double **D_mat;
    double d_i;

    D_mat = create_empty_mat(n);
    if (D_mat == NULL){
        return NULL;
    }
    for (int i = 0; i < n; i++){
        d_i = 0;
        for (int j = 0; j < d; j++){
            d_i += W_mat[i][j];
        }
        D_mat[i][i] = sqrt(d_i);
    }
    return D_mat;
}

static double** create_I_mat(int n){
    /* return an identity matrix of size n */
    double ** I_mat;
    I_mat = create_empty_mat(n);
    if (I_mat == NULL)
    {
        return NULL;
    }
    for (int i = 0; i < n; i++){
        I_mat[i][i] = 1;
    }
    return I_mat;
}

tatic double** calc_mat_product(double ** mat1, double** mat2, int n){
    /* return a mat which is the product of mat1, mat2 */
    double** mat_product;
    double sum;

    mat_product = create_empty_mat(n);
    if (mat_product == NULL)
    {
        return NULL;
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            sum = 0;
            for (int l = 0; l < n; l++){
                sum += mat1[i][l] * mat2[l][j];
            }
            mat_product[i][j] = sum;
        }
    }
    return mat_product;
}

static double** calc_DWD_product(double ** D_mat, double** W_mat, int n){
    /* return a mat which is the product of D_mat^(1/2),W_mat,D_mat^(1/2)  */
    double ** DW_product, **DWD_product;
    DW_product = calc_mat_product(D_mat, W_mat, n);
    if (DW_product == NULL)
    {
        return NULL;
    }
    
    DWD_product = calc_mat_product(DW_product, D_mat, n);
    if (DWD_product == NULL)
    {
        return NULL;
    }
    free_mat(DW_product, n);
    return DWD_product;
}

static double** calc_mat_sub(double ** mat1, double** mat2, int n){
    /* return a mat which is the subtraction mat1, mat2  */
    double ** sub_mat;
    sub_mat = create_empty_mat(n);
    for(int i = 0; i < n; i++ ){
        for(int j = 0; j < n; j++){
            sub_mat[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return sub_mat;
}


static double** create_lnorm_mat(double ** D_mat, double** W_mat, int n){
    /* return a mat Lnorm s.t Lnorm_mat[i][j] = Normalized Graph Laplacian */
    double ** I_mat, ** lnorm, **DWD_product;

    I_mat = create_I_mat(n);
    if (I_mat == NULL)
    {
        return NULL;
    }
    DWD_product = calc_DWD_product(D_mat, W_mat, n);
    if (DWD_product == NULL)
    {
        return NULL;
    }
    lnorm = calc_mat_sub(I_mat, DWD_product, n);

    free_mat(I_mat, n);
    free_mat(DWD_product, n);

    return lnorm;
    
}

static double calc_Delta(double *current_centroid, double *point, int d)
/* a function the calculates the distance between a given centroid (mu) of a cluster
and a given point (x) */
{
    double Delta, current_coord_dist;
    int i;
    Delta = (double) 0;

    for(i = 0; i < d ;i++)
    {
        /*there's the pow() function that will raise a number to any power
        But it would be much more efficient to just multiply the number with itself*/
        current_coord_dist = (double) (current_centroid[i] - point[i])*(current_centroid[i] - point[i]);
        /* instead of returning a Delta dector and calculating the rooted and summing the squared distances,
        we'll just return the summation of the squared distances of all the coordinates*/
        Delta += current_coord_dist;
    }
    return Delta;
}

static int find_closest_cluster(double* point, cluster* clusters, int k, int d) {
    /* a function that finds the cluster with closest centroid to a given point */
    double min_Delta = DBL_MAX;
    double Delta;
    int cluster_with_min_dist = -1;
    double *current_centroid;
    int i;

    for (i = 0 ; i < k ; i++) {
        current_centroid = clusters[i].centroid;
        Delta = calc_Delta(current_centroid, point, d);
        if (Delta < min_Delta) {
            min_Delta = Delta;
            cluster_with_min_dist = i;
        }
    }
    return cluster_with_min_dist;
}

static int find_best_cluster(cluster* clusters, double* points, int k, int points_num, int d) {
    /* a function that gets the clusters and chooses the one with the shortest distance from the point */
    int i,j;
    double *point;
    int best_cluster;

    point = (double*) calloc (d, sizeof(double)); /* create a point of dimension d */

    if (point == NULL){
        printf("An Error Has Occurred");
        return 0;
    }
    for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++){
            clusters[i].points_sum[j] = (double) 0;
        }
        clusters[i].size = 0;
    }
    for (i = 0; i < points_num; i++){
        for (j = 0; j < d; j++){
            point[j] = points[i * d + j];
        }
        best_cluster = find_closest_cluster(point, clusters, k, d);

        /* adding the point to the best cluster */
        clusters[best_cluster].size += 1 ;
        for (j = 0; j < d; j++){
            clusters[best_cluster].points_sum[j] += point[j] ;
        }
    }
    free(point);
    return 1;
}


/* a function that returns 0 if the centroids convegre, returns 1 elsewise */
static int does_converge(cluster *clusters, double *centroids, int k, int d, double epsilon){
    int i,j;
    double Delta;

    for (i = 0; i < k; i++){
        for (j = 0; j < d; j ++){
            Delta = fabs(centroids[i * d + j] - clusters[i].centroid[j]); /* fabs returns the absolute value */
            if (Delta == epsilon || Delta > epsilon){
                return 0;
            }
        }
    }
    return 1;
}

static void update_centroids(cluster* clusters, int k, int d) {
    /* the function updates the centroid to the mean value of each centroid */
    int i,j;

    for (i = 0; i < k; i++){
        for (j = 0 ; j < d; j++){
            if(clusters[i].size==0){
                clusters[i].centroid[j] = 0;
            }else{
                clusters[i].centroid[j] = (double) (clusters[i].points_sum[j] / clusters[i].size);
            }
        }
    }

    /* no need to return anything, it updates the clusters array :) */
}

static double **kmeanspp(double *points, double **centroids, int n, int d, int k, int max_iter, double epsilon)
{
    int converge_flag, best_cluster, iteration_counter;
    int i, j;
    double *current_centroids;
    cluster *clusters;

    converge_flag = 0;
    iteration_counter = 0;

    /* allocate space for updated centroids */
    current_centroids = (double*) calloc(k * d, sizeof(double));
    if(current_centroids == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    clusters = (cluster*) calloc(k, sizeof(cluster)); /* each cluster is an object that contains the points in the cluster, their number and their sum */
    if(clusters == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < k; i++) {
        clusters[i].centroid = (double*) calloc(d, sizeof(double)); /* create a centroid point of dimension d */
        clusters[i].points_sum = (double*) calloc(d, sizeof(double)); /* initialize points sum in each coordinate */
        if((clusters[i].points_sum == NULL) || (clusters[i].centroid == NULL)){
            printf("An Error Has Occurred");
            return NULL;
        }
        for (j = 0; j < d; j ++){
            clusters[i].centroid[j] = centroids[i][j];
            clusters[i].points_sum[j] = centroids[i][j];
        }
        clusters[i].size = 1;
    }
    while (converge_flag == 0 && iteration_counter < max_iter)
    {
        best_cluster = find_best_cluster(clusters, points, k, n, d);
        if (best_cluster == 0){
            return NULL;
        }
        for (i = 0; i < k; i ++){
            for (j = 0; j < d; j++){
                current_centroids[i * d + j] = clusters[i].centroid[j];
            }
        }
        /* if the algorithm doesn't converge,
            we want to update the centroid until it will return all the centroids
            with an error smaller than epsilon */
        update_centroids(clusters, k, d);
        converge_flag = does_converge(clusters, current_centroids, k, d, epsilon);
        /* converge_flag = 0 if the centroids convegre, returns 1 elsewise */
        iteration_counter++;
    }
    for(i = 0; i < k; i++)
    {
        for(j = 0; j < d; j++)
        {
            centroids[i][j] = clusters[i].centroid[j];
        }
    }
    for(i = 0; i < k; i++){
        free(clusters[i].centroid);
        free(clusters[i].points_sum);
    }
    free(clusters);
    free(current_centroids);
    return centroids;
}
static PyObject* calc_lnorm(PyObject* self, PyObject *args){
    PyObject *py_points, *py_points2return, *point, *index, row;
    double* C_points;
    double ** DD_mat, ** W_mat, **lnorm_mat; /* The Weighted Adjacency Matrix */

    if (!PyArg_ParseTuple(args, "Oii", &py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(py_points)){
        return NULL;
    }
    
    C_points = convert_Py_to_C_mat(py_points, n, d);
    if (C_points == NULL)
    {
        return NULL;
    }
    
    W_mat = create_w_mat(C_points, n);
    if (W_mat == NULL){
        return NULL;
    }

    DD_mat = create_DD_mat(W_mat, n);
    if (DD_mat == NULL){
        return NULL;
    }

    lnorm_mat = create_lnorm_mat(W_mat, n);
    if (lnorm_mat == NULL){
        return NULL;
    }
    
    py_points2return = convert_C_to_Py_mat(lnorm_mat, n);
 
    free(C_points);
    free_mat(W_mat, n);
    free_mat(DD_mat, n);
    free_mat(lnorm_mat, n);

    return py_points2return;
}
double** calc_off(double** A, int n){
    double off;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < j; i++){
            off = off + A[i+j*n] * A[i+j*n];
        }
    }
    return off;
}
int* get_max_not_in_diag(double** sym_mat, int n){
    /* calculates the Aij entery for pivoting */
    int* indexes_to_return;
    int i = 0;
    int j = 1;
    double current = 0;

    indexes_to_return = (int*) calloc(2, sizeof(int));

    for(int k = 0; k < vicCnt; k++){
        for(l = k; l < vicCnt; l++){
            if ((k != l) && (fabs(sym_mat[k][l]) > fabs(current))){
                tmp = mat[k][l];
                i = k;
                j = l;
            }
        }
    }

    indexes_to_return[0] = i;
    indexes_to_return[1] = j;

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
    transpose = create_empty_mat(n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            transpose[j][i] = mat[i][j];
        }
    }
    return transpose;
}

static PyObject* apply_Jacobi(PyObject* self, PyObject *args){
    PyObject *py_points, *py_points2return, *point;
    double **P, **P_T, **A, **PTA, **PTAP, **last_V, **V; 
    double theta,t,c,s;
    int i,j; /* indices for Aij */
    int does_converge = 0;
    int off_A, off_Anew;
    double epsilon = 0.00001;
    int iteration_counter = 0;
    int max_iter = 100;
    double ***return_arr;


    return_arr = (double ***) calloc(2, sizeof(double **)); /* will contain A and V to return */
    if (return_arr == NULL){
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "Oii", &py_points, &n, &d)){
    return NULL;
    }
    if (!PyList_Check(py_points)){
        return NULL;
    }
    A = convert_Py_to_C_mat(py_points, n, d);
    if (A == NULL){
        return NULL;
    }

    last_V = create_I_mat(n);

    while (iteration_counter < max_iter && does_converge != 1){
        /* calculate P mat */
        indices = get_max_not_in_diag(A, n); /* returns an array with 2 indices */
        i = indices[0];
        j = indices[1];

        theta = (A[j][j] - A[i][i])/ (2*A[i][j]);
        t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
        c = 1 / (sqrt(pow(t, 2) + 1));
        s = t * c;

        P = create_I_mat(n);
        P[i][j] = s;
        P[i][i] = c;
        P[j][j] = c;
        P[j][i] = (-1) * s;

        free(indices);

        /* calculate P transpose */
        P_T = transpose(P ,n);

        PTA = calc_mat_product(P_T, A, n);
        if (PTA == NULL){
            return NULL;
        }

        PTAP = calc_mat_product(PTA, P, n);
        if (PTAP == NULL){
            return NULL;
        }
        last_V = V;
        V = calc_mat_product(last_V, P, n);
        if (V == NULL){
            return NULL;
        }

        /* check if the algorithm converges */
        off_A = calc_off(A);
        off_Anew = calc_off(PTAP); 
        if ((off_A - off_Anew) < epsilon){
            does_converge = 1;
        }

        A = PTAP; /* update A = A' */
        iteration_counter++;
    }

    /* free all matrixes */
    free_mat(last_V);
    free_mat(P);
    free_mat(P_T);
    free_mat(PTA);
    free_mat(PTAP);

    return_arr[0] = A; /* The diagonal of the final A' is the eigenvalues of the original A */
    return_arr[1] = V; /* The eigenvectors of A */
    return return_arr;
}


static PyObject* calc_DDM(PyObject* self, PyObject *args){
    PyObject *py_points, *py_points2return, *point;
    double* C_points;
    double ** DD_mat, ** W_mat; /* The Weighted Adjacency Matrix */

    if (!PyArg_ParseTuple(args, "Oii", &py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(py_points)){
        return NULL;
    }
    
    C_points = convert_Py_to_C_mat(py_points, n, d);
    if (C_points == NULL)
    {
        return NULL;
    }
    
    W_mat = create_w_mat(C_points, n);
    if (W_mat == NULL){
        return NULL;
    }

    DD_mat = create_DD_mat(W_mat, n);
    if (DD_mat == NULL){
        return NULL;
    }
    
    py_points2return = convert_C_to_Py_mat(DD_mat, n);
 
    free(C_points);
    free_mat(W_mat, n);
    free_mat(DD_mat, n);

    return py_points2return;
}

static PyObject* calc_WAM(PyObject* self, PyObject *args){
    PyObject *py_points, *py_points2return, *point;
    double* C_points;
    double ** W_mat; /* The Weighted Adjacency Matrix */

    if (!PyArg_ParseTuple(args, "Oii", &py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(py_points)){
        return NULL;
    }
    
    C_points = convert_Py_to_C_mat(py_points, n, d);
    if (C_points == NULL)
    {
        return NULL;
    }

    W_mat = create_w_mat(C_points, n);
    if (W_mat == NULL){
        return NULL;
    }
    
    py_points2return = convert_C_to_Py_mat(W_mat, n);
 
    free(C_points);
    free_mat(W_mat, n);

    return py_points2return;
}
    

static PyObject* kmeans(PyObject* self, PyObject *args){
    PyObject *py_points, *py_centroids, *point, *centroid, *index, *py_centroinds2return, *py_centroid;
    int d; /* the dimention of the given points */
    int i, j; /* indices for iterating centroids */
    int k; /* number of centroids to return */
    int n; /* number of points to cluster; */
    int max_iter; /* number of max iteration givenas input to python */
    double *points; /* points to cluster */
    double **centroids, **centroids_after_kmeans;
    int casting2double_flag; /* flag to check if casting failed */



    casting2double_flag = 0;

    /* checking if the input is valid */
    /*
    The PyArg_ParseTuple function allows you to cast directly to a Python object subtype using the format string "O!"
     (notice-this is different than just plain "O").

     "OOiiiisi" -> the types of objects we expect: 2xobjects, 4xint, 1xstring

     If the argument does not match the specified PyObject type, it will throw a TypeError
    */
    if (!PyArg_ParseTuple(args, "OOiiiisi", &py_points, &py_centroids, &max_iter, &n, &d, &k, &goal, &goal_len))
        return NULL;
    if (!PyList_Check(py_centroids)) /* Return true if py_centroids is a list object or an instance of a subtype of the list type. This function always succeeds. */
        return NULL;
    if (!PyList_Check(py_points))
        return NULL;




    /* allocate space for points */
    points = (double*) calloc(d*n, sizeof(double));
    if(points == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }

    /* get input points from python */
    for (i = 0; i < n; i++) {
        point = PyList_GetItem(py_points, i); /* Return the object at position index in the list pointed to by list */
        if (!PyList_Check(point)){
            continue;
        }
        /* save py_points in points */
        for (j = 0; j < d; j++) {
            index = PyList_GetItem(point, j);
            points[i*d + j] = PyFloat_AsDouble(index);
            /* Return a C double representation of the contents of pyfloat. If pyfloat is not a Python floating point object but has a __float__() method,
             this method will first be called to convert pyfloat into a float. If __float__() is not defined then it falls back to __index__().
              This method returns -1.0 upon failure, so one should call PyErr_Occurred() to check for errors. */
            if (PyErr_Occurred() && points[d*i + j]  == -1.0){
                casting2double_flag = 1;
                break;
            }
        }
        if (casting2double_flag == 1){
            break;
        }
    }

    /* allocate 2-D array space */
    centroids = (double **)calloc(k, sizeof(double *));
    if(centroids==NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    /* initialize */
    for (i = 0; i < k; i++){
        centroids[i] = (double *)calloc(d, sizeof(double));
        if(centroids[i] == NULL){
            printf("An Error Has Occurred");
            return NULL;
        }
    }

    /* get first k centroids from python */
    for (i = 0; i < k; i++) {
        centroid = PyList_GetItem(py_centroids, i);
        if (!PyList_Check(centroid)){
            continue;
        }
        for (j = 0; j < d; j++) {
            index = PyList_GetItem(centroid, j);
            centroids[i][j] = PyFloat_AsDouble(index);
            /* Return a C double representation of the contents of pyfloat. If pyfloat is not a Python floating point object but has a __float__() method,
             this method will first be called to convert pyfloat into a float. If __float__() is not defined then it falls back to __index__().
              This method returns -1.0 upon failure, so one should call PyErr_Occurred() to check for errors. */
            if (PyErr_Occurred() && centroids[i][j]  == -1.0){
                casting2double_flag = 1; /*casting failed */
                break;
            }
        }
        if (casting2double_flag == 1){
            break;
        }
    }


    /* apply kmeans on centroids */
    centroids_after_kmeans = NULL;
    if (casting2double_flag == 0){
        centroids_after_kmeans = kmeanspp(points, centroids, n, d, k, max_iter, epsilon);
        if (centroids_after_kmeans == NULL)
        {
            return NULL;
        }
    }

    py_centroinds2return = PyList_New(k); /* Return a new list of length k on success, or NULL on failure. */
    if (py_centroinds2return == NULL)
        return NULL;
    for (i = 0; i < k; i++)
    {
        py_centroid = PyList_New(d); /* Return a new list of length d on success, or NULL on failure. */
        if (py_centroid == NULL)
            return NULL;
        for (j = 0; j < d; j++)
        {
            PyList_SET_ITEM(py_centroid, j, Py_BuildValue("d", centroids_after_kmeans[i][j]));
            /* Macro form of PyList_SetItem() without error checking. This is normally only used to fill in new lists where there is no previous content.
            Py_BuildValue - Create a new value based on a format string similar to those accepted by the PyArg_Parse*.
            this line saves the centroid as a double value as python object */
        }
        PyList_SetItem(py_centroinds2return, i, Py_BuildValue("O", py_centroid)); /* "O" -> make sure it's an object */
    }
    free(points);
    free_mat(centroids, n)

    return py_centroinds2return;
}


/* define the module in order to be able to import from python */

/*
 * The module definition struct, which holds all information needed to create a module object.
 * There is usually only one statically initialized variable of this type for each module.
 */

static PyMethodDef kmeansppMethods[] = {
        {"apply_Jacobi", (PyCFunction) apply_Jacobi, METH_VARARGS, PyDoc_STR("aplying Jacobis algorithms and returns eigenvalues and eigenvectors of a given matrix")}
        {"calc_WAM", (PyCFunction) calc_WAM, METH_VARARGS, PyDoc_STR("calculate weighted adjacency matrix")},
        {"calc_DDM", (PyCFunction) calc_DDM, METH_VARARGS, PyDoc_STR("calculate diagonal degree matrix")},
        {"calc_lnorm", (PyCFunction) calc_lnorm, METH_VARARGS, PyDoc_STR("calculate normalized graph laplacian")},
        {"kmeans", (PyCFunction) kmeans, METH_VARARGS, PyDoc_STR("apply kmeans algorithm")}
        {NULL, NULL, 0, NULL}};
/* This function must be registered with the interpreter using the METH_VARARGS flag;
 this is described in section The Module's Method Table and Initialization Function.
 PyDoc_STR - The docstring for the function */


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT, "myspkmeans", NULL, -1, kmeansppMethods};


/* When the Python program imports module "myspkmeans" for the first time, PyInit_myspkmeans() is called.
    The initialization function must be named PyInit_name(), where name is the name of the module,
     and should be the only non-static item defined in the module file. */
PyMODINIT_FUNC
PyInit_myspkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}