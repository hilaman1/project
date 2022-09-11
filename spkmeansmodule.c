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


typedef struct
{
    double *centroid; /* the centroid of this cluster */
    int size; /* number of points in this cluster */
    double *points_sum; /* the sum of the points in this cluster - it's an array because we sum over each dimension */
} cluster;



double* convert_Py_to_C_array(PyObject *py_points, int n, int d){
    /* converts numpy python array to 1-D C array */
    PyObject *point, *index;
    double *C_points;
    int i,j;
    C_points = (double*) calloc(n * d, sizeof(double));
    /* get input points from python */
    for ( i = 0; i < n; i++) {
        point = PyList_GetItem(py_points, i); /* Return the object at position index in the list pointed to by list */
        if (!PyList_Check(point)){
            continue;
        }
        /* save py_points in C_points */
        for ( j = 0; j < d; j++) {
            index = PyList_GetItem(point, j);
            C_points[i*d + j] = PyFloat_AsDouble(index);
            if (PyErr_Occurred() && C_points[i*d + j]  == -1.0){
                return NULL;
            }
        }
    }
    return C_points;
}

PyObject* convert_C_to_Py_mat(double** c_mat, int mat_size){
    /* converts C mat to numpy python array*/
    PyObject *py_mat, *row;
    int i,j;
    py_mat =  PyList_New(mat_size);
    if (py_mat == NULL)
        return NULL;
    
    for ( i = 0; i < mat_size; i++){
        row = PyList_New(mat_size);
        if (row == NULL){
            return NULL;
        }
        for ( j = 0; j < mat_size; j++){
            PyList_SET_ITEM(row, j, Py_BuildValue("d", c_mat[i][j]));
        }
        PyList_SetItem(py_mat, i, Py_BuildValue("O", row));
    }
    return py_mat;
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

static double **kmeanspp(double *points, double **centroids, int n, int d, int k, int max_iter, double epsilon){
    /* implementation of kmeans */
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


PyObject* create_mat_tuple_py(double ***mat, int n){
    /* creates a pointer to 2 matrices to return to python */
    PyObject *tuple_mat, *mat1, *mat2;
    tuple_mat = PyList_New(2);
    mat1 = convert_C_to_Py_mat(mat[0], n);
    mat2 = convert_C_to_Py_mat(mat[1], n);
    PyList_SetItem(tuple_mat, 0, Py_BuildValue("O", mat1));
    PyList_SetItem(tuple_mat, 1, Py_BuildValue("O", mat2));
    return tuple_mat;
}

static PyObject* apply_Jacobi_py(PyObject* self, PyObject *args){
    /* implementation on Jacobi's algorithm */
    PyObject *py_points;
    PyObject *eigen_vals_and_vec_py;
    double *C_points;
    double ***eigen;
    double **points;
    int n, d;


    if (!PyArg_ParseTuple(args, "Oii", &py_points, &n, &d)){
    return NULL;
    }
    if (!PyList_Check(py_points)){
        return NULL;
    }
    C_points = convert_Py_to_C_array(py_points, n, d);
    if (C_points == NULL){
        return NULL;
    }
    points = convert_1D_arr_to_mat(C_points, n, d);
    if (points == NULL)
    {
        return NULL;
    }
    eigen = apply_Jacobi(points, n, d);
    if (eigen == NULL){
        return NULL;
    }
    eigen_vals_and_vec_py = create_mat_tuple_py(eigen, n);

    free(eigen);
    free(C_points);

    return eigen_vals_and_vec_py;

}


static PyObject* calc_DDM(PyObject* self, PyObject *args){
    PyObject *py_points, *py_D_mat_2return;
    double* C_points;
    double ** DD_mat, ** W_mat, **points; 
    int n, d;

    if (!PyArg_ParseTuple(args, "Oii", &py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(py_points)){
        return NULL;
    }
    
    C_points = convert_Py_to_C_array(py_points, n, d);
    if (C_points == NULL)
    {
        return NULL;
    }
    points = convert_1D_arr_to_mat(C_points, n, d);
    if (points == NULL)
    {
        return NULL;
    }

    W_mat = create_DD_mat(points, n,d);
    if (W_mat == NULL)
    {
        return NULL;
    }
    DD_mat = create_DD_mat(W_mat, n,d);
    if (DD_mat == NULL)
    {
        return NULL;
    }
    py_D_mat_2return = convert_C_to_Py_mat(DD_mat, n);
 
    free(C_points);
    free_mat(W_mat, n);
    free_mat(DD_mat, n);

    return py_D_mat_2return;
}


static PyObject* calc_lnorm(PyObject* self, PyObject *args){
    PyObject *py_points, *py_lnorm_2return;
    double *C_points;
    double ** lnorm_mat, **points; 
    int n, d;

    if (!PyArg_ParseTuple(args, "Oii", &py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(py_points)){
        return NULL;
    }
    
    C_points = convert_Py_to_C_array(py_points, n, d);
    if (C_points == NULL)
    {
        return NULL;
    }
    points = convert_1D_arr_to_mat(C_points, n, d);
    if (points == NULL)
    {
        return NULL;
    }

    lnorm_mat = lnorm_calc(points, n,d);
    if (lnorm_mat == NULL)
    {
        return NULL;
    }
    
    py_lnorm_2return = convert_C_to_Py_mat(lnorm_mat, n);
 
    free(C_points);
    free_mat(lnorm_mat, n);
    free_mat(C_points, n);
    return py_lnorm_2return;
}

static PyObject* calc_WAM(PyObject* self, PyObject *args){
    PyObject *py_points, *py_W_mat_2return;
    double *C_points;
    double **points;
    double **W_mat; /* The Weighted Adjacency Matrix */
    int n ,d;

    if (!PyArg_ParseTuple(args, "Oii", &py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(py_points)){
        return NULL;
    }
    
    C_points = convert_Py_to_C_array(py_points, n, d);
    if (C_points == NULL)
    {
        return NULL;
    }
    points = convert_1D_arr_to_mat(C_points, n, d);
    if (points == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    W_mat = create_W_mat(points, n, d);
    if (W_mat == NULL){
        return NULL;
    }
    
    py_W_mat_2return = convert_C_to_Py_mat(W_mat, n);
 
    free(C_points);
    free_mat(points, n);
    free_mat(W_mat, n);

    return py_W_mat_2return;
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
    double epsilon;


    casting2double_flag = 0;

    /* checking if the input is valid
        If the argument does not match the specified PyObject type, it will throw a TypeError
    */
    if (!PyArg_ParseTuple(args, "OOiiiisi", &py_points, &py_centroids, &max_iter, &n, &d, &k, &epsilon))
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
    free_mat(centroids, n);

    return py_centroinds2return;
}

static PyObject* calc_T(PyObject* self, PyObject *args)
{
    PyObject *py_points, *py_T_mat_2return;

    int n, d, k, i;
    double *C_points;
    double **W_mat, **DD_mat, **L_norm, **eigenvalus_mat, **eigenvectors_mat, **T_mat, **points;
    double*** eigen;
    eigenvector_object* eigen_pairs;
   
    if (!PyArg_ParseTuple(args, "Oiii", &py_points, &n, &d, &k))
        return NULL;
    if (!PyList_Check(py_points))
        return NULL;

    C_points= convert_Py_to_C_array(py_points, n, d);
    if (C_points== NULL)
    {
        return NULL;
    }
    points = convert_1D_arr_to_mat(C_points, n, d);
    if (points == NULL)
    {
        return NULL;
    }
    W_mat= create_W_mat(points, n, d);
    if (W_mat== NULL)
    {
        return NULL;

    }
    DD_mat= create_DD_mat(W_mat, n, n);
    if (DD_mat== NULL)
    {
        return NULL;

    }
    L_norm = lnorm_calc(DD_mat, W_mat, n);
    if (L_norm == NULL)
    {
        return NULL;

    }
    eigen = apply_Jacobi(L_norm, n, d);
    if (eigen == NULL)
    {
        return NULL;
    }
    eigenvalus_mat = eigen[0];
    eigenvectors_mat = eigen[1];

    if (k == 0){
        k = get_k(eigenvalus_mat, eigenvectors_mat, n);
    }

    T_mat = get_T_mat(eigen_pairs, n, k);
    if (T_mat == NULL)
    {
        return NULL;
    }

    py_T_mat_2return = convert_C_to_Py_mat(T_mat, n);

    free(C_points);
    free_mat(W_mat, n);
    free_mat(DD_mat, n);
    free_mat(L_norm, n);
    free(eigen);
    free_mat(eigenvalus_mat, n);
    free_mat(eigenvectors_mat, n);  
    free_mat(T_mat, n);
    free_mat(points, n);

    for (i = 0; i < n; i++){
        free(eigen_pairs[i].eigenvector);
    }
    free(eigen_pairs);
    return py_T_mat_2return;
}





/* define the module in order to be able to import from python */

/*
 * The module definition struct, which holds all information needed to create a module object.
 * There is usually only one statically initialized variable of this type for each module.
 */

static PyMethodDef kmeansppMethods[] = {
        {"apply_Jacobi_py", (PyCFunction) apply_Jacobi_py, METH_VARARGS, PyDoc_STR("aplying Jacobis algorithms and returns eigenvalues and eigenvectors of a given matrix")},
        {"calc_T", (PyCFunction) calc_T, METH_VARARGS, PyDoc_STR("calc T matrix for spkmeans algorithm implementation")},
        {"calc_WAM", (PyCFunction) calc_WAM, METH_VARARGS, PyDoc_STR("calculate weighted adjacency matrix")},
        {"calc_DDM", (PyCFunction) calc_DDM, METH_VARARGS, PyDoc_STR("calculate diagonal degree matrix")},
        {"calc_lnorm", (PyCFunction) calc_lnorm, METH_VARARGS, PyDoc_STR("calculate normalized graph laplacian")},
        {"kmeans", (PyCFunction) kmeans, METH_VARARGS, PyDoc_STR("apply kmeans algorithm")},
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