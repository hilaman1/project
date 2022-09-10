import sys
import numpy as np
import pandas as pd
import spkmeans as spk

np.random.seed(0)
max_iter = 300


def get_CMD_input():
    if (len(sys.argv) < 3) or (len(sys.argv) > 4):
        print("Invalid Input!")
        sys.exit(1)
    if len(sys.argv) == 4:
        temp_k = sys.argv[1]
        goal = sys.argv[2]
        file_name = sys.argv[3]
        if (sys.argv[1]).isdigit():
            k = int(temp_k)
            if k < 1:
                print("Invalid input!")
                sys.exit(1)
    elif len(sys.argv) == 3:
        goal = sys.argv[1]
        file_name = sys.argv[3]
        k = 0
    return goal, k, file_name


def kmeans_pp(points, n, d, k):
    start_ind = np.random.choice(n)
    indices_of_chosen_centroids = [0 for i in range(k)]
    indices_of_chosen_centroids[0] = start_ind
    initialized_centroids = np.zeros((0, d))
    initialized_centroids = np.append(initialized_centroids, [points[start_ind]], axis=0)
    current_ind = 0
    # Dl[i] will hold the closest centroid for each x_i
    temp_vec = np.power(np.subtract(points, initialized_centroids[current_ind]), 2)
    dist_mat = np.sum(temp_vec, axis=1).reshape(1, n)
    Dl = np.min(dist_mat, axis=0)
    P_xl = Dl / np.sum(Dl)
    # choose a random centroid
    current_ind += 1
    rand_ind = np.random.choice(n, p=P_xl)
    indices_of_chosen_centroids[current_ind] = rand_ind
    initialized_centroids = np.append(initialized_centroids, [points[rand_ind]], axis=0)
    while current_ind < k - 1:
        # row[i] is the distances between all points and centroid[i]
        temp_vec = np.power(np.subtract(points, initialized_centroids[current_ind]), 2)
        temp_mat = np.sum(temp_vec, axis=1).reshape(1, n)
        dist_mat = np.concatenate((dist_mat, temp_mat), axis=0)
        Dl = np.min(dist_mat, axis=0)
        P_xl = Dl / np.sum(Dl)
        # choose a random centroid
        current_ind += 1
        rand_ind = np.random.choice(n, p=P_xl)
        indices_of_chosen_centroids[current_ind] = rand_ind
        initialized_centroids = np.append(initialized_centroids, [points[rand_ind]], axis=0)

    return initialized_centroids, indices_of_chosen_centroids


def main():
    goal, k, file_name = get_CMD_input()
    df = pd.read_csv(file_name, header=None)
    df = df.sort_values([0])
    points = df.to_numpy()
    n = df.shape[0]
    d = df.shape[1]
    if k > n:
        print("Invalid Input!")
        sys.exit(1)
    initialized_centroids, indices_of_chosen_centroids = kmeans_pp(points, n, d, k)
    points2cluster = df.values.tolist()
    # All implementations of the different goals
    # must be performed by calling the C extension
    spk_module_output = spk.kmeans(points2cluster, initialized_centroids.tolist(), max_iter, n, d, k)  # call spkmeans C module
    if spk_module_output is None:
        print("An Error Has Occurred")
        sys.exit(1)

    if goal == 'spk':
        final_centroids = np.array(spk_module_output)
        final_centroids = np.round(final_centroids, 4)
        for i in range(len(indices_of_chosen_centroids)):
            if i != (len(indices_of_chosen_centroids) - 1):
                print(indices_of_chosen_centroids[i], end=",")
            else:
                print(indices_of_chosen_centroids[i])
        for centroid in final_centroids.tolist():
            print(','.join([format(centroid[j], ".4f") for j in range(len(centroid))]))
    elif goal == 'Jacobi':
        eigenvalues = np.array(spk_module_output)[0]
        eigenvectors = np.array(spk_module_output)[1]
        eigenvectors_mat = np.zeros(eigenvectors[0].shape)
        print(eigenvalues)
        for eigenvector in eigenvectors:  # TODO change that when write C module
            np.append(eigenvectors_mat, eigenvector)
        print(eigenvectors_mat)
    elif goal == 'wam':
        wam_mat = spk.calc_WAM(points2cluster, n, d)  # call spkmeans C module
    elif goal == 'ddg':
        dd_mat = spk.calc_DDM(points2cluster, n, d)  # call spkmeans C module
    elif goal == 'lnorm':
        lnorm_mat = spk.calc_lnorm(points2cluster, n, d)  # call spkmeans C module
    elif goal == 'jacobi':
        # input is a symmetric matrix
        Jacobi_output = spk.apply_Jacobi(points2cluster, n, d)  # call spkmeans C module




main()
