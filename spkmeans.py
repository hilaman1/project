import sys
import numpy as np
import pandas as pd
import myspkmeans as spk



max_iter = 300


def get_CMD_input():
    # get output from user from cmd
    if (len(sys.argv) != 4):
        print("Invalid Input!")
        sys.exit(1)
    temp_k = sys.argv[1]
    goal = sys.argv[2]
    file_name = sys.argv[3]
    if (sys.argv[1]).isdigit():
        k = int(temp_k)
    else:
        print("Invalid Input!")
        sys.exit(1)
    return goal, k, file_name


def kmeans_pp(points, n, d, k):
    np.random.seed(0)
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
        if (np.sum(Dl) == 0):
            continue
        P_xl = Dl / np.sum(Dl)
        # choose a random centroid
        current_ind += 1
        rand_ind = np.random.choice(n, p=P_xl)
        indices_of_chosen_centroids[current_ind] = rand_ind
        initialized_centroids = np.append(initialized_centroids, [points[rand_ind]], axis=0)
    return initialized_centroids, indices_of_chosen_centroids

def print_mat(mat_as_list):
    # print matrix in the wanted format
    for row in mat_as_list:
         print(','.join([format(row[i], ".4f") for i in range(len(row))]))

def print_eigenval(mat_as_list):
    # print eigenvalues in the wanted format
    print(','.join([format(mat_as_list[i][i], ".4f") for i in range(len(mat_as_list))]))

def print_centroids(indices_of_chosen_centroids, final_centroids):
    for i in range(len(indices_of_chosen_centroids)):
        if i != (len(indices_of_chosen_centroids) - 1):
            print(indices_of_chosen_centroids[i], end=",")
        else:
            print(indices_of_chosen_centroids[i])
    for centroid in final_centroids.tolist():
        print(','.join([format(centroid[j], ".4f") for j in range(len(centroid))]))

goal, k, file_name = get_CMD_input()
if goal not in {"spk", "lnorm", "wam", "ddg", "jacobi"}:
    print("Invalid Input!")
    sys.exit(1)
if k < 0:
    print("Invalid Input!")
    sys.exit(1)
df = pd.read_csv(file_name, header=None)
n = df.shape[0]
d = df.shape[1]
points = df.to_numpy().tolist()
# implementations of the different goals
if goal == "spk":
    if k == 1 or k < 0:
        print("Invalid Input!")
        sys.exit(1)
    if k == 0:
        spk_output = spk.calc_T(points, n, d, k)  # call spkmeans C module
        if (spk_output == None):
            print("An Error Has Occurred")
            sys.exit(1)
        updated_k = len(spk_output[0])

    else:
        updated_k = k
        spk_output = spk.calc_T(points, n, d, updated_k)
        if (spk_output == None):
            print("An Error Has Occurred")
            sys.exit(1)

    initialized_centroids, indices_of_chosen_centroids = kmeans_pp(spk_output, n, updated_k, updated_k)
    centroid_after_kmeans = spk.kmeans(spk_output, initialized_centroids.tolist(), max_iter, n, updated_k,
                                   updated_k, 0)  # call spkmeans C module
    if centroid_after_kmeans is None:
        print("An Error Has Occurred")
        sys.exit(1)

    final_centroids = np.array(centroid_after_kmeans)
    final_centroids = np.round(final_centroids, 4)
    # print output
    print_centroids(indices_of_chosen_centroids, final_centroids)


if goal == "wam":
    wam_mat = spk.calc_W(points, n, d)  # call C module
    if wam_mat == None:
        print("An Error Has Occurred")
        sys.exit(1)
    print_mat(wam_mat)
if goal == "ddg":
    dd_mat = spk.calc_DDM(points, n, d)  # call C module
    if dd_mat == None:
        print("An Error Has Occurred")
        sys.exit(1)
    print_mat(dd_mat)
if goal == "lnorm":
    lnorm_mat = spk.calc_lnorm(points, n, d)  # call C module
    if lnorm_mat == None:
        print("An Error Has Occurred")
        sys.exit(1)
    print_mat(lnorm_mat)
if goal == "jacobi":
    # input is a symmetric matrix
    Jacobi_output = spk.apply_Jacobi_py(points, n, d)  # call C module
    if Jacobi_output == None:
        print("An Error Has Occurred")
        sys.exit(1)
    print_eigenval(Jacobi_output[0])
    print_mat(Jacobi_output[1])
