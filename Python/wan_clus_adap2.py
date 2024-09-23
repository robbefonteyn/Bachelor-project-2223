import math

from mean_misval import mean_misval
from dist_misval import dist_misval
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
import plotly.express as px
import sys

if not sys.warnoptions:
    import warnings

    warnings.simplefilter("error")


# The plot generation does not work all that well

def wan_clus_adap2(filename, min_number_genes_per_sphere=2, rk_prelim=None, max_iterations=500, div=30, transpose=False,
                   cluster_assignment_array=np.empty(shape=(0, 0)), dimension=None,
                   dataset_original=np.empty(shape=(0, 0)),
                   quality_criterion=2, plot=False):
    dataset = np.genfromtxt(filename, delimiter=",")

    if transpose:
        dataset = dataset.transpose()

    # samples = Number of columns in the dataset, variable is used to set up the cluster_assignment_array if this is
    # not provided by the user
    samples = dataset.shape[1]

    # Calculate the rk_prelim value if the user has not specified one. Formula can be found in the paper
    if not rk_prelim:
        vector_components = dataset.shape[0]
        rk_prelim = math.sqrt(vector_components - 1) / 2

    if rk_prelim <= 0:
        return print("rk_prelim should be greater than 0!")

    # Variables that dictate the maximum iterations the program will go through (standard = 500) and div is a
    # variable used to calculate the delta_rad variable later on. This will determine when the program will stop
    # looping (default = 30)

    # Make cluster_assignment_array an array of zeros with with dimensions (1, samples) if this is not a given variable
    if not cluster_assignment_array.any():
        cluster_assignment_array = np.zeros(samples)

    # Changes the dimension (or the number of rows) in the dataset if this is given (default is None)
    if dimension:
        dataset = dataset[0:dimension, :]

    # cluster should start at max value of cluster_assignment_array
    # cluster_means is an empty array where the means of each cluster will be appended to
    # convergence and toofewpoints variables are just used for the while and for loops that follow
    cluster = int(np.max(cluster_assignment_array))
    cluster_means = np.empty(shape=(dataset.shape[0], 1))
    convergence = 1
    toofewpoints = 0

    while convergence == 1 and toofewpoints < 5:
        # index_of_zeros_dataset is an array with the indexes of each 0 present in the cluster_assignment_array
        # data_not_in_cluster is a submatrix of the dataset that only has the data of the genes that are still
        # not assigned to a cluster
        data_not_in_cluster = dataset[:, cluster_assignment_array == 0]
        index_of_zeros_dataset = np.where(cluster_assignment_array == 0)
        if not data_not_in_cluster.any():
            break

        # Calculates mean of the array data_not_in_cluster with the MODE set to the Quality criterion
        mean1 = mean_misval(data_not_in_cluster, quality_criterion)

        # calculation of starting radius the algorithm uses by looking at the maximum distance from mean1
        # the radius will be lowered by the delta_radius variable each time it iterates without finding convergence
        radius = max(dist_misval(data_not_in_cluster, mean1, quality_criterion))
        # print(radius)
        delta_radius = (radius - rk_prelim) / div
        radius -= delta_radius

        # genes_in_sphere is the amount of genes in the sphere if mean1 is used as center
        genes_in_sphere = dist_misval(data_not_in_cluster, mean1, quality_criterion) < radius
        convergence = 0
        tel = 0

        for iteration in range(max_iterations):
            if iteration == 0 and not genes_in_sphere.any():
                break
            iteration += 1
            # Calculate the new mean (mean2) using the new genes in sphere
            mean2 = mean_misval(data_not_in_cluster[:, genes_in_sphere], quality_criterion)

            # Statement used to determine wether the algorithm is should still lower the radius variable or should
            # just assign rk_prelim (given variable) to radius
            if iteration < div + tel:
                radius -= delta_radius
            else:
                radius = rk_prelim

            # genes_in_sphere is calculated again, but now is the amount of genes in the sphere if mean2 is used as
            # center
            genes_in_sphere = dist_misval(data_not_in_cluster, mean2, quality_criterion) < radius

            # if this is empty (meaning the radius is too small, tel is increased and the delta_radius is again added
            # to the radius and the genes_in_sphere is calculated again using mean2 as center
            if not genes_in_sphere.any() and iteration < div + tel:
                tel += 1
                radius += delta_radius
                genes_in_sphere = dist_misval(data_not_in_cluster, mean2, quality_criterion) < radius

            # When the iteration is bigger or equal to div + tel and if mean1 == mean2 then convergence is reached
            # (set to 1) and we break out of the for loop, moving on to the assigning of the genes to a cluster
            if iteration >= div + tel:
                if (mean1 == mean2).all():  # and mean1.sum() == np.inf and mean2.sum() == np.inf:
                    convergence = 1
                    break

            # If no convergence is reached the mean1 will be set to mean2 and the for loop starts again, calculating
            # a new mean2 by lowering the radius
            mean1 = mean2

        # if convergence and if genes_in_sphere is bigger or equal to the minimum number of genes per sphere
        # (given variable) we increase the cluster number and append the means to cluster_means and we assign the
        # new cluster number to the genes in the cluster_assignment_array
        # if the genes_in_sphere variable is lower then the minimum required, the cluster is discarded and set to -1
        # if no convergence is found the algorithm will warn the user about not fining all clusters
        if convergence == 1:
            if genes_in_sphere.sum() >= min_number_genes_per_sphere:
                cluster += 1
                np.append(cluster_means, mean2, axis=1)
                genes_in_sphere_true = np.where(genes_in_sphere)
                cluster_assignment_array[np.take(index_of_zeros_dataset, genes_in_sphere_true)] = cluster
                output = f"Found cluster with {genes_in_sphere.sum()}"
                if genes_in_sphere.sum() == 1:
                    print(output + f" gene  / Convergence after {iteration} loops")
                else:
                    print(output + f" genes / Convergence after {iteration} loops")
                toofewpoints = 0

            else:
                genes_in_sphere_true = np.where(genes_in_sphere)
                cluster_assignment_array[np.take(index_of_zeros_dataset, genes_in_sphere_true)] = -1
                output = f"Discarded cluster with {genes_in_sphere.sum()}"
                if genes_in_sphere.sum() == 1:
                    print(output + " gene")
                else:
                    print(output + " genes")
                toofewpoints += 1

        else:
            print("Warning! No convergence: Not all clusters may have been found!")

    # From here on this is code to generate a plot of all the clustered data
    if dataset_original.any():
        at = np.transpose(dataset_original)
    else:
        at = np.transpose(dataset)

    cluster_assignment_array = cluster_assignment_array - (cluster_assignment_array > -1)

    if plot:
        cluster_expression_dict = {}
        for loop in range(cluster):
            mv = at[cluster_assignment_array == loop, :]
            # cluster_expression_dict[loop] = mv.flatten()
            for gene_idx in range(len(mv)):
                mv_gene = mv[gene_idx]
                cluster_expression_dict[f'{loop} {gene_idx}'] = mv_gene
        fig, ax = plt.subplots()
        ax.boxplot(cluster_expression_dict.values())
        ax.set_xticklabels(cluster_expression_dict.keys())
        plt.show()

        df = dataset
        pca = PCA()

        cluster_names_list = []
        for cluster in cluster_assignment_array:
            cluster_name = "Cluster " + str(cluster)
            cluster_names_list.append(cluster_name)

        df = dataset

        components = pca.fit_transform(df)

        labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }
        fig = px.scatter_matrix(
            components,
            labels=labels,
            dimensions=range(4),
            # color=cluster_names_list
        )
        fig.update_traces(diagonal_visible=False)
        fig.show()

    print(f"Cluster assignment = {cluster_assignment_array}")
    return cluster_assignment_array

