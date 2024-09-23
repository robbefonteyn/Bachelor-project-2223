import numpy as np


def dist_misval(array1, array2, mode):
    # variable array2 here is a single column vector in wan_clus_adap2, with as many rows as there are in array1
    # array2 is actually the means of the rows of arra1 (in wan_clus_adap2)
    # nr_mis_values = np.count_nonzero(array1 == float('inf')  + np.count_nonzero(array2 == float('inf') or np.nan)
    nr_mis_values = np.count_nonzero(np.isnan(array1)) + np.count_nonzero(np.isinf(array1)) + \
                    np.count_nonzero(np.isnan(array2)) + np.count_nonzero(np.isinf(array2))

    dif = array1 - (array2 * np.ones(array1.shape[1]))
    if nr_mis_values == 0:
        if mode == 1:
            # dist here is a row vector containing the highest value of each column when the absolute values are used
            dist = np.amax(abs(dif), 0)

        elif mode == 2:
            dist = np.sqrt(sum(np.square(dif)))

        else:
            dist = sum(abs(dif))

    else:

        nr_columns_dif = dif.shape[0]

        sum_of_finite_values = sum(np.isfinite(dif))
        sum_of_finite_values = sum_of_finite_values.astype('float')
        sum_of_finite_values[sum_of_finite_values == 0] = np.NaN

        dif[np.isfinite(dif) == False] = 0

        if mode == 1:
            # Same if number of missing values is 0
            dist = np.amax(abs(dif), 0)

        elif mode == 2:
            # dist = sqrt((Nrdim*Q.^-1).*sum(dif.^2,1))
            dist = np.sqrt((nr_columns_dif * (sum_of_finite_values ** -1)) * sum(np.square(dif)))

        else:
            dist = ((nr_columns_dif * sum_of_finite_values ** -1) * sum(abs(dif)))

    return dist
