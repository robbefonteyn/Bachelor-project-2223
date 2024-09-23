import numpy as np


def mean_misval(array, mode):
    nr_missing_values = np.count_nonzero(np.isnan(array)) + np.count_nonzero(np.isinf(array))

    if nr_missing_values == 0:
        if mode == 1:
            mean = (array.min(1) + array.max(1)) / 2
            mean = mean.reshape(mean.shape[0], 1)

        elif mode == 2:
            # Need to make it so it only calculates the mean per row
            mean = np.mean(array, axis=1)
            mean = mean.reshape(mean.shape[0], 1)

        else:
            # Need to make it so it only calculates the median per row
            mean = np.median(array, axis=1)
            mean = mean.reshape(mean.shape[0], 1)

    else:
        if mode == 1:
            mean = (array.min(1) + array.max(1)) / 2
            mean = mean.reshape(mean.shape[0], 1)

        elif mode == 2:
            # sum_of_finite_values needs to be a column vector that only has the amount of finite values per row
            sum_of_finite_values = np.sum(np.isfinite(array), axis=1)
            sum_of_finite_values = sum_of_finite_values.reshape(sum_of_finite_values.shape[0], 1)
            sum_of_finite_values = sum_of_finite_values.astype('float')
            sum_of_finite_values[sum_of_finite_values == 0] = np.NaN

            array[np.isfinite(array) == False] = 0

            # mn = sum(A, 2)./Q
            sum_per_row = np.sum(array, axis=1)
            sum_per_row = sum_per_row.reshape(sum_per_row.shape[0], 1)
            mean = sum_per_row / sum_of_finite_values

        else:
            # Get the number of rows of the array
            nr_rows = array.shape[0]
            mean = np.empty((0, 1))
            for row in range(nr_rows):
                # M=median(A(row,find(finite(A(row,:)))))
                to_calc = array[row][np.isfinite(array[row])]
                median = np.median(to_calc)
                if median:
                    # mean.append(m)
                    mean = np.append(mean, np.array([[median]]))
                else:
                    # mean.append(None)
                    mean = np.append(mean, np.nan)
            mean = mean.reshape(mean.shape[0], 1)

    return mean
