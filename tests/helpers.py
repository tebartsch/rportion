import numpy as np
from sortedcontainers import SortedList


def max_pos(ind: int, vector: np.ndarray[(None,), int]):
    """
    Get lower and upper indices of the largest zero region in vector containing ind.

                         0  1  2  3  4  5  6
    max_pos(3, np.array([0, 1, 0, 0, 0, 0, 1]) == (2, 5)
    max_pos(3, np.array([0, 1, 0, 0, 0, 0, 0]) == (2, 6)
    """
    used_inds = np.nonzero(vector != 0)[0]
    higher_inds = used_inds[used_inds >= ind]
    lower_inds = used_inds[used_inds <= ind]
    if higher_inds.size > 0:
        highest_ind = higher_inds.min() - 1
    else:
        highest_ind = vector.size - 1
    if lower_inds.size > 0:
        lowest_ind = lower_inds.max() + 1
    else:
        lowest_ind = 0
    if lowest_ind <= highest_ind:
        return lowest_ind, highest_ind
    else:
        return None


def get_maximal_rectangles_from_numpy(occupation_mat: np.ndarray):
    """
    Return a subset of maximal rectangles than can be fitted into the free spaces specified by matrix. Zero value
    means that the location is free. I.e. for

                         0  1  2  3  4
                    ---+---------------
                     0 | 0  0  0  0  0
                     1 | 0  1  1  0  0
        usage_mat =  2 | 0  1  1  1  0
                     3 | 0  0  1  1  0
                     4 | 0  0  0  0  0

    the return will contain six rectangles.

        x0 x1   y0 y1
        -------------
         0  5    0  1
         0  5    4  5
         0  2    3  5
         0  1    0  5
         3  5    0  2
         4  5    0  5

    Parameters
    ----------
    occupation_mat : ndarray
        2d numpy matrix containing zero values for free pixels and non zero values for non free pixels
    mode : str
        One of "all", "greedy-1" and "greedy-2".
          - "all": find all rectangle
          - "greedy-1": faster than "all" but does not find all rectangles
          - "greedy-2": faster than "greedy-1" but finds even less rectangles

    Returns
    -------
    free_sectors: SortedList[(int, int, int, int)]
    """
    free_rectangles = SortedList([])

    for x in range(occupation_mat.shape[1]):
        for y in range(occupation_mat.shape[0]):
            # For every pixel start by growing in y direction and after that in x direction.
            max_x_pos = max_pos(y, occupation_mat[:, x])
            if max_x_pos is None:
                continue
            for search_x in range(max_x_pos[0], max_x_pos[1]+1):
                pos_y = (min(search_x, y), max(search_x, y))
                pos_x = max_pos(x, np.max(occupation_mat[pos_y[0]:pos_y[1] + 1, :] != 0, axis=0))
                pos_y = max_pos(y, occupation_mat[:, pos_x[0]:pos_x[1]+1].max(axis=1))
                if pos_x is not None:
                    new_coords = pos_x[0], pos_x[1] + 1, pos_y[0], pos_y[1] + 1
                    if new_coords not in free_rectangles:
                        free_rectangles.add(new_coords)
            # For every pixel start by growing in x direction and after that in y direction.
            max_y_pos = max_pos(x, occupation_mat[y, :])
            if max_y_pos is None:
                continue
            for search_y in range(max_y_pos[0], max_y_pos[1]+1):
                pos_x = (min(search_y, x), max(search_y, x))
                pos_y = max_pos(y, np.max(occupation_mat[:, pos_x[0]:pos_x[1] + 1] != 0, axis=1))
                pos_x = max_pos(x, occupation_mat[pos_y[0]:pos_y[1]+1, :].max(axis=0))
                if pos_x is not None:
                    new_coords = pos_x[0], pos_x[1] + 1, pos_y[0], pos_y[1] + 1
                    if new_coords not in free_rectangles:
                        free_rectangles.add(new_coords)

    return free_rectangles