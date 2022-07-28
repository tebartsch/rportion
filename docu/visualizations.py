import os
import random
import time

import numpy as np
import portion as P
from matplotlib.gridspec import GridSpec
from numpy import ndarray
from sortedcontainers import SortedList
from tqdm import tqdm
from matplotlib.axes import Axes
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle as MplRectangle
import imageio

from rportion import rclosed, rclosedopen, RectBisection
from matplotlib import pyplot as plt

from rportion.rportion import rempty, RectPolygon, ropen

algo_tree_polygon = "tree (polygon)"
algo_tree_bisection = "tree (bisection)"
algo_array_all = "array all"
algo_array_greedy_1 = "array greedy 1"
algo_array_greedy_2 = "array greedy 2"
algorithms = [
    algo_tree_polygon,
    algo_array_all,
    algo_array_greedy_1,
    algo_array_greedy_2
]


def max_pos(ind: int, vector: ndarray[(None,), int]):
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


def get_maximal_rectangles_from_numpy(occupation_mat: ndarray, mode: str):
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

    if mode == algo_array_all:
        for y in range(occupation_mat.shape[0]):
            for x in range(occupation_mat.shape[1]):
                # For every pixel start by growing in y direction and after that in x direction.
                max_x_pos = max_pos(x, occupation_mat[:, y])
                if max_x_pos is None:
                    continue
                for search_x in range(max_x_pos[0], max_x_pos[1] + 1):
                    pos_x = (min(search_x, x), max(search_x, x))
                    pos_y = max_pos(y, np.max(occupation_mat[pos_x[0]:pos_x[1] + 1, :] != 0, axis=0))
                    pos_x = max_pos(x, occupation_mat[:, pos_y[0]:pos_y[1] + 1].max(axis=1))
                    if pos_y is not None:
                        new_coords = pos_x[0], pos_x[1] + 1, pos_y[0], pos_y[1] + 1
                        if new_coords not in free_rectangles:
                            free_rectangles.add(new_coords)
                # For every pixel start by growing in x direction and after that in y direction.
                max_y_pos = max_pos(y, occupation_mat[x, :])
                if max_y_pos is None:
                    continue
                for search_y in range(max_y_pos[0], max_y_pos[1] + 1):
                    pos_y = (min(search_y, y), max(search_y, y))
                    pos_x = max_pos(x, np.max(occupation_mat[:, pos_y[0]:pos_y[1] + 1] != 0, axis=1))
                    pos_y = max_pos(y, occupation_mat[pos_x[0]:pos_x[1] + 1, :].max(axis=0))
                    if pos_y is not None:
                        new_coords = pos_x[0], pos_x[1] + 1, pos_y[0], pos_y[1] + 1
                        if new_coords not in free_rectangles:
                            free_rectangles.add(new_coords)

    if mode == algo_array_greedy_1:
        for y in range(occupation_mat.shape[0]):
            for x in range(occupation_mat.shape[1]):
                # For every pixel start by growing in y direction and after that in x direction.
                pos_x = max_pos(x, occupation_mat[:, y])
                if pos_x is not None:
                    pos_y = max_pos(y, np.max(occupation_mat[pos_x[0]:pos_x[1] + 1, :] != 0, axis=0))
                    if pos_y is not None:
                        new_coords = pos_x[0], pos_x[1] + 1, pos_y[0], pos_y[1] + 1
                        if new_coords not in free_rectangles:
                            free_rectangles.add(new_coords)
                # For every pixel start by growing in x direction and after that in y direction.
                pos_y = max_pos(y, occupation_mat[x, :])
                if pos_y is not None:
                    pos_x = max_pos(x, np.max(occupation_mat[:, pos_y[0]:pos_y[1] + 1] != 0, axis=1))
                    if pos_y is not None:
                        new_coords = pos_x[0], pos_x[1] + 1, pos_y[0], pos_y[1] + 1
                        if new_coords not in free_rectangles:
                            free_rectangles.add(new_coords)

    if mode == algo_array_greedy_2:
        # For every pixel start by growing in y direction and after that in x direction.
        finished_pixels = occupation_mat.copy()
        while True:
            if np.all(finished_pixels != 0):
                break
            argmax = np.argmax((finished_pixels == 0).flatten())
            x, y = np.unravel_index(argmax, occupation_mat.shape)
            pos_x = max_pos(x, occupation_mat[:, y])
            if pos_x is not None:
                pos_y = max_pos(y, np.max(occupation_mat[pos_x[0]:pos_x[1] + 1, :] != 0, axis=0))
                if pos_y is not None:
                    new_coords = pos_x[0], pos_x[1] + 1, pos_y[0], pos_y[1] + 1
                    finished_pixels[pos_x[0]:pos_x[1] + 1, pos_y[0]:pos_y[1] + 1] = 1
                    if new_coords not in free_rectangles:
                        free_rectangles.add(new_coords)
        # For every pixel start by growing in x direction and after that in y direction.
        finished_pixels = occupation_mat.copy()
        while True:
            if np.all(finished_pixels != 0):
                break
            argmax = np.argmax((finished_pixels == 0).flatten())
            x, y = np.unravel_index(argmax, occupation_mat.shape)
            pos_y = max_pos(y, occupation_mat[x, :])
            if pos_y is not None:
                pos_x = max_pos(x, np.max(occupation_mat[:, pos_y[0]:pos_y[1] + 1] != 0, axis=1))
                if pos_y is not None:
                    new_coords = pos_x[0], pos_x[1] + 1, pos_y[0], pos_y[1] + 1
                    finished_pixels[pos_x[0]:pos_x[1] + 1, pos_y[0]:pos_y[1] + 1] = 1
                    if new_coords not in free_rectangles:
                        free_rectangles.add(new_coords)

    return free_rectangles


def plot_rectangles(ax: Axes,
                    rectangles: list[(int, int, int, int)],
                    color_ind: int,
                    label: str,
                    alpha: float = 0.35):
    cmap = get_cmap('Set3')
    flag = True
    for x0, x1, y0, y1 in rectangles:
        if flag:
            kwargs = {'label': label}
            flag = False
        else:
            kwargs = {}
        ax.add_patch(MplRectangle((x0, y0), x1 - x0, y1 - y0,
                                  facecolor=cmap(color_ind),
                                  edgecolor="black",
                                  linewidth=0,
                                  alpha=alpha,
                                  **kwargs))


def bounding_coords(rectangles: RectPolygon):
    x_int = rectangles.enclosing_x_interval
    y_int = rectangles.enclosing_y_interval
    if x_int.empty or y_int.empty:
        return None
    else:
        return x_int.lower, x_int.upper, y_int.lower, y_int.upper


def plot_rect_bisection(ax: Axes, bisec: RectBisection,
                        box: tuple[tuple[int, int], tuple[int, int]] | None = None):
    if box is None:
        x_int, y_int = bisec.used_polygon.enclosing_intervals
        assert (x_int.lower != -P.inf and x_int.upper != P.inf
                and y_int.lower != -P.inf and y_int.upper != P.inf), (
            f"If parameter box == None the rectangle must be bounded. "
            f"Boundaries are x = {x_int}, y = {y_int}")
        x_extra = (x_int.upper - x_int.lower) // 10 + 2
        y_extra = (y_int.upper - y_int.lower) // 10 + 2
        box = ((x_int.lower - x_extra, x_int.upper + x_extra),
               (y_int.lower - y_extra, y_int.upper + y_extra))

    ax.set_xlim(*box[0])
    ax.set_ylim(*box[1])
    enclosing_rec = rclosed(box[0][0], box[0][1], box[1][0], box[1][1])
    used_coords = [bounding_coords(enclosing_rec & rec)
                   for rec in bisec.maximal_rectangles()]
    plot_rectangles(ax, [e for e in used_coords if e is not None], 3, "used")
    free_coords = [bounding_coords(enclosing_rec & rec)
                   for rec in bisec.maximal_free_rectangles()]
    plot_rectangles(ax, [e for e in free_coords if e is not None], 0, "free")


def create_gif():
    r_list = [
        ("add", rempty()),
        ("add", rclosedopen(0, 3, 0, 6)),
        ("sub", rclosedopen(2, 5, 2, 5)),
        ("add", rclosedopen(-1, 1, -1, 1)),
        ("add", rclosedopen(5, 8, 5, 8)),
        ("add", rclosedopen(5, 6, 5, 6)),
        ("add", rclosedopen(6, 7, 0, 6)),
        ("add", rclosedopen(-1, 8, 3, 4)),
        ("sub", rclosedopen(1, 8, 2, 3)),
    ]

    filename_base = "image"
    output_folder = "images"
    os.makedirs(output_folder, exist_ok=True)

    poly = RectBisection()
    for i, (operation, r) in enumerate(r_list):
        fig, ax = plt.subplots(1, 1)
        ax.set_xlim([-2 - 0.5, 10 - 0.5])
        ax.set_ylim([-2 - 0.5, 10 - 0.5])
        if operation == "add":
            poly |= r
        elif operation == "sub":
            poly -= r
        fig, ax = plt.subplots(1, 1)
        x_lim = (-2, 9)
        y_lim = (-2, 9)
        plot_rect_bisection(ax, poly, (x_lim, y_lim))
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim, )
        ax.legend(loc="upper right")
        fig.savefig(os.path.join(output_folder, f"{filename_base}-{i}.png"))

    with imageio.get_writer(os.path.join(f"{filename_base}.gif"),
                            mode='I', duration=1) as writer:
        for i in range(len(r_list)):
            image = imageio.v2.imread(os.path.join(output_folder, f"{filename_base}-{i}.png"))
            writer.append_data(image)


def create_random_polygon(n: int,
                          x_max: int, max_x_len: int,
                          y_max: int, max_y_len: int,
                          algorithms: list[str],
                          show_progress: bool = False):
    assert set(algorithms).issubset(set(algorithms))

    r_list = []
    for i in range(n):
        x_left = random.randint(0, x_max - max_x_len)
        x_right = x_left + random.randint(1, max_x_len)
        y_left = random.randint(0, y_max - max_y_len)
        y_right = y_left + random.randint(1, max_y_len)
        r_list.append((x_left, x_right, y_left, y_right))

    if show_progress:
        iterator = tqdm(list(zip(r_list[::2], r_list[1::2])))
    else:
        iterator = list(zip(r_list[::2], r_list[1::2]))

    polygons = {}
    times = {}

    for algo in algorithms:
        if algo == algo_tree_polygon:
            polygons[algo] = ropen(-P.inf, P.inf, -P.inf, P.inf)
        elif algo == algo_tree_bisection:
            polygons[algo] = RectBisection.from_rect_polygon(ropen(-P.inf, P.inf, -P.inf, P.inf))
        else:
            polygons[algo] = np.zeros((x_max, y_max))
        times[algo] = np.zeros(n)

    i = 0
    for r1, r2 in iterator:
        # Add Rectangle
        for algo in algorithms:
            start_time = time.time()
            if algo in {algo_tree_polygon, algo_tree_bisection}:
                polygons[algo] -= rclosedopen(*r1)
            else:
                polygons[algo][r1[0]:r1[1], r1[2]:r1[3]] = 1
                get_maximal_rectangles_from_numpy(polygons[algo], mode=algo)
            end_time = time.time()
            times[algo][i] = (end_time - start_time) * 1000.0

        i += 1

        # Subtract Rectangle
        for algo in algorithms:
            start_time = time.time()
            if algo in {algo_tree_polygon, algo_tree_bisection}:
                polygons[algo] |= rclosedopen(*r2)
            else:
                polygons[algo][r2[0]:r2[1], r2[2]:r2[3]] = 1
                get_maximal_rectangles_from_numpy(polygons[algo], mode=algo)
            end_time = time.time()
            times[algo][i] = (end_time - start_time) * 1000.0

        i += 1

    enclosing_rec = rclosed(0, x_max, 0, y_max)
    max_rectangles = {}
    for algo in algorithms:
        if algo == algo_tree_polygon:
            max_rectangles[algo] = (
                [],
                [coords
                 for rec in polygons[algo].maximal_rectangles()
                 if (coords := bounding_coords(rec & enclosing_rec)) is not None]
            )
        elif algo == algo_tree_bisection:
            max_rectangles[algo] = (
                [coords
                 for rec in polygons[algo].maximal_free_rectangles()
                 if (coords := bounding_coords(rec & enclosing_rec)) is not None],
                [coords
                 for rec in polygons[algo].maximal_rectangles()
                 if (coords := bounding_coords(rec & enclosing_rec)) is not None]
            )
        else:
            max_rectangles[algo] = (
                get_maximal_rectangles_from_numpy(polygons[algo] == 0, mode=algo),
                get_maximal_rectangles_from_numpy(polygons[algo], mode=algo)
            )

    return max_rectangles, times


def evaluate_random_polygon(save_image=True, show_progress=True):
    n = 100
    x_max = 512
    max_x_len = 75
    y_max = 512
    max_y_len = 75

    algos = [
        algo_tree_polygon,
        algo_tree_bisection,
        # algo_array_all,
        # algo_array_greedy_1,
        algo_array_greedy_2
    ]

    repetitions = 1
    max_rectangles_dict, times_dict = create_random_polygon(n=n,
                                                            x_max=x_max, max_x_len=max_x_len,
                                                            y_max=y_max, max_y_len=max_y_len,
                                                            algorithms=algos,
                                                            show_progress=show_progress)
    # Perform more iterations to average time consumption
    for _ in range(repetitions - 1):

        _, curr_times_dict = create_random_polygon(
            n=n,
            x_max=x_max, max_x_len=max_x_len,
            y_max=y_max, max_y_len=max_y_len,
            algorithms=algos,
            show_progress=True
        )
        for name, curr_times in curr_times_dict.items():
            times_dict[name] += curr_times

    for name in times_dict:
        times_dict[name] /= repetitions

    if algo_tree_bisection in algos:
        for algo in list(set(algos) - {algo_tree_bisection}):
            print(f"ALGORITHM {algo}")
            common_used_rectangles = (set(max_rectangles_dict[algo][0])
                                      & set(max_rectangles_dict[algo_tree_bisection][0]))
            common_free_rectangles = (set(max_rectangles_dict[algo][1])
                                      & set(max_rectangles_dict[algo_tree_bisection][1]))
            only_tree_used_recs = (set(max_rectangles_dict[algo_tree_bisection][0])
                                   - set(max_rectangles_dict[algo][0]))
            only_array_used_recs = (set(max_rectangles_dict[algo][0])
                                    - set(max_rectangles_dict[algo_tree_bisection][0]))
            only_tree_free_recs = (set(max_rectangles_dict[algo_tree_bisection][1])
                                   - set(max_rectangles_dict[algo][1]))
            only_array_free_recs = (set(max_rectangles_dict[algo][1])
                                    - set(max_rectangles_dict[algo_tree_bisection][1]))

            print(f"Used rectangles")
            print(f"  array & tree: {len(common_used_rectangles)} rectangles")
            print(f"  tree - array: {only_tree_used_recs}")
            print(f"  array - tree: {only_array_used_recs}")

            print(f"Free rectangles")
            print(f"  array & tree: {len(common_free_rectangles)} rectangles")
            print(f"  tree - array: {only_tree_free_recs}")
            print(f"  array - tree: {only_array_free_recs}")
            print()

    if save_image:
        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, len(algos), figure=fig, hspace=0.4)
        fig.show()
        algo_axes = [fig.add_subplot(gs[0, i]) for i in range(len(algos))]
        time_ax = fig.add_subplot(gs[1, :])

        x_lim = (0, x_max)
        y_lim = (0, y_max)

        for i, (name, max_rectangles) in enumerate(max_rectangles_dict.items()):
            algo_axes[i].set_title(f"{name}")
            plot_rectangles(algo_axes[i], max_rectangles[0], 3, "used")
            plot_rectangles(algo_axes[i], max_rectangles[1], 0, "free")
            algo_axes[i].set_title(f"{name}")
            algo_axes[i].set_xlim(x_lim)
            algo_axes[i].set_ylim(y_lim)
        algo_axes[-1].legend(loc="upper right")

        for name, times in times_dict.items():
            time_ax.plot(times, label=name)
        time_ax.set_xlabel("index of union/subtraction")
        time_ax.set_ylabel("ms")
        time_ax.set_title("time consumption per step")
        time_ax.set_yscale('log')
        time_ax.legend(loc="upper left")

        fig.savefig(f"random.png")


def benchmark():
    x_max = 12
    max_x_len = 4
    y_max = 12
    max_y_len = 4
    N = np.arange(10, 100)

    algos = [
        algo_tree_polygon,
        algo_tree_bisection,
        # algo_array_all,
        algo_array_greedy_1,
        algo_array_greedy_2,
    ]

    times = {algo: np.zeros(N.size) for algo in algos}

    repetitions = 1
    for i, n in tqdm(list(enumerate(N))):
        for _ in range(repetitions):
            _, times_dict = create_random_polygon(n=n,
                                                  x_max=x_max, max_x_len=max_x_len,
                                                  y_max=y_max, max_y_len=max_y_len,
                                                  algorithms=algos)
            for algo in algos:
                times[algo][i] += times_dict[algo].sum() / repetitions

    fig, ax = plt.subplots(1, 1)
    for algo, arr in times.items():
        ax.plot(N, arr, label=algo)
    ax.set_title("time consumption")
    ax.set_xlabel("number of rectangles")
    ax.set_ylabel("ms")
    ax.legend(loc="upper right")
    fig.savefig(f"benchmark.png")


def main():
    random.seed(452351)
    print("creating gif")
    create_gif()
    print("evaluate a random sequence of unions and subtractions")
    evaluate_random_polygon(show_progress=True)
    print("benchmark random sequence of unions and subtractions of increasing length")
    benchmark()


if __name__ == "__main__":
    main()
