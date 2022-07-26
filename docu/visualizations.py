import os
import random
import time

import numpy as np
from numpy import ndarray
from sortedcontainers import SortedList
from tqdm import tqdm
from matplotlib.axes import Axes
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle as MplRectangle
import imageio

from rportion import rclosed, rclosedopen, RPolygon
from matplotlib import pyplot as plt

from rportion.rportion import rempty


def max_pos(ind: int, vector: ndarray[(None,), int]):
    """
    Get lower and upper indices of the largest zero region in vector containing ind.

                   0  1  2  3  4  5  6
    f(3, np.array([0, 1, 0, 0, 0, 0, 1]) == (2, 5)
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


def get_maximal_rectangles_from_numpy(occupation_mat: ndarray):
    """
    Return all maximal rectangles than can be fitted into the free spaces specified by matrix. Zero value
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

    Returns
    -------
    free_sectors: SortedList[(int, int, int, int)]
    """
    free_rectangles = SortedList([])
    # For every pixel start by growing in word line direction and after that in bit line direction.
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
    # For every pixel start by growing in bit line direction and after that in word line direction.
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
                    label: str):
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
                                  alpha=0.35,
                                  **kwargs))


def bounding_coords(rectangles: RPolygon):
    x_int = rectangles.enclosing_x_interval
    y_int = rectangles.enclosing_y_interval
    if x_int.empty or y_int.empty:
        return None
    else:
        return x_int.lower, x_int.upper, y_int.lower, y_int.upper


def plot_rpolygon(ax: Axes, poly: RPolygon, box=((int, int), (int, int))):
    enclosing_rec = rclosed(box[0][0], box[0][1], box[1][0], box[1][1])
    used_coords = [bounding_coords(rec & enclosing_rec)
                   for rec in poly.maximal_used_rectangles()]
    plot_rectangles(ax, [e for e in used_coords if e is not None], 3, "used")
    free_coords = [bounding_coords(rec & enclosing_rec)
                   for rec in poly.maximal_free_rectangles()]
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

    poly = RPolygon()
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
        plot_rpolygon(ax, poly, (x_lim, y_lim))
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim, )
        ax.legend(loc="upper right")
        fig.savefig(os.path.join(output_folder, f"{filename_base}-{i}.png"))

    with imageio.get_writer(os.path.join(f"{filename_base}.gif"),
                            mode='I', duration=1) as writer:
        for i in range(len(r_list)):
            image = imageio.v2.imread(os.path.join(output_folder, f"{filename_base}-{i}.png"))
            writer.append_data(image)


def create_random_polygon(n: int, x_max: int, max_x_len: int, y_max: int, max_y_len: int,
                          show_progress: bool = False):
    r_list = []
    for i in range(n):
        x_left = random.randint(0, x_max - max_x_len)
        x_right = x_left + random.randint(0, max_x_len)
        y_left = random.randint(0, y_max - max_y_len)
        y_right = y_left + random.randint(0, max_y_len)
        r_list.append((x_left, x_right, y_left, y_right))

    if show_progress:
        iterator = tqdm(list(zip(r_list[::2], r_list[1::2])))
    else:
        iterator = list(zip(r_list[::2], r_list[1::2]))

    polygons = {
        "numpy array": np.zeros((x_max + 1, y_max + 1)),
        "interval tree": RPolygon(),
    }
    times = {
        "numpy array": np.zeros(n),
        "interval tree": np.zeros(n),
    }
    i = 0
    for r1, r2 in iterator:
        # Rectangle 1
        start_time = time.time()
        polygons["numpy array"][r1[0]:r1[1], r1[2]:r1[3]] = 1
        get_maximal_rectangles_from_numpy(polygons["numpy array"])
        end_time = time.time()
        times["numpy array"][i] = (end_time - start_time) * 1000.0

        start_time = time.time()
        polygons["interval tree"] |= rclosedopen(*r1)
        end_time = time.time()
        times["interval tree"][i] = (end_time - start_time) * 1000.0

        i += 1

        # Rectangle 2
        start_time = time.time()
        polygons["numpy array"][r2[0]:r2[1], r2[2]:r2[3]] = 0
        get_maximal_rectangles_from_numpy(polygons["numpy array"])
        end_time = time.time()
        times["numpy array"][i] = (end_time - start_time) * 1000.0

        start_time = time.time()
        polygons["interval tree"] -= rclosedopen(*r2)
        end_time = time.time()
        times["interval tree"][i] = (end_time - start_time) * 1000.0

        i += 1

    enclosing_rec = rclosed(0, x_max, 0, y_max)
    max_rectangles = {
        "numpy array": (
            get_maximal_rectangles_from_numpy(polygons["numpy array"] == 0),
            get_maximal_rectangles_from_numpy(polygons["numpy array"]),
        ),
        "interval tree": (
            [coords
             for rec in polygons["interval tree"].maximal_used_rectangles()
             if (coords := bounding_coords(rec & enclosing_rec)) is not None],
            [coords
             for rec in polygons["interval tree"].maximal_free_rectangles()
             if (coords := bounding_coords(rec & enclosing_rec)) is not None],
        )
    }

    return max_rectangles, times


def evaluate_random_polygon(save_image=True):
    n = 30
    x_max = 512
    max_x_len = 128
    y_max = 512
    max_y_len = 128

    repetitions = 10
    max_rectangles_dict, times_dict = create_random_polygon(n=n,
                                                            x_max=x_max, max_x_len=max_x_len,
                                                            y_max=y_max, max_y_len=max_y_len,
                                                            show_progress=True)
    for _ in range(repetitions - 1):

        _, curr_times_dict = create_random_polygon(
            n=n,
            x_max=x_max, max_x_len=max_x_len,
            y_max=y_max, max_y_len=max_y_len,
            show_progress=True
        )
        for name, curr_times in curr_times_dict.items():
            times_dict[name] += curr_times

    for name in times_dict:
        times_dict[name] /= repetitions

    if save_image:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        x_lim = (0, x_max)
        y_lim = (0, y_max)

        for i, (name, max_rectangles) in enumerate(max_rectangles_dict.items()):
            plot_rectangles(axes[i], max_rectangles[0], 3, "used")
            plot_rectangles(axes[i], max_rectangles[1], 0, "free")
            axes[i].set_xlim(x_lim)
            axes[i].set_ylim(y_lim, )
            axes[i].legend(loc="upper right")
            axes[i].set_title(f"{name}")

        for name, times in times_dict.items():
            axes[-1].plot(times, label=name)
            axes[-1].set_xlabel("index of union/subtraction")
            axes[-1].set_ylabel("ms")
            axes[-1].set_title("time consumption per step")
        axes[-1].legend(loc="upper right")

        fig.savefig(f"random.png")


def benchmark():
    x_max = 512
    max_x_len = 128
    y_max = 512
    max_y_len = 128
    N = np.arange(10, 30)
    times = np.zeros((2, N.shape[0]))

    repetitions = 3

    for i, n in tqdm(list(enumerate(N))):
        for _ in range(repetitions):
            _, times_dict = create_random_polygon(n=n,
                                                  x_max=x_max, max_x_len=max_x_len,
                                                  y_max=y_max, max_y_len=max_y_len)
        times[0, i] = times_dict["numpy array"].sum() / repetitions
        times[1, i] = times_dict["interval tree"].sum() / repetitions

    fig, ax = plt.subplots(1, 1)
    ax.plot(N, times[0], label=f"numpy array")
    ax.plot(N, times[1], label=f"interval tree")
    ax.set_title("time consumption")
    ax.set_xlabel("number of rectangles")
    ax.set_ylabel("ms")
    ax.legend(loc="upper right")
    fig.savefig(f"benchmark.png")


def main():
    print("creating gif")
    create_gif()
    print("evaluate a random sequence of unions and subtractions")
    evaluate_random_polygon()
    print("benchmark random sequence of unions and subtractions of increasing length")
    benchmark()


if __name__ == "__main__":
    main()
