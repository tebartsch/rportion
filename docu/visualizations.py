import os
import random
import time

import numpy as np
from tqdm import tqdm
from matplotlib.axes import Axes
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle as MplRectangle
import imageio

from rportion import rclosed, rclosedopen, RPolygon
from matplotlib import pyplot as plt

from rportion.rportion import rempty


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
        fig.close()

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
        r_list.append(rclosedopen(x_left, x_right, y_left, y_right))

    if show_progress:
        iterator = tqdm(list(zip(r_list[::2], r_list[1::2])))
    else:
        iterator = list(zip(r_list[::2], r_list[1::2]))

    poly = RPolygon()
    i = 0
    times = np.zeros(n)
    for r1, r2 in iterator:
        start_time = time.time()
        poly |= r1
        end_time = time.time()
        times[i] = (end_time - start_time) * 1000.0
        i += 1

        start_time = time.time()
        poly -= r2
        end_time = time.time()
        times[i] = (end_time - start_time) * 1000.0
        i += 1

    return poly, times


def evaluate_random_polygon(save_image=True):
    N = 100
    x_max = N
    y_max = N

    poly, times = create_random_polygon(n=N,
                                        x_max=x_max, max_x_len=20,
                                        y_max=y_max, max_y_len=20,
                                        show_progress=True)

    if save_image:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        x_lim = (0, x_max)
        y_lim = (0, y_max)

        plot_rpolygon(axes[0], poly, (x_lim, y_lim))
        axes[0].set_xlim(x_lim)
        axes[0].set_ylim(y_lim, )
        axes[0].legend(loc="upper right")
        axes[0].set_title("random sequence of unions/subtractions")

        axes[1].plot(times, label="")
        axes[1].set_ylabel("ms")
        axes[1].set_ylabel("index of union/subtraction")
        axes[1].set_title("time consumption of every union/subtraction")

        fig.savefig(f"random.png")


def benchmark():
    x_max = 100
    y_max = 100
    N = np.arange(10, 30)
    times = np.zeros(N.shape)

    repetitions = 10

    for i, n in tqdm(list(enumerate(N))):
        start_time = time.time()
        for _ in range(repetitions):
            create_random_polygon(n=n,
                                  x_max=x_max, max_x_len=10,
                                  y_max=y_max, max_y_len=10)
        end_time = time.time()
        times[i] = (end_time - start_time) * 1000.0 / repetitions

    fig, ax = plt.subplots(1, 1)
    ax.plot(N, times, label="time consumption")
    ax.set_xlabel("number of rectangles")
    ax.set_ylabel("ms")
    ax.legend(loc="upper right")
    fig.savefig(f"benchmark.png")


def main():
    print("creating gif")
    # create_gif()
    print("evaluate a random sequence of unions and subtractions")
    evaluate_random_polygon()
    print("benchmark random sequence of unions and subtractions of increasing length")
    benchmark()


if __name__ == "__main__":
    main()
