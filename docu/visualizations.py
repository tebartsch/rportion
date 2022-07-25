from matplotlib.axes import Axes
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle as MplRectangle
from portion import closed

from rportion import rclosed, rclosedopen, RPolygon
from matplotlib import pyplot as plt


def plot_free_rectangles(ax: Axes,
                         rectangles: list[(int, int, int, int)]):
    cmap = get_cmap('Set3')
    for x0, x1, y0, y1 in rectangles:
        ax.add_patch(MplRectangle((x0 - 0.5, y0 - 0.5), x1 - x0 + 1, y1 - y0 + 1,
                                  facecolor=cmap(-1),
                                  edgecolor="black",
                                  linewidth=0,
                                  alpha=0.35))


def create_gif():
    rec = RPolygon()

    enclosing_rec = rclosed(-2, 10, -2, 10)
    r_list = [
        rclosedopen(0, 4, 0, 4),
        rclosedopen(2, 6, 2, 6),
    ]

    for r in r_list:
        fig, ax = plt.subplots(1, 1)
        ax.set_xlim([-2-0.5, 10-0.5])
        ax.set_ylim([-2-0.5, 10-0.5])
        rec -= r
        ints = [(r.x_interval & closed(-2, 10), r.y_interval & closed(-2, 10))
                for r in rec]
        print(ints)
        coords = [(x_int.lower, x_int.upper-1, y_int.lower, y_int.upper-1) for x_int, y_int in ints]
        plot_free_rectangles(ax, coords)
        plt.show()


if __name__ == "__main__":
    create_gif()
