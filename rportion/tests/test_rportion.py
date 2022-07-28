import unittest
from itertools import permutations, combinations
import random

import numpy as np
import portion as P
from portion.interval import open, closedopen, Atomic, empty, Interval

from rportion import RPolygon
from rportion.rportion import ropen, rclosedopen, rclosed, rempty, RBoundary

from rportion.tests.algorithms import get_maximal_rectangles_from_numpy


def data_tree_to_string(x_boundaries: list[RBoundary],
                        y_intervals: list[list[Interval]],
                        spacing: int):
    n = len(y_intervals)
    msg = " "*spacing + "|"
    for x_b in x_boundaries[-1:0:-1]:
        msg += f"{str(x_b):>{spacing}}"
    msg += "\n" + f"-"*spacing + "+"
    for i in range(n):
        msg += f"-"*spacing
    msg += "\n"
    for i, row in enumerate(y_intervals):
        x_b = x_boundaries[i]
        msg += f"{str(~x_b):>{spacing}}|"
        for val in row:
            msg += f"{str(val):>{spacing}}"
        msg += "\n"
    return msg

def print_mat(mat: list[list[Interval]]):
    for row in mat:
        for val in row:
            print(f"{str(val):<25}", end="")
        print()
    print()


def print_rpolygon(rpoly: RPolygon, show_trees=False, spacing=20):
    if show_trees:
        print("USED")
        msg = data_tree_to_string(list(rpoly._x_boundaries), rpoly._used_y_ranges, spacing)
        print(msg)
        for i, e in enumerate(rpoly.maximal_used_rectangles(), start=1):
            print(i, "-->", e.enclosing_x_interval, e.enclosing_y_interval)
        print("FREE")
        msg = data_tree_to_string(list(rpoly._x_boundaries), rpoly._free_y_ranges, spacing)
        print(msg)
        for i, e in enumerate(rpoly.maximal_free_rectangles(), start=1):
            print(i, "-->", e.enclosing_x_interval, e.enclosing_y_interval)



class TestRPortion(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestRPortion, self).__init__(*args, **kwargs)

    def test_RPolygon_maximal_atomic_x_rectangles(self):
        # Add two finite (in x dimension) polygons.
        #        (1)          (2)         (3)         (4)           (5)
        #          +--+       +---+       +--+       +---+       +--+
        #          |  |       |   |       |  |       |   |       |  |
        #     +--+ |  |     +---+ |     +------+     | +---+     |  | +--+
        #     |  | |  |     | | | |     | |  | |     | | | |     |  | |  |
        #     |  | +--+     | +-|-+     | +--+ |     +-|-+ |     +--+ |  |
        #     |  |          |   |       |      |       |   |          |  |
        #     +--+          +---+       +------+       +---+          +--+
        x1, x2, x3, x4 = (1, 2, 3, 4)
        y_interval_1 = closedopen(-1, 1)
        y_interval_2 = closedopen(0, 2)
        # (1)
        poly = rempty()
        poly._add_atomic(closedopen(x1, x2), y_interval_1)
        poly._add_atomic(closedopen(x3, x4), y_interval_2)
        self.assertListEqual(list(poly._maximal_used_atomic_x_rectangles()), [
            (closedopen(x1, x2), y_interval_1),
            (closedopen(x3, x4), y_interval_2)
        ])
        self.assertListEqual(list(poly._maximal_free_atomic_x_rectangles()), [
            (open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2)),
            (open(-P.inf, x3),
             Interval.from_atomic(~y_interval_1.right, y_interval_1.upper, P.inf, P.OPEN)),
            (closedopen(x2, P.inf),
             Interval.from_atomic(P.OPEN, -P.inf, y_interval_2.lower, ~y_interval_2.left)),
            (open(-P.inf, x1), open(-P.inf, P.inf)),
            (closedopen(x2, x3), open(-P.inf, P.inf)),
            (closedopen(x4, P.inf), open(-P.inf, P.inf)),
        ])
        # (2)
        poly = rempty()
        poly._add_atomic(closedopen(x1, x3), y_interval_1)
        poly._add_atomic(closedopen(x2, x4), y_interval_2)
        self.assertListEqual(list(poly._maximal_used_atomic_x_rectangles()), [
            (closedopen(x1, x4), y_interval_1 & y_interval_2),
            (closedopen(x1, x3), y_interval_1),
            (closedopen(x2, x4), y_interval_2),
            (closedopen(x2, x3), y_interval_1 | y_interval_2)
        ])
        self.assertListEqual(list(poly._maximal_free_atomic_x_rectangles()), [
            (open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2)),
            (open(-P.inf, x2),
             Interval.from_atomic(~y_interval_1.right, y_interval_1.upper, P.inf, P.OPEN)),
            (closedopen(x3, P.inf),
             Interval.from_atomic(P.OPEN, -P.inf, y_interval_2.lower, ~y_interval_2.left)),
            (open(-P.inf, x1), open(-P.inf, P.inf)),
            (closedopen(x4, P.inf), open(-P.inf, P.inf)),
        ])
        # (3)
        poly = rempty()
        poly._add_atomic(closedopen(x1, x4), y_interval_1)
        poly._add_atomic(closedopen(x2, x3), y_interval_2)
        self.assertListEqual(list(poly._maximal_used_atomic_x_rectangles()), [
            (closedopen(x1, x4), y_interval_1),
            (closedopen(x2, x3), y_interval_1 | y_interval_2),
        ])
        self.assertListEqual(list(poly._maximal_free_atomic_x_rectangles()), [
            (open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2)),
            (open(-P.inf, x2),
             Interval.from_atomic(~y_interval_1.right, y_interval_1.upper, P.inf, P.OPEN)),
            (closedopen(x3, P.inf),
             Interval.from_atomic(~y_interval_2.right, y_interval_1.upper, P.inf, P.OPEN)),
            (open(-P.inf, x1), open(-P.inf, P.inf)),
            (closedopen(x4, P.inf), open(-P.inf, P.inf)),
        ])
        # skip (4) & (5)

        # Test the following polygon addition sequence.
        #
        #   x1    x3        x1    x3  x4   x6        x1 x2 x3  x4 x5 x6
        #    +----+          +----+                   +----+
        #    |    |          |    |                   |  +-|------+
        #    +----+          +----+                   +--|-+      |
        #              ->                       ->       |        |
        #                              +----+            |      +-|--+
        #                              |    |            +------|-+  |
        #                              +----+                   +----+
        y_interval_1 = closedopen(4, 6)
        y_interval_2 = closedopen(1, 3)
        y_interval_3 = closedopen(2, 5)
        x1, x2, x3, x4, x5, x6 = (1, 2, 3, 4, 5, 6)
        # (1)
        poly = rempty()
        poly._add_atomic(closedopen(x1, x3), y_interval_1)
        poly._add_atomic(closedopen(x4, x6), y_interval_2)
        poly._add_atomic(closedopen(x2, x5), y_interval_3)
        self.assertListEqual(list(poly._maximal_used_atomic_x_rectangles()), [
            (closedopen(x1, x5), closedopen(x4, x5)),
            (closedopen(x2, x6), closedopen(x2, x3)),
            (closedopen(x2, x5), closedopen(x2, x5)),
            (closedopen(x1, x3), closedopen(x4, x6)),
            (closedopen(x4, x6), closedopen(x1, x3)),
            (closedopen(x2, x3), closedopen(x2, x6)),
            (closedopen(x4, x5), closedopen(x1, x5)),
        ])
        self.assertListEqual(list(poly._maximal_free_atomic_x_rectangles()), [
            (open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2 | y_interval_3)),
            (open(-P.inf, x4),
             Interval.from_atomic(P.OPEN, -P.inf, y_interval_3.lower, ~y_interval_3.left)),
            (closedopen(x3, P.inf),
             Interval.from_atomic(~y_interval_3.right, y_interval_3.upper, P.inf, P.OPEN)),
            (open(-P.inf, x2), Interval.from_atomic(P.OPEN, -P.inf, y_interval_1.lower, ~y_interval_1.left)),
            (closedopen(x5, P.inf),
             Interval.from_atomic(~y_interval_2.right, y_interval_2.upper, P.inf, P.OPEN)),
            (open(-P.inf, x1), open(-P.inf, P.inf)),
            (closedopen(x6, P.inf), open(-P.inf, P.inf)),
        ])

        # Test the following polygon addition sequence, where the third rectangle directly touches the second one.
        #
        #    x1   x2         x1   x2    x3   x4       x1   x2   x3   x4   x5
        #    +----+          +----+     +----+        +----+    +---------+
        #    |    |    ->    |    |     |    |    ->  |    |    |         |
        #    +----+          +----+     +----+        +----+    +---------+
        y_interval = closedopen(0, 1)
        x1, x2, x3, x4, x5 = 1, 2, 3, 4, 5
        poly = rempty()
        # poly._add_atomic(closedopen(x1, x2), y_interval)
        poly._add_atomic(closedopen(x3, x4), y_interval)
        poly._add_atomic(closedopen(x4, x5), y_interval)
        print_rpolygon(poly, show_trees=True)

        # Add three sectors in each six possible different orders and test if the result underlying
        # data structure of RPolygon is the same for all.
        #
        #      x1 x2 x3 x4 x5
        #      +--+--+--+--+
        #      |  |  |  |  |
        #      +--+--+--+--+
        x1, x2, x3, x4, x5 = (1, 2, 3, 4, 5)
        y_interval = closedopen(0, 1)
        x_lims = [(x1, x3), (x2, x4), (x3, x5)]
        y_intervals = [y_interval, y_interval, y_interval]
        poly_list = []
        for arrangement in permutations(zip(x_lims, y_intervals)):
            poly = rempty()
            for (x_a, x_b), y_int in arrangement:
                poly._add_atomic(closedopen(x_a, x_b), y_int)
            poly_list.append(poly)
        for poly_1, poly_2 in combinations(poly_list, 2):
            self.assertListEqual(poly_1._used_y_ranges, poly_2._used_y_ranges)
            self.assertListEqual(poly_1._free_y_ranges, poly_2._free_y_ranges)


    def test_RPolygon___sub__(self):
        poly = rclosedopen(-P.inf, 0, -P.inf, P.inf)
        poly._sub_interval_product(open(-P.inf, P.inf), open(-P.inf, 0) | open(1, P.inf))
        # poly._sub_atomic(closedopen(-P.inf, 0), open(-P.inf, P.inf))
        print_rpolygon(poly, show_trees=True)
        print("----")
        print_rpolygon(ropen(-P.inf, 0, -P.inf, P.inf) & rclosed(0, 3, 0, 3), show_trees=True)

    def test_RPolygon__add_atomic(self):
        # empty x_interval
        poly = rempty()
        x_atom = closedopen(1, 0)
        y_interval = open(1, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[empty()]])
        self.assertListEqual(poly._free_y_ranges, [[open(-P.inf, P.inf)]])

        # empty y_interval
        poly = rempty()
        x_atom = closedopen(-P.inf, 1)
        y_interval = empty()
        poly._add_atomic(x_atom, y_interval)
        print_rpolygon(poly, show_trees=True, spacing=15)
        self.assertListEqual(poly._used_y_ranges, [[empty()]])
        self.assertListEqual(poly._free_y_ranges, [[open(-P.inf, P.inf)]])

        # add the whole area
        poly = rempty()
        x_atom = closedopen(-P.inf, P.inf)
        y_interval = open(-P.inf, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[open(-P.inf, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[empty()]])

        # add half planes
        # # (a) Right half space
        poly = rempty()
        x_atom = closedopen(0, P.inf)
        y_interval = open(-P.inf, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [empty(), empty()],
            [open(-P.inf, P.inf)]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [empty(), open(-P.inf, P.inf)],
            [empty()]
        ])
        # # (b) Left half space
        poly = rempty()
        x_atom = closedopen(-P.inf, 0)
        y_interval = open(-P.inf, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [empty(), open(-P.inf, P.inf)],
            [empty()]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [empty(), empty()],
            [open(-P.inf, P.inf)]
        ])
        # # (c) Upper half space
        poly = rempty()
        x_atom = closedopen(-P.inf, P.inf)
        y_interval = open(-P.inf, 0)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[open(-P.inf, 0)]])
        self.assertListEqual(poly._free_y_ranges, [[closedopen(0, P.inf)]])
        # # (d) Lower half space
        poly = rempty()
        x_atom = closedopen(-P.inf, P.inf)
        y_interval = closedopen(0, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[closedopen(0, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[open(-P.inf, 0)]])

        # add a single bounded polygon
        poly = rempty()
        x_atom = closedopen(1, 3)
        y_interval = closedopen(2, 4)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [empty(), empty(), empty()],
            [empty(), closedopen(2, 4)],
            [empty()]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [open(-P.inf, 2) | closedopen(4, P.inf),
             open(-P.inf, 2) | closedopen(4, P.inf),
             open(-P.inf, P.inf)],
            [open(-P.inf, 2) | closedopen(4, P.inf),
             open(-P.inf, 2) | closedopen(4, P.inf)],
            [open(-P.inf, P.inf)]
        ])

    def test_RPolygon__sub_atomic(self):
        # Removing an arbitrary rectangle from an empty polygon is effect less
        # TODO

        # Remove rectangle with empty x_interval from the whole plane
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = closedopen(1, 0)
        y_interval = open(1, P.inf)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[open(-P.inf, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[empty()]])

        # Remove rectangle with empty y_interval from the whole plane
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = closedopen(-P.inf, 1)
        y_interval = empty()
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[open(-P.inf, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[empty()]])

        # Remove half planes from the whole plane
        # # (a) Right half space
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = closedopen(0, P.inf)
        y_interval = open(-P.inf, P.inf)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [empty(), open(-P.inf, P.inf)],
            [empty()]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [empty(), empty()],
            [open(-P.inf, P.inf)]
        ])
        # # (b) Left half space
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = closedopen(-P.inf, 0)
        y_interval = open(-P.inf, P.inf)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [empty(), empty()],
            [open(-P.inf, P.inf)]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [empty(), open(-P.inf, P.inf)],
            [empty()]
        ])
        # # (c) Upper half space
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = closedopen(-P.inf, P.inf)
        y_interval = open(-P.inf, 0)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[closedopen(0, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[open(-P.inf, 0)]])
        # # (d) Lower half space
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = closedopen(-P.inf, P.inf)
        y_interval = closedopen(0, P.inf)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[open(-P.inf, 0)]])
        self.assertListEqual(poly._free_y_ranges, [[closedopen(0, P.inf)]])

        # remove a single bounded rectangle from the plane
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = closedopen(1, 3)
        y_interval = closedopen(2, 4)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [open(-P.inf, 2) | closedopen(4, P.inf),
             open(-P.inf, 2) | closedopen(4, P.inf),
             open(-P.inf, P.inf)],
            [open(-P.inf, 2) | closedopen(4, P.inf),
             open(-P.inf, 2) | closedopen(4, P.inf)],
            [open(-P.inf, P.inf)]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [empty(), empty(), empty()],
            [empty(), closedopen(2, 4)],
            [empty()]
        ])

    def test_maximal_rectangles_extra(self):
        n = 15  # Number of rectangles

        # repeat test which generates random polygons multiply times.
        for _ in range(3):

            x_max = 16
            max_x_len = 2
            y_max = 4
            max_y_len = 2

            rec_list = []
            for i in range(n):
                x_left = random.randint(0, x_max - max_x_len)
                x_right = x_left + random.randint(1, max_x_len)
                y_left = random.randint(0, y_max - max_y_len)
                y_right = y_left + random.randint(1, max_y_len)
                rec_list.append((x_left, x_right, y_left, y_right))

            arr = np.zeros((y_max, x_max))
            poly = rempty()

            for i, r in enumerate(rec_list):
                arr[r[2]:r[3], r[0]:r[1]] = 1
                poly |= rclosedopen(*r)

            arr_used_rectangles = set(get_maximal_rectangles_from_numpy(arr == 0))
            arr_free_rectangles = set(get_maximal_rectangles_from_numpy(arr))

            poly_used_rectangles = set([
                (r.enclosing_x_interval.lower, r.enclosing_x_interval.upper,
                 r.enclosing_y_interval.lower, r.enclosing_y_interval.upper)
                for r in poly.maximal_used_rectangles()
            ])
            poly_free_rectangles = set([
                (r.enclosing_x_interval.lower, r.enclosing_x_interval.upper,
                 r.enclosing_y_interval.lower, r.enclosing_y_interval.upper)
                for r in (poly | (~ropen(0, x_max, 0, y_max))).maximal_free_rectangles()])

            def matrix_to_str(usage_arr: np.array):
                msg = ""
                msg += "   "
                for i in range(0, usage_arr.shape[1], 3):
                    msg += f"{i:<3}"
                msg += "\n"
                msg += "  +" + ("-" * usage_arr.shape[1]) + "+\n"
                for i in range(usage_arr.shape[0]):
                    msg += f"{i:>2}|"
                    for j in range(usage_arr.shape[1]):
                        msg += "x" if usage_arr[i, j] else " "
                    msg += "|\n"
                msg += "  +" + ("-" * usage_arr.shape[1]) + "+"
                msg += "\n"
                return msg

            def difference_string(usage_arr: np.ndarray, rectangles: set[tuple[int, int, int, int]]):
                msg = "Used areas\n"
                msg += matrix_to_str(usage_arr)
                msg += "\n"
                msg += "Rectangles\n"
                for x0, x1, y0, y1 in rectangles:
                    r_arr = np.zeros(usage_arr.shape)
                    r_arr[y0:y1, x0:x1] = 1
                    msg += matrix_to_str(r_arr) + "\n"
                    msg += "\n"
                return msg

            diff_used_1 = poly_used_rectangles - arr_used_rectangles
            if diff_used_1 != set():
                msg = ("The class RPolygon provided a used rectangle which has not been found by "
                       "the array algorithm.\n")
                msg += difference_string(arr, diff_used_1)
                msg += f"The following rectangles have been added in the given order: {rec_list}."
                self.fail(msg)

            diff_used_2 = arr_used_rectangles - poly_used_rectangles
            if diff_used_2 != set():
                msg = ("The class RPolygon failed to provide a used rectangle which has been found by "
                       "the array algorithm.\n")
                msg += difference_string(arr, diff_used_2)
                msg += f"The following rectangles have been added in the given order: {rec_list}."
                self.fail(msg)

            diff_free_1 = poly_free_rectangles - arr_free_rectangles
            if diff_free_1 != set():
                msg = "The class RPolygon provided a used rectangle which has not been found by" \
                      "the array algorithm.\n"
                msg += difference_string(arr, diff_free_1)
                msg += f"The following rectangles have been added in the given order: {rec_list}."
                self.fail(msg)

            diff_free_2 = arr_free_rectangles - poly_free_rectangles
            if diff_free_2 != set():
                msg = "The class RPolygon failed to provide a used rectangle which has been found by" \
                      "the array algorithm.\n"
                msg += f"The following rectangles have been added in the given order: {rec_list}."
                msg += difference_string(arr, diff_free_2)
                self.fail(msg)

            print("Tested maximal rectangle calculation for the following poylgon.")
            print(matrix_to_str(arr))
