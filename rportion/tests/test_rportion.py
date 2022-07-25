import unittest
from itertools import permutations, combinations

import portion as P
from portion.interval import open, closedopen, Atomic, empty, Interval

from rportion import RPolygon


def print_mat(mat: list[list[Interval]]):
    for row in mat:
        for val in row:
            print(f"{str(val):<25}", end="")
        print()
    print()


def print_rpolygon(rpoly: RPolygon):
    print(rpoly._x_boundaries)
    print_mat(rpoly._used_y_ranges)
    print_mat(rpoly._free_y_ranges)
    print()
    for i, e in enumerate(rpoly.maximal_used_rectangles(), start=1):
        print(i, "-->", e)
    print()
    for i, e in enumerate(rpoly.maximal_free_rectangles(), start=1):
        print(i, "-->", e)


class TestValidation(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestValidation, self).__init__(*args, **kwargs)

    def test_RPolygon_maximum_rectangles(self):
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
        poly = RPolygon()
        poly._add_atomic(Atomic(P.CLOSED, x1, x2, P.OPEN), y_interval_1)
        poly._add_atomic(Atomic(P.CLOSED, x3, x4, P.OPEN), y_interval_2)
        self.assertListEqual(list(poly.maximal_used_rectangles()), [
            (closedopen(x1, x2), y_interval_1),
            (closedopen(x3, x4), y_interval_2)
        ])
        self.assertListEqual(list(poly.maximal_free_rectangles()), [
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
        poly = RPolygon()
        poly._add_atomic(Atomic(P.CLOSED, x1, x3, P.OPEN), y_interval_1)
        poly._add_atomic(Atomic(P.CLOSED, x2, x4, P.OPEN), y_interval_2)
        self.assertListEqual(list(poly.maximal_used_rectangles()), [
            (closedopen(x1, x4), y_interval_1 & y_interval_2),
            (closedopen(x1, x3), y_interval_1),
            (closedopen(x2, x4), y_interval_2),
            (closedopen(x2, x3), y_interval_1 | y_interval_2)
        ])
        self.assertListEqual(list(poly.maximal_free_rectangles()), [
            (open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2)),
            (open(-P.inf, x2),
             Interval.from_atomic(~y_interval_1.right, y_interval_1.upper, P.inf, P.OPEN)),
            (closedopen(x3, P.inf),
             Interval.from_atomic(P.OPEN, -P.inf, y_interval_2.lower, ~y_interval_2.left)),
            (open(-P.inf, x1), open(-P.inf, P.inf)),
            (closedopen(x4, P.inf), open(-P.inf, P.inf)),
        ])
        # (3)
        poly = RPolygon()
        poly._add_atomic(Atomic(P.CLOSED, x1, x4, P.OPEN), y_interval_1)
        poly._add_atomic(Atomic(P.CLOSED, x2, x3, P.OPEN), y_interval_2)
        self.assertListEqual(list(poly.maximal_used_rectangles()), [
            (closedopen(x1, x4), y_interval_1),
            (closedopen(x2, x3), y_interval_1 | y_interval_2),
        ])
        self.assertListEqual(list(poly.maximal_free_rectangles()), [
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
        poly = RPolygon()
        poly._add_atomic(Atomic(P.CLOSED, x2, x5, P.OPEN), y_interval_3)
        poly._add_atomic(Atomic(P.CLOSED, x1, x3, P.OPEN), y_interval_1)
        poly._add_atomic(Atomic(P.CLOSED, x4, x6, P.OPEN), y_interval_2)
        self.assertListEqual(list(poly.maximal_used_rectangles()), [
            (closedopen(x1, x5), closedopen(x4, x5)),
            (closedopen(x2, x6), closedopen(x2, x3)),
            (closedopen(x2, x5), closedopen(x2, x5)),
            (closedopen(x1, x3), closedopen(x4, x6)),
            (closedopen(x4, x6), closedopen(x1, x3)),
            (closedopen(x2, x3), closedopen(x2, x6)),
            (closedopen(x4, x5), closedopen(x1, x5)),
        ])
        self.assertListEqual(list(poly.maximal_free_rectangles()), [
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

        # Add three sectors in each six possible different orders and test if the result underlying
        # data structure of RPolygon stays the same.
        #
        #     x1 x2 x3  x4 x5 x6
        #      +----+
        #      |  +-|------+
        #      +--|-+      |
        #         |        |
        #         |      +-|--+
        #         +------|-+  |
        #                +----+
        x1, x2, x3, x4 = (1, 2, 3, 4)
        y_interval_1 = closedopen(-1, 1)
        y_interval_2 = closedopen(0, 2)
        x_lims = [(x1, x3), (x4, x6), (x2, x5)]
        y_intervals = [y_interval_1, y_interval_2, y_interval_3]
        poly_list = []
        for arrangement in permutations(zip(x_lims, y_intervals)):
            poly = RPolygon()
            for (x_a, x_b), y_int in arrangement:
                poly._add_atomic(Atomic(P.CLOSED, x_a, x_b, P.OPEN), y_int)
            poly_list.append(poly)
        for poly_1, poly_2 in combinations(poly_list, 2):
            self.assertListEqual(poly_1._used_y_ranges, poly_2._used_y_ranges)
            self.assertListEqual(poly_1._free_y_ranges, poly_2._free_y_ranges)

    def test_RPolygon__add_atomic(self):
        # empty x_interval
        poly = RPolygon()
        x_atom = Atomic(P.CLOSED, 1, 0, P.OPEN)
        y_interval = open(1, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[empty()]])
        self.assertListEqual(poly._free_y_ranges, [[open(-P.inf, P.inf)]])

        # empty y_interval
        poly = RPolygon()
        x_atom = Atomic(P.OPEN, -P.inf, 1, P.OPEN)
        y_interval = empty()
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[empty()]])
        self.assertListEqual(poly._free_y_ranges, [[open(-P.inf, P.inf)]])

        # add the whole area
        poly = RPolygon()
        x_atom = Atomic(P.OPEN, -P.inf, P.inf, P.OPEN)
        y_interval = open(-P.inf, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[open(-P.inf, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[empty()]])

        # add half planes
        # # (a) Right half space
        poly = RPolygon()
        x_atom = Atomic(P.CLOSED, 0, P.inf, P.OPEN)
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
        poly = RPolygon()
        x_atom = Atomic(P.OPEN, -P.inf, 0, P.CLOSED)
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
        poly = RPolygon()
        x_atom = Atomic(P.OPEN, -P.inf, P.inf, P.OPEN)
        y_interval = open(-P.inf, 0)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[open(-P.inf, 0)]])
        self.assertListEqual(poly._free_y_ranges, [[closedopen(0, P.inf)]])
        # # (d) Lower half space
        poly = RPolygon()
        x_atom = Atomic(P.OPEN, -P.inf, P.inf, P.OPEN)
        y_interval = closedopen(0, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[closedopen(0, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[open(-P.inf, 0)]])

        # add a single bounded polygon
        poly = RPolygon()
        x_atom = Atomic(P.CLOSED, 1, 3, P.OPEN)
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
