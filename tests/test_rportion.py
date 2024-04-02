import unittest
from copy import copy
from itertools import permutations, combinations
import random

import numpy as np
import portion as P
from typing import List, Set, Tuple
from rportion import RPolygon
from rportion.rportion import ropen, rclosedopen, rclosed, rempty, RBoundary, rsingleton, ropenclosed, \
    _traverse_diagonally, rproduct

from tests.helpers import get_maximal_rectangles_from_numpy


def data_tree_to_string(x_boundaries: List[RBoundary],
                        y_intervals: List[List[P.Interval]],
                        spacing: int):
    n = len(y_intervals)
    msg = " " * spacing + "|"
    for x_b in x_boundaries[-1:0:-1]:
        msg += f"{str(x_b):>{spacing}}"
    msg += "\n" + f"-" * spacing + "+"
    for i in range(n):
        msg += f"-" * spacing
    msg += "\n"
    for i, row in enumerate(y_intervals):
        x_b = x_boundaries[i]
        msg += f"{str(~x_b):>{spacing}}|"
        for val in row:
            msg += f"{str(val):>{spacing}}"
        msg += "\n"
    return msg


def print_rpolygon(rpoly: RPolygon, show_trees=False, spacing=20):
    print("USED")
    if show_trees:
        msg = data_tree_to_string(list(rpoly._x_boundaries), rpoly._used_y_ranges, spacing)
        print(msg)
    for i, e in enumerate(rpoly.maximal_rectangles(), start=1):
        print(i, "-->", e)
    print("FREE")
    if show_trees:
        msg = data_tree_to_string(list(rpoly._x_boundaries), rpoly._free_y_ranges, spacing)
        print(msg)
    for i, e in enumerate((~rpoly).maximal_rectangles(), start=1):
        print(i, "-->", e)


class TestRBoundary(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestRBoundary, self).__init__(*args, **kwargs)

    def test_init(self):
        r = RBoundary(0, P.OPEN)
        self.assertEqual(r.val, 0)
        self.assertEqual(r.btype, P.OPEN)
        r = RBoundary(1, P.CLOSED)
        self.assertEqual(r.val, 1)
        self.assertEqual(r.btype, P.CLOSED)
        r = RBoundary(-P.inf, P.CLOSED)
        self.assertEqual(r.val, -P.inf)
        self.assertEqual(r.btype, P.OPEN)
        r = RBoundary(P.inf, P.CLOSED)
        self.assertEqual(r.val, P.inf)
        self.assertEqual(r.btype, P.OPEN)

    def test_bool(self):
        with self.assertRaises(TypeError):
            r = RBoundary(0, P.OPEN)
            bool(r)

    def test_invert(self):
        r = RBoundary(0, P.OPEN)
        self.assertEqual(~r, RBoundary(0, P.CLOSED))
        r = RBoundary(1, P.CLOSED)
        self.assertEqual(~r, RBoundary(1, P.OPEN))

    def test_equal(self):
        self.assertEqual(RBoundary(0, P.OPEN), RBoundary(0, P.OPEN))
        self.assertNotEqual(RBoundary(0, P.OPEN), RBoundary(0, P.CLOSED))
        self.assertNotEqual(RBoundary(0, P.CLOSED), RBoundary(0, P.OPEN))
        self.assertEqual(RBoundary(0, P.CLOSED), RBoundary(0, P.CLOSED))

        self.assertNotEqual(RBoundary(0, P.OPEN), RBoundary(1, P.OPEN))
        self.assertNotEqual(RBoundary(0, P.OPEN), RBoundary(1, P.CLOSED))
        self.assertNotEqual(RBoundary(0, P.CLOSED), RBoundary(1, P.OPEN))
        self.assertNotEqual(RBoundary(0, P.CLOSED), RBoundary(1, P.CLOSED))

        self.assertNotEqual(RBoundary(1, P.OPEN), RBoundary(0, P.OPEN))
        self.assertNotEqual(RBoundary(1, P.OPEN), RBoundary(0, P.CLOSED))
        self.assertNotEqual(RBoundary(1, P.CLOSED), RBoundary(0, P.OPEN))
        self.assertNotEqual(RBoundary(1, P.CLOSED), RBoundary(0, P.CLOSED))

        self.assertEqual(RBoundary(1, P.OPEN), RBoundary(1, P.OPEN))
        self.assertNotEqual(RBoundary(1, P.OPEN), RBoundary(1, P.CLOSED))
        self.assertNotEqual(RBoundary(1, P.CLOSED), RBoundary(1, P.OPEN))
        self.assertEqual(RBoundary(1, P.CLOSED), RBoundary(1, P.CLOSED))

    def test_le(self):
        self.assertTrue(RBoundary(0, P.OPEN) <= RBoundary(0, P.OPEN))
        self.assertTrue(RBoundary(0, P.OPEN) <= RBoundary(0, P.CLOSED))
        self.assertFalse(RBoundary(0, P.CLOSED) <= RBoundary(0, P.OPEN))
        self.assertTrue(RBoundary(0, P.CLOSED) <= RBoundary(0, P.CLOSED))

        self.assertTrue(RBoundary(0, P.OPEN) <= RBoundary(1, P.OPEN))
        self.assertTrue(RBoundary(0, P.OPEN) <= RBoundary(1, P.CLOSED))
        self.assertTrue(RBoundary(0, P.CLOSED) <= RBoundary(1, P.OPEN))
        self.assertTrue(RBoundary(0, P.CLOSED) <= RBoundary(1, P.CLOSED))

        self.assertFalse(RBoundary(1, P.OPEN) <= RBoundary(0, P.OPEN))
        self.assertFalse(RBoundary(1, P.OPEN) <= RBoundary(0, P.CLOSED))
        self.assertFalse(RBoundary(1, P.CLOSED) <= RBoundary(0, P.OPEN))
        self.assertFalse(RBoundary(1, P.CLOSED) <= RBoundary(0, P.CLOSED))

        self.assertTrue(RBoundary(1, P.OPEN) <= RBoundary(1, P.OPEN))
        self.assertTrue(RBoundary(1, P.OPEN) <= RBoundary(1, P.CLOSED))
        self.assertFalse(RBoundary(1, P.CLOSED) <= RBoundary(1, P.OPEN))
        self.assertTrue(RBoundary(1, P.CLOSED) <= RBoundary(1, P.CLOSED))

    def test_lt(self):
        self.assertFalse(RBoundary(0, P.OPEN) < RBoundary(0, P.OPEN))
        self.assertTrue(RBoundary(0, P.OPEN) < RBoundary(0, P.CLOSED))
        self.assertFalse(RBoundary(0, P.CLOSED) < RBoundary(0, P.OPEN))
        self.assertFalse(RBoundary(0, P.CLOSED) < RBoundary(0, P.CLOSED))

        self.assertTrue(RBoundary(0, P.OPEN) < RBoundary(1, P.OPEN))
        self.assertTrue(RBoundary(0, P.OPEN) < RBoundary(1, P.CLOSED))
        self.assertTrue(RBoundary(0, P.CLOSED) < RBoundary(1, P.OPEN))
        self.assertTrue(RBoundary(0, P.CLOSED) < RBoundary(1, P.CLOSED))

        self.assertFalse(RBoundary(1, P.OPEN) < RBoundary(0, P.OPEN))
        self.assertFalse(RBoundary(1, P.OPEN) < RBoundary(0, P.CLOSED))
        self.assertFalse(RBoundary(1, P.CLOSED) < RBoundary(0, P.OPEN))
        self.assertFalse(RBoundary(1, P.CLOSED) < RBoundary(0, P.CLOSED))

        self.assertFalse(RBoundary(1, P.OPEN) < RBoundary(1, P.OPEN))
        self.assertTrue(RBoundary(1, P.OPEN) < RBoundary(1, P.CLOSED))
        self.assertFalse(RBoundary(1, P.CLOSED) < RBoundary(1, P.OPEN))
        self.assertFalse(RBoundary(1, P.CLOSED) < RBoundary(1, P.CLOSED))

    def test_ge(self):
        self.assertTrue(RBoundary(0, P.OPEN) >= RBoundary(0, P.OPEN))
        self.assertFalse(RBoundary(0, P.OPEN) >= RBoundary(0, P.CLOSED))
        self.assertTrue(RBoundary(0, P.CLOSED) >= RBoundary(0, P.OPEN))
        self.assertTrue(RBoundary(0, P.CLOSED) >= RBoundary(0, P.CLOSED))

        self.assertFalse(RBoundary(0, P.OPEN) >= RBoundary(1, P.OPEN))
        self.assertFalse(RBoundary(0, P.OPEN) >= RBoundary(1, P.CLOSED))
        self.assertFalse(RBoundary(0, P.CLOSED) >= RBoundary(1, P.OPEN))
        self.assertFalse(RBoundary(0, P.CLOSED) >= RBoundary(1, P.CLOSED))

        self.assertTrue(RBoundary(1, P.OPEN) >= RBoundary(0, P.OPEN))
        self.assertTrue(RBoundary(1, P.OPEN) >= RBoundary(0, P.CLOSED))
        self.assertTrue(RBoundary(1, P.CLOSED) >= RBoundary(0, P.OPEN))
        self.assertTrue(RBoundary(1, P.CLOSED) >= RBoundary(0, P.CLOSED))

        self.assertTrue(RBoundary(1, P.OPEN) >= RBoundary(1, P.OPEN))
        self.assertFalse(RBoundary(1, P.OPEN) >= RBoundary(1, P.CLOSED))
        self.assertTrue(RBoundary(1, P.CLOSED) >= RBoundary(1, P.OPEN))
        self.assertTrue(RBoundary(1, P.CLOSED) >= RBoundary(1, P.CLOSED))

    def test_gt(self):
        self.assertFalse(RBoundary(0, P.OPEN) > RBoundary(0, P.OPEN))
        self.assertFalse(RBoundary(0, P.OPEN) > RBoundary(0, P.CLOSED))
        self.assertTrue(RBoundary(0, P.CLOSED) > RBoundary(0, P.OPEN))
        self.assertFalse(RBoundary(0, P.CLOSED) > RBoundary(0, P.CLOSED))

        self.assertFalse(RBoundary(0, P.OPEN) > RBoundary(1, P.OPEN))
        self.assertFalse(RBoundary(0, P.OPEN) > RBoundary(1, P.CLOSED))
        self.assertFalse(RBoundary(0, P.CLOSED) > RBoundary(1, P.OPEN))
        self.assertFalse(RBoundary(0, P.CLOSED) > RBoundary(1, P.CLOSED))

        self.assertTrue(RBoundary(1, P.OPEN) > RBoundary(0, P.OPEN))
        self.assertTrue(RBoundary(1, P.OPEN) > RBoundary(0, P.CLOSED))
        self.assertTrue(RBoundary(1, P.CLOSED) > RBoundary(0, P.OPEN))
        self.assertTrue(RBoundary(1, P.CLOSED) > RBoundary(0, P.CLOSED))

        self.assertFalse(RBoundary(1, P.OPEN) > RBoundary(1, P.OPEN))
        self.assertFalse(RBoundary(1, P.OPEN) > RBoundary(1, P.CLOSED))
        self.assertTrue(RBoundary(1, P.CLOSED) > RBoundary(1, P.OPEN))
        self.assertFalse(RBoundary(1, P.CLOSED) > RBoundary(1, P.CLOSED))

    def test_copy(self):
        r = RBoundary(0, P.OPEN)
        print(copy(r))
        self.assertEqual(r, copy(r))
        self.assertNotEqual(id(r), id(copy(r)))


class TestRPolygonConstructors(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestRPolygonConstructors, self).__init__(*args, **kwargs)

    def test_rsingleton(self):
        x = 0
        y = 3
        p = rsingleton(x, y)
        self.assertListEqual(list(p._x_boundaries), [
            RBoundary(-P.inf, P.OPEN),
            RBoundary(x, P.OPEN),
            RBoundary(x, P.CLOSED),
            RBoundary(P.inf, P.OPEN),
        ])
        self.assertListEqual(p._used_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.singleton(y)],
            [P.empty()],
        ])
        self.assertListEqual(p._free_y_ranges, [
            [P.open(-P.inf, y) | P.open(y, P.inf), P.open(-P.inf, y) | P.open(y, P.inf), P.open(-P.inf, P.inf)],
            [P.open(-P.inf, y) | P.open(y, P.inf), P.open(-P.inf, y) | P.open(y, P.inf)],
            [P.open(-P.inf, P.inf)],
        ])

    def test_rempty(self):
        p = rempty()
        self.assertListEqual(list(p._x_boundaries), [
            RBoundary(-P.inf, P.OPEN),
            RBoundary(P.inf, P.OPEN),
        ])
        self.assertListEqual(p._used_y_ranges, [[P.empty()]])
        self.assertListEqual(p._free_y_ranges, [[P.open(-P.inf, P.inf)]])

    def test_ropen(self):
        x1, x2 = 0, 3
        y1, y2 = 1, 2
        p = ropen(x1, x2, y1, y2)
        self.assertListEqual(list(p._x_boundaries), [
            RBoundary(-P.inf, P.OPEN),
            RBoundary(x1, P.CLOSED),
            RBoundary(x2, P.OPEN),
            RBoundary(P.inf, P.OPEN),
        ])
        self.assertListEqual(p._used_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.open(y1, y2)],
            [P.empty()],
        ])
        self.assertListEqual(p._free_y_ranges, [
            [P.openclosed(-P.inf, y1) | P.closedopen(y2, P.inf),
             P.openclosed(-P.inf, y1) | P.closedopen(y2, P.inf),
             P.open(-P.inf, P.inf)],
            [P.openclosed(-P.inf, y1) | P.closedopen(y2, P.inf), P.openclosed(-P.inf, y1) | P.closedopen(y2, P.inf)],
            [P.open(-P.inf, P.inf)],
        ])

    def test_rclosed(self):
        x1, x2 = 0, 3
        y1, y2 = 1, 2
        p = rclosed(x1, x2, y1, y2)
        self.assertListEqual(list(p._x_boundaries), [
            RBoundary(-P.inf, P.OPEN),
            RBoundary(x1, P.OPEN),
            RBoundary(x2, P.CLOSED),
            RBoundary(P.inf, P.OPEN),
        ])
        self.assertListEqual(p._used_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.closed(y1, y2)],
            [P.empty()],
        ])
        self.assertListEqual(p._free_y_ranges, [
            [P.open(-P.inf, y1) | P.open(y2, P.inf), P.open(-P.inf, y1) | P.open(y2, P.inf), P.open(-P.inf, P.inf)],
            [P.open(-P.inf, y1) | P.open(y2, P.inf), P.open(-P.inf, y1) | P.open(y2, P.inf)],
            [P.open(-P.inf, P.inf)],
        ])

    def test_ropenclosed(self):
        x1, x2 = 0, 3
        y1, y2 = 1, 2
        p = ropenclosed(x1, x2, y1, y2)
        self.assertListEqual(list(p._x_boundaries), [
            RBoundary(-P.inf, P.OPEN),
            RBoundary(x1, P.CLOSED),
            RBoundary(x2, P.CLOSED),
            RBoundary(P.inf, P.OPEN),
        ])
        self.assertListEqual(p._used_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.openclosed(y1, y2)],
            [P.empty()],
        ])
        self.assertListEqual(p._free_y_ranges, [
            [P.openclosed(-P.inf, y1) | P.open(y2, P.inf), P.openclosed(-P.inf, y1) | P.open(y2, P.inf), P.open(-P.inf, P.inf)],
            [P.openclosed(-P.inf, y1) | P.open(y2, P.inf), P.openclosed(-P.inf, y1) | P.open(y2, P.inf)],
            [P.open(-P.inf, P.inf)],
        ])

    def test_rclosedopen(self):
        x1, x2 = 0, 3
        y1, y2 = 1, 2
        p = rclosedopen(x1, x2, y1, y2)
        self.assertListEqual(list(p._x_boundaries), [
            RBoundary(-P.inf, P.OPEN),
            RBoundary(x1, P.OPEN),
            RBoundary(x2, P.OPEN),
            RBoundary(P.inf, P.OPEN),
        ])
        self.assertListEqual(p._used_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.closedopen(y1, y2)],
            [P.empty()],
        ])
        self.assertListEqual(p._free_y_ranges, [
            [P.open(-P.inf, y1) | P.closedopen(y2, P.inf), P.open(-P.inf, y1) | P.closedopen(y2, P.inf), P.open(-P.inf, P.inf)],
            [P.open(-P.inf, y1) | P.closedopen(y2, P.inf), P.open(-P.inf, y1) | P.closedopen(y2, P.inf)],
            [P.open(-P.inf, P.inf)],
        ])


class TestRPolygonProperties(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestRPolygonProperties, self).__init__(*args, **kwargs)

    def test_empty(self):
        p1 = rclosed(0, 3, 0, 3)
        p2 = ropen(1, 2, 1, 2)
        self.assertTrue(rempty().empty)
        self.assertFalse(p1.empty)
        self.assertFalse(p2.empty)
        self.assertTrue((p2 - p1).empty)
        self.assertFalse((p1 - p2).empty)

    def test_atomic(self):
        p1 = rclosed(0, 3, 0, 3)
        p2 = ropen(1, 2, 1, 2)
        self.assertTrue(rempty().atomic)
        self.assertTrue(p1.atomic)
        self.assertTrue(p2.atomic)
        self.assertFalse((p1 - p2).atomic)

    def test_enclosure(self):
        x1, x2, x3, x4 = 0, 1, 2, 3
        y1, y2, y3 = 0, 1, 2
        self.assertEqual(
            (rclosed(x1, x3, y1, y3) | rclosed(x2, x4, y1, y2)).enclosure,
            rclosed(x1, x4, y1, y3)
        )

    def test_x_enclosure_interval(self):
        x1, x2, x3, x4 = 0, 1, 2, 3
        y1, y2, y3 = 0, 1, 2
        self.assertEqual(
            (rclosed(x1, x3, y1, y3) | rclosed(x2, x4, y1, y2)).x_enclosure_interval,
            P.closed(x1, x4)
        )
        self.assertEqual(
            (ropen(-P.inf, x1, y1, y2) | ropen(x2, P.inf, y1, y2)).x_enclosure_interval,
            P.open(-P.inf, P.inf)
        )

    def test_y_enclosure_interval(self):
        x1, x2, x3, x4 = 0, 1, 2, 3
        y1, y2, y3 = 0, 1, 2
        self.assertEqual(
            (rclosed(x1, x3, y1, y3) | rclosed(x2, x4, y1, y2)).y_enclosure_interval,
            P.closed(y1, y3)
        )

    def test_boundaries(self):
        x1, x2 = 0, 2
        y1, y2 = 0, 3
        p = rclosedopen(x1, x2, y1, y2)
        self.assertEqual(p.x_lower, x1)
        self.assertEqual(p.x_upper, x2)
        self.assertEqual(p.y_lower, y1)
        self.assertEqual(p.y_upper, y2)
        self.assertEqual(p.x_left, P.CLOSED)
        self.assertEqual(p.x_right, P.OPEN)
        self.assertEqual(p.y_left, P.CLOSED)
        self.assertEqual(p.y_right, P.OPEN)


class TestRPolygonClassMethods(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestRPolygonClassMethods, self).__init__(*args, **kwargs)

    def test_from_interval_product(self):
        x1, x2 = 0, 1
        y1, y2, y3, y4 = 0, 2, 3, 4
        p = RPolygon.from_interval_product(P.open(x1, x2), P.openclosed(y1, y2))
        self.assertListEqual(list(p._x_boundaries), [
            RBoundary(-P.inf, P.OPEN),
            RBoundary(x1, P.CLOSED),
            RBoundary(x2, P.OPEN),
            RBoundary(P.inf, P.OPEN),
        ])
        self.assertListEqual(p._used_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.openclosed(y1, y2)],
            [P.empty()],
        ])
        self.assertListEqual(p._free_y_ranges, [
            [P.openclosed(-P.inf, y1) | P.open(y2, P.inf), P.openclosed(-P.inf, y1) | P.open(y2, P.inf), P.open(-P.inf, P.inf)],
            [P.openclosed(-P.inf, y1) | P.open(y2, P.inf), P.openclosed(-P.inf, y1) | P.open(y2, P.inf)],
            [P.open(-P.inf, P.inf)],
        ])
        p = RPolygon.from_interval_product(P.open(x1, x2), P.openclosed(y1, y2) | P.closedopen(y3, y4))
        self.assertListEqual(list(p._x_boundaries), [
            RBoundary(-P.inf, P.OPEN),
            RBoundary(x1, P.CLOSED),
            RBoundary(x2, P.OPEN),
            RBoundary(P.inf, P.OPEN),
        ])
        self.assertListEqual(p._used_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.openclosed(y1, y2) | P.closedopen(y3, y4)],
            [P.empty()],
        ])
        self.assertListEqual(p._free_y_ranges, [
            [P.openclosed(-P.inf, y1) | P.open(y2, y3) | P.closedopen(y4, P.inf),
             P.openclosed(-P.inf, y1) | P.open(y2, y3) | P.closedopen(y4, P.inf),
             P.open(-P.inf, P.inf)],
            [P.openclosed(-P.inf, y1) | P.open(y2, y3) | P.closedopen(y4, P.inf),
             P.openclosed(-P.inf, y1) | P.open(y2, y3) | P.closedopen(y4, P.inf)],
            [P.open(-P.inf, P.inf)],
        ])


class TestRPolygonOperations(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestRPolygonOperations, self).__init__(*args, **kwargs)

    def test_repr(self):
        self.assertEqual(repr(rempty()), "(x=(), y=())")
        self.assertEqual(repr(ropen(0, 2, 0, 1)), "(x=(0,2), y=(0,1))")
        self.assertEqual(repr(rclosed(0, 2, 0, 1)), "(x=[0,2], y=[0,1])")
        self.assertEqual(repr(ropenclosed(0, 2, 0, 1)), "(x=(0,2], y=(0,1])")
        self.assertEqual(repr(rclosedopen(0, 2, 0, 1)), "(x=[0,2), y=[0,1))")
        self.assertEqual(repr(rclosed(0, 2, 0, 2) | rclosed(1, 3, 0, 1)),
                         "(x=[0,3], y=[0,1]) | (x=[0,2], y=[0,2])")

    def test_eq(self):
        self.assertTrue(rempty() == rempty())
        self.assertFalse(rempty() == ropen(0, 1, 0, 1))
        self.assertFalse(ropen(0, 1, 0, 1) == ropen(0, 1, 0, 2))
        self.assertFalse(ropen(0, 1, 0, 1) == ropen(0, 2, 0, 1))

    def test__add_atomic(self):
        # empty x_interval
        poly = rempty()
        x_atom = P.closedopen(1, 0)
        y_interval = P.open(1, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.empty()]])
        self.assertListEqual(poly._free_y_ranges, [[P.open(-P.inf, P.inf)]])

        # empty y_interval
        poly = rempty()
        x_atom = P.closedopen(-P.inf, 1)
        y_interval = P.empty()
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.empty()]])
        self.assertListEqual(poly._free_y_ranges, [[P.open(-P.inf, P.inf)]])

        # add the whole area
        poly = rempty()
        x_atom = P.closedopen(-P.inf, P.inf)
        y_interval = P.open(-P.inf, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.open(-P.inf, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[P.empty()]])

        # add half planes
        # # (a) Right half space
        poly = rempty()
        x_atom = P.closedopen(0, P.inf)
        y_interval = P.open(-P.inf, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [P.empty(), P.empty()],
            [P.open(-P.inf, P.inf)]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [P.empty(), P.open(-P.inf, P.inf)],
            [P.empty()]
        ])
        # # (b) Left half space
        poly = rempty()
        x_atom = P.closedopen(-P.inf, 0)
        y_interval = P.open(-P.inf, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [P.empty(), P.open(-P.inf, P.inf)],
            [P.empty()]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [P.empty(), P.empty()],
            [P.open(-P.inf, P.inf)]
        ])
        # # (c) Upper half space
        poly = rempty()
        x_atom = P.closedopen(-P.inf, P.inf)
        y_interval = P.open(-P.inf, 0)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.open(-P.inf, 0)]])
        self.assertListEqual(poly._free_y_ranges, [[P.closedopen(0, P.inf)]])
        # # (d) Lower half space
        poly = rempty()
        x_atom = P.closedopen(-P.inf, P.inf)
        y_interval = P.closedopen(0, P.inf)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.closedopen(0, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[P.open(-P.inf, 0)]])

        # add a single bounded polygon
        poly = rempty()
        x_atom = P.closedopen(1, 3)
        y_interval = P.closedopen(2, 4)
        poly._add_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.closedopen(2, 4)],
            [P.empty()]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [P.open(-P.inf, 2) | P.closedopen(4, P.inf),
             P.open(-P.inf, 2) | P.closedopen(4, P.inf),
             P.open(-P.inf, P.inf)],
            [P.open(-P.inf, 2) | P.closedopen(4, P.inf),
             P.open(-P.inf, 2) | P.closedopen(4, P.inf)],
            [P.open(-P.inf, P.inf)]
        ])

    def test__sub_atomic(self):
        # Remove rectangle with empty x_interval from the whole plane
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = P.closedopen(1, 0)
        y_interval = P.open(1, P.inf)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.open(-P.inf, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[P.empty()]])

        # Remove rectangle with empty y_interval from the whole plane
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = P.closedopen(-P.inf, 1)
        y_interval = P.empty()
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.open(-P.inf, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[P.empty()]])

        # Remove half planes from the whole plane
        # # (a) Right half space
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = P.closedopen(0, P.inf)
        y_interval = P.open(-P.inf, P.inf)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [P.empty(), P.open(-P.inf, P.inf)],
            [P.empty()]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [P.empty(), P.empty()],
            [P.open(-P.inf, P.inf)]
        ])
        # # (b) Left half space
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = P.closedopen(-P.inf, 0)
        y_interval = P.open(-P.inf, P.inf)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [P.empty(), P.empty()],
            [P.open(-P.inf, P.inf)]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [P.empty(), P.open(-P.inf, P.inf)],
            [P.empty()]
        ])
        # # (c) Upper half space
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = P.closedopen(-P.inf, P.inf)
        y_interval = P.open(-P.inf, 0)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.closedopen(0, P.inf)]])
        self.assertListEqual(poly._free_y_ranges, [[P.open(-P.inf, 0)]])
        # # (d) Lower half space
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = P.closedopen(-P.inf, P.inf)
        y_interval = P.closedopen(0, P.inf)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [[P.open(-P.inf, 0)]])
        self.assertListEqual(poly._free_y_ranges, [[P.closedopen(0, P.inf)]])

        # remove a single bounded rectangle from the plane
        poly = ropen(-P.inf, P.inf, -P.inf, P.inf)
        x_atom = P.closedopen(1, 3)
        y_interval = P.closedopen(2, 4)
        poly._sub_atomic(x_atom, y_interval)
        self.assertListEqual(poly._used_y_ranges, [
            [P.open(-P.inf, 2) | P.closedopen(4, P.inf),
             P.open(-P.inf, 2) | P.closedopen(4, P.inf),
             P.open(-P.inf, P.inf)],
            [P.open(-P.inf, 2) | P.closedopen(4, P.inf),
             P.open(-P.inf, 2) | P.closedopen(4, P.inf)],
            [P.open(-P.inf, P.inf)]
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.closedopen(2, 4)],
            [P.empty()]
        ])

        # Another example
        poly = rclosedopen(-P.inf, 0, -P.inf, P.inf)
        poly._sub_atomic(P.open(-P.inf, P.inf), P.open(-P.inf, 0))
        poly._sub_atomic(P.open(-P.inf, P.inf), P.open(1, P.inf))
        self.assertListEqual(poly._used_y_ranges, [
            [P.empty(), P.closed(0, 1)],
            [P.empty()],
        ])
        self.assertListEqual(poly._free_y_ranges, [
            [P.open(-P.inf, 0) | P.open(1, P.inf), P.open(-P.inf, 0) | P.open(1, P.inf)],
            [P.open(-P.inf, P.inf)],
        ])

    def test_intersection(self):
        p1 = rclosed(0, 2, 0, 2)
        p2 = rclosed(1, 3, 1, 3)
        p3 = rclosed(3, 4, 3, 4)
        self.assertEqual(
            p1.intersection(p2),
            p1 & p2
        )
        self.assertEqual(
            p1 & p2,
            rclosed(1, 2, 1, 2)
        )
        self.assertEqual(
            p1.intersection(p3),
            p1 & p3
        )
        self.assertEqual(
            p1.intersection(p3),
            rempty()
        )

    def test_union(self):
        p1 = rclosed(0, 2, 0, 2)
        p2 = rclosed(1, 3, 1, 3)
        self.assertEqual(
            p1.union(rempty()),
            p1 | rempty()
        )
        self.assertEqual(
            p1 | rempty(),
            p1
        )
        self.assertEqual(
            p1.union(p2),
            p1 | p2
        )
        self.assertListEqual(
            list((p1 | p2).maximal_rectangles()),
            [rclosed(0, 3, 1, 2), rclosed(0, 2, 0, 2), rclosed(1, 3, 1, 3), rclosed(1, 2, 0, 3)]
        )

    def test_complement(self):
        p = rclosed(0, 2, 0, 2)
        self.assertEqual(
            p.complement(),
            ~p
        )
        self.assertListEqual(
            list(p._used_y_ranges),
            list((~p)._free_y_ranges)
        )
        self.assertListEqual(
            list(p._free_y_ranges),
            list((~p)._used_y_ranges)
        )
        self.assertEqual(p, ~(~p))

    def test_difference(self):
        p1 = rclosed(0, 2, 0, 2)
        p2 = rclosed(1, 3, 1, 3)
        self.assertEqual(
            p1.difference(p2),
            p1 - p2
        )
        self.assertListEqual(
            list((p1 - p2).maximal_rectangles()),
            [RPolygon.from_interval_product(P.closed(0, 2), P.closedopen(0, 1)),
             RPolygon.from_interval_product(P.closedopen(0, 1), P.closed(0, 2))]
        )

    def test_rectangle_partition(self):
        x0, x1, x2, x3 = 0, 1, 2, 3
        y0, y1, y2, y3 = 0, 1, 2, 3

        poly = ropen(x0, x1, y0, y1)
        self.assertListEqual(list(poly.rectangle_partitioning()),
                             [ropen(x0, x1, y0, y1)])

        poly = ropen(x0, x2, y0, y2) | ropen(x1, x3, y1, y3)
        self.assertListEqual(list(poly.rectangle_partitioning()),
                             [ropen(x0, x3, y1, y2),
                              RPolygon.from_interval_product(P.open(x0, x2), P.openclosed(y0, y1)),
                              RPolygon.from_interval_product(P.open(x1, x3), P.closedopen(y2, y3))])

    def test_maximal_atomic_x_rectangles(self):
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
        y_interval_1 = P.closedopen(-1, 1)
        y_interval_2 = P.closedopen(0, 2)
        # (1)
        poly = rempty()
        poly._add_atomic(P.closedopen(x1, x2), y_interval_1)
        poly._add_atomic(P.closedopen(x3, x4), y_interval_2)
        self.assertListEqual(list(poly._maximal_atomic_x_rectangles()), [
            (P.closedopen(x1, x2), y_interval_1),
            (P.closedopen(x3, x4), y_interval_2)
        ])
        self.assertListEqual(list((~poly)._maximal_atomic_x_rectangles()), [
            (P.open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2)),
            (P.open(-P.inf, x3),
             P.Interval.from_atomic(~y_interval_1.right, y_interval_1.upper, P.inf, P.OPEN)),
            (P.closedopen(x2, P.inf),
             P.Interval.from_atomic(P.OPEN, -P.inf, y_interval_2.lower, ~y_interval_2.left)),
            (P.open(-P.inf, x1), P.open(-P.inf, P.inf)),
            (P.closedopen(x2, x3), P.open(-P.inf, P.inf)),
            (P.closedopen(x4, P.inf), P.open(-P.inf, P.inf)),
        ])
        # (2)
        poly = rempty()
        poly._add_atomic(P.closedopen(x1, x3), y_interval_1)
        poly._add_atomic(P.closedopen(x2, x4), y_interval_2)
        self.assertListEqual(list(poly._maximal_atomic_x_rectangles()), [
            (P.closedopen(x1, x4), y_interval_1 & y_interval_2),
            (P.closedopen(x1, x3), y_interval_1),
            (P.closedopen(x2, x4), y_interval_2),
            (P.closedopen(x2, x3), y_interval_1 | y_interval_2)
        ])
        self.assertListEqual(list((~poly)._maximal_atomic_x_rectangles()), [
            (P.open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2)),
            (P.open(-P.inf, x2),
             P.Interval.from_atomic(~y_interval_1.right, y_interval_1.upper, P.inf, P.OPEN)),
            (P.closedopen(x3, P.inf),
             P.Interval.from_atomic(P.OPEN, -P.inf, y_interval_2.lower, ~y_interval_2.left)),
            (P.open(-P.inf, x1), P.open(-P.inf, P.inf)),
            (P.closedopen(x4, P.inf), P.open(-P.inf, P.inf)),
        ])
        # (3)
        poly = rempty()
        poly._add_atomic(P.closedopen(x1, x4), y_interval_1)
        poly._add_atomic(P.closedopen(x2, x3), y_interval_2)
        self.assertListEqual(list(poly._maximal_atomic_x_rectangles()), [
            (P.closedopen(x1, x4), y_interval_1),
            (P.closedopen(x2, x3), y_interval_1 | y_interval_2),
        ])
        self.assertListEqual(list((~poly)._maximal_atomic_x_rectangles()), [
            (P.open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2)),
            (P.open(-P.inf, x2),
             P.Interval.from_atomic(~y_interval_1.right, y_interval_1.upper, P.inf, P.OPEN)),
            (P.closedopen(x3, P.inf),
             P.Interval.from_atomic(~y_interval_2.right, y_interval_1.upper, P.inf, P.OPEN)),
            (P.open(-P.inf, x1), P.open(-P.inf, P.inf)),
            (P.closedopen(x4, P.inf), P.open(-P.inf, P.inf)),
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
        y_interval_1 = P.closedopen(4, 6)
        y_interval_2 = P.closedopen(1, 3)
        y_interval_3 = P.closedopen(2, 5)
        x1, x2, x3, x4, x5, x6 = (1, 2, 3, 4, 5, 6)
        # (1)
        poly = rempty()
        poly._add_atomic(P.closedopen(x1, x3), y_interval_1)
        poly._add_atomic(P.closedopen(x4, x6), y_interval_2)
        poly._add_atomic(P.closedopen(x2, x5), y_interval_3)
        self.assertListEqual(list(poly._maximal_atomic_x_rectangles()), [
            (P.closedopen(x1, x5), P.closedopen(x4, x5)),
            (P.closedopen(x2, x6), P.closedopen(x2, x3)),
            (P.closedopen(x2, x5), P.closedopen(x2, x5)),
            (P.closedopen(x1, x3), P.closedopen(x4, x6)),
            (P.closedopen(x4, x6), P.closedopen(x1, x3)),
            (P.closedopen(x2, x3), P.closedopen(x2, x6)),
            (P.closedopen(x4, x5), P.closedopen(x1, x5)),
        ])
        self.assertListEqual(list((~poly)._maximal_atomic_x_rectangles()), [
            (P.open(-P.inf, P.inf), ~(y_interval_1 | y_interval_2 | y_interval_3)),
            (P.open(-P.inf, x4),
             P.Interval.from_atomic(P.OPEN, -P.inf, y_interval_3.lower, ~y_interval_3.left)),
            (P.closedopen(x3, P.inf),
             P.Interval.from_atomic(~y_interval_3.right, y_interval_3.upper, P.inf, P.OPEN)),
            (P.open(-P.inf, x2), P.Interval.from_atomic(P.OPEN, -P.inf, y_interval_1.lower, ~y_interval_1.left)),
            (P.closedopen(x5, P.inf),
             P.Interval.from_atomic(~y_interval_2.right, y_interval_2.upper, P.inf, P.OPEN)),
            (P.open(-P.inf, x1), P.open(-P.inf, P.inf)),
            (P.closedopen(x6, P.inf), P.open(-P.inf, P.inf)),
        ])

        # Test the following polygon addition sequence, where the third rectangle directly touches the second one.
        #
        #    x1   x2         x1   x2    x3   x4       x1   x2   x3   x4   x5
        #    +----+          +----+     +----+        +----+    +---------+
        #    |    |    ->    |    |     |    |    ->  |    |    |         |
        #    +----+          +----+     +----+        +----+    +---------+
        y_interval = P.closedopen(0, 1)
        x1, x2, x3, x4, x5 = 1, 2, 3, 4, 5
        poly = rempty()
        # poly._add_atomic(P.closedopen(x1, x2), y_interval)
        poly._add_atomic(P.closedopen(x3, x4), y_interval)
        poly._add_atomic(P.closedopen(x4, x5), y_interval)

        # Add three sectors in each six possible different orders and test if the result underlying
        # data structure of RPolygon is the same for all.
        #
        #      x1 x2 x3 x4 x5
        #      +--+--+--+--+
        #      |  |  |  |  |
        #      +--+--+--+--+
        x1, x2, x3, x4, x5 = (1, 2, 3, 4, 5)
        y_interval = P.closedopen(0, 1)
        x_lims = [(x1, x3), (x2, x4), (x3, x5)]
        y_intervals = [y_interval, y_interval, y_interval]
        poly_list = []
        for arrangement in permutations(zip(x_lims, y_intervals)):
            poly = rempty()
            for (x_a, x_b), y_int in arrangement:
                poly._add_atomic(P.closedopen(x_a, x_b), y_int)
            poly_list.append(poly)
        for poly_1, poly_2 in combinations(poly_list, 2):
            self.assertListEqual(poly_1._used_y_ranges, poly_2._used_y_ranges)
            self.assertListEqual(poly_1._free_y_ranges, poly_2._free_y_ranges)

    def test_maximal_rectangles_validation(self):
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
                (r.x_enclosure_interval.lower, r.x_enclosure_interval.upper,
                 r.y_enclosure_interval.lower, r.y_enclosure_interval.upper)
                for r in poly.maximal_rectangles()
            ])
            poly_free_rectangles = set([
                (r.x_enclosure_interval.lower, r.x_enclosure_interval.upper,
                 r.y_enclosure_interval.lower, r.y_enclosure_interval.upper)
                for r in (~(poly | (~ropen(0, x_max, 0, y_max)))).maximal_rectangles()])

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

            def difference_string(usage_arr: np.ndarray, rectangles: Set[Tuple[int, int, int, int]]):
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

            print("Tested maximal rectangle calculation for the following polygon.")
            print(matrix_to_str(arr))

    def test_boundary(self):
        # Single rectangle
        #     +--+
        #     |  |
        #     +--+
        x0, x1 = 0, 2
        y0, y1 = 1, 3
        for constructor in [rclosed, ropen, rclosedopen, ropenclosed]:
            poly = constructor(x0, x1, y0, y1)
            self.assertEqual(poly.boundary(),
                             rproduct(P.closed(x0, x1), P.singleton(y0) | P.singleton(y1))
                             | rproduct(P.singleton(x0) | P.singleton(x1), P.closed(y0, y1)))

        # Two overlapping rectangles
        #        +--+
        #        |+-|+
        #        +|-+|
        #         +--+
        x0, x1, x2, x3 = 0, 1, 2, 3
        y0, y1, y2, y3 = 4, 5, 6, 7
        for constructor in [rclosed, ropen, rclosedopen, ropenclosed]:
            poly = constructor(x0, x2, y0, y2) | constructor(x1, x3, y1, y3)
            self.assertEqual(poly.boundary(),
                             rproduct(P.singleton(x0), P.closed(y0, y2))
                             | rproduct(P.closed(x0, x1), P.singleton(y2))
                             | rproduct(P.singleton(x1), P.closed(y2, y3))
                             | rproduct(P.closed(x1, x3), P.singleton(y3))
                             | rproduct(P.singleton(x3), P.closed(y1, y3))
                             | rproduct(P.closed(x2, x3), P.singleton(y1))
                             | rproduct(P.singleton(x2), P.closed(y0, y1))
                             | rproduct(P.closed(x0, x2), P.singleton(y0)))

        # Rectangle with a hole
        #       +------+
        #       | +--+ |
        #       | |  | |
        #       | +--+ |
        #       +------+
        x0, x1, x2, x3 = 0, 1, 2, 3
        y0, y1, y2, y3 = 4, 5, 6, 7
        for constructor in [rclosed, ropen, rclosedopen, ropenclosed]:
            poly = constructor(x0, x3, y0, y3)
            poly -= constructor(x1, x2, y1, y2)
            self.assertEqual(poly.boundary(),
                             rproduct(P.closed(x0, x3), P.singleton(y0) | P.singleton(y3))
                             | rproduct(P.singleton(x0) | P.singleton(x3), P.closed(y0, y3))
                             | rproduct(P.closed(x1, x2), P.singleton(y1) | P.singleton(y2))
                             | rproduct(P.singleton(x1) | P.singleton(x2), P.closed(y1, y2)))


class TestIntervalTreeFunctions(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestIntervalTreeFunctions, self).__init__(*args, **kwargs)

    def test__traverse_diagonally(self):
        next_accumulator = lambda curr, l_parent, r_parent: curr | l_parent | r_parent
        adj_y_interval = lambda curr, l_parent, r_parent: curr

        boundaries = [RBoundary(-P.inf, P.OPEN), RBoundary(P.inf, P.OPEN)]
        interval_tree = [[P.empty()]]
        with self.assertRaises(StopIteration):
            next(_traverse_diagonally(boundaries, interval_tree, next_accumulator, adj_y_interval))

        boundaries = [RBoundary(-P.inf, P.OPEN), RBoundary(0, P.CLOSED),
                      RBoundary(1, P.OPEN), RBoundary(P.inf, P.OPEN)]
        interval_tree = [
            [P.empty(), P.empty(), P.empty()],
            [P.empty(), P.open(0, 1)],
            [P.empty()]
        ]
        iterator = _traverse_diagonally(boundaries, interval_tree, next_accumulator, adj_y_interval)
        self.assertEqual(next(iterator), (P.open(0, 1), P.open(0, 1)))
        with self.assertRaises(StopIteration):
            next(iterator)
