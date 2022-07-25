from copy import deepcopy
from typing import Iterator

import portion as P
from portion.interval import Interval, open, empty, closedopen, closed, Atomic
from sortedcontainers import SortedList


def _sub_contained_atomics(int_1: Interval, int_2: Interval) -> Interval:
    """
    Remove all atomic sectors from int_1 which are contained completely in int_2, i.e.


         int_1:     ---     ----
         int_2: ----------    ---------
        return:             ----
    """
    ret_int = empty()
    for y_atomic in int_1:
        if not int_2.contains(y_atomic):
            ret_int |= y_atomic
    return ret_int


class RBoundary:
    """
    This class represents a right boundary and is used internally by the class RPolygon. A comparison operator
    is defined, such that

        (1, CLOSED) < (2, OPEN)
         (2, OPEN)  < (2, CLOSED)

    etc. holds.
    """

    __slots__ = ("val", "btype")

    def __init__(self, val, btype):
        self.val = val
        self.btype = btype

    def __bool__(self):
        raise ValueError("The truth value of a bound is ambiguous.")

    def __invert__(self):
        return RBoundary(self.val, ~self.btype)

    def __eq__(self, other):
        return self.val == other.val and self.btype == other.btype

    def __le__(self, other):
        if self.val != other.val:
            return self.val <= other.val
        elif self.btype == P.OPEN:
            return True
        else:
            return other.btype == P.CLOSED

    def __lt__(self, other):
        if self.val != other.val:
            return self.val < other.val
        elif self.btype == P.OPEN:
            return other.btype == P.CLOSED
        else:
            return False

    def __repr__(self):
        return f"({self.val}, {self.btype})"


def traverse_diagonally(boundaries: list[RBoundary],
                        interval_triangle: list[list[Interval]]) -> Iterator[tuple[Interval, Interval]]:
    """
    Traverse self._free_y_ranges diagonally from the top left to obtain all maximal free rectangles.

    I.e. if we have boundaries

        [-inf, b1, b2, ..., bn, inf]

    then we traverse interval_triangle like shown below.

               | inf  bn  b(n-1)   ...     b1
        -------+------------------------------
          -inf |  1    3     6     ...    _16_/
            b1 |  2    5     :       _15_/
            b2 |  4    :     :  _14_/
             : |  :    :   _13_/
             : |  :   _12_/
            bn |_11__/
    """
    n = len(interval_triangle)
    next_parent_y_intervals = [empty(), empty()]
    for start_row in range(n):
        parent_y_intervals = next_parent_y_intervals
        next_parent_y_intervals = [empty()]
        for row, col in enumerate(range(start_row, -1, -1)):
            curr_y_int = interval_triangle[row][col]

            adj_curr_y_int0 = _sub_contained_atomics(curr_y_int, parent_y_intervals[row])
            adj_curr_y_int = _sub_contained_atomics(adj_curr_y_int0, parent_y_intervals[row + 1])

            l_bound = ~boundaries[row]
            r_bound = boundaries[-col - 1]
            curr_x_int = Interval.from_atomic(l_bound.btype, l_bound.val,
                                              r_bound.val, r_bound.btype)
            if not curr_x_int.empty and not adj_curr_y_int.empty:
                yield curr_x_int, adj_curr_y_int

            next_parent_y_intervals.append(adj_curr_y_int
                                           | parent_y_intervals[row]
                                           | parent_y_intervals[row + 1])
        next_parent_y_intervals.append(empty())


class RPolygon:
    """
    This class represents an orthogonal Polygon.

    An orthogonal polygon is a union of rectangles. A rectangle is represented by two
    intervals of the `portion`-library.
    """
    __slots__ = ("_x_boundaries", "_used_y_ranges", "_free_y_ranges")

    def __init__(self):
        # [-inf, b1, b2, ..., bn, inf]
        self._x_boundaries: SortedList[RBoundary] = SortedList([
            RBoundary(-P.inf, P.OPEN),
            RBoundary(P.inf, P.OPEN)
        ])

        # For every pair of boundaries in x-dimension [-inf, b1, b2, ..., bn, inf] we store the used/free interval
        # of in y-dimension. This results in two array-like data structures in triangular for.
        #
        #         | inf  bn  ... b2  b1
        #    -----+---------------------
        #    -inf |                 ___/
        #      b1 |             ___/
        #      b2 |         ___/
        #       : |      __/
        #      bn |_____/
        self._used_y_ranges: list[list[Interval]] = [[empty()]]
        self._free_y_ranges: list[list[Interval]] = [[open(-P.inf, P.inf)]]

    @classmethod
    def from_interval_product(cls, x_interval: Interval, y_interval: Interval) -> 'RPolygon':
        """
        Create a (simple) rectangular polygon as the product of two intervals.

        :param x_interval:
        :param y_interval:
        :return RPolygon:
        """
        instance = cls()
        instance._add_interval_product(x_interval, y_interval)
        return instance

    def maximal_used_rectangles(self) -> Iterator[tuple[Interval, Interval]]:
        """
        Yield all maximal rectangle which are contained in this polygon.

        A rectangle is called maximal if it is not contained in any larger rectangle of this polygon.

        :return: Iterator[(Interval, Interval)]:
        """
        return traverse_diagonally(list(self._x_boundaries), self._used_y_ranges)

    def maximal_free_rectangles(self) -> (Interval, Interval):
        """
        Yield all maximal rectangle which do not intersect the polygon uniquely.

        A rectangle is called maximal if it is not contained in any larger rectangle not intersecting this polygon.

        :return: Iterator[(Interval, Interval)]:
        """
        return traverse_diagonally(list(self._x_boundaries), self._free_y_ranges)

    def __add__(self, other) -> 'RPolygon':
        rec_copy = deepcopy(self)
        for x_interval, y_interval in other.maximal_used_rectangles(self):
            rec_copy._add_interval_product(x_interval, y_interval)
        return rec_copy

    def _add_interval_product(self, x_interval: Interval, y_interval: Interval):
        if not x_interval.empty and not y_interval.empty:
            for x_atom in x_interval:
                x_atom = Atomic(x_atom.left, x_atom.lower, x_atom.upper, x_atom.right)
                self._add_atomic(x_atom, y_interval)

    def _add_atomic(self, add_x_atom: Atomic, add_y_interval: Interval):
        """
        Update self._used_y_ranges and self._free_y_ranges such that they represent the union of the old rectangular
        polygon and the provided interval product.
        """
        add_x_interval = Interval.from_atomic(add_x_atom.left, add_x_atom.lower,
                                              add_x_atom.upper, add_x_atom.right)

        if add_x_interval.empty or add_y_interval.empty:
            return

        # Left and right boundary of the x_interval to be added
        l_bound_type = ~add_x_atom.left if add_x_atom.lower != -P.inf else P.OPEN
        x_b_left = RBoundary(add_x_atom.lower, l_bound_type)
        r_bound_type = add_x_atom.right if add_x_atom.upper != P.inf else P.OPEN
        x_b_right = RBoundary(add_x_atom.upper, r_bound_type)

        # (a) Add new rows/columns to self._used_y_ranges and self._free_y_ranges
        def extend_ranges_mat(mat: list[list[Interval]], lbound: RBoundary):
            n = len(mat)
            index = self._x_boundaries.bisect_right(lbound)
            prev_y_int_row = mat[index - 1]
            mat.insert(index - 1, list(prev_y_int_row))
            for row in range(0, index):
                prev_y_int = mat[row][n - index]
                mat[row].insert(n - index, prev_y_int)

        if x_b_left not in self._x_boundaries:
            extend_ranges_mat(self._used_y_ranges, x_b_left)
            extend_ranges_mat(self._free_y_ranges, x_b_left)
            self._x_boundaries.add(x_b_left)
        if x_b_right not in self._x_boundaries:
            extend_ranges_mat(self._used_y_ranges, x_b_right)
            extend_ranges_mat(self._free_y_ranges, x_b_right)
            self._x_boundaries.add(x_b_right)

        # (b) Update the y interval in self._used_y_ranges and self._free_y_ranges
        n = len(self._used_y_ranges)
        for i in range(n):
            for j in range(n - i):
                l_bound = self._x_boundaries[i]
                r_bound = self._x_boundaries[-j - 1]
                x_interval = Interval.from_atomic(~l_bound.btype, l_bound.val, r_bound.val, r_bound.btype)

                if add_x_interval.contains(x_interval):
                    self._used_y_ranges[i][j] |= add_y_interval
                else:
                    adj_add_x_interval = add_x_interval

                    add_x_int_l_bound = RBoundary(add_x_interval.lower, ~add_x_interval.left)
                    add_y_int_r_bound = RBoundary(add_x_interval.upper, add_x_interval.right)

                    left_ind = n - 1 - self._x_boundaries.bisect_left(add_x_int_l_bound)
                    if l_bound < add_x_int_l_bound and (0 < left_ind < n - i):
                        # left = self._x_boundaries[left_ind]
                        adj_add_x_interval |= Interval.from_atomic(~l_bound.btype, l_bound.val,
                                                                   adj_add_x_interval.upper, adj_add_x_interval.right)
                        y_int_left = self._used_y_ranges[i][left_ind]
                    else:
                        y_int_left = open(-P.inf, P.inf)

                    right_ind = self._x_boundaries.bisect_left(add_y_int_r_bound)
                    if add_y_int_r_bound < r_bound and right_ind < n - j:
                        adj_add_x_interval |= Interval.from_atomic(adj_add_x_interval.left, adj_add_x_interval.lower,
                                                                   r_bound.val, r_bound.btype)
                        y_int_right = self._used_y_ranges[right_ind][j]
                    else:
                        y_int_right = open(-P.inf, P.inf)

                    if adj_add_x_interval.contains(x_interval):
                        self._used_y_ranges[i][j] |= y_int_left & add_y_interval & y_int_right

        n = len(self._free_y_ranges)
        for i in range(n):
            for j in range(n - i):
                l_bound = ~self._x_boundaries[i]
                r_bound = self._x_boundaries[-j - 1]
                x_interval = Interval.from_atomic(l_bound.btype, l_bound.val, r_bound.val, r_bound.btype)
                if add_x_interval.overlaps(x_interval):
                    self._free_y_ranges[i][j] -= add_y_interval


def rempty() -> RPolygon:
    """
    Create an empty rectangular polygon.

    :return RPolygon:
    """
    return RPolygon.from_interval_product(empty(), empty())


def rclosed(x_lower, x_upper, y_lower, y_upper) -> RPolygon:
    """
    Create a closed rectangular polygon.

    :param x_lower: value of the lower bound in x-dimension
    :param x_upper: value of the upper bound in x-dimension
    :param y_lower: value of the lower bound in y-dimension
    :param y_upper: value of the upper bound in y-dimension
    :return RPolygon:
    """
    return RPolygon.from_interval_product(closed(x_lower, x_upper), closed(y_lower, y_upper))


def rclosedopen(x_lower, x_upper, y_lower, y_upper) -> RPolygon:
    """
    Create a rectangular polygon which is closed on both lower bounds
    and open on both upper bounds.

    :param x_lower: value of the lower bound in x-dimension
    :param x_upper: value of the upper bound in x-dimension
    :param y_lower: value of the lower bound in y-dimension
    :param y_upper: value of the upper bound in y-dimension
    :return RPolygon:
    """
    return RPolygon.from_interval_product(closedopen(x_lower, x_upper), closedopen(y_lower, y_upper))
