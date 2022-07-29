from copy import copy
from typing import Iterator, List, Tuple

import portion as P
from portion.interval import Interval, open, empty, closedopen, closed, openclosed
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
        if self.val == P.inf or self.val == -P.inf:
            self.btype = P.OPEN
        else:
            self.btype = btype

    def __bool__(self):
        raise ValueError("The truth value of a bound is ambiguous.")

    def __invert__(self):
        if self.val == P.inf or self.val == -P.inf:
            return RBoundary(self.val, P.OPEN)
        else:
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

    def __copy__(self):
        return self.__class__(self.val, self.btype)


def _extend_ranges_mat(mat: List[List[Interval]],
                       x_boundaries: 'SortedList[RBoundary]',
                       rbound: RBoundary):
    n = len(mat)
    index = x_boundaries.bisect_right(rbound)
    prev_y_int_row = mat[index - 1]
    mat.insert(index - 1, list(prev_y_int_row))
    for row in range(0, index):
        prev_y_int = mat[row][n - index]
        mat[row].insert(n - index, prev_y_int)


def _y_interval_triangle_prunable(bound_index: int, y_interval_triangle: List[List[Interval]]) -> bool:
    """
    Test if the given bound can be removed from y_interval_triangle.

    :bound_index int:
    :y_interval_triangle list[list[Interval]]:
    :return bool:

    For the following triangle the boundary 4 can be pruned, because
      - the y-intervals of row (Y) are contained in the corresponding
        y-intervals of row (X), and
      - the y-intervals of column (N) are contained in the corresponding
        y-intervals of column (M).

                              (M)    (N)
                   |  +inf     5      4     3
            -------+-------------------------------
             -inf  |   ()     ()     ()    ()
        (X)     3  |   ()    [0,1)  [0,1)
        (Y)     4  |   ()    [0,1)
                5  |   ()

    In contrast for this triangle the boundary 4 cannot be pruned.

                                  (M)    (N)
                   |  +inf     3      4     5
            -------+-------------------------------
             -inf  |   ()     ()     ()    ()
        (X)     3  |   ()    [0,1)  [0,1)
        (Y)     4  |   ()    [0,2)
                5  |   ()
    """
    n = len(y_interval_triangle)
    if bound_index == 0 or bound_index == n:
        return False

    tr = y_interval_triangle

    i = bound_index
    cols = range(0, n - i)
    for j in cols:
        if not tr[i - 1][j].contains(tr[i][j]):
            return False

    j = n - bound_index
    rows = range(0, i)
    for i in rows:
        if not tr[i][j - 1].contains(tr[i][j]):
            return False

    return True


def _interval_triangle_prunable_indices(y_interval_triangle: List[List[Interval]]) -> List[int]:
    """Return a list of indices which can be removed from y_interval_triangle."""
    n = len(y_interval_triangle)
    return [b_index for b_index in range(1, n+1)
            if _y_interval_triangle_prunable(b_index, y_interval_triangle)]


def _remove_boundaries(bound_indices: List[int], y_interval_triangle: List[List[Interval]]):
    """
    Remove rows and columnss from y_interval_triangle specified by bound_indices.

    :bound_indices list[int]:
    :y_interval_triangle list[list[Interval]]: Data

    This function is meant to be performed if the boundary is prunable, see function
    `_y_interval_triangle_prunable`.

    I.e. for bound_indices == [1] we get

             old y-interval-triangle                         new y-interval-triangle

               |  +inf     5      4     3                       |  +inf     5      3
        -------+----------------------------             -------+----------------------
         -inf  |   ()     ()     ()    ()                  -inf |   ()     ()      ()
            3  |   ()    [0,1)  [0,1)           ==>          3  |   ()    [0,1)
            4  |   ()    [0,1)                               5  |   ()
            5  |   ()
    """
    n = len(y_interval_triangle)
    n_removed = 0
    for bound_index in sorted(bound_indices):
        for i in range(bound_index - n_removed):
            y_interval_triangle[i].pop(n - bound_index)
        y_interval_triangle.pop(bound_index - n_removed)
        n_removed += 1


def _update_x_boundaries_and_y_interval_triangles(
        x_boundaries: 'SortedList[RBoundary]',
        y_interval_triangle_add: List[List[Interval]],
        y_interval_triangle_sub: List[List[Interval]],
        other_x_atom: Interval, other_y_interval: Interval):
    """
    Update
       - y_interval_triangle_add
                such that it represents the union of the old polygon and the one provided as parameters
       - y_interval_triangle_sub
                such that it represents the difference of the old polygon and the one provided as paramters

    :x_boundaries SortedList[RBoundary]:
    :y_interval_triangle_add list[list[Interval]]:
    :y_interval_triangle_sub list[list[Interval]]:
    :new_x_atom Atomic, new_y_interval: Interval:
    """
    assert other_x_atom.atomic

    # Left and right boundary of other_x_atom
    l_bound_type = ~other_x_atom.left
    x_b_left = RBoundary(other_x_atom.lower, l_bound_type)
    r_bound_type = other_x_atom.right
    x_b_right = RBoundary(other_x_atom.upper, r_bound_type)

    # (a) Add new rows/columns to the y-interval-triangles and the x-boundaries
    if x_b_left not in x_boundaries:
        _extend_ranges_mat(y_interval_triangle_add, x_boundaries, x_b_left)
        _extend_ranges_mat(y_interval_triangle_sub, x_boundaries, x_b_left)
        x_boundaries.add(x_b_left)
    if x_b_right not in x_boundaries:
        _extend_ranges_mat(y_interval_triangle_add, x_boundaries, x_b_right)
        _extend_ranges_mat(y_interval_triangle_sub, x_boundaries, x_b_right)
        x_boundaries.add(x_b_right)

    # (b) Update y_interval_triangle_add
    n = len(y_interval_triangle_add)
    for i in range(n):
        for j in range(n - i):
            l_bound = x_boundaries[i]
            r_bound = x_boundaries[-j - 1]
            x_interval = Interval.from_atomic(~l_bound.btype, l_bound.val, r_bound.val, r_bound.btype)

            if other_x_atom.contains(x_interval):
                y_interval_triangle_add[i][j] |= other_y_interval
            else:
                adj_other_x_interval = other_x_atom

                other_x_int_l_bound = ~RBoundary(other_x_atom.lower, other_x_atom.left)
                other_y_int_r_bound = RBoundary(other_x_atom.upper, other_x_atom.right)

                left_ind = x_boundaries.bisect_left(other_x_int_l_bound)
                col = n - left_ind
                if l_bound < other_x_int_l_bound and (0 < col < n - i):
                    adj_other_x_interval |= Interval.from_atomic(~l_bound.btype, l_bound.val,
                                                                 adj_other_x_interval.upper,
                                                                 adj_other_x_interval.right)
                    y_int_left = y_interval_triangle_add[i][col]
                else:
                    y_int_left = open(-P.inf, P.inf)

                right_ind = x_boundaries.bisect_left(other_y_int_r_bound)
                if other_y_int_r_bound < r_bound and right_ind < n - j:
                    adj_other_x_interval |= Interval.from_atomic(adj_other_x_interval.left,
                                                                 adj_other_x_interval.lower,
                                                                 r_bound.val, r_bound.btype)
                    y_int_right = y_interval_triangle_add[right_ind][j]
                else:
                    y_int_right = open(-P.inf, P.inf)

                if adj_other_x_interval.contains(x_interval):
                    y_interval_triangle_add[i][j] |= y_int_left & other_y_interval & y_int_right

    # (c) Update y_interval_triangle_sub
    n = len(y_interval_triangle_sub)
    for i in range(n):
        for j in range(n - i):
            l_bound = ~x_boundaries[i]
            r_bound = x_boundaries[-j - 1]
            x_interval = Interval.from_atomic(l_bound.btype, l_bound.val, r_bound.val, r_bound.btype)
            if other_x_atom.overlaps(x_interval):
                y_interval_triangle_sub[i][j] -= other_y_interval

    # (d) Prune the y-interval-triangles if possible
    indices_add = _interval_triangle_prunable_indices(y_interval_triangle_add)
    indices_sub = _interval_triangle_prunable_indices(y_interval_triangle_sub)
    assert indices_add == indices_sub, (
        "If one y-interval-triangle is prunable then the other should be too"
    )
    _remove_boundaries(indices_add, y_interval_triangle_add)
    _remove_boundaries(indices_add, y_interval_triangle_sub)
    for ind in sorted(indices_add, reverse=True):
        x_boundaries.pop(ind)


def _traverse_diagonally(boundaries: List[RBoundary],
                         interval_triangle: List[List[Interval]]) -> Iterator[Tuple[Interval, Interval]]:
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
        self._used_y_ranges: List[List[Interval]] = [[empty()]]
        self._free_y_ranges: List[List[Interval]] = [[open(-P.inf, P.inf)]]

    @property
    def empty(self):
        try:
            next(self.maximal_used_rectangles())
            return False
        except StopIteration:
            return True

    @property
    def atomic(self):
        try:
            max_rec_iter = self.maximal_used_rectangles()
            next(max_rec_iter)
            next(max_rec_iter)
            return False
        except StopIteration:
            return True

    @property
    def x_enclosure_interval(self) -> Interval:
        """
        smallest y_interval enclosing the y-dimension of the polygon.
        """
        n = len(self._free_y_ranges[0])
        for l_int in range(0, n):
            if self._free_y_ranges[0][l_int] == open(-P.inf, P.inf):
                break
        else:
            l_int += 1
        n = len(self._free_y_ranges)
        for r_int in range(0, n):
            if self._free_y_ranges[r_int][0] == open(-P.inf, P.inf):
                break
        else:
            r_int += 1
        l_bound = self._x_boundaries[n - l_int]
        r_bound = self._x_boundaries[r_int]
        return Interval.from_atomic(~l_bound.btype, l_bound.val,
                                    r_bound.val, r_bound.btype)

    @property
    def y_enclosure_interval(self) -> Interval:
        """
        smallest y_interval enclosing the y-dimension of the polygon.
        """
        return ~self._free_y_ranges[0][0]

    @property
    def enclosure(self) -> 'RPolygon':
        return self.__class__.from_interval_product(self.x_enclosure_interval, self.y_enclosure_interval)

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

    def maximal_used_rectangles(self) -> Iterator['RPolygon']:
        """
        Yield all maximal rectangle which are contained in this polygon.

        A rectangle is called maximal if it is not contained in any larger rectangle of this polygon.

        :return: Iterator[RPolygon]:
        """
        for x_atom, y_interval in self._maximal_used_atomic_x_rectangles():
            for y_atom in y_interval:
                if y_atom.empty:
                    continue
                yield self.__class__.from_interval_product(x_atom, y_atom)

    def __copy__(self):
        poly_copy = self.__class__()
        poly_copy._x_boundaries = SortedList(copy(b) for b in self._x_boundaries)
        poly_copy._used_y_ranges = [
            [y_int for y_int in row]
            for row in self._used_y_ranges
        ]
        poly_copy._free_y_ranges = [
            [y_int for y_int in row]
            for row in self._free_y_ranges
        ]
        return poly_copy

    def _maximal_used_atomic_x_rectangles(self) -> Iterator[Tuple[Interval, Interval]]:
        return _traverse_diagonally(list(self._x_boundaries), self._used_y_ranges)

    def maximal_free_rectangles(self) -> Iterator['RPolygon']:
        """
        Yield all maximal rectangle which do not intersect the polygon uniquely.

        A rectangle is called maximal if it is not contained in any larger rectangle not intersecting this polygon.

        :return: Iterator['RPolygon']:
        """
        for x_atom, y_interval in self._maximal_free_atomic_x_rectangles():
            for y_atom in y_interval:
                if y_atom.empty:
                    continue
                yield self.__class__.from_interval_product(x_atom, y_atom)

    def _maximal_free_atomic_x_rectangles(self) -> Iterator[Tuple[Interval, Interval]]:
        return _traverse_diagonally(list(self._x_boundaries), self._free_y_ranges)

    def __or__(self, other) -> 'RPolygon':
        rec_copy = copy(self)
        for x_interval, y_interval in other._maximal_used_atomic_x_rectangles():
            rec_copy._add_interval_product(x_interval, y_interval)
        return rec_copy

    def __sub__(self, other) -> 'RPolygon':
        poly_copy = copy(self)
        for x_interval, y_interval in other._maximal_used_atomic_x_rectangles():
            poly_copy._sub_interval_product(x_interval, y_interval)
        return poly_copy

    def __invert__(self) -> 'RPolygon':
        poly_copy = copy(self)
        poly_copy._invert()
        return poly_copy

    def __and__(self, other) -> 'RPolygon':
        return self - (~other)

    def __repr__(self):
        if self.empty:
            return "((),())"

        string = []
        for atomic_poly in self.maximal_used_rectangles():
            x_int = atomic_poly.x_enclosure_interval
            y_int = atomic_poly.y_enclosure_interval
            string.append(f"(x={repr(x_int)}, y={repr(y_int)})")
        return " | ".join(string)

    def _invert(self):
        temp = self._used_y_ranges
        self._used_y_ranges = self._free_y_ranges
        self._free_y_ranges = temp

    def _add_interval_product(self, x_interval: Interval, y_interval: Interval):
        if not x_interval.empty and not y_interval.empty:
            for x_atom in x_interval:
                self._add_atomic(x_atom, y_interval)

    def _sub_interval_product(self, x_interval: Interval, y_interval: Interval):
        if not x_interval.empty and not y_interval.empty:
            for x_atom in x_interval:
                self._sub_atomic(x_atom, y_interval)

    def _add_atomic(self, add_x_atom: Interval, add_y_interval: Interval):
        """
        Update self._used_y_ranges and self._free_y_ranges such that they represent the union of the old rectangular
        polygon and the provided interval product.
        """
        _update_x_boundaries_and_y_interval_triangles(
            self._x_boundaries,
            self._used_y_ranges,
            self._free_y_ranges,
            add_x_atom,
            add_y_interval,
        )

    def _sub_atomic(self, sub_x_atom: Interval, sub_y_interval: Interval):
        """
        Update self._used_y_ranges and self._free_y_ranges such that they represent the set difference of the old
        rectangular polygon and the provided interval product.
        """
        _update_x_boundaries_and_y_interval_triangles(
            self._x_boundaries,
            self._free_y_ranges,
            self._used_y_ranges,
            sub_x_atom,
            sub_y_interval,
        )


def rempty() -> RPolygon:
    """
    Create an empty rectangular polygon.

    :return RPolygon:
    """
    return RPolygon.from_interval_product(empty(), empty())


def ropen(x_lower, x_upper, y_lower, y_upper) -> RPolygon:
    """
    Create an open rectangular polygon.

    :param x_lower: value of the lower bound in x-dimension
    :param x_upper: value of the upper bound in x-dimension
    :param y_lower: value of the lower bound in y-dimension
    :param y_upper: value of the upper bound in y-dimension
    :return RPolygon:
    """
    return RPolygon.from_interval_product(open(x_lower, x_upper), open(y_lower, y_upper))


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


def ropenclosed(x_lower, x_upper, y_lower, y_upper) -> RPolygon:
    """
    Create a rectangular polygon which is open on both lower bounds
    and closed on both upper bounds.

    :param x_lower: value of the lower bound in x-dimension
    :param x_upper: value of the upper bound in x-dimension
    :param y_lower: value of the lower bound in y-dimension
    :param y_upper: value of the upper bound in y-dimension
    :return RPolygon:
    """
    return RPolygon.from_interval_product(openclosed(x_lower, x_upper), openclosed(y_lower, y_upper))


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
