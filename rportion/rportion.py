from copy import copy
from typing import Iterator

import portion as P
from portion.interval import Interval, open, empty, closedopen, closed
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


def _extend_ranges_mat(mat: list[list[Interval]],
                       x_boundaries: SortedList[RBoundary],
                       rbound: RBoundary):
    n = len(mat)
    index = x_boundaries.bisect_right(rbound)
    prev_y_int_row = mat[index - 1]
    mat.insert(index - 1, list(prev_y_int_row))
    for row in range(0, index):
        prev_y_int = mat[row][n - index]
        mat[row].insert(n - index, prev_y_int)


def _y_interval_triangle_prunable(bound_index: int, y_interval_triangle: list[list[Interval]]) -> bool:
    """
    Test if the given bound can be removed from y_interval_triangle.

    :bound_index int:
    :y_interval_triangle list[list[Interval]]:

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


def _remove_boundary(bound_index: int, y_interval_triangle: list[list[Interval]]):
    """
    Remove row and columns from y_interval_triangle specified by bound_index.

    :bound_index int:
    :y_interval_triangle list[list[Interval]]: Data

    This function is meant to be performed if the boundary is prunable, see function
    `_y_interval_triangle_prunable`.

    I.e. for bound_index == 1 we get

             old y-interval-triangle                         new y-interval-triangle

               |  +inf     5      4     3                       |  +inf     5      3
        -------+----------------------------             -------+----------------------
         -inf  |   ()     ()     ()    ()                  -inf |   ()     ()      ()
            3  |   ()    [0,1)  [0,1)           ==>          3  |   ()    [0,1)
            4  |   ()    [0,1)                               5  |   ()
            5  |   ()
    """
    n = len(y_interval_triangle)
    for i in range(bound_index):
        y_interval_triangle[i].pop(n - bound_index)
    y_interval_triangle.pop(bound_index)


def _add_rectangle_to_interval_tree(
        x_boundaries: SortedList[RBoundary],
        y_interval_triangle: list[list[Interval]],
        other_x_atom: Interval, other_y_interval: Interval):
    """
    Update y_interval_triangle such that it represents the union of the old polygon and the one
    provided in the parameters.

    :x_boundaries SortedList[RBoundary]:
    :y_interval_triangle list[list[Interval]]:
    :new_x_atom Atomic, new_y_interval: Interval:
    """
    assert other_x_atom.atomic

    # Left and right boundary of other_x_atom
    x_b_left = RBoundary(other_x_atom.lower, ~other_x_atom.left)
    x_b_right = RBoundary(other_x_atom.upper, other_x_atom.right)

    # (a) Add new rows/columns to the y-interval-triangles and the x-boundaries
    if x_b_left not in x_boundaries:
        _extend_ranges_mat(y_interval_triangle, x_boundaries, x_b_left)
        x_boundaries.add(x_b_left)
    if x_b_right not in x_boundaries:
        _extend_ranges_mat(y_interval_triangle, x_boundaries, x_b_right)
        x_boundaries.add(x_b_right)

    # (b) Update y_interval_triangle_add
    n = len(y_interval_triangle)
    for i in range(n):
        for j in range(n - i):
            l_bound = x_boundaries[i]
            r_bound = x_boundaries[-j - 1]
            x_interval = Interval.from_atomic(~l_bound.btype, l_bound.val, r_bound.val, r_bound.btype)

            if other_x_atom.contains(x_interval):
                y_interval_triangle[i][j] |= other_y_interval
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
                    y_int_left = y_interval_triangle[i][col]
                else:
                    y_int_left = open(-P.inf, P.inf)

                right_ind = x_boundaries.bisect_left(other_y_int_r_bound)
                if other_y_int_r_bound < r_bound and right_ind < n - j:
                    adj_other_x_interval |= Interval.from_atomic(adj_other_x_interval.left,
                                                                 adj_other_x_interval.lower,
                                                                 r_bound.val, r_bound.btype)
                    y_int_right = y_interval_triangle[right_ind][j]
                else:
                    y_int_right = open(-P.inf, P.inf)

                if adj_other_x_interval.contains(x_interval):
                    y_interval_triangle[i][j] |= y_int_left & other_y_interval & y_int_right

    # (c) Prune the y-interval-triangles if possible
    l_ind = x_boundaries.bisect_left(x_b_left)
    l_prunable = _y_interval_triangle_prunable(l_ind, y_interval_triangle)
    if l_prunable:
        _remove_boundary(l_ind, y_interval_triangle)
        x_boundaries.remove(x_b_left)

    r_ind = x_boundaries.bisect_left(x_b_right)
    r_prunable = _y_interval_triangle_prunable(r_ind, y_interval_triangle)
    if r_prunable:
        _remove_boundary(r_ind, y_interval_triangle)
        x_boundaries.remove(x_b_right)


def _sub_rectangle_to_interval_tree(
        x_boundaries: SortedList[RBoundary],
        y_interval_triangle: list[list[Interval]],
        other_x_atom: Interval, other_y_interval: Interval):
    """
    Update y_interval_triangle such that it represents the difference of the old polygon and the one provided
    in the parameters.

    :x_boundaries SortedList[RBoundary]:
    :y_interval_triangle list[list[Interval]]:
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
        _extend_ranges_mat(y_interval_triangle, x_boundaries, x_b_left)
        x_boundaries.add(x_b_left)
    if x_b_right not in x_boundaries:
        _extend_ranges_mat(y_interval_triangle, x_boundaries, x_b_right)
        x_boundaries.add(x_b_right)

    # (b) Update y_interval_triangle
    n = len(y_interval_triangle)
    for i in range(n):
        for j in range(n - i):
            l_bound = ~x_boundaries[i]
            r_bound = x_boundaries[-j - 1]
            x_interval = Interval.from_atomic(l_bound.btype, l_bound.val, r_bound.val, r_bound.btype)
            if other_x_atom.overlaps(x_interval):
                y_interval_triangle[i][j] -= other_y_interval

    # (c) Prune the y-interval-triangles if possible
    l_ind = x_boundaries.bisect_left(x_b_left)
    l_prunable = _y_interval_triangle_prunable(l_ind, y_interval_triangle)
    if l_prunable:
        _remove_boundary(l_ind, y_interval_triangle)
        x_boundaries.remove(x_b_left)

    r_ind = x_boundaries.bisect_left(x_b_right)
    r_prunable = _y_interval_triangle_prunable(r_ind, y_interval_triangle)
    if r_prunable:
        _remove_boundary(r_ind, y_interval_triangle)
        x_boundaries.remove(x_b_right)


def _traverse_diagonally(boundaries: list[RBoundary],
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


class RectPolygon:
    """
    This class represents an orthogonal Polygon.
    """
    __slots__ = ("_x_boundaries", "_y_intervals")

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
        self._y_intervals: list[list[Interval]] = [[empty()]]

    @property
    def enclosing_intervals(self) -> tuple[Interval, Interval]:
        return self.enclosing_x_interval, self.enclosing_y_interval

    @property
    def enclosing_x_interval(self) -> Interval:
        """
        smallest y_interval enclosing the y-dimension of the polygon.
        """
        n = len(self._y_intervals)

        for i in range(n):
            j = n - 1 - i
            if not self._y_intervals[i][j].empty:
                break
        else:
            i += 1

        l_ind = i


        for j in range(n):
            i = n - 1 - j
            if not self._y_intervals[i][j].empty:
                break
        else:
            j += 1

        r_ind = j

        l_bound = self._x_boundaries[l_ind]
        r_bound = self._x_boundaries[-r_ind-1]
        return Interval.from_atomic(~l_bound.btype, l_bound.val,
                                    r_bound.val, r_bound.btype)

    @property
    def enclosing_y_interval(self) -> Interval:
        """
        smallest y_interval enclosing the y-dimension of the polygon.
        """
        n = len(self._y_intervals)
        interval = empty()
        for i in range(n):
            for j in range(n-i):
                interval |= self._y_intervals[i][j]
        return interval

    @classmethod
    def from_intervals(cls, x_interval: Interval, y_interval: Interval) -> 'RectPolygon':
        """
        Create a (simple) rectangular polygon as the product of two intervals.

        :param x_interval:
        :param y_interval:
        :return RPolygon:
        """
        instance = cls()
        instance._add_interval_product(x_interval, y_interval)
        return instance

    def maximal_rectangles(self) -> Iterator['RectPolygon']:
        """
        Yield all maximal rectangle which are contained in this polygon.

        A rectangle is called maximal if it is not contained in any larger rectangle of this polygon.

        :return: Iterator[RPolygon]:
        """
        for x_atom, y_interval in self._maximal_atomic_x_rectangles():
            for y_atom in y_interval:
                if y_atom.empty:
                    continue
                yield self.__class__.from_intervals(x_atom, y_atom)

    def __copy__(self):
        poly_copy = self.__class__()
        poly_copy._x_boundaries = SortedList(copy(b) for b in self._x_boundaries)
        poly_copy._y_intervals = [
            [y_int for y_int in row]
            for row in self._y_intervals
        ]
        return poly_copy

    def __or__(self, other) -> 'RectPolygon':
        rec_copy = copy(self)
        for x_interval, y_interval in other._maximal_atomic_x_rectangles():
            rec_copy._add_interval_product(x_interval, y_interval)
        return rec_copy

    def __sub__(self, other) -> 'RectPolygon':
        rec_copy = copy(self)
        for rec in other.maximal_rectangles():
            rec_copy._sub_interval_product(
                rec.enclosing_x_interval,
                rec.enclosing_y_interval)
        return rec_copy

    def __invert__(self) -> 'RectPolygon':
        inverted = ropen(-P.inf, P.inf, -P.inf, P.inf)
        for x_atom, y_interval in self._maximal_atomic_x_rectangles():
            if x_atom.empty:
                continue
            inverted._sub_interval_product(x_atom, y_interval)
        return inverted

    def __and__(self, other) -> 'RectPolygon':
        return self - (~other)


    def _maximal_atomic_x_rectangles(self) -> Iterator[tuple[Interval, Interval]]:
        return _traverse_diagonally(list(self._x_boundaries), self._y_intervals)

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
        Update self._y_intervals such that they represent the union of the old rectangular
        polygon and the provided interval product.
        """
        _add_rectangle_to_interval_tree(
            self._x_boundaries,
            self._y_intervals,
            add_x_atom,
            add_y_interval,
        )

    def _sub_atomic(self, sub_x_atom: Interval, sub_y_interval: Interval):
        """
        Update self._y_intervals such that they represent the set difference of the old
        rectangular polygon and the provided interval product.
        """
        _sub_rectangle_to_interval_tree(
            self._x_boundaries,
            self._y_intervals,
            sub_x_atom,
            sub_y_interval,
        )


class RectBisection:
    """
    This class represents the bisection of the plane as it is induced by a rectilinear polygon.
    """

    __slots__ = ("used_polygon", "free_polygon")

    def __init__(self):
        self.used_polygon = RectPolygon()
        self.free_polygon = RectPolygon.from_intervals(
            open(-P.inf, P.inf),
            open(-P.inf, P.inf)
        )

    @classmethod
    def from_interval_product(cls, x_interval: Interval, y_interval: Interval) -> 'RectBisection':
        """
        Create a (simple) rectangular polygon as the product of two intervals.

        :param x_interval:
        :param y_interval:
        :return RPolygon:
        """
        instance = cls()
        instance.used_polygon |= RectPolygon.from_intervals(x_interval, y_interval)
        instance.free_polygon -= RectPolygon.from_intervals(x_interval, y_interval)
        return instance

    @classmethod
    def from_rect_polygon(cls, rect_polygon: RectPolygon) -> 'RectBisection':
        instance = cls()
        for rec in rect_polygon.maximal_rectangles():
            instance |= rec
        return instance

    def maximal_rectangles(self) -> Iterator['RectPolygon']:
        """
        Yield all maximal rectangle which are contained in this polygon.

        A rectangle is called maximal if it is not contained in any larger rectangle of this polygon.

        :return: Iterator[RPolygon]:
        """
        return self.used_polygon.maximal_rectangles()

    def maximal_free_rectangles(self) -> Iterator['RectPolygon']:
        """
        Yield all maximal rectangle which do not intersect the polygon uniquely.

        A rectangle is called maximal if it is not contained in any larger rectangle not intersecting this polygon.

        :return: Iterator['RPolygon']:
        """
        return self.free_polygon.maximal_rectangles()

    def __copy__(self):
        rect_bisec_copy = self.__class__()
        rect_bisec_copy.used_polygon = copy(self.used_polygon)
        rect_bisec_copy.free_polygon = copy(self.free_polygon)
        return rect_bisec_copy

    def __or__(self, other) -> 'RectBisection':
        if isinstance(other, RectBisection) or isinstance(other, RectPolygon):
            poly_copy = copy(self)
            for rectangle in other.maximal_rectangles():
                poly_copy.used_polygon |= rectangle
                poly_copy.free_polygon -= rectangle
            return poly_copy
        else:
            return NotImplemented

    def __sub__(self, other) -> 'RectBisection':
        if isinstance(other, RectBisection) or isinstance(other, RectPolygon):
            poly_copy = copy(self)
            for rectangle in other.maximal_rectangles():
                poly_copy.used_polygon -= rectangle
                poly_copy.free_polygon |= rectangle
            return poly_copy
        else:
            return NotImplemented

    def __invert__(self) -> 'RectBisection':
        poly_copy = copy(self)
        poly_copy._invert()
        return poly_copy

    def __and__(self, other) -> 'RectBisection':
        return self - (~other)

    def _invert(self):
        temp = self.used_polygon
        self.used_polygon = self.free_polygon
        self.free_polygon = temp


def rempty() -> RectPolygon:
    """
    Create an empty rectangular polygon.

    :return RPolygon:
    """
    return RectPolygon.from_intervals(empty(), empty())


def ropen(x_lower, x_upper, y_lower, y_upper) -> RectPolygon:
    """
    Create an open rectangular polygon.

    :param x_lower: value of the lower bound in x-dimension
    :param x_upper: value of the upper bound in x-dimension
    :param y_lower: value of the lower bound in y-dimension
    :param y_upper: value of the upper bound in y-dimension
    :return RPolygon:
    """
    return RectPolygon.from_intervals(open(x_lower, x_upper), open(y_lower, y_upper))


def rclosed(x_lower, x_upper, y_lower, y_upper) -> RectPolygon:
    """
    Create a closed rectangular polygon.

    :param x_lower: value of the lower bound in x-dimension
    :param x_upper: value of the upper bound in x-dimension
    :param y_lower: value of the lower bound in y-dimension
    :param y_upper: value of the upper bound in y-dimension
    :return RPolygon:
    """
    return RectPolygon.from_intervals(closed(x_lower, x_upper), closed(y_lower, y_upper))


def rclosedopen(x_lower, x_upper, y_lower, y_upper) -> RectPolygon:
    """
    Create a rectangular polygon which is closed on both lower bounds
    and open on both upper bounds.

    :param x_lower: value of the lower bound in x-dimension
    :param x_upper: value of the upper bound in x-dimension
    :param y_lower: value of the lower bound in y-dimension
    :param y_upper: value of the upper bound in y-dimension
    :return RPolygon:
    """
    return RectPolygon.from_intervals(closedopen(x_lower, x_upper), closedopen(y_lower, y_upper))
