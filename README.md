# rportion - data structure and operations for rectilinear polygons


[![PyPI pyversions](https://img.shields.io/pypi/pyversions/pytest-codeblocks.svg?branch=master)](https://test.pypi.org/project/rportion/)
[![Tests](https://github.com/tilmann-bartsch/rportion/actions/workflows/test.yaml/badge.svg?branch=master)](https://github.com/tilmann-bartsch/portion/actions/workflows/test.yaml)
[![Coverage Status](https://coveralls.io/repos/github/tilmann-bartsch/rportion/badge.svg?branch=master)](https://coveralls.io/github/tilmann-bartsch/rportion?branch=master)
[![LGTM](https://img.shields.io/lgtm/grade/python/github/tilmann-bartsch/rportion.svg?branch=master)](https://lgtm.com/projects/g/tilmann-bartsch/rportion)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Commits](https://img.shields.io/github/last-commit/tilmann-bartsch/rportion/master)](https://github.com/tilmann-bartsch/rportion/commits/master)

The `rportion` library provides data structure to represent rectilinear 2D polygons (unions of 2D-intervals) in Python 3.7+.
It is built upon the library [`portion`](https://github.com/AlexandreDecan/portion) and follows its concepts.

 - 2D-Intervals (rectangles) which can be open/closed and finite/infinite at every boundary
 - support for intersection, union, complement and difference
 - obtain all maximum rectangles inside and outside the given polygon

In the case of integers/floats it can be used keep track of the area resulting 
from the union/difference of rectangles:

<p align="center">
  <img src="https://github.com/tilmann-bartsch/rportion/raw/master/docu/simple-example_solid.gif">
</p>

Internally the library uses an [interval tree](https://en.wikipedia.org/wiki/Interval_tree) to represent a polygon.

## Table of contents

  * [Installation](#installation)
  * [Documentation & usage](#documentation--usage)
      * [Polygon creation](#polygon-creation)
      * [Polygon bounds & attributes](#polygon-bounds--attributes)
      * [Polygon operations](#polygon-operations)
      * [Maximum rectangle iterators](#maximum-rectangle-iterators)
      * [Internal data structure](#internal-data-structure)
  * [Changelog](#changelog)
  * [Contributions](#contributions)
  * [License](#license)

## Installation

Install `rportion` from [PyPi-test](https://test.pypi.org/project/rportion/) with `pip` using 
`pip install -i https://test.pypi.org/simple/ rportion`.

Install `rportion` with the development environment using `pip install -e ".[test]"`.

## Documentation & usage

### Polygon creation

Atomic polygons (rectangles) can be created by one of the following:
```python
>>> import rportion as rp
>>> rp.ropen(0, 2, 0, 1)
(x=(0,2), y=(0,1))
>>> rp.rclosed(0, 2, 0, 1)
(x=[0,2], y=[0,1])
>>> rp.ropenclosed(0, 2, 0, 1)
(x=(0,2], y=(0,1])
>>> rp.rclosedopen(0, 2, 0, 1)
(x=[0,2), y=[0,1))
>>> rp.rsingleton(0, 1)
(x=[0], y=[1])
>>> rp.rempty()
(x=(), y=())
```

Polygons can also be created by using two intervals of the underlying library 
[`portion`](https://github.com/AlexandreDecan/portion):
```python
>>> import portion as P
>>> import rportion as rp
>>> rp.RPolygon.from_interval_product(P.openclosed(0, 2), P.closedopen(0, 1))
(x=(0,2], y=[0,1))
```

[&uparrow; back to top](#table-of-contents)
### Polygon bounds & attributes

An `RPolygon` defines the following properties
 - `empty` is true if the polygon is empty.
   ```python
   >>> rp.rclosed(0, 2, 1, 2).empty
   False
   >>> rp.rempty().empty
   True
   ```
 - `atomic` is true if the polygon can be expressed by a single rectangle.
   ```python
   >>> rp.rempty().atomic
   True
   >>> rp.rclosedopen(0, 2, 1, 2).atomic
   True
   >>> (rp.rclosed(0, 2, 1, 2) | rp.rclosed(0, 2, 1, 3)).atomic
   True
   >>> (rp.rclosed(0, 2, 1, 2) | rp.rclosed(1, 2, 1, 3)).atomic
   False
   ```
 - `enclosure` is the smallest rectangle containing the polygon.
   ```python
   >>> (rp.rclosed(0, 2, 0, 2) | rp.rclosed(1, 3, 0, 1)).enclosure
   (x=[0,3], y=[0,2])
   >>> (rp.rclosed(0, 1, -3, 3) | rp.rclosed(-P.inf, P.inf, -1, 1)).enclosure
   (x=(-inf,+inf), y=[-3,3])
   ```
 - `enclosure_x_interval` is the smallest rectangle containing the polygon's extension in x-dimension.
   ```python
   >>> (rp.rclosed(0, 2, 0, 2) | rp.rclosed(1, 3, 0, 1)).x_enclosure_interval
   x=[0,3]
   >>> (rp.rclosed(0, 1, -3, 3) | rp.rclosed(-P.inf, P.inf, -1, 1)).x_enclosure_interval
   (-inf,+inf)
   ```
 - `enclosure_y_interval` is the smallest interval containing the polygon's extension in y-dimension.
   ```python
   >>> (rp.rclosed(0, 2, 0, 2) | rp.rclosed(1, 3, 0, 1)).y_enclosure_interval
   [0,2]
   >>> (rp.rclosed(0, 1, -3, 3) | rp.rclosed(-P.inf, P.inf, -1, 1)).y_enclosure_interval
   [-3,3]
   ```
 - `x_lower`, `x_upper`, `y_lower` and `y_upper` yield the boundaries of the rectangle enclosing
   the polygon.
   ```python
   >>> p = rp.rclosedopen(0, 2, 1, 3)
   >>> p.x_lower, p.x_upper, p.y_lower, p.y_upper
   (0, 2, 1, 3)
   ```
 - `x_left`, `x_right`, `y_left` and `y_right` yield the type of the boundaries of the rectangle enclosing
    the polygon.
    ```python
    >>> p = rp.rclosedopen(0, 2, 1, 3)
    >>> p.x_left, p.x_right, p.y_left, p.y_right
    (CLOSED, OPEN, CLOSED, OPEN)
    ```

### Polygon operations

`RPolygon` instances support the following operations:
 - `p.intersection(other)` and `p & other` return the intersection of two rectilinear polygons.
   ```python
   >>> rp.rclosed(0, 2, 0, 2) & rp.rclosed(1, 3, 0, 1)
   (x=[1,2], y=[0,1])
   ```
 - `p.union(other)` and `p | other` return the union of two rectilinear polygons.
   ```python
   >>> rp.rclosed(0, 2, 0, 2) | rp.rclosed(1, 3, 0, 1)
   (x=[0,3], y=[0,1]) | (x=[0,2], y=[0,2])
   ```
   Note that the resulting polygon is represented by the union of all maximal rectangles contained in
   in the polygon, see [Maximum rectangle iterators](#maximum-rectangle-iterators).
 - `p.complement()` and `~p` return the complement of the rectilinear polygon.
   ```python
   >>> ~rp.ropen(-P.inf, 0, -P.inf, P.inf)
   ((x=[0,+inf), y=(-inf,+inf))
   ```
 - `p.difference(other)` and `p - other` return the difference of two rectilinear polygons.
   ```python
   rp.rclosed(0, 3, 0, 2) - rp.ropen(2, 4, 1, 3)
   (x=[0,3], y=[0,1]) | (x=[0,2], y=[0,2])
   ```
   Note that the resulting polygon is represented by the union of all maximal rectangles contained in
   in the polygon, see [Maximum rectangle iterators](#maximum-rectangle-iterators).

[&uparrow; back to top](#table-of-contents)
### Maximum rectangle iterators

The method `maximal_used_rectangles` of a `RPolygon` instance returns an iterator over all maximal rectangles contained
in the rectilinear polygon. Similarly, the method `maximal_free_rectangles` of a `RPolygon` instance return an
iterator over all maximal rectangles outside the rectilinear polygon.

A maximal rectangle is rectangle in the polygon which is not a real subset of any other rectangle contained in
the rectilinear polygon.

I.e. for the polygon
```python
>>> poly = rp.rclosedopen(2, 5, 1, 4) | rp.rclosedopen(1, 8, 2, 3) | rp.rclosedopen(6, 8, 1, 3)
>>> poly = poly - rp.rclosedopen(4, 7, 2, 4)
>>> list(poly.maximal_used_rectangles())
[(x=[1,4), y=[2,3)), (x=[2,5), y=[1,2)), (x=[6,8), y=[1,2)), (x=[2,4), y=[1,4)), (x=[7,8), y=[1,3))]
>>> list(poly.maximal_free_rectangles())
[(x=(-inf,+inf), y=(-inf,1)), (x=(-inf,+inf), y=[4,+inf)), ...]
```
which can be visualized as follows:
<p align="center">
  <img width="95%" src="https://github.com/tilmann-bartsch/rportion/raw/master/docu/simple-example.png">
</p>

**Left:** Simple Rectilinear polygon. The red areas are part of the polygon.<br>
**Right:** Maximal contained rectangles are drawn above each other transparently. The maximum
              rectangles lying outside the polygon are represented analogously in green.

[&uparrow; back to top](#table-of-contents)
## Internal data structure

The polygon is internally stored using an [interval tree](https://en.wikipedia.org/wiki/Interval_tree). Every
node of the tree corresponds to an interval in x-dimension which is representable by boundaries (in x-dimension) 
present in the polygon. Each node contains an 1D-interval (by using the library
[`portion`](https://github.com/AlexandreDecan/portion)) in y-dimension. Combining those 1D-intervals
yields a rectangle contained in the polygon.

I.e. for the rectangle `(x=[0, 2), y=[1, 3))` this can be visualized as follows.
```
  interval tree with      x-interval corresponding       y-interval stored in
 a lattice-like shape             to each node                each node
       ┌─x─┐                      ┌─(-∞,+∞)─┐                  ┌─()──┐
       │   │                      │         │                  │     │
     ┌─x─┬─x─┐               ┌─(-∞,2)──┬──[0,+∞)─┐          ┌─()──┬──()─┐
     │   │   │               │         │         │          │     │     │
     x   x   x            (-∞,0]     [0,2)     [2,+∞)      ()   [1,3)   ()
```
The class `RPolygon` holds three data values to store its data.
  - `_x_boundaries`: Sorted list of necessary boundaries in x-dimension with type (`OPEN` or `CLOSED`)
  - `_used_y_ranges`: List of lists in a triangular shape representing the interval tree for the
                      space occupied by the rectilinear polygon.
  - `_free_y_ranges`: List of list in a triangular shape representing the interval tree of
                      for the space not occupied by the rectilinear polygon.

I.e.
```python
>>> poly = rp.rclosedopen(0, 2, 1, 3)
>>> poly._x_boundaries
SortedList([(-inf, OPEN), (0, OPEN), (2, OPEN), (+inf, OPEN)])
>>> poly._used_y_ranges
[[(), (), ()], 
 [(), [1,3)], 
 [()]]
>>> poly._free_y_ranges
[[(-inf,1) | [3,+inf), (-inf,1) | [3,+inf), (-inf,+inf)], 
 [(-inf,1) | [3,+inf), (-inf,1) | [3,+inf)], 
 [(-inf,+inf)]]
```

You can use the function `data_tree_to_string` as noted below to print the internal data structure:

```python
>>> poly = rp.rclosedopen(0, 2, 1, 3)
>>> print(data_tree_to_string(poly._x_boundaries, poly._used_y_ranges, 6))
                |  +inf     2     0
----------------+------------------
     -inf (OPEN)|    ()    ()    ()
      0 (CLOSED)|    () [1,3)
      2 (CLOSED)|    ()
```

```python
>>> poly = rp.rclosedopen(2, 5, 1, 4) | rp.rclosedopen(1, 8, 2, 3) | rp.rclosedopen(6, 8, 1, 3)
>>> poly = poly - rp.rclosedopen(4, 7, 2, 4)
>>> print(data_tree_to_string(poly._x_boundaries, poly._used_y_ranges, 6))
                |  +inf     8     7     6     5     4     2     1
----------------+------------------------------------------------
     -inf (OPEN)|    ()    ()    ()    ()    ()    ()    ()    ()
      1 (CLOSED)|    ()    ()    ()    ()    () [2,3) [2,3)
      2 (CLOSED)|    ()    ()    ()    () [1,2) [1,4)
      4 (CLOSED)|    ()    ()    ()    () [1,2)
      5 (CLOSED)|    ()    ()    ()    ()
      6 (CLOSED)|    () [1,2) [1,2)
      7 (CLOSED)|    () [1,3)
```

```python
def data_tree_to_string(x_boundaries,
                        y_intervals,
                        spacing: int):
    col_space = 10
    n = len(y_intervals)
    msg = " " * (spacing + col_space) + "|"
    for x_b in x_boundaries[-1:0:-1]:
        msg += f"{str(x_b.val):>{spacing}}"
    msg += "\n" + f"-" * (spacing+col_space) + "+"
    for i in range(n):
        msg += f"-" * spacing
    msg += "\n"
    for i, row in enumerate(y_intervals):
        x_b = x_boundaries[i]
        msg += f"{str((~x_b).val) + ' (' + str((~x_b).btype) + ')':>{spacing+ col_space}}|"
        for val in row:
            msg += f"{str(val):>{spacing}}"
        msg += "\n"
    return msg
```

[&uparrow; back to top](#table-of-contents)
## Changelog
This library adheres to a [semantic versioning](https://semver.org/) scheme.
See [CHANGELOG.md](https://github.com/tilmann-bartsch/rportion/blob/master/CHANGELOG.md) for the list of changes.

## Contributions
Contributions are very welcome! Feel free to report bugs or suggest new features using GitHub issues and/or pull requests.

## License
Distributed under [MIT License](https://github.com/tilmann-bartsch/rportion/blob/master/LICENSE).