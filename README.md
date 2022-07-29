# rportion - data structure and operations for rectilinear polygons

[![Tests](https://github.com/tilmann-bartsch/rportion/actions/workflows/test.yaml/badge.svg?branch=master)](https://github.com/tilmann-bartsch/portion/actions/workflows/test.yaml)
[![Coverage Status](https://coveralls.io/repos/github/tilmann-bartsch/rportion/badge.svg?branch=master)](https://coveralls.io/github/tilmann-bartsch/rportion?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Commits](https://badgen.net/github/last-commit/tilmann-bartsch/rportion)](https://github.com/tilmann-bartsch/rportion/commits/)

![](https://github.com/tilmann-bartsch/rportion/blob/master/docu/simple-example.gif)

The `rportion` represent for rectilinear 2D polygons (unions of 2D-intervals) in Python 3.7+.
It is built upon the library [`portion`](https://github.com/AlexandreDecan/portion).

 - 2D-Intervals (rectangles) which can be open/closed and finite/infinite at every boundary
 - support for intersection, union, complement and difference
 - obtain all maximum rectangles inside and outside the given polygon

## Table of contents

  * [Installation](#installation)
  * [Documentation & usage](#documentation--usage)
      * [Polygon creation](#polygon-creation)
      * [Polygon bounds & attributes](#polygon-bounds--attributes)
      * [Polygon operations](#polygon-operations)
      * [Maximum rectangle iterators](#maximum-rectangle-iterators)
  * [Changelog](#changelog)
  * [Contributions](#contributions)
  * [License](#license)

## Installation

TODO

## Documentation & usage

### Polygon creation

Atomic polygons (rectangles) can be created by one of the following:
```python
>>> import rportion as rp
>>> rp.ropen(0, 2, 0, 1)
(x=(0,2), y=(0,1))
>>> rp.rclosed(1, 2)
(x=[0,2], y=[0,1])
>>> rp.ropenclosed(0, 2, 0, 1)
(x=(0,2], y=(0,1])
>>> rp.rclosedopen(0, 2, 0, 1)
(x=[0,2), y=[0,1))
>>> rp.rempty()
((), ())
```

Polygons can also created by using two intervals of the underlying library 
[`portion`](https://github.com/AlexandreDecan/portion):
```python
>>> import portion as P
>>> import rportion as rp
>>> rp.RPolygon.from_interval_product(P.openclosed(0, 2), P.closedopen(0, 1))
(x=(0,2], y=[0,1))
```

### Polygon bounds & attributes

TODO

### Polygon operations

TODO

### Maximum rectangle iterators

TODO

## Changelog
This library adheres to a [semantic versioning](https://semver.org/) scheme.
See [CHANGELOG.md](https://github.com/tilmann-bartsch/rportion/blob/master/CHANGELOG.md) for the list of changes.

## Contributions
Contributions are very welcome! Feel free to report bugs or suggest new features using GitHub issues and/or pull requests.

## License
Distributed under [MIT License](https://github.com/tilmann-bartsch/rportion/blob/master/LICENSE).