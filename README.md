
# Parallel Minimum Cycle Basis Library

This project implements algorithms to compute exact and approximate *Minimum Cycle Bases* of weighted graphs. 
The library is written in C++-14 using the [Boost Graph](https://www.boost.org/) libraries for the underlying
graph implementation. It supports parallel and distributed execution using the
[TBB (Intel Threading Building Blocks)](https://software.intel.com/en-us/tbb) library
and [MPI](https://www.mpi-forum.org/).

An older package containing several sequential minimum cycle basis algorithms, using the LEDA library, can 
be found at https://github.com/d-michail/mcb.

## Cite

If you use this package please cite the following paper:

- K. Mehlhorn and D. Michail.
   **Implementing Minimum Cycle Basis Algorithms.**
   ACM Journal of Experimental Algorithmics, 11(2):1-14, 2006.
   <i class="far fa-file-pdf"></i> [pdf](https://d-michail.github.io/assets/papers/implMCBjournal.pdf),
   <i class="fas fa-link"></i> [web](https://portal.acm.org/citation.cfm?id=1187436.1216582)

Citing software is just as important as citing any other important sources in your research.
If you’re not sure whether or not to cite something, [Shouldacite](http://bit.ly/shouldacite) can help
you decide if you should.

## Algorithms

### Undirected Graphs

The following functions are available which implement different algorithmic variants. All of them use a technique
called _support vector approach_ in order to establish linear independence. Their main differences are how they
compute the actual cycles. 

- signed graph 
   * `mcb_sva_signed`
   * `mcb_sva_signed_tbb`  
   * `mcb_sva_signed_mpi`
- cycles collection from a feedback vertex set 
   * `mcb_sva_fvs_trees`
   * `mcb_sva_fvs_trees_tbb`
   * `mcb_sva_fvs_trees_mpi`
   * `mcb_sva_fvs_trees_tbb_mpi`
- isometric cycles collection
   * `mcb_sva_iso_trees`
   * `mcb_sva_iso_trees_tbb`
   * `mcb_sva_iso_trees_mpi`
   * `mcb_sva_iso_trees_tbb_mpi`

### Undirected Graphs (approximation algorithms)

A (2k-1)-approximate algorithm for any integer k >= 1. 

- signed graph
  * `approx_mcb_sva_signed`
- cycles collection from a feedback vertex set 
  * `approx_mcb_sva_fvs_trees`
- isometric cycles collection
  * `approx_mcb_sva_iso_trees`

The signed graph variation `approx_mcb_sva_signed`  has time complexity equal to  
O( m n^{1+1/k} + min(m^3 + m n^2 \log n, n^{3+3/k}) ).


All these different implementations support weighted undirected graphs which (a) do not contain self-loops, 
(b) do not contain multiple edges, and (c) all edge weights are positive. If your graph has self-loops, 
multiple edges or zero weight edges, you need to handle those in a preprocessing step.

Two demos programs are included which read a graph in DIMACS format and compute a minimum cycle basis.

## Examples

See [examples/README](examples/README.md) for a few basic examples to get you started.

## Develop

In order to build the library you need to have CMake installed. A C++-14 compiler is required, such as 
GCC or Clang.

### Using Eclipse

Assume the source file is downloaded in a folder called `parmcb`. Create a folder called `parmcb-build`
parallel to the `parmcb` folder. This works best when using Eclipse. 

Switch to your new folder `parmcb-build` and issue 

```
cmake ../parmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
```

in order to build without debugging symbols. Use Debug if debugging symbols are required.
Open up Eclipse and import an existing project into the workspace.

### Using Clang

Use 

```
CC=/usr/bin/clang CXX=/usr/bin/clang++ cmake ../parmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
```

### Custom TBB location 

See the following [guide](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo) 
on how to install TBB in Ubuntu based systems.

If cmake fails to locate TBB, try something like 

```
TBBROOT=/apps/compilers/intel/19.0.1/tbb cmake ../parmcb/
```

or 

```
TBBROOT=/opt/intel/oneapi/tbb/2021.9.0 cmake ../parmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
```

### Logging

The library has some log statements at various locations that can be helpful. To compile with these statements 
use the parameter `PARMCB_LOGGING` as shown in the following:

```
cmake ../parmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DPARMCB_LOGGING=ON
```

### Linker errors

If you are experiencing any linker errors with the boost libraries, then you boost libraries
might have been compiled with an older ABI. In that case a possible workaround is to use 

```
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
```

in the `CMakeLists.txt` file.

Happy coding!

## License 

The library may be used under the terms of the [Boost Software License, Version 1.0](https://www.boost.org/LICENSE_1_0.txt). 
Please note that the library is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.

Please refer to the license for details.

SPDX-License-Identifier: BSL-1.0

## Author

(C) Copyright 2019-2023, by Dimitrios Michail



