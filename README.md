
# Parallel Minimum Cycle Basis Library

This project implements algorithms to compute exact and approximate *Minimum Cycle Bases* of weighted graphs. 
The library is written in C++-14 using the [Boost Graph](https://www.boost.org/) libraries for the underlying
graph implementation. It supports parallel and distributed execution using the
[TBB (Intel Threading Building Blocks)](https://software.intel.com/en-us/tbb) library
and [MPI](https://www.mpi-forum.org/).

An older package containing several sequential minimum cycle basis algorithms, using the LEDA library, can 
be found at https://github.com/d-michail/mcb.

## Cite

If you use this package please consider citing the following paper:

- K. Mehlhorn and D. Michail.
   **Implementing Minimum Cycle Basis Algorithms.**
   ACM Journal of Experimental Algorithmics, 11(2):1-14, 2006.
   <i class="far fa-file-pdf"></i> [pdf](https://d-michail.github.io/assets/papers/implMCBjournal.pdf),
   <i class="fas fa-link"></i> [web](https://portal.acm.org/citation.cfm?id=1187436.1216582)

Citing software is just as important as citing any other important sources in your research.
If you’re not sure whether or not to cite something, [Shouldacite](http://bit.ly/shouldacite) can help
you decide if you should.

## Algorithms

The following functions are available which implement different algorithmic variants. All of them use a technique
called _support vector approach_ in order to establish 
linear independence. There main differences are how they
compute the actual cycles.

- signed graph 
   * `mcb_sva_signed`, 
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

Two demos programs are included which read a graph in DIMACS format and compute a minimum cycle basis.

## Example Execution

TODO

## Develop

In order to build the library you need to have CMake installed. A C++-14 compiler is required, such as 
GCC or Clang.

### Using Eclipse

Assume the source file is downloaded in a folder called `parmcb`. Create a 
folder called `parmcb-build` parallel to the `parmcb` folder. This works best 
when using Eclipse. 

Switch to your new folder `parmcb-build` and issue 

```
cmake ../parmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug
```

in order to build using debugging symbols. 

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

### Linker errors

If you are experiencing any linker errors with the boost libraries, then you boost libraries
might have been compiled with an older ABI. In that case a possible workaround is to use 

```
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
```

in the `CMakeLists.txt` file.




Happy coding!
