
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

All these different implementations support weighted undirected graphs which (a) do not contain self-loops, 
(b) do not contain multiple edges, and (c) all edge weights are positive. If your graph has self-loops, 
multiple edges or zero weight edges, you need to handle those in a preprocessing step.

Two demos programs are included which read a graph in DIMACS format and compute a minimum cycle basis.

## Example Usage

```
#include <iostream>
#include <fstream>
#include <map>
#include <list>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <parmcb/parmcb.hpp>

using namespace boost;

int main(int argc, char *argv[]) {

    typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double> > Graph;
    typedef graph_traits<Graph>::edge_descriptor Edge;

    //
    //   0 -- 1 -- 2 -- 6 -- 9
    //   |         |         |
    //   5 -- 4 -- 3 -- 7 -- 8
    //

    Graph graph;

    auto v0 = add_vertex(graph);
    auto v1 = add_vertex(graph);
    auto v2 = add_vertex(graph);
    auto v3 = add_vertex(graph);
    auto v4 = add_vertex(graph);
    auto v5 = add_vertex(graph);
    auto v6 = add_vertex(graph);
    auto v7 = add_vertex(graph);
    auto v8 = add_vertex(graph);
    auto v9 = add_vertex(graph);

    auto e01 = add_edge(v0, v1, graph).first;
    auto e12 = add_edge(v1, v2, graph).first;
    auto e23 = add_edge(v2, v3, graph).first;
    auto e34 = add_edge(v3, v4, graph).first;
    auto e45 = add_edge(v4, v5, graph).first;
    auto e05 = add_edge(v0, v5, graph).first;
    auto e26 = add_edge(v2, v6, graph).first;
    auto e37 = add_edge(v3, v7, graph).first;
    auto e69 = add_edge(v6, v9, graph).first;
    auto e78 = add_edge(v7, v8, graph).first;
    auto e89 = add_edge(v8, v9, graph).first;

    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    weight[e01] = 1.0;
    weight[e12] = 1.0;
    weight[e23] = 100.0;
    weight[e34] = 1.0;
    weight[e45] = 1.0;
    weight[e05] = 1.0;
    weight[e26] = 1.0;
    weight[e37] = 2.0;
    weight[e69] = 1.0;
    weight[e78] = 5.0;
    weight[e89] = 1.0;

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_signed(graph, weight, std::back_inserter(cycles));

    std::cout << "MCB weight = " << mcb_weight << std::endl;
    std::cout << "MCB cycles" << std::endl;
    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;

        for (auto eit = cycle.begin(); eit != cycle.end(); eit++) {
            std::cout << *eit << " ";
        }
        std::cout << std::endl;

    }

    return EXIT_SUCCESS;
}

```

## Develop

In order to build the library you need to have CMake installed. A C++-14 compiler is required, such as 
GCC or Clang.

### Using Eclipse

Assume the source file is downloaded in a folder called `parmcb`. Create a folder called `parmcb-build`
parallel to the `parmcb` folder. This works best when using Eclipse. 

Switch to your new folder `parmcb-build` and issue 

```
cmake ../parmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug
```

in order to build using debugging symbols. Use Release if no debugging symbols are required.
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
