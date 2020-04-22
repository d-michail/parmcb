
# Parallel and Distributed Minimum Cycle Basis Library

This project implements algorithms to compute exact and approximate Minimum Cycle Bases of weighted graphs. 
The library is written in C++-14 using the [Boost Graph](https://www.boost.org/) libraries for the underlying
graph implementation. Multicore support is provided using the
[TBB (Intel Threading Building Blocks)](https://software.intel.com/en-us/tbb) library.

TODO

# Develop

In order to build the library you need to have CMake installed. A C++-14 compiler is required, such as 
GCC or Clang.

## Using Eclipse

Assume the source file is downloaded in a folder called `parmcb`. Create a 
folder called `parmcb-build` parallel to the `parmcb` folder. This works best 
when using Eclipse. 

Switch to your new folder `parmcb-build` and issue 

```
cmake ../parmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug
```

in order to build using debugging symbols. 

Open up Eclipse and import an existing project into the workspace.

## Using Clang

Use 

```
CC=/usr/bin/clang CXX=/usr/bin/clang++ cmake ../parmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
```

# Custom TBB location 

See the following [guide](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo) 
on how to install TBB in Ubuntu based systems.

If cmake fails to locate TBB, try something like 

```
TBBROOT=/apps/compilers/intel/19.0.1/tbb cmake ../parmcb/
```

## Linker errors

If you are experiencing any linker errors with the boost libraries, then you boost libraries
might have been compiled with an older ABI. In that case a possible workaround is to use 

```
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
```

in the `CMakeLists.txt` file.




Happy coding!
