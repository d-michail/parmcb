

# Develop using Eclipse 

Assume the source file is downloaded in a folder called `libmcb`. Create a 
folder called `libmcb-build` parallel to the `libmcb` folder. This works best 
when using Eclipse. 

Switch to your new folder `libmcb-build` and issue 

```
cmake ../libmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug
```

in order to build using debugging symbols. 

Open up Eclipse and import an existing project into the workspace.


# Compile using Clang

Use 

``
CC=/usr/bin/clang CXX=/usr/bin/clang++ cmake ../libmcb/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
```

# Requirements

The library requires boost and TBB (Intel Threading Building Blocks).

See the following [guide](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo) 
on how to install TBB in Ubuntu based systems.

Happy coding!
