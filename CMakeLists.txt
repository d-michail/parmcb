cmake_minimum_required (VERSION 3.10 FATAL_ERROR)

project (mcb)

set(VERSION "0.1")
set(PACKAGE_VERSION "0.1")
set(DESCRIPTION "Minimum Cycle Basis Library")

set(LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message(STATUS "MCB library version: ${VERSION}")

## Add custom path with cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/etc/cmake/")

## ensure out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "Requires an out-of-source build. Please create a separate build directory and run 'cmake path/to/parmcb [options]' there.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

## additional user-defined library directories
foreach(library ${CMAKE_EXTRA_LIBRARIES})
    link_directories(${library})
endforeach(library)

# Boost Libraries
set(Boost_DEBUG OFF)
set(Boost_USE_MULTITHREADED ON)
find_package (Boost REQUIRED COMPONENTS timer program_options thread mpi)
if(Boost_FOUND)
    message(STATUS "Will link against the following boost libraries: ${Boost_LIBRARIES}")
    message(STATUS "Boost includes: ${Boost_INCLUDE_DIRS}")
endif()

# Intel TBB
find_package(TBB)

# MPI
find_package(MPI COMPONENTS CXX)

#configure_file(include/mcb/config.h.in config.h @ONLY IMMEDIATE)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

#add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

add_subdirectory(include/mcb)

add_library(mcb INTERFACE)
target_include_directories(mcb INTERFACE ${PROJECT_SOURCE_DIR}/include/)

add_executable(dimacs src/dimacs.cpp)
target_include_directories(dimacs PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(dimacs PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(dimacs Boost::timer Boost::program_options Boost::thread)
target_include_directories(dimacs PUBLIC ${TBB_INCLUDE_DIRS})
target_compile_definitions(dimacs PUBLIC ${TBB_DEFINITIONS})
target_link_libraries(dimacs ${TBB_LIBRARIES})

add_executable(dimacs_mpi src/dimacs_mpi.cpp)
target_include_directories(dimacs_mpi PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(dimacs_mpi PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(dimacs_mpi Boost::timer Boost::program_options Boost::thread Boost::mpi)
target_include_directories(dimacs_mpi PUBLIC ${TBB_INCLUDE_DIRS})
target_compile_definitions(dimacs_mpi PUBLIC ${TBB_DEFINITIONS})
target_link_libraries(dimacs_mpi ${TBB_LIBRARIES})
target_include_directories(dimacs_mpi PUBLIC ${MPI_C_INCLUDE_PATH})
target_link_libraries(dimacs_mpi MPI::MPI_CXX)


enable_testing()
include(CTest)
add_executable(test1 ${PROJECT_SOURCE_DIR}/test/test1.cpp)
target_include_directories(test1 PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test1 PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(test1 Boost::timer)
target_include_directories(test1 PUBLIC ${TBB_INCLUDE_DIRS})
target_compile_definitions(test1 PUBLIC ${TBB_DEFINITIONS})
target_link_libraries(test1 ${TBB_LIBRARIES})
add_test(NAME test1 COMMAND test1)

add_executable(test2 ${PROJECT_SOURCE_DIR}/test/test2.cpp)
target_include_directories(test2 PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test2 PUBLIC ${Boost_INCLUDE_DIR})
target_link_libraries(test2 Boost::timer)
target_include_directories(test2 PUBLIC ${TBB_INCLUDE_DIRS})
target_compile_definitions(test2 PUBLIC ${TBB_DEFINITIONS})
target_link_libraries(test2 ${TBB_LIBRARIES})
add_test(NAME test2 COMMAND test2)

add_executable(test_fvs ${PROJECT_SOURCE_DIR}/test/test_fvs.cpp)
target_include_directories(test_fvs PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_fvs PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(test_fvs PUBLIC ${TBB_INCLUDE_DIRS})
target_compile_definitions(test_fvs PUBLIC ${TBB_DEFINITIONS})
target_link_libraries(test_fvs ${TBB_LIBRARIES})
add_test(NAME test_fvs COMMAND test_fvs)

install(TARGETS mcb dimacs
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )

# Packaging support
set(CPACK_PACKAGE_VENDOR "Dimitrios Michail")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Minimum Cycle Basis Library")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/LICENSE)
set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/README.md)
set(CPACK_GENERATOR "DEB;TGZ")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "DM")
include(CPack)

