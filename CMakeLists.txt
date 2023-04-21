cmake_minimum_required (VERSION 3.10 FATAL_ERROR)

project (parmcb CXX)

set(VERSION "0.1")
set(PACKAGE_VERSION "0.1")
set(DESCRIPTION "Parallel Minimum Cycle Basis Library")

set(LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message(STATUS "Building PARMCB library version: ${VERSION}")

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

# Use classic install dirs
include(GNUInstallDirs)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

add_subdirectory(include/parmcb)
add_subdirectory(include/parmcb/detail)
add_subdirectory(include/parmcb/mpi)
set(GENERATED_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include/parmcb)
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)

option(PARMCB_LOGGING "PARMCB_LOGGING" OFF)

set(PARMCB_HAVE_BOOST ${Boost_FOUND})
set(PARMCB_HAVE_TBB ${TBB_FOUND})
set(PARMCB_HAVE_MPI ${MPI_FOUND})

# If you have ABI problems, try this one
#add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

add_library(parmcb INTERFACE)
target_include_directories(parmcb INTERFACE ${PROJECT_SOURCE_DIR}/include/)
configure_file(${PROJECT_SOURCE_DIR}/include/parmcb/config.hpp.in ${GENERATED_INCLUDE_DIR}/config.hpp)

install(FILES ${GENERATED_INCLUDE_DIR}/config.hpp DESTINATION include/parmcb)

set(
    DEMO_SOURCES
    "mcb-dimacs.cpp"
    "collection-stats-dimacs.cpp"
)
foreach(demosourcefile ${DEMO_SOURCES})
    string(REPLACE ".cpp" "" demoname ${demosourcefile})
    add_executable(${demoname} ${PROJECT_SOURCE_DIR}/src/${demosourcefile})
    target_include_directories(${demoname} PUBLIC ${PROJECT_SOURCE_DIR}/include)
    target_include_directories(${demoname} PUBLIC ${Boost_INCLUDE_DIRS})
    target_link_libraries(${demoname} Boost::timer Boost::program_options Boost::thread)
    target_include_directories(${demoname} PUBLIC ${TBB_INCLUDE_DIRS})
    target_compile_definitions(${demoname} PUBLIC ${TBB_DEFINITIONS})
    target_link_libraries(${demoname} ${TBB_LIBRARIES})
endforeach(demosourcefile ${DEMO_SOURCES})

if (PARMCB_HAVE_TBB AND PARMCB_HAVE_MPI)
    set(
        DEMO_MPI_SOURCES
        "mcb-dimacs-mpi.cpp"
    )
    foreach(demosourcefile ${DEMO_MPI_SOURCES})
        string(REPLACE ".cpp" "" demoname ${demosourcefile})
        add_executable(${demoname} ${PROJECT_SOURCE_DIR}/src/${demosourcefile})
        target_include_directories(${demoname} PUBLIC ${PROJECT_SOURCE_DIR}/include)
        target_include_directories(${demoname} PUBLIC ${Boost_INCLUDE_DIRS})
        target_link_libraries(${demoname} Boost::timer Boost::program_options Boost::thread Boost::mpi)
        target_include_directories(${demoname} PUBLIC ${TBB_INCLUDE_DIRS})
        target_compile_definitions(${demoname} PUBLIC ${TBB_DEFINITIONS})
        target_link_libraries(${demoname} ${TBB_LIBRARIES})
        target_include_directories(${demoname} PUBLIC ${MPI_C_INCLUDE_PATH})
        target_link_libraries(${demoname} MPI::MPI_CXX)
    endforeach(demosourcefile ${DEMO_MPI_SOURCES})
endif()


enable_testing()
include(CTest)
set(
    TEST_SOURCES
    "test_forest_index.cpp"
    "test_fvs.cpp"
    "test_mcb.cpp"
    "test_approx_mcb.cpp"
)
foreach(testsourcefile ${TEST_SOURCES})
    string(REPLACE ".cpp" "" testname ${testsourcefile})
    add_executable(${testname} ${PROJECT_SOURCE_DIR}/test/${testsourcefile})
    target_include_directories(${testname} PUBLIC ${PROJECT_SOURCE_DIR}/include)
    target_include_directories(${testname} PUBLIC ${Boost_INCLUDE_DIRS})
    target_link_libraries(${testname} Boost::timer)
    target_include_directories(${testname} PUBLIC ${TBB_INCLUDE_DIRS})
    target_compile_definitions(${testname} PUBLIC ${TBB_DEFINITIONS})
    target_link_libraries(${testname} ${TBB_LIBRARIES})
    add_test(NAME ${testname} COMMAND ${testname})
endforeach(testsourcefile ${TEST_SOURCES})

if(PARMCB_HAVE_TBB AND PARMCB_HAVE_MPI)
   set(INSTALL_TARGETS "parmcb" "mcb-dimacs" "mcb-dimacs-mpi")
else()
   set(INSTALL_TARGETS "parmcb" "mcb-dimacs")
endif()

install(TARGETS ${INSTALL_TARGETS}
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


