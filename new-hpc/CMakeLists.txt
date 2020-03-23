 #cmake verson
cmake_minimum_required(VERSION 3.2)

#project name
 PROJECT(hpc)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE"Debug")

#head file path 
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(Boost COMPONENTS system iostreams iostreams system program_options REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

find_package(BLAS REQUIRED)
INCLUDE_DIRECTORIES(${BLAS_INCLUDE_DIR})

find_package(LAPACK REQUIRED)
INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIR})

#source directory
AUX_SOURCE_DIRECTORY(. DIR_SRCS)

#add excutable file
ADD_EXECUTABLE(LidSolver ${DIR_SRCS})
#ADD_EXECUTABLE(LidSolver src/LidDrivenCavity.cpp  src/LidDrivenCavitySolver.cpp  src/Poisson.cpp)

#add link library
TARGET_LINK_LIBRARIES(LidSolver ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

