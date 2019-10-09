cmake_minimum_required(VERSION 3.10)
project(p4est-test)

set(CMAKE_CXX_STANDARD 17)

#p4est library
SET(P4EST_DIR "/Users/jleen/code/p4est/src")

#SC library
SET(SC_DIR "/Users/jleen/code/p4est/sc/src")

SET(ANACONDA_DIR "/Users/jleen/anaconda")

SET(HDF5_DIR "/usr/local/hdf5")

SET(OPTIMIZATION_FLAGS "-O0 -g")


# compiler flags and settings
SET(CXX_COMPILE_FLAGS "${OPTIMIZATION_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -I${P4EST_DIR} -I${SC_DIR} -I${HDF5_DIR}/include")

SET(C_COMPILE_FLAGS "${OPTIMIZATION_FLAGS} -I${P4EST_DIR} -I${SC_DIR} -I$(BZ_DIR)")

SET(LINK_FLAGS "-std=c++${CMAKE_CXX_STANDARD}  ${OPTIMIZATION_FLAGS} -L${SC_DIR}/.libs -L${P4EST_DIR}/.libs -L${HDF5_DIR}/lib -lsc -lp4est -lm -lhdf5_hl -lhdf5")

SET(CMAKE_C_COMPILER /usr/local/bin/mpicc)
SET(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)

# set up compiler commands
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CXX_COMPILE_FLAGS}" )
SET( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${C_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${LINK_FLAGS}" )

include_directories(${PROJECT_BINARY_DIR}/include)
file(GLOB INCLUDES "./include/*.hpp")

file(GLOB SOURCES "./src/*.cpp")

message(${CXX_COMPILE_FLAGS})

add_executable(p4est-test.x ${SOURCES} ${INCLUDES})

set_target_properties(p4est-test.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run )
