############################################################################
##
## BSD 3-Clause License
##
## [License text]
##
############################################################################

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include("openroad")

option(LOAD_CPLEX "Load CPLEX" ON)
if (LOAD_CPLEX)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -m64 -Wall -fPIC -std=c++11")
  set(LINKER_OPTIONS "-Wl,--export-dynamic")
  find_package(Cplex)
  include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})
endif()

option(LOAD_GUROBI "Load Gurobi" OFF)
if (LOAD_GUROBI)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -m64 -Wall -fPIC -std=c++11")

  find_package(GUROBI)
  #set(GUROBI_ROOT_DIR $ENV{GUROBI_HOME})
  #find_library(GUROBI_LIB gurobi PATHS ${GUROBI_ROOT_DIR}/lib)
  #find_path(GUROBI_INCLUDE_DIR gurobi_c++.h PATHS ${GUROBI_ROOT_DIR}/include)

  message(STATUS "Gurobi include directory: ${GUROBI_INCLUDE_DIRS}")
  message(STATUS "Gurobi library: ${GUROBI_LIBRARY}")
  message(STATUS "Gurobi C++ library: ${GUROBI_CXX_LIBRARY}")
endif()

find_package(Threads REQUIRED)
find_package(ortools REQUIRED)

add_library(par_lib
  src/PartitionMgr.cpp
  src/TritonPart.cpp
  src/Hypergraph.cpp
  src/Utilities.cpp
  src/Coarsener.cpp
  src/Multilevel.cpp
  src/Refiner.cpp
  src/Partitioner.cpp
  src/Evaluator.cpp
  src/GreedyRefine.cpp
  src/ILPRefine.cpp
  src/KWayFMRefine.cpp
  src/KWayPMRefine.cpp
  src/PriorityQueue.cpp
)

target_include_directories(par_lib
  PUBLIC
    include
)

if (LOAD_GUROBI)
  target_include_directories(par_lib
    PRIVATE
      ${GUROBI_INCLUDE_DIRS}
  )
endif()

target_link_libraries(par_lib
  PRIVATE
    odb
    OpenSTA
    utl_lib
    dbSta_lib
    ortools::ortools
)

if (LOAD_CPLEX)
  target_compile_definitions(par_lib
    PRIVATE
      LOAD_CPLEX
  )

  target_link_libraries(par_lib
    PRIVATE
      ${CPLEX_LIBRARIES}
      ${CMAKE_THREAD_LIBS_INIT}
      dl
  )
endif()

if (LOAD_GUROBI)
  target_compile_definitions(par_lib
    PRIVATE
      LOAD_GUROBI
  )

  # Link Gurobi C++ library first, followed by the main Gurobi library
  target_link_libraries(par_lib
    PRIVATE
      ${Concert_LIBRARIES}
      ${GUROBI_LIBRARY}
      ${GUROBI_CXX_LIBRARY}
      dl
  )
endif()

swig_lib(NAME par NAMESPACE par I_FILE src/partitionmgr.i SCRIPTS src/partitionmgr.tcl)

target_sources(par PRIVATE src/MakePartitionMgr.cpp)
messages(TARGET par_lib)
target_include_directories(par PUBLIC include)
target_link_libraries(par PRIVATE odb OpenSTA par_lib)

if (LOAD_CPLEX)
  target_compile_definitions(par PRIVATE LOAD_CPLEX)
  target_link_libraries(par PRIVATE ${CPLEX_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} dl)
endif()

if (LOAD_GUROBI)
  target_compile_definitions(par PRIVATE LOAD_GUROBI)
  target_link_libraries(par PRIVATE ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} dl)
endif()

if (Python3_FOUND AND BUILD_PYTHON)
  swig_lib(NAME par_py NAMESPACE par LANGUAGE python I_FILE src/partitionmgr-py.i SWIG_INCLUDES ${PROJECT_SOURCE_DIR}/include/par SCRIPTS ${CMAKE_CURRENT_BINARY_DIR}/par_py.py)
  target_include_directories(par_py PUBLIC include)
  target_link_libraries(par_py PUBLIC par_lib)
endif()

# Uncomment if required
# add_subdirectory(test)
