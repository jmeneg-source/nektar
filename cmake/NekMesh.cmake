#
# CMake configuration file for NekMesh installation.
#

SET(NEKTAR_BUILD_DEMOS       OFF CACHE BOOL "Build demos.")
SET(NEKTAR_BUILD_SOLVERS     OFF CACHE BOOL "Build example solvers.")
SET(NEKTAR_BUILD_UNIT_TESTS  OFF CACHE BOOL "Build unit tests.")
SET(NEKTAR_BUILD_TESTS       OFF CACHE BOOL "Build regression tests.")
SET(NEKTAR_BUILD_PYTHON      ON  CACHE BOOL "Build Nektar++ Python bindings")
SET(NEKTAR_USE_CCM           ON  CACHE BOOL "Use CCMIO library for binary Star-CCM+ in NekMesh.")
SET(NEKTAR_USE_HDF5          ON  CACHE BOOL "Build HDF5 libraries.")
SET(NEKTAR_USE_MESHGEN       ON  CACHE BOOL "Build mesh generation utilities.")
SET(NEKTAR_USE_MPI           ON  CACHE BOOL "Enable MPI.")
SET(NEKTAR_USE_THREAD_SAFETY ON  CACHE BOOL "Guarantee thread safety in certain core Nektar++ classes.")
SET(NEKTAR_USE_VTK           ON  CACHE BOOL "Use VTK library for utilities.")