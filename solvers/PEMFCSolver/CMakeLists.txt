CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_PEMFC
    "Build the PEM Fuel Cell Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_PEMFC OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_PEMFC )
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)

    ADD_SOLVER_EXECUTABLE(PEMFCSolver SOURCES
        ./PEMFCSolver.cpp
        ./EquationSystems/PEMFCSpecies.cpp
        ./EquationSystems/PEMFCSteady.cpp
        ./EquationSystems/PEMFCMultilayer.cpp
        ./EquationSystems/PEMFCCathodeMS.cpp)

    ADD_NEKTAR_TEST(PEMFC_CathodeMS_Coarse)
ENDIF( NEKTAR_SOLVER_PEMFC )
