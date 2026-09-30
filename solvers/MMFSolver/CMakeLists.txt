CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_MMF
    "Build the ADR Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_MMF )
    ADD_SOLVER_EXECUTABLE(MMFSolver SOURCES
       ./MMFSolver.cpp
       ./EquationSystems/MMFMaxwell.cpp)

    IF(NEKTAR_BUILD_SOLVER_LIBS)
        MESSAGE(WARNING, "Solver library build not set up for MMFSolver - skipping.")
    ENDIF()

    ADD_NEKTAR_TEST        (TestMaxwellSphere)
ENDIF( NEKTAR_SOLVER_MMF )
