CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_PULSEWAVE
    "Build the Pulse-wave solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_PULSEWAVE OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_PULSEWAVE )
    ADD_SOLVER_EXECUTABLE(PulseWaveSolver
        SOURCES PulseWaveSolver.cpp
        LIBRARY_SOURCES
        EquationSystems/PulseWaveSystem.h
        EquationSystems/PulseWaveSystem.cpp
        EquationSystems/PulseWavePropagation.h
        EquationSystems/PulseWavePropagation.cpp
        EquationSystems/PulseWaveBoundary.h
        EquationSystems/PulseWaveBoundary.cpp
        EquationSystems/PulseWavePressureArea.h
        EquationSystems/PulseWavePressureArea.cpp
        EquationSystems/TimeDependentInflow.h
        EquationSystems/TimeDependentInflow.cpp
        EquationSystems/QInflow.h
        EquationSystems/QInflow.cpp
        EquationSystems/TerminalOutflow.h
        EquationSystems/TerminalOutflow.cpp
        EquationSystems/ROutflow.h
        EquationSystems/ROutflow.cpp
        EquationSystems/RCROutflow.h
        EquationSystems/RCROutflow.cpp
        EquationSystems/UndefinedInOutflow.h
        EquationSystems/UndefinedInOutflow.cpp
        EquationSystems/BetaPressureArea.h
        EquationSystems/BetaPressureArea.cpp
        EquationSystems/EmpiricalPressureArea.h
        EquationSystems/EmpiricalPressureArea.cpp
        EquationSystems/PowerPressureArea.h
        EquationSystems/PowerPressureArea.cpp
        EquationSystems/AInflow.h
        EquationSystems/AInflow.cpp
        EquationSystems/UInflow.h
        EquationSystems/UInflow.cpp
        RiemannSolvers/UpwindPulseSolver.cpp
        )

    IF(NEKTAR_BUILD_SOLVER_LIBS)
        SET(PulseWaveSolverHeaders
            EquationSystems/AInflow.h
            EquationSystems/BetaPressureArea.h
            EquationSystems/EmpiricalPressureArea.h
            EquationSystems/PowerPressureArea.h
            EquationSystems/PulseWaveBoundary.h
            EquationSystems/PulseWavePressureArea.h
            EquationSystems/PulseWavePropagation.h
            EquationSystems/PulseWaveSystem.h
            EquationSystems/QInflow.h
            EquationSystems/RCROutflow.h
            EquationSystems/ROutflow.h
            EquationSystems/TemplatePressureArea.h
            EquationSystems/TerminalOutflow.h
            EquationSystems/TimeDependentInflow.h
            EquationSystems/UInflow.h
            EquationSystems/UndefinedInOutflow.h
            RiemannSolvers/UpwindPulseSolver.h
            Utilities/PulseWaveSystemOutput.h)
        # Create separate library. No need to pass source files, as the executable's
        # LIBRARY_SOURCES objects will be reused
        ADD_SOLVER_LIBRARY(PulseWaveSolver HEADERS ${PulseWaveSolverHeaders}
                          DEPENDS SolverUtils)
    ENDIF()

    SUBDIRS(Utilities)

    ADD_NEKTAR_PERFORMANCE_TEST(Perf_TwoBifurcations)

    ADD_NEKTAR_TEST(55_Artery_Network)
    ADD_NEKTAR_TEST(Bifurcation)
    ADD_NEKTAR_TEST(Junction)
    ADD_NEKTAR_TEST(TwoBifurcations)
    ADD_NEKTAR_TEST(Merging)
    ADD_NEKTAR_TEST(VariableAreaTest)
    ADD_NEKTAR_TEST(VariableMatPropTest)
    ADD_NEKTAR_TEST(QInflow)
    ADD_NEKTAR_TEST(RCR_boundary)
    ADD_NEKTAR_TEST(AInflow)
    ADD_NEKTAR_TEST(UInflow)
    ADD_NEKTAR_TEST(Viscoelasticity)
    ADD_NEKTAR_TEST(EmpiricalPressure)
    ADD_NEKTAR_TEST(PowerPressure)


    IF (NEKTAR_USE_MPI)
        ADD_NEKTAR_TEST(55_Artery_Network_par)
        ADD_NEKTAR_TEST(TwoBifurcations_par)
        TARGET_LINK_LIBRARIES(PulseWaveSolver LINK_PRIVATE ${GSMPI_LIBRARY} ${XXT_LIBRARY})
    ENDIF (NEKTAR_USE_MPI)

ENDIF( NEKTAR_SOLVER_PULSEWAVE )
