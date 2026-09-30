CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_CARDIAC_EP
    "Build the Cardiac electrophysiology solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_CARDIAC_EP OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_CARDIAC_EP )

    # The CellMLToNektar script converts CellML repository models into
    # Nektar++ code, which can be compiled into the CardiacEPSolver
    # To allow this script to be tested, a test model is automatically
    # converted and compiled into the solver, so it can be executed as part
    # of the test suite. To do this Python3 and the myokit package is required.
    # This is ignored, if these pre-requisites are not satisfied.

    # Assume support is not available.
    SET(NEKTAR_HAS_MYOKIT FALSE)

    # Find the Python3 interpreter
    find_package(
        Python3
        COMPONENTS Interpreter
    )

    if ( Python3_FOUND )
        # Check if 'myokit' is installed
        execute_process(
            COMMAND pip show myokit
            RESULT_VARIABLE EXIT_CODE
            OUTPUT_QUIET
        )

        if (NOT ${EXIT_CODE} EQUAL 0)
            message(STATUS "Found Python3, but without myokit package. "
                "Skipping testing of CellMLToNektar script.")
        else()
            message(STATUS "Found Python3 with myokit package")

            SET(MYOKIT_TEST_MODEL courtemanche_ramirez_nattel_1998.cellml)
            MESSAGE(STATUS "Building test model: ${MYOKIT_TEST_MODEL}")

            # Run the conversion script to convert the test model and put the
            # resulting code into the CardiacEPSolver/CellModels directory
            execute_process(
                COMMAND python3 convert.py -n CellMLTest
                    -o ${CMAKE_CURRENT_SOURCE_DIR}/CellModels
                    ${MYOKIT_TEST_MODEL}
                WORKING_DIRECTORY
                ${CMAKE_CURRENT_SOURCE_DIR}/Utilities/CellMLToNektar
    		)

            # Enable the compilation of the test code and subsequent test.
            SET(NEKTAR_HAS_MYOKIT TRUE)	
        endif()
    else()
        message(STATUS "Python3 not found. "
                 "Skipping testing of CellMLToNektar script.")

    endif ()

    SET(CardiacEPSolverSources
        EquationSystems/Monodomain.cpp
        EquationSystems/Bidomain.cpp
        EquationSystems/BidomainRoth.cpp
        CellModels/CellModel.cpp
        CellModels/FitzhughNagumo.cpp
        CellModels/AlievPanfilov.cpp
        CellModels/CourtemancheRamirezNattel98.cpp
        CellModels/FentonKarma.cpp
        CellModels/LuoRudy91.cpp
        CellModels/Fox02.cpp
        CellModels/PanditGilesDemir03.cpp
        CellModels/Winslow99.cpp
        CellModels/TenTusscher06.cpp
        Filters/FilterCheckpointCellModel.cpp
        Filters/FilterElectrogram.cpp
        Filters/FilterBenchmark.cpp
        Filters/FilterCellHistoryPoints.cpp
        Filters/FilterOffsetPhase.cpp
        Stimuli/Stimulus.cpp
        Stimuli/StimulusCircle.cpp
        Stimuli/StimulusRect.cpp
        Stimuli/StimulusPoint.cpp
        Stimuli/Protocol.cpp
        Stimuli/ProtocolSingle.cpp
        Stimuli/ProtocolS1.cpp
        Stimuli/ProtocolS1S2.cpp
        )

	IF( NEKTAR_HAS_MYOKIT )
        SET(CardiacEPSolverSources ${CardiacEPSolverSources}
            CellModels/CellMLTest.cpp)
	ENDIF()

    IF( NEKTAR_USE_FFTW )
        SET(CardiacEPSolverSources ${CardiacEPSolverSources}
            Filters/FilterHilbertFFTPhase.cpp
        )
    ENDIF( NEKTAR_USE_FFTW )

    ADD_SOLVER_EXECUTABLE(CardiacEPSolver
        SOURCES CardiacEPSolver.cpp
        LIBRARY_SOURCES
        ${CardiacEPSolverSources})

    IF(NEKTAR_BUILD_SOLVER_LIBS)
        SET(CardiacEPSolverHeaders
            CellModels/AlievPanfilov.h
            CellModels/CellModel.h
            CellModels/CourtemancheRamirezNattel98.h
            CellModels/FentonKarma.h
            CellModels/FitzhughNagumo.h
            CellModels/Fox02.h
            CellModels/LuoRudy91.h
            CellModels/PanditGilesDemir03.h
            CellModels/TenTusscher06.h
            CellModels/Winslow99.h
            EquationSystems/Bidomain.h
            EquationSystems/BidomainRoth.h
            EquationSystems/Monodomain.h
            Filters/FilterBenchmark.h
            Filters/FilterCellHistoryPoints.h
            Filters/FilterCheckpointCellModel.h
            Filters/FilterElectrogram.h
            Filters/FilterHilbertFFTPhase.h
            Filters/FilterOffsetPhase.h
            Stimuli/Protocol.h
            Stimuli/ProtocolS1.h
            Stimuli/ProtocolS1S2.h
            Stimuli/ProtocolSingle.h
            Stimuli/StimulusCircle.h
            Stimuli/Stimulus.h
            Stimuli/StimulusPoint.h
            Stimuli/StimulusRect.h)
        IF (NEKTAR_HAS_MYOKIT)
            SET(CardiacEPSolverHeaders ${CardiacEPSolverHeaders}
				CellModels/CellMLTest.h)
        ENDIF()

        # Create separate library. No need to pass source files, as the executable's
        # LIBRARY_SOURCES objects will be reused
        IF( NEKTAR_USE_FFTW )
            ADD_SOLVER_LIBRARY(CardiacEPSolver HEADERS ${CardiacEPSolverHeaders}
                                DEPENDS SolverUtils ${FFTW_LIBRARY})
        ELSE()
            ADD_SOLVER_LIBRARY(CardiacEPSolver HEADERS ${CardiacEPSolverHeaders}
                                DEPENDS SolverUtils)
        ENDIF()
    ENDIF()

    ADD_SUBDIRECTORY(Utilities)
    IF( NEKTAR_USE_FFTW )
        TARGET_LINK_LIBRARIES(CardiacEPSolver LINK_PRIVATE ${FFTW_LIBRARY})
    ENDIF()
ENDIF( NEKTAR_SOLVER_CARDIAC_EP )
