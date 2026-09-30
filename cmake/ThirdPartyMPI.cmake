########################################################################
#
# ThirdParty configuration for Nektar++
#
# MPI
#
########################################################################

OPTION(NEKTAR_USE_MPI "Use MPI for parallelisation." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_GSMPI
    "Build GSMPI if needed" ON
    "NEKTAR_USE_MPI" OFF)

IF( NEKTAR_USE_MPI )
    # First check to see if our compiler has MPI built in to avoid linking
    # libraries etc.
    INCLUDE (CheckIncludeFiles)
    INCLUDE (CheckFunctionExists)
    CHECK_INCLUDE_FILES  (mpi.h    HAVE_MPI_H)
    CHECK_FUNCTION_EXISTS(MPI_Send HAVE_MPI_SEND)

    SET(MPI_BUILTIN OFF CACHE INTERNAL
        "Determines whether MPI is built into the compiler")
    IF (NOT "${HAVE_MPI_H}" OR NOT "${HAVE_MPI_SEND}")
        # Enable version detection to identify if MPICH or OPENMPI
        SET(MPI_DETERMINE_LIBRARY_VERSION YES)
        FIND_PACKAGE(MPI REQUIRED)

        INCLUDE_DIRECTORIES(SYSTEM ${MPI_CXX_INCLUDE_PATH} )
        MESSAGE(STATUS "Found MPI: ${MPI_CXX_LIBRARIES}")
    ELSE()
        SET(MPI_BUILTIN ON)
        MESSAGE(STATUS "Found MPI: built in")
	FIND_PROGRAM(HAVE_APRUN aprun)
	IF (HAVE_APRUN)
	    # Probably on Cray
            SET(MPIEXEC "aprun" CACHE STRING "MPI job launching command")
	    SET(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING
                "MPI job launcher flag to specify number of processes")
	ELSE()
            SET(MPIEXEC "mpirun" CACHE STRING "MPI job launching command")
	    SET(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING
                "MPI job launcher flag to specify number of processes")
	ENDIF()
	MARK_AS_ADVANCED(MPIEXEC)
	MARK_AS_ADVANCED(MPIEXEC_NUMPROC_FLAG)
	UNSET(HAVE_APRUN CACHE)
    ENDIF()

    # Add compiler definition for MPI type. This is needed by the Tester
    # when defining the hostfile to ensure it is formatted correctly.
    IF ("${MPI_CXX_LIBRARY_VERSION_STRING}" MATCHES "MPICH")
	    ADD_DEFINITIONS(-DNEKTAR_MPI_TYPE=1)
    ELSE()
	    ADD_DEFINITIONS(-DNEKTAR_MPI_TYPE=0)
    ENDIF()

    ADD_DEFINITIONS(-DNEKTAR_USE_MPI)

    IF (THIRDPARTY_BUILD_GSMPI)
        INCLUDE(ExternalProject)
        THIRDPARTY_LIBRARY(GSMPI_LIBRARY STATIC gsmpi DESCRIPTION "GSMPI Library")
        THIRDPARTY_LIBRARY(XXT_LIBRARY STATIC xxt DESCRIPTION "XXT Library")
        MARK_AS_ADVANCED(GSMPI_LIBRARY)
        MARK_AS_ADVANCED(XXT_LIBRARY)

    	IF(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0")
            UNSET(PATCH CACHE)
            FIND_PROGRAM(PATCH patch)
            IF(NOT PATCH)
	        MESSAGE(FATAL_ERROR
	            "'patch' tool for modifying files not found. Cannot build gsmpi.")
            ENDIF()
	    SET(GSMPI_PATCH_COMMAND ${PATCH} -p1 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/gsmpi-1.2.1_2-cmake4.0.patch)
        ENDIF()

        EXTERNALPROJECT_ADD(
            gsmpi-1.2.1_2
            URL ${TPURL}/gsmpi-1.2.1_2.tar.bz2
            URL_MD5 3690b4a658324cd5ad17cf8f6115fb5d
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/gsmpi-1.2.1_2
            BINARY_DIR ${TPBUILD}/gsmpi-1.2.1_2
            TMP_DIR ${TPBUILD}/gsmpi-1.2.1_2-tmp
            INSTALL_DIR ${TPDIST}
            BUILD_BYPRODUCTS ${GSMPI_LIBRARY} ${XXT_LIBRARY}
	    PATCH_COMMAND ${GSMPI_PATCH_COMMAND}
            CONFIGURE_COMMAND
                ${CMAKE_COMMAND}
                ${NEKTAR_EXTERNAL_PROJECT_CMAKE_GENERATOR_ARGS}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_BUILD_TYPE:STRING=Debug
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                ${TPSRC}/gsmpi-1.2.1_2
        )
        MESSAGE(STATUS "Build GSMPI: ${GSMPI_LIBRARY}")
        MESSAGE(STATUS "Build XXT: ${XXT_LIBRARY}")

        ADD_DEPENDENCIES(thirdparty gsmpi-1.2.1_2)
    ELSE (THIRDPARTY_BUILD_GSMPI)
        MESSAGE(FATAL_ERROR "Must build GSMPI and XXT")
    ENDIF (THIRDPARTY_BUILD_GSMPI)
ENDIF( NEKTAR_USE_MPI )
