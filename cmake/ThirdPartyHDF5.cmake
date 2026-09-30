########################################################################
#
# ThirdParty configuration for Nektar++
#
# HDF5
#
########################################################################

# Default HDF5 on for non-Windows platforms
IF(WIN32)
    SET(NEKTAR_DEFAULT_HDF5 OFF)
ELSE()
    SET(NEKTAR_DEFAULT_HDF5 ON)
ENDIF()

OPTION(NEKTAR_USE_HDF5
    "Enable HDF5 I/O support." ${NEKTAR_DEFAULT_HDF5})

UNSET(NEKTAR_DEFAULT_HDF5)

IF (NEKTAR_USE_HDF5)
    # Try to find parallel system HDF5 first.
    IF (NEKTAR_USE_MPI)
        SET(HDF5_PREFER_PARALLEL ON)
    ENDIF()

    FIND_PACKAGE(HDF5 QUIET)

    IF (HDF5_IS_PARALLEL)
        ADD_DEFINITIONS(-DNEKTAR_HDF5_PARALLEL)
    ENDIF()

    IF (HDF5_FOUND)
        SET(BUILD_HDF5 OFF)

        IF (NOT HDF5_IS_PARALLEL AND NEKTAR_USE_MPI)
            MESSAGE(WARNING "Using non-parallel HDF5 with MPI enabled: possible performance issues.")
        ELSEIF (NOT NEKTAR_USE_MPI)
            # Parallel built HDF5 will not compile against non-MPI enabled code.
            SET(BUILD_HDF5 ON)
        ENDIF()
    ELSE()
        SET(BUILD_HDF5 ON)
    ENDIF()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_HDF5
        "Build HDF5 from ThirdParty" ${BUILD_HDF5}
        "NEKTAR_USE_HDF5" OFF)

    IF(THIRDPARTY_BUILD_HDF5)
        IF (NEKTAR_USE_MPI)
            SET(HDF5_MPI_CONFIG -DHDF5_ENABLE_PARALLEL=ON)
        ELSE()
            SET(HDF5_MPI_CONFIG -DHDF5_ENABLE_PARALLEL=OFF)
        ENDIF()

        INCLUDE(ExternalProject)

        # In debug mode HDF5 is built either as hdf5_D or hdf5_debug
        IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
            IF (WIN32)
                set(HDF5_LIB_NAME hdf5_D)
            ELSE()
                set(HDF5_LIB_NAME hdf5_debug)
            ENDIF()
        ELSE()
            SET(HDF5_LIB_NAME hdf5)
        ENDIF()

        THIRDPARTY_LIBRARY(HDF5_LIBRARIES SHARED ${HDF5_LIB_NAME}
            DESCRIPTION "HDF5 library")

        EXTERNALPROJECT_ADD(
            hdf5-1.12.3
            PREFIX ${TPSRC}
            URL ${TPURL}/hdf5-1.12.3.tar.bz2
            URL_MD5 5d609bf2a74f980aa42dbe61de452185
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/hdf5-1.12.3
            BINARY_DIR ${TPBUILD}/hdf5-1.12.3 
            TMP_DIR ${TPBUILD}/hdf5-1.12.3-tmp
            INSTALL_DIR ${TPDIST}
            BUILD_BYPRODUCTS ${HDF5_LIBRARIES}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                ${NEKTAR_EXTERNAL_PROJECT_CMAKE_GENERATOR_ARGS}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_C_FLAGS="-w"
                -DCMAKE_CXX_FLAGS="-w"
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
		-DDEFAULT_API_VERSION=v110
                ${HDF5_MPI_CONFIG}
                -DHDF5_BUILD_CPP_LIB=OFF
                -DBUILD_TESTING=OFF
                -DHDF5_BUILD_TOOLS=OFF
                ${TPSRC}/hdf5-1.12.3
            )

        # Add definition to enable dllexport flag.
        IF(WIN32)
            ADD_DEFINITIONS(-DH5_BUILT_AS_DYNAMIC_LIB)

            # Make sure dlls are installed
            INSTALL(CODE "FILE(GLOB hdf5dlls \"${TPDIST}/bin/hdf5*.dll\")
                IF (NOT hdf5dlls STREQUAL \"\")
                    FILE(INSTALL \${hdf5dlls} DESTINATION \${CMAKE_INSTALL_PREFIX}/${NEKTAR_BIN_DIR})
                ENDIF()
            ")
        ENDIF()

        SET(HDF5_INCLUDE_DIRS ${TPDIST}/include CACHE FILEPATH
            "HDF5 include directory" FORCE)

        MESSAGE(STATUS "Build HDF5: ${HDF5_LIBRARIES}")

        SET(HDF5_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE()
        MESSAGE(STATUS "Found HDF5: ${HDF5_LIBRARIES}")
        SET(HDF5_CONFIG_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})
        ADD_CUSTOM_TARGET(hdf5-1.12.3 ALL)

        # Newer HDF5 versions have changed the API
        # We compile deprecated symbols using the old API if newer than 1.10.0
        IF(HDF5_VERSION VERSION_GREATER_EQUAL 1.10.0)
            ADD_DEFINITIONS(-DH5_USE_110_API)
        ENDIF()
    ENDIF()

    ADD_DEPENDENCIES(thirdparty hdf5-1.12.3)

    MARK_AS_ADVANCED(HDF5_LIBRARIES)
    MARK_AS_ADVANCED(HDF5_INCLUDE_DIRS)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
ENDIF()
