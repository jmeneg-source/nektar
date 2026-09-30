########################################################################
#
# ThirdParty configuration for Nektar++
#
# ARPACK
#
########################################################################

OPTION(NEKTAR_USE_ARPACK
    "Use Arpack routines for evaluating eigenvalues and eigenvectors" OFF)

SET(BUILD_ARPACK OFF)

IF (NEKTAR_USE_ARPACK)
    FIND_LIBRARY(ARPACK_LIBRARY NAMES "arpack.1" "arpack")

    IF (ARPACK_LIBRARY)
        MESSAGE(STATUS "Found Arpack: ${ARPACK_LIBRARY}")
        MARK_AS_ADVANCED(ARPACK_LIBRARY)
    ELSE()
    	IF(CMAKE_Fortran_COMPILER)
    	    SET(BUILD_ARPACK ON)
    	ELSE()
            MESSAGE(FATAL_ERROR "Could not find or build Arpack")
        ENDIF()
    ENDIF()

    OPTION(THIRDPARTY_BUILD_ARPACK "Build arpack libraries from ThirdParty."
        ${BUILD_ARPACK})

    IF(THIRDPARTY_BUILD_ARPACK)
    	INCLUDE(ExternalProject)

        THIRDPARTY_LIBRARY(ARPACK_LIBRARY SHARED arpack DESCRIPTION "ARPACK library")

        IF(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0")
            UNSET(PATCH CACHE)
            FIND_PROGRAM(PATCH patch)
            IF(NOT PATCH)
	        MESSAGE(FATAL_ERROR
	            "'patch' tool for modifying files not found. Cannot build lst.")
            ENDIF()
	    SET(ARPACK_PATCH_COMMAND ${PATCH} -p1 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/arpack-ng-1.0-cmake4.0.patch)
        ENDIF()

    	EXTERNALPROJECT_ADD(
            arpack-ng-1.0
            PREFIX ${TPSRC}
            URL ${TPURL}/arpack-ng.tar.gz
            URL_MD5 "26cb30275d24eb79c207ed403e794736"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/arpack-ng-1.0
            BINARY_DIR ${TPBUILD}/arpack-ng-1.0
            TMP_DIR ${TPBUILD}/arpack-ng-1.0-tmp
            INSTALL_DIR ${TPDIST}
            BUILD_BYPRODUCTS ${ARPACK_LIBRARY}
	    PATCH_COMMAND ${ARPACK_PATCH_COMMAND}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                ${NEKTAR_EXTERNAL_PROJECT_CMAKE_GENERATOR_ARGS}
                "-DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}"
                "-DCMAKE_Fortran_FLAGS=-w -O3 -fallow-argument-mismatch"
                "-DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}"
                "-DCMAKE_INSTALL_LIBDIR:PATH=${TPDIST}/lib"
                "-DBUILD_SHARED_LIBS:STRING=ON"
                ${TPSRC}/arpack-ng-1.0
            )

        INCLUDE_DIRECTORIES(${TPDIST}/include)

        MESSAGE(STATUS "Build arpack: ${ARPACK_LIBRARY}")
    ELSE()
    	ADD_CUSTOM_TARGET(arpack-ng-1.0 ALL)
    ENDIF()

    ADD_DEPENDENCIES(thirdparty arpack-ng-1.0)
ENDIF()

