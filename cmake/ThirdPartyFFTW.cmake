########################################################################
#
# ThirdParty configuration for Nektar++
#
# FFTW
#
########################################################################

OPTION(NEKTAR_USE_FFTW
    "Use FFTW routines for performing the Fast Fourier Transform." OFF)

IF (NEKTAR_USE_FFTW)
    # Use static linkage to avoid issue with MKL
    SET(ORIG_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

    # Define static suffixes based on the platform
    IF(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
    ELSE()
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    ENDIF()

    # Set some common FFTW search paths for the library.
    SET(FFTW_SEARCH_PATHS $ENV{LD_LIBRARY_PATH} $ENV{FFTW_HOME}/lib)
    FIND_LIBRARY(FFTW_LIBRARY NAMES fftw3 fftw3f PATHS ${FFTW_SEARCH_PATHS})

    FIND_PATH(FFTW_INCLUDE_DIR NAMES fftw3.h CACHE FILEPATH 
        "FFTW include directory.")

    SET(CMAKE_FIND_LIBRARY_SUFFIXES ${ORIG_SUFFIXES})

    IF (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
        SET(BUILD_FFTW OFF)
    ELSE()
        SET(BUILD_FFTW ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_FFTW
        "Build FFTW from ThirdParty" ${BUILD_FFTW}
        "NEKTAR_USE_FFTW" OFF)

    IF (THIRDPARTY_BUILD_FFTW)
        INCLUDE(ExternalProject)
        THIRDPARTY_LIBRARY(FFTW_LIBRARY STATIC fftw3 DESCRIPTION "FFTW library")

        IF(APPLE)
            UNSET(PATCH CACHE)
            FIND_PROGRAM(PATCH patch)
            IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                    "'patch' tool for modifying files not found. Cannot build fftw3.")
            ENDIF()
            SET(FFTW_PATCH_COMMAND 
	        cp ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/config.guess . &&
	        cp ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/config.sub . &&
                ${PATCH} -p0 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/fftw-3.2.2.patch)
            IF(NOT CMAKE_OSX_SYSROOT)
                EXECUTE_PROCESS(
                    COMMAND xcrun --sdk macosx --show-sdk-path
                    OUTPUT_VARIABLE MACOSX_SDK_PATH
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                )
            ELSE()
                SET(MACOSX_SDK_PATH ${CMAKE_OSX_SYSROOT})
            ENDIF()

            MESSAGE(STATUS "macOS SDK for FFTW: ${MACOSX_SDK_PATH}")
            SET(FFTW_CONFIG SDKROOT=${CMAKE_OSX_SYSROOT}; 
                "CFLAGS=-O3 -isysroot ${MACOSX_SDK_PATH} -arch arm64";
                "CPPFLAGS=-isysroot ${MACOSX_SDK_PATH}";
                "LDFLAGS=-isysroot ${MACOSX_SDK_PATH} -arch arm64";)
        ELSE()
            SET(FFTW_PATCH_COMMAND 
	        cp ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/config.guess . &&
	        cp ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/config.sub .)
            SET(FFTW_CONFIG "CFLAGS=-w"; "CPPFLAGS=-w")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)
        EXTERNALPROJECT_ADD(
            fftw-3.2.2
            URL ${TPURL}/fftw-3.2.2.tar.gz
            URL_MD5 "b616e5c91218cc778b5aa735fefb61ae"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/fftw-3.2.2
            BINARY_DIR ${TPBUILD}/fftw-3.2.2
            TMP_DIR ${TPBUILD}/fftw-3.2.2-tmp
            INSTALL_DIR ${TPDIST}
            BUILD_BYPRODUCTS ${FFTW_LIBRARY}
	    PATCH_COMMAND ${FFTW_PATCH_COMMAND}
            CONFIGURE_COMMAND
                CC=${CMAKE_C_COMPILER}
                ${TPSRC}/fftw-3.2.2/configure
                ${FFTW_CONFIG}
                --prefix=${TPDIST}
                --libdir=${TPDIST}/lib
                --quiet
                --enable-static
                --enable-shared
                --disable-dependency-tracking
        )

        SET(FFTW_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "FFTW include" FORCE)

        MESSAGE(STATUS "Build FFTW: ${TPDIST}/lib/${FFTW_LIBRARY}")
        SET(FFTW_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE ()
        ADD_CUSTOM_TARGET(fftw-3.2.2 ALL)
        MESSAGE(STATUS "Found FFTW: ${FFTW_LIBRARY}")
        SET(FFTW_CONFIG_INCLUDE_DIR ${FFTW_INCLUDE_DIR})
    ENDIF()

    ADD_DEPENDENCIES(thirdparty fftw-3.2.2)

    # Test if FFTW path is a system path. Only add to include path if not an
    # implicitly defined CXX include path (due to GCC 6.x now providing its own
    # version of some C header files and -isystem reorders include paths).
    GET_FILENAME_COMPONENT(X "${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES}"  ABSOLUTE)
    GET_FILENAME_COMPONENT(Y ${FFTW_INCLUDE_DIR} ABSOLUTE)
    STRING(FIND "${X}" "${Y}" X_FIND)

    IF (X_FIND EQUAL -1)
        INCLUDE_DIRECTORIES(SYSTEM ${FFTW_INCLUDE_DIR})
    ENDIF()

    MARK_AS_ADVANCED(FFTW_LIBRARY)
    MARK_AS_ADVANCED(FFTW_INCLUDE_DIR)
ENDIF( NEKTAR_USE_FFTW )
