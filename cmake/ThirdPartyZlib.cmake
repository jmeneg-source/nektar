########################################################################
#
# ThirdParty configuration for Nektar++
#
# ZLib
#
########################################################################

# Find a system ZLIB library. If not found enable the THIRDPARTY_BUILD_ZLIB
# option.
# On Windows, we want to force the use of third party zlib
# since this will be used with the boost build if boost is being
# built as a third party lib
IF(WIN32)
    MESSAGE(STATUS "On a WIN32 platform, zlib will be built as a third party library...")
    SET(BUILD_ZLIB ON)
ELSE()
    FIND_PACKAGE(ZLIB QUIET)
    IF (ZLIB_FOUND AND ZLIB_VERSION_STRING VERSION_GREATER 1.2.8)
        SET(BUILD_ZLIB OFF)
    ELSE ()
        SET(BUILD_ZLIB ON)
    ENDIF()
ENDIF()

OPTION(THIRDPARTY_BUILD_ZLIB "Build ZLib library" ${BUILD_ZLIB})

# If we or the user
IF (THIRDPARTY_BUILD_ZLIB)
    INCLUDE(ExternalProject)

    IF (WIN32)
        SET(ZLIB_NAME zlib)
        SET(ZLIB_NAME_DEBUG zlibd)
    ELSE ()
        SET(ZLIB_NAME z)
        SET(ZLIB_NAME_DEBUG z)

        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
            MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build zlib.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)
	IF(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0")
	    SET(ZLIB_PATCH_COMMAND ${PATCH} -p1 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/zlib-1.2.9.patch &&
	                           ${PATCH} -p1 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/zlib-1.2.9-cmake4.0.patch)
        ELSE()
	    SET(ZLIB_PATCH_COMMAND ${PATCH} -p1 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/zlib-1.2.9.patch)
        ENDIF()
    ENDIF ()

    THIRDPARTY_LIBRARY(ZLIB_LIBRARIES SHARED ${ZLIB_NAME} DESCRIPTION "Zlib library")
    THIRDPARTY_LIBRARY(ZLIB_LIBRARIES_DEBUG SHARED ${ZLIB_NAME_DEBUG} DESCRIPTION "Zlib library")

    EXTERNALPROJECT_ADD(
        zlib-1.2.9
        URL ${TPURL}/zlib-1.2.9.tar.gz
        URL_MD5 "e453644539a07783aa525e834491134e"
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/zlib-1.2.9
        BINARY_DIR ${TPBUILD}/zlib-1.2.9
        TMP_DIR ${TPBUILD}/zlib-1.2.9-tmp
        INSTALL_DIR ${TPDIST}
	PATCH_COMMAND ${ZLIB_PATCH_COMMAND}
        BUILD_BYPRODUCTS ${ZLIB_LIBRARIES}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            ${NEKTAR_EXTERNAL_PROJECT_CMAKE_GENERATOR_ARGS}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            -DCMAKE_C_FLAGS:STRING=-fPIC
            -DCMAKE_MACOSX_RPATH=1
            ${TPSRC}/zlib-1.2.9
        )

    IF (APPLE)
        EXTERNALPROJECT_ADD_STEP(zlib-1.2.9 patch-install-path
            COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id  @rpath/libz.1.2.9.dylib ${TPDIST}/lib/libz.1.2.9.dylib
            DEPENDEES install)
    ENDIF ()

    MESSAGE(STATUS "Build Zlib: ")
    MESSAGE(STATUS " -- Optimized: ${ZLIB_LIBRARIES}")
    MESSAGE(STATUS " -- Debug:     ${ZLIB_LIBRARIES_DEBUG}")
    SET(ZLIB_INCLUDE_DIR ${TPDIST}/include CACHE PATH "Zlib include" FORCE)
    SET(ZLIB_CONFIG_INCLUDE_DIR ${TPINC})
ELSE (THIRDPARTY_BUILD_ZLIB)
    ADD_CUSTOM_TARGET(zlib-1.2.9 ALL)
    MESSAGE(STATUS "Found Zlib: ${ZLIB_LIBRARIES} (version ${ZLIB_VERSION_STRING})")

    # We use the found library also for debug builds.
    SET(ZLIB_LIBRARIES_DEBUG ${ZLIB_LIBRARIES})

    SET(ZLIB_CONFIG_INCLUDE_DIR ${ZLIB_INCLUDE_DIRS})
ENDIF (THIRDPARTY_BUILD_ZLIB)

ADD_DEPENDENCIES(thirdparty zlib-1.2.9)

INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
