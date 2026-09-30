########################################################################
#
# ThirdParty configuration for Nektar++
#
# Scotch partitioner
#
########################################################################

IF (NOT WIN32)
    OPTION(NEKTAR_USE_SCOTCH
        "Use Scotch library for performing mesh partitioning." ON)
ENDIF(NOT WIN32)

IF (NEKTAR_USE_SCOTCH)
    IF (NEKTAR_USE_MPI)
        FIND_PACKAGE(Scotch 5 COMPONENTS ptscotch)
    ELSE()
        FIND_PACKAGE(Scotch 5)
    ENDIF()

    IF (SCOTCH_FOUND)
        SET(BUILD_SCOTCH OFF)
    ELSE()
        SET(BUILD_SCOTCH ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_SCOTCH
        "Build Scotch library from ThirdParty" ${BUILD_SCOTCH}
        "NEKTAR_USE_SCOTCH" OFF)

    ADD_DEFINITIONS(-DNEKTAR_USE_SCOTCH)

    IF (THIRDPARTY_BUILD_SCOTCH)
        INCLUDE(ExternalProject)

        UNSET(FLEX CACHE)
        FIND_PROGRAM(FLEX flex)
        IF(NOT FLEX)
            MESSAGE(FATAL_ERROR
                "'flex' lexical parser not found. Cannot build scotch.")
        ENDIF(NOT FLEX)
        MARK_AS_ADVANCED(FLEX)

        # Note that scotch is compiled in the source-tree, so we unpack the
        # source code in the ThirdParty builds directory.
        SET(SCOTCH_SRC ${TPBUILD}/scotch-6.0.4/src)

        IF (APPLE)
            SET(SCOTCH_MAKE Makefile.inc.i686_mac_darwin10)
            SET(SCOTCH_LDFLAGS "")
            SET(SCOTCH_CFLAGS "-w -O3 -Drestrict=__restrict -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DCOMMON_TIMING_OLD -DSCOTCH_RENAME -DCOMMON_PTHREAD_BARRIER")
        ELSE ()
            IF (CMAKE_SIZEOF_VOID_P EQUAL 8)
                SET(SCOTCH_MAKE Makefile.inc.x86-64_pc_linux2)
                SET(SCOTCH_CFLAGS "-w -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict -DIDXSIZE64")
            ELSE ()
                SET(SCOTCH_MAKE Makefile.inc.i686_pc_linux2)
                SET(SCOTCH_CFLAGS "-w -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict")
            ENDIF ()
            SET(SCOTCH_LDFLAGS "-lz -lm -lrt -lpthread")
        ENDIF ()

        # Determine the build target and compiler to use for scotch.
        # We use the normal C compiler by default. If MPI is being used and is
        # not built into the normal compiler, we use the MPI-specified compiler.
        IF (NEKTAR_USE_MPI)
            SET(SCOTCH_BUILD_TARGET "ptscotch")
        ELSE ()
            SET(SCOTCH_BUILD_TARGET "scotch")
        ENDIF ()
        IF (NEKTAR_USE_MPI AND NOT MPI_BUILTIN)
            SET(SCOTCH_C_COMPILER ${MPI_C_COMPILER})
        ELSE ()
            SET(SCOTCH_C_COMPILER ${CMAKE_C_COMPILER})
        ENDIF ()

        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
            MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build scotch.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        FIND_PROGRAM(SCOTCH_MAKE_EXECUTABLE NAMES gmake make mingw32-make REQUIRED)
        MARK_AS_ADVANCED(SCOTCH_MAKE_EXECUTABLE)

        SET(SCOTCH_BUILD_BYPRODUCTS
            ${TPDIST}/include/scotch.h
            ${TPDIST}/include/scotchf.h
            ${TPDIST}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}scotch${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${TPDIST}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}scotcherr${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
        IF (NEKTAR_USE_MPI)
            LIST(APPEND SCOTCH_BUILD_BYPRODUCTS
                ${TPDIST}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ptscotch${CMAKE_STATIC_LIBRARY_SUFFIX}
                ${TPDIST}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ptscotcherr${CMAKE_STATIC_LIBRARY_SUFFIX}
            )
        ENDIF()

        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            scotch-6.0.4
            PREFIX ${TPSRC}
            URL ${TPURL}/scotch_6.0.4.tar.gz
            URL_MD5 "d58b825eb95e1db77efe8c6ff42d329f"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/scotch-6.0.4
            BINARY_DIR ${TPBUILD}/scotch-6.0.4
            TMP_DIR ${TPBUILD}/scotch-6.0.4-tmp
            INSTALL_DIR ${TPDIST}
            BUILD_BYPRODUCTS ${SCOTCH_BUILD_BYPRODUCTS}
            CONFIGURE_COMMAND rm -f ${SCOTCH_SRC}/Makefile.inc
            COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${SCOTCH_SRC}/Make.inc/${SCOTCH_MAKE}
                ${SCOTCH_SRC}/Makefile.inc
            COMMAND ${PATCH} -p0 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/scotch-6_0_4-implicit-function.patch
            BUILD_COMMAND ${SCOTCH_MAKE_EXECUTABLE} -C ${SCOTCH_SRC}
                "CFLAGS=-I${TPDIST}/include ${SCOTCH_CFLAGS}"
                "LDFLAGS=-L${TPDIST}/lib ${SCOTCH_LDFLAGS}"
                "CLIBFLAGS=-fPIC"
                "CCP=${SCOTCH_C_COMPILER}"
                "CCD=${SCOTCH_C_COMPILER}"
                "YACC=bison -pscotchyy -y -b y -Wno-yacc"
                ${SCOTCH_BUILD_TARGET}
            INSTALL_COMMAND ${SCOTCH_MAKE_EXECUTABLE} -C ${SCOTCH_SRC}
                prefix=${TPDIST} install
        )

        EXECUTE_PROCESS(
            COMMAND ${FLEX} --version
            OUTPUT_VARIABLE FLEX_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        # PATCH USED TO SOLVE COMPILATION ERROR (undefined reference to `scotchyywrap') 
        # WHEN THIRD PARTY SCOTCH IS COMPILED WITH FLEX 2.6.3. THE PROBLEM HAS BEEN SOLVED
        # WITH FLEX 2.6.4
        IF (FLEX_VERSION STREQUAL "flex 2.6.3")
            EXTERNALPROJECT_ADD_STEP(scotch-6.0.4 patch-flex
                WORKING_DIRECTORY ${TPBUILD}/scotch-6.0.4
                COMMAND ${PATCH} -p0 -f < ${CMAKE_SOURCE_DIR}/cmake/thirdparty-patches/scotch-6_0_4-flex-2_6_3-yy-compatibility.patch
                DEPENDERS configure
                DEPENDEES patch)
        ENDIF()

        THIRDPARTY_LIBRARY(SCOTCH_LIBRARY STATIC scotch
            DESCRIPTION "Scotch library")
        THIRDPARTY_LIBRARY(SCOTCHERR_LIBRARY STATIC scotcherr
            DESCRIPTION "Scotch error library")
        THIRDPARTY_LIBRARY(PTSCOTCH_LIBRARY STATIC ptscotch;scotch
            DESCRIPTION "PT-Scotch library")
        THIRDPARTY_LIBRARY(PTSCOTCHERR_LIBRARY STATIC ptscotcherr
            DESCRIPTION "PT-Scotch error library")
        SET(SCOTCH_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "Scotch include directory" FORCE)
        SET(SCOTCH_LIBRARY_DIR ${TPDIST}/lib CACHE FILEPATH
            "Scotch library directory" FORCE)
        IF(NEKTAR_USE_MPI)
            MESSAGE(STATUS "Build PT-Scotch: ${PTSCOTCH_LIBRARY}")
        ELSE()
            MESSAGE(STATUS "Build Scotch: ${SCOTCH_LIBRARY}")
        ENDIF()
        SET(SCOTCH_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE (THIRDPARTY_BUILD_SCOTCH)
        ADD_CUSTOM_TARGET(scotch-6.0.4 ALL)
        SET(SCOTCH_CONFIG_INCLUDE_DIR ${SCOTCH_INCLUDE_DIR})
    ENDIF (THIRDPARTY_BUILD_SCOTCH)

    ADD_DEPENDENCIES(thirdparty scotch-6.0.4)

    INCLUDE_DIRECTORIES(SYSTEM ${SCOTCH_INCLUDE_DIR} ${PTSCOTCH_INCLUDE_DIR})

    MARK_AS_ADVANCED(SCOTCH_LIBRARY)
    MARK_AS_ADVANCED(SCOTCHERR_LIBRARY)
    MARK_AS_ADVANCED(SCOTCH_LIBRARY_DIR)
    MARK_AS_ADVANCED(SCOTCH_INCLUDE_DIR)
    MARK_AS_ADVANCED(PTSCOTCH_LIBRARY)
    MARK_AS_ADVANCED(PTSCOTCHERR_LIBRARY)
ENDIF()
