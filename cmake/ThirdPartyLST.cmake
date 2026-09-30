########################################################################
#
# ThirdParty confiuration for Nektar++
#
# Linear stability analysis (LST) module in FieldConvert
#
########################################################################

OPTION(NEKTAR_USE_LST "Use LST module in FieldConvert" OFF)

IF(NEKTAR_USE_LST)

    INCLUDE(ExternalProject)
    
    IF(NOT CMAKE_Fortran_COMPILER)
        MESSAGE(FATAL_ERROR
                "Could not find a Fortran compiler to build LST library")
    ENDIF()

    THIRDPARTY_LIBRARY(LST_LIBRARY
        SHARED lst DESCRIPTION "Linear stability analysis library")

    IF(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0")
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
	    MESSAGE(FATAL_ERROR
	        "'patch' tool for modifying files not found. Cannot build lst.")
        ENDIF()
	SET(LST_PATCH_COMMAND ${PATCH} -p1 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/lst-1.6-cmake4.0.patch)
    ENDIF()

    EXTERNALPROJECT_ADD(
        lst-1.6
        URL                ${TPURL}/lst_v1.6.zip
        URL_MD5            5820c5c37016f036b0ecffe289b2f761
        PREFIX             ${TPSRC}
        STAMP_DIR          ${TPBUILD}/stamp
        DOWNLOAD_DIR       ${TPSRC}
        SOURCE_DIR         ${TPSRC}/lst-1.6
        BINARY_DIR         ${TPBUILD}/lst-1.6
        TMP_DIR            ${TPBUILD}/lst-1.6-tmp
        INSTALL_DIR        ${TPDIST}
        BUILD_BYPRODUCTS   ${LST_LIBRARY}
	PATCH_COMMAND      ${LST_PATCH_COMMAND}
        CONFIGURE_COMMAND  ${CMAKE_COMMAND} 
        ${NEKTAR_EXTERNAL_PROJECT_CMAKE_GENERATOR_ARGS}
        "-DCMAKE_MACOSX_RPATH=1"
        "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
        "-DCMAKE_Fortran_FLAGS=-w -O3"
        "-DCMAKE_INSTALL_LIBDIR:PATH=${TPDIST}/lib"
        "-DCMAKE_INSTALL_INCDIR:PATH=${TPDIST}/include"
        ${TPSRC}/lst-1.6)

    INCLUDE_DIRECTORIES(${TPDIST}/include)
    ADD_DEPENDENCIES(thirdparty lst-1.6)

    MARK_AS_ADVANCED(LST_LIBRARY)
    MESSAGE(STATUS "Build LST: ${LST_LIBRARY}")

ENDIF(NEKTAR_USE_LST)



