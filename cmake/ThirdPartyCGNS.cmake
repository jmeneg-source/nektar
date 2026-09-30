########################################################################
#
# ThirdParty configuration for Nektar++
#
# CNGS
#
########################################################################

OPTION(NEKTAR_USE_CGNS "Use CGNS library for NekMesh" OFF)

SET(BUILD_CGNS OFF)

IF( NEKTAR_USE_CGNS )  
    #First search for system cgns install 
    FIND_LIBRARY(CGNS_LIBRARY NAMES "cgns")
    FIND_PATH(CGNS_INCLUDE_DIR NAMES cgnslib.h CACHE FILEPATH 
            			   "CGNS include directory.")

    IF (CGNS_LIBRARY AND CGNS_INCLUDE_DIR)
        MESSAGE(STATUS "Found CGNS: ${CGNS_LIBRARY}")
        MESSAGE(STATUS "Found CGNS include dir: ${CGNS_INCLUDE_DIR}")
	SET(BUILD_CGNS OFF)
    ELSE()
	MESSAGE(STATUS "System  CGNS not found")
        SET(BUILD_CGNS ON)
    ENDIF()

    MESSAGE(STATUS "Build CGNS ${BUILD_CGNS}")
	
    OPTION(THIRDPARTY_BUILD_CGNS
        "Build CGNS from ThirdParty" ${BUILD_CGNS})
     
    MESSAGE(STATUS "CGNS Thirdparty build: ${THIRDPARTY_BUILD_CGNS}")
    IF(THIRDPARTY_BUILD_CGNS )
        INCLUDE( ExternalProject )

        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
            MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build Tetgen.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

	#disable gcc warning
    IF (NOT WIN32)
	   SET(WARNING_FLAGS "-w")
	ELSE()
	   SET(WARNING_FLAGS " ")
    ENDIF()


	#setup hdf5 options
    IF(NEKTAR_USE_HDF5)
	    IF(HDF5_IS_PARALLEL)
	   	    SET(ENABLE_HDF5 "ON")
	   	    SET(NEED_MPI "ON")
            MESSAGE(STATUS "Build CGNS with HDF5 with MPI")
	    ELSE()		
	   	    SET(ENABLE_HDF5 "ON")
	   	    SET(NEED_MPI "OFF")
            MESSAGE(STATUS "Build CGNS with HDF5")
	    ENDIF()
     ELSE()
	   	SET(ENABLE_HDF5 "OFF")
	   	SET(NEED_MPI "OFF")
        MESSAGE(STATUS "Build CGNS")
    ENDIF()
	THIRDPARTY_LIBRARY(CGNS_LIBRARY SHARED cgns
            DESCRIPTION "CGNS library")
        
	EXTERNALPROJECT_ADD(
	    libcgns-4.4
	    PREFIX ${TPSRC}
	    URL ${TPURL}/CGNS-4.4.0.tar.gz
       	    URL_MD5 "7d82f6834c11ee873232cf131fadfba6"
       	    STAMP_DIR ${TPBUILD}/stamp
       	    DOWNLOAD_DIR ${TPSRC}
       	    SOURCE_DIR ${TPSRC}/cgns-4.4
       	    BINARY_DIR ${TPBUILD}/cgns-4.4
       	    TMP_DIR ${TPBUILD}/cgns-4.4-tmp
       	    INSTALL_DIR ${TPDIST}
            BUILD_BYPRODUCTS ${CGNS_LIBRARY}
            PATCH_COMMAND ${PATCH} -p1 < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/cgns-hdf5-prefer-parallel.patch
       	    CONFIGURE_COMMAND ${CMAKE_COMMAND} 
            ${NEKTAR_EXTERNAL_PROJECT_CMAKE_GENERATOR_ARGS}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST} 
            -DCMAKE_C_FLAGS=${WARNING_FLAGS}
            -DBUILD_SHARED_LIBS:BOOL=ON
	    -DCGNS_ENABLE_HDF5=${ENABLE_HDF5}
	    -DHDF5_NEED_MPI=${NEED_MPI}
	    -DCMAKE_BUILD_TYPE:STRING=Release 
             ${TPSRC}/cgns-4.4
             )

        INCLUDE_DIRECTORIES(SYSTEM NekMesh ${TPDIST}/include)
        ADD_DEPENDENCIES(thirdparty libcgns-4.4)
    ELSE()
        ADD_CUSTOM_TARGET(libcgns ALL)
        INCLUDE_DIRECTORIES(SYSTEM NekMesh ${CNGS_INCLUDE_DIR})
        ADD_DEPENDENCIES(thirdparty libcgns)
	SET(CNGS_CONFIG_INCLUDE_DIR ${CGNS_INCLUDE_DIR})
    ENDIF()

    MARK_AS_ADVANCED(CGNS_LIBRARY)
    MARK_AS_ADVANCED(CGNS_DIR)
    MARK_AS_ADVANCED(CGNS_INCLUDE_DIR)
ENDIF( NEKTAR_USE_CGNS )

