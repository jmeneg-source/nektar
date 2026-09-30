OPTION(NEKTAR_USE_VTK "Use VTK library for utilities." OFF)

IF( NEKTAR_USE_VTK )
    # VTK9 uses modified component names - VTK9 can still be discovered with
    # the old names but these are deprecated and produce a number of warnings.
    MESSAGE(STATUS "Looking for VTK >= 9...")
    FIND_PACKAGE(VTK 9 QUIET COMPONENTS
        FiltersCore IOLegacy IOXML IOImage RenderingCore)
    IF(NOT VTK_FOUND)
        MESSAGE(STATUS "VTK 9+ not found, looking for earlier VTK versions...")
        # If we didn't find VTK9+, search VTK<9 using the older component names
        FIND_PACKAGE(VTK COMPONENTS
            vtkFiltersGeometry vtkIOLegacy vtkIOXML vtkIOImage vtkRenderingCore)
    ENDIF()

    IF (VTK_FOUND)
        MESSAGE(STATUS "Found VTK: ${VTK_USE_FILE}")
        IF (VTK_MAJOR_VERSION EQUAL 6 AND VTK_MINOR_VERSION EQUAL 0 AND VTK_BUILD_VERSION EQUAL 0)
            ADD_DEFINITIONS(-DNEKTAR_HAS_VTK_6_0_0)
        ENDIF()

        IF (VTK_MAJOR_VERSION LESS 9)
            # Not required for VTK9+, see https://vtk.org/doc/nightly/html/
            # md__builds_gitlab-kitware-sciviz-ci_Documentation_Doxygen_ModuleMigration.html
            INCLUDE(${VTK_USE_FILE})
        ENDIF()
        SET(BUILD_VTK OFF)
    ELSE (VTK_FOUND)
        SET(BUILD_VTK ON)
    ENDIF (VTK_FOUND)

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_VTK
        "Build VTK library from ThirdParty" ${BUILD_VTK}
        "NEKTAR_USE_VTK" OFF)

    IF( THIRDPARTY_BUILD_VTK )
        INCLUDE( ExternalProject )
        SET(VTK_LIB_LIST vtkFiltersGeometry-9.3 vtkIOLegacy-9.3 vtkIOXML-9.3 vtkIOImage-9.3 vtkRenderingCore-9.3 vtkIOCore-9.3 vtkCommonCore-9.3 vtkFiltersCore-9.3 vtkCommonDataModel-9.3 vtkCommonExecutionModel-9.3 vtksys-9.3)

        THIRDPARTY_LIBRARY(VTK_LIBRARIES SHARED ${VTK_LIB_LIST} DESCRIPTION "VTK libs")

        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
            MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build VTK.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        # The cmake package has been modified due to a bug in the CMake files
        # which causes it to produce an error when the path includes a '+'.
        # Obviously this is inconvenient for us.
        IF(APPLE)
           SET(VTK_USE_X OFF)
           SET(VTK_USE_COCOA ON)
        ELSE()
           SET(VTK_USE_X ON)
           SET(VTK_USE_COCOA OFF)
        ENDIF()
        EXTERNALPROJECT_ADD(
            vtk-9.3.0
            URL "https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz"
            URL_MD5 "8b4dbb0ec85a6c0cf39803b6f891a8f2"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/vtk-9.3.0
            BINARY_DIR ${TPBUILD}/vtk-9.3.0
            TMP_DIR ${TPBUILD}/vtk-9.3.0-tmp
            INSTALL_DIR ${TPDIST}
            PATCH_COMMAND ${PATCH} -p0 -f < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/vtk-9.3.0.patch
            BUILD_BYPRODUCTS ${VTK_LIBRARIES}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} 
                ${NEKTAR_EXTERNAL_PROJECT_CMAKE_GENERATOR_ARGS}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST} 
                -DBUILD_SHARED_LIBS:BOOL=ON 
                -DCMAKE_BUILD_TYPE:STRING=Release 
                -DCMAKE_C_FLAGS="-w"
                -DCMAKE_CXX_FLAGS="-w"
                -DVTK_USE_X=${VTK_USE_X}
                -DVTK_USE_COCOA=${VTK_USE_COCOA}
                -DVTK_MODULE_USE_EXTERNAL_VTK_zlib=ON
                -DZLIB_INCLUDE_DIR=${ZLIB_INCLUDE_DIR}
                -DZLIB_LIBRARY=${ZLIB_LIBRARIES}
                ${TPSRC}/vtk-9.3.0
        )
        SET(VTK_USE_FILE ${VTK_DIR}/UseVTK.cmake)
        SET(VTK_INCLUDE_DIRS ${TPDIST}/include/vtk-9.3 CACHE FILEPATH
            "VTK include directory" FORCE)
        ADD_DEPENDENCIES(thirdparty vtk-9.3.0)
        ADD_DEPENDENCIES(vtk-9.3.0 zlib-1.2.9)
    ENDIF()

    # Force VTK headers to be treated as system headers.
    INCLUDE_DIRECTORIES(SYSTEM ${VTK_INCLUDE_DIRS})

    ADD_DEFINITIONS(-DNEKTAR_USING_VTK)
ENDIF( NEKTAR_USE_VTK )
