########################################################################
#
# ThirdParty configuration for Nektar++
#
# Python interfaces
#
########################################################################

IF (NEKTAR_BUILD_PYTHON)
    FIND_PACKAGE(Python3 COMPONENTS Interpreter Development REQUIRED)

    INCLUDE_DIRECTORIES(${Python3_INCLUDE_DIRS})
    ADD_DEFINITIONS(-DWITH_PYTHON)

    MESSAGE(STATUS "Searching for Python:")
    MESSAGE(STATUS "-- Found interpreter: ${Python3_EXECUTABLE}")
    MESSAGE(STATUS "-- Found development library: ${Python3_LIBRARIES}")

    # Build directory for binary installations
    SET(NEKPY_BASE_DIR "${CMAKE_BINARY_DIR}/python" CACHE INTERNAL "")
    FILE(MAKE_DIRECTORY ${NEKPY_BASE_DIR})
    FILE(MAKE_DIRECTORY ${NEKPY_BASE_DIR}/NekPy)

    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/python/setup.py.in ${NEKPY_BASE_DIR}/setup.py)
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/python/setup.cfg.in ${NEKPY_BASE_DIR}/setup.cfg)
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/python/pyproject.toml.in ${NEKPY_BASE_DIR}/pyproject.toml)

    ADD_CUSTOM_TARGET(nekpy-install-user
        DEPENDS _MultiRegions
        COMMAND ${Python3_EXECUTABLE} -m pip install --user .
        WORKING_DIRECTORY ${NEKPY_BASE_DIR})

    ADD_CUSTOM_TARGET(nekpy-install-dev
        DEPENDS _MultiRegions
        COMMAND ${Python3_EXECUTABLE} -m pip install --user -e .
        WORKING_DIRECTORY ${NEKPY_BASE_DIR})

    ADD_CUSTOM_TARGET(nekpy-install-system
        DEPENDS _MultiRegions
        COMMAND ${Python3_EXECUTABLE} -m pip install .
        WORKING_DIRECTORY ${NEKPY_BASE_DIR})

    # Require pybind11 v3 or later for unique_ptr transfer
    FIND_PACKAGE(pybind11 3 CONFIG QUIET)

    IF (pybind11_FOUND)
        MESSAGE(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR}")
        ADD_CUSTOM_TARGET(pybind11 ALL)
        INCLUDE_DIRECTORIES(${pybind11_INCLUDE_DIR})
    ELSE()
        MESSAGE(STATUS "Build pybind11: ${TPDIST}/include/pybind11")
        EXTERNALPROJECT_ADD(
            pybind11
            PREFIX ${TPSRC}
            URL ${TPURL}/pybind11-3.0.1.zip
            URL_MD5 53f015a45ffaeeec2ad605ac436526ba
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/pybind11
            BINARY_DIR ${TPBUILD}/pybind11
            TMP_DIR ${TPBUILD}/pybind11-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
            ${NEKTAR_EXTERNAL_PROJECT_CMAKE_GENERATOR_ARGS}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            -DPYBIND11_TEST=OFF
            ${TPSRC}/pybind11)

        # Add third-party include to include path.
        INCLUDE_DIRECTORIES(${TPDIST}/include)

        ADD_DEPENDENCIES(thirdparty pybind11)
    ENDIF()

    FILE(WRITE ${NEKPY_BASE_DIR}/NekPy/__init__.py "# Adjust dlopen flags to avoid OpenMPI issues
try:
    import DLFCN as dl
    import sys
    sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
except ImportError:
    pass")
ENDIF()
