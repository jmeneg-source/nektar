CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Operators/SwitchLimits.h.in
                  Operators/SwitchLimits.h)

SET(STD_REGIONS_SOURCES
StdExpansion.cpp
StdExpansion0D.cpp
StdExpansion1D.cpp
StdExpansion2D.cpp
StdExpansion3D.cpp
StdHexExp.cpp
StdMatrixKey.cpp
StdNodalPrismExp.cpp
StdNodalTetExp.cpp
StdNodalTriExp.cpp
StdPrismExp.cpp
StdPyrExp.cpp
StdQuadExp.cpp
StdSegExp.cpp
StdPointExp.cpp
StdTetExp.cpp
StdTriExp.cpp
)

SET(STD_REGIONS_HEADERS
StdExpansion0D.h
StdExpansion1D.h
StdExpansion2D.h
StdExpansion3D.h
StdExpansion.h
StdHexExp.h
StdMatrixKey.h
StdNodalPrismExp.h
StdNodalTetExp.h
StdNodalTriExp.h
StdPrismExp.h
StdPyrExp.h
StdQuadExp.h
StdRegionsDeclspec.h
StdSegExp.h
StdPointExp.h
StdTetExp.h
StdTriExp.h
StdRegions.hpp
Operators/BwdTransSumFacStdKernels.hpp
Operators/IProductWRTBaseSumFacStdKernels.hpp
Operators/PhysDerivSumFacStdKernels.hpp
#${CMAKE_BINARY_DIR}/library/StdRegions/Operators/SwitchLimits.h
)

ADD_DEFINITIONS(-DSTD_REGIONS_EXPORTS)


ADD_NEKTAR_LIBRARY(StdRegions
    SOURCES ${STD_REGIONS_SOURCES}
    HEADERS ${STD_REGIONS_HEADERS}
    DEPENDS LibUtilities
    SUMMARY "Nektar++ StdRegions library"
    DESCRIPTION "This library provides construction of the reference expansions for the various 1D, 2D and 3D regions.")

TARGET_INCLUDE_DIRECTORIES(StdRegions PRIVATE ${CMAKE_BINARY_DIR}/library/StdRegions)

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()


INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/StdRegions COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

