SET(COLLECTIONS_SOURCES
  CoalescedGeomData.cpp
  Collection.cpp
  CollectionOptimisation.cpp
  Operator.cpp
  BwdTrans.cpp
  Helmholtz.cpp
  LinearAdvectionDiffusionReaction.cpp
  IProductWRTBase.cpp
  IProductWRTDerivBase.cpp
  IProduct.cpp
  PhysDeriv.cpp
  PhysInterp1DScaled.cpp
  )

SET(COLLECTIONS_HEADERS
  CoalescedGeomData.h
  Collection.h
  CollectionOptimisation.h
  IProduct.h
  MatrixFreeBase.h
  Operator.h      
)

# The current source dir is needed to get the header files.
ADD_DEFINITIONS(-DCOLLECTIONS_EXPORTS -I${CMAKE_CURRENT_SOURCE_DIR})

ADD_NEKTAR_LIBRARY(Collections
    SOURCES ${COLLECTIONS_SOURCES}
    HEADERS ${COLLECTIONS_HEADERS}
    DEPENDS LocalRegions MatrixFreeOps
    SUMMARY "Nektar++ Collections library"
    DESCRIPTION "This library provides efficient multi-expansion implementations of operators on the various supported regions.")


INSTALL(DIRECTORY ./
        DESTINATION ${NEKTAR_INCLUDE_DIR}/Collections
        COMPONENT dev
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
