ADD_NEKTAR_EXECUTABLE(PartitionAnalyse
    COMPONENT demos DEPENDS SpatialDomains SOURCES PartitionAnalyse.cpp)
ADD_NEKTAR_EXECUTABLE(LinMeshGraph
    COMPONENT demos DEPENDS SpatialDomains SOURCES LinMeshGraph.cpp)