SET(NekMeshUnitTestsSources
    main.cpp
    TestTetrahedron.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(NekMeshUnitTests
    COMPONENT unit-test DEPENDS libNekMesh SOURCES ${NekMeshUnitTestsSources})
ADD_TEST(NAME NekMeshUnitTests COMMAND NekMeshUnitTests --detect_memory_leaks=0)
