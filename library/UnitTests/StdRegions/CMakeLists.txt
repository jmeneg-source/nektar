SET(StdRegionsUnitTestsSources
    main.cpp
    TestLocCollapsedCoords.cpp
    TestEquiSpaced.cpp
    TestGLLSpaced.cpp
    TestInvBwdTrans.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(StdRegionsUnitTests
    COMPONENT unit-test DEPENDS StdRegions SOURCES ${StdRegionsUnitTestsSources})
ADD_TEST(NAME StdRegionsUnitTests COMMAND StdRegionsUnitTests --detect_memory_leaks=0)
