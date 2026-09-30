SET(SpatialDomainsUnitTestsSources
    main.cpp
    TestMovement.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(SpatialDomainsUnitTests
    COMPONENT unit-test DEPENDS SpatialDomains SOURCES ${SpatialDomainsUnitTestsSources})
ADD_TEST(NAME SpatialDomainsUnitTests COMMAND SpatialDomainsUnitTests --detect_memory_leaks=0)
