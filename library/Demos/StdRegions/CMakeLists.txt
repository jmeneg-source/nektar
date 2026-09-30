ADD_NEKTAR_EXECUTABLE(StdProject
    COMPONENT demos DEPENDS StdRegions SOURCES StdProject.cpp)
ADD_NEKTAR_EXECUTABLE(StdInterp
    COMPONENT demos DEPENDS StdRegions SOURCES StdInterp.cpp)
ADD_NEKTAR_EXECUTABLE(StdInterpDeriv
    COMPONENT demos DEPENDS StdRegions SOURCES StdInterpDeriv.cpp)
ADD_NEKTAR_EXECUTABLE(StdInterpBasis
    COMPONENT demos DEPENDS StdRegions SOURCES StdInterpBasis.cpp)
ADD_NEKTAR_EXECUTABLE(StdEquiToCoeff2D
    COMPONENT demos DEPENDS StdRegions SOURCES StdEquiToCoeff2D.cpp)

# Interpolation tests.
ADD_NEKTAR_TEST(StdInterp)

# Basis interpolation tests.
ADD_NEKTAR_TEST(StdInterpBasis)

# interpolation tests for derivatives.
ADD_NEKTAR_TEST(StdInterpDeriv)

# Projection tests.
ADD_NEKTAR_TEST(StdProject1D)
ADD_NEKTAR_TEST(StdProject1D_Diff)
ADD_NEKTAR_TEST(StdProject2D)
ADD_NEKTAR_TEST(StdProject2D_Diff)
ADD_NEKTAR_TEST(StdProject3D)
ADD_NEKTAR_TEST(StdProject3D_Diff)
