ADD_NEKPY_EXECUTABLE(UnsteadyDiffusion.py UnsteadyDiffusion.py)
ADD_NEKPY_EXECUTABLE(Helmholtz.py Helmholtz.py)

ADD_NEKPY_TEST(UnsteadyDiffusion)

IF (NEKTAR_BUILD_UNIT_TESTS)
    ADD_NEKPY_EXECUTABLE(Filter_UnitTest.py Filter_UnitTest.py)
    ADD_NEKPY_EXECUTABLE(EquationSystem_UnitTest.py EquationSystem_UnitTest.py)
    ADD_NEKPY_EXECUTABLE(UnitTestUtils.py UnitTestUtils.py)
    ADD_NEKPY_TEST(Filter_UnitTest)
    ADD_NEKPY_TEST(EquationSystem_UnitTest)
ENDIF()
