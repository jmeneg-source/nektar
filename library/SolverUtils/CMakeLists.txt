SET(SOLVER_UTILS_SOURCES
  Core/Coupling.cpp
  Core/CouplingFile.cpp
  Core/Misc.cpp
  Core/SessionFunction.cpp
  AdvectionSystem.cpp
  ALEHelper.cpp
  Advection/Advection.cpp
  Advection/Advection3DHomogeneous1D.cpp
  Advection/AdvectionNonConservative.cpp
  Advection/AdvectionWeakDG.cpp
  Diffusion/Diffusion.cpp
  Diffusion/Diffusion3DHomogeneous1D.cpp
  Diffusion/DiffusionLDG.cpp
  Diffusion/DiffusionIP.cpp
  Driver.cpp
  DriverAdaptive.cpp
  DriverArnoldi.cpp
  DriverModifiedArnoldi.cpp
  DriverParallelInTime.cpp
  DriverParareal.cpp
  DriverPFASST.cpp
  DriverStandard.cpp
  DriverSteadyState.cpp	
  EquationSystem.cpp
  Filters/Filter.cpp
  Filters/FilterAeroForces.cpp
  Filters/FilterAverageFields.cpp
  Filters/EvaluatePoints.cpp
  Filters/FilterLagrangianPoints.cpp
  Filters/FilterMaxMinFields.cpp
  Filters/FilterCheckpoint.cpp
  Filters/FilterEnergy1D.cpp
  Filters/FilterEnergy.cpp
  Filters/FilterError.cpp
  Filters/FilterHistoryPoints.cpp
  Filters/FilterIntegral.cpp
  Filters/FilterMean.cpp
  Filters/FilterModalEnergy.cpp
  Filters/FilterMovingAverage.cpp
  Filters/FilterFieldConvert.cpp
  Filters/FilterThresholdMax.cpp
  Filters/FilterThresholdMin.cpp
  Filters/FilterBodyFittedVelocity.cpp
  RiemannSolvers/RiemannSolver.cpp
  RiemannSolvers/UpwindSolver.cpp
  UnsteadySystem.cpp
  MMFSystem.cpp
  Forcing/Forcing.cpp
  Forcing/ForcingAbsorption.cpp
  Forcing/ForcingBody.cpp
  Forcing/ForcingMovingReferenceFrame.cpp
  Forcing/ForcingNoise.cpp
  Forcing/ForcingProgrammatic.cpp
  Forcing/ForcingSyntheticEddy.cpp
)

SET(SOLVER_UTILS_HEADERS
  Core/Coupling.h
  Core/CouplingFile.h
  Core/Misc.h
  Core/SessionFunction.h
  AdvectionSystem.h
  ALEHelper.h
  Advection/Advection.h
  Advection/Advection3DHomogeneous1D.h
  Advection/AdvectionNonConservative.h
  Advection/AdvectionWeakDG.h
  Diffusion/Diffusion.h
  Diffusion/Diffusion3DHomogeneous1D.h
  Diffusion/DiffusionLDG.h
  Diffusion/DiffusionIP.h
  Driver.h
  DriverAdaptive.h
  DriverArnoldi.h
  DriverModifiedArnoldi.h
  DriverParallelInTime.h
  DriverParareal.h
  DriverPFASST.h
  DriverStandard.h
  DriverSteadyState.h
  EquationSystem.h
  Filters/Filter.h
  Filters/FilterAeroForces.h
  Filters/FilterAverageFields.h
  Filters/EvaluatePoints.h
  Filters/FilterLagrangianPoints.h
  Filters/FilterMaxMinFields.h
  Filters/FilterCheckpoint.h
  Filters/FilterEnergy1D.h
  Filters/FilterEnergy.h
  Filters/FilterError.h
  Filters/FilterHistoryPoints.h
  Filters/FilterInterfaces.hpp
  Filters/FilterIntegral.h
  Filters/FilterMean.h
  Filters/FilterModalEnergy.h
  Filters/FilterMovingAverage.h
  Filters/FilterFieldConvert.h
  Filters/FilterThresholdMax.h
  Filters/FilterThresholdMin.h
  Filters/FilterBodyFittedVelocity.h
  RiemannSolvers/RiemannSolver.h
  RiemannSolvers/UpwindSolver.h
  SolverUtils.hpp
  SolverUtilsDeclspec.h
  UnsteadySystem.h
  Forcing/Forcing.h
  Forcing/ForcingAbsorption.h
  Forcing/ForcingBody.h
  Forcing/ForcingMovingReferenceFrame.h
  Forcing/ForcingNoise.h
  Forcing/ForcingProgrammatic.h
  Forcing/ForcingSyntheticEddy.h
)

IF (NEKTAR_USE_ARPACK)
  SET(SOLVER_UTILS_SOURCES ${SOLVER_UTILS_SOURCES}
    DriverArpack.cpp)
  SET(SOLVER_UTILS_HEADERS ${SOLVER_UTILS_HEADERS}
    DriverArpack.h)
ENDIF (NEKTAR_USE_ARPACK)

IF (NEKTAR_USE_CWIPI)
  SET(SOLVER_UTILS_SOURCES ${SOLVER_UTILS_SOURCES}
    Core/CouplingCwipi.cpp)
  SET(SOLVER_UTILS_HEADERS ${SOLVER_UTILS_HEADERS}
    Core/CouplingCwipi.h)
ENDIF (NEKTAR_USE_CWIPI)

ADD_DEFINITIONS(-DSOLVER_UTILS_EXPORTS)

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
    SET(SOLVER_UTILS_SOURCES ${SOLVER_UTILS_SOURCES}
        Filters/FilterPython.cpp)
    SET(SOLVER_UTILS_HEADERS ${SOLVER_UTILS_HEADERS}
        Filters/FilterPython.h)
ENDIF()

ADD_NEKTAR_LIBRARY(SolverUtils
    SOURCES ${SOLVER_UTILS_SOURCES}
    HEADERS ${SOLVER_UTILS_HEADERS}
    DEPENDS FieldUtils
    SUMMARY "Nektar++ SolverUtils library"
    DESCRIPTION "This library provides support classes and routines for constructing complete spectral/hp element solvers.")

IF (NEKTAR_BUILD_PYTHON)
    TARGET_LINK_LIBRARIES(SolverUtils LINK_PRIVATE ${Python3_LIBRARIES})
ENDIF()

IF (NEKTAR_USE_VTK)
    IF (THIRDPARTY_BUILD_VTK)
        TARGET_LINK_LIBRARIES(SolverUtils LINK_PUBLIC ${VTK_LIBRARIES})
        ADD_DEPENDENCIES(SolverUtils vtk-9.3.0)
    ELSE()
        IF (VTK_MAJOR_VERSION LESS 6)
            TARGET_LINK_LIBRARIES(SolverUtils LINK_PUBLIC vtkCommon vtkIO vtkGraphics vtkFiltersCore)
        ELSEIF (VTK_MAJOR_VERSION  GREATER_EQUAL 9)
            TARGET_LINK_LIBRARIES(SolverUtils LINK_PUBLIC VTK::CommonCore VTK::IOLegacy VTK::IOXML VTK::FiltersCore)
            vtk_module_autoinit(
                    TARGETS SolverUtils
                    MODULES VTK::CommonCore VTK::IOLegacy VTK::IOXML VTK::FiltersCore)
        ELSE ()
            TARGET_LINK_LIBRARIES(SolverUtils LINK_PUBLIC vtkCommonCore vtkIOLegacy vtkIOXML vtkFiltersCore)
        ENDIF ()
    ENDIF()
ENDIF (NEKTAR_USE_VTK)
