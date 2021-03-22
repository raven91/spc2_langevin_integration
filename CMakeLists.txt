cmake_minimum_required(VERSION 3.13)
project(spc2OdeIntegration)

set(CMAKE_CXX_STANDARD 17)

add_executable(spc2OdeIntegration main.cpp
        DynamicalSystems/ChimeraSystem.cpp DynamicalSystems/ChimeraSystem.hpp
        DynamicalSystems/ChimeraSystemPtr.cpp DynamicalSystems/ChimeraSystemPtr.hpp
        ControlEngines/SimulationEngine.cpp ControlEngines/SimulationEngine.hpp
        ControlEngines/SimulationEnginePtr.cpp ControlEngines/SimulationEnginePtr.hpp
        Observers/BinaryObserver.cpp Observers/BinaryObserver.hpp
        Observers/BinaryObserverPtr.cpp Observers/BinaryObserverPtr.hpp
        Definitions.cpp Definitions.hpp
        BoundaryConfiguration/PeriodicBoundaryConditions.cpp BoundaryConfiguration/PeriodicBoundaryConditions.hpp
        Steppers/RungeKutta4Stepper.cpp Steppers/RungeKutta4Stepper.hpp
        Steppers/RungeKutta4StepperPtr.cpp Steppers/RungeKutta4StepperPtr.hpp
        Steppers/StochasticRungeKuttaStepper.hpp Steppers/StochasticEulerStepper.hpp
        Parallelization/Parallelization.cpp Parallelization/Parallelization.hpp
        Parallelization/Thread.cpp Parallelization/Thread.hpp
        Parallelization/ThreadSharedMemory.cpp Parallelization/ThreadSharedMemory.hpp
        Parallelization/ThreadTwoSided.cpp Parallelization/ThreadTwoSided.hpp
        ParametersConfiguration/ChimeraParameters.cpp ParametersConfiguration/ChimeraParameters.hpp
        DynamicalSystems/VortexArraysSystem.cpp DynamicalSystems/VortexArraysSystem.hpp
        DynamicalSystems/VortexArraysSystemPtr.cpp DynamicalSystems/VortexArraysSystemPtr.hpp
        Steppers/StochasticRungeKuttaStepperPtr.hpp Steppers/StochasticEulerStepperPtr.hpp ForceSummation/LinkedList.cpp ForceSummation/LinkedList.hpp ForceSummation/VerletList.cpp ForceSummation/VerletList.hpp)
#add_executable(spc2OdeIntegration main.cpp DynamicalSystems/ChimeraSystem.cpp DynamicalSystems/ChimeraSystem.hpp ControlEngines/SimulationEngine.cpp ControlEngines/SimulationEngine.hpp Observers/BinaryObserver.cpp Observers/BinaryObserver.hpp Definitions.cpp Definitions.hpp BoundaryConfiguration/PeriodicBoundaryConditions.cpp BoundaryConfiguration/PeriodicBoundaryConditions.hpp Steppers/RungeKutta4Stepper.cpp Steppers/RungeKutta4Stepper.hpp Steppers/StochasticRungeKuttaStepper.hpp Steppers/StochasticEulerStepper.hpp Parallelization/Parallelization.cpp Parallelization/Parallelization.hpp Parallelization/Thread.cpp Parallelization/Thread.hpp ParametersConfiguration/ChimeraParameters.cpp ParametersConfiguration/ChimeraParameters.hpp DynamicalSystems/VortexArraysSystem.cpp DynamicalSystems/VortexArraysSystem.hpp)

include_directories(/usr/local/mpich-install/include/ /usr/local/include/)
link_directories(/usr/local/mpich-install/lib/ /usr/local/lib/)
set(CMAKE_OSX_SYSROOT /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk)