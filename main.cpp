#include "Definitions.hpp"
#include "ControlEngines/SimulationEngine.hpp"
#include "ControlEngines/SimulationEnginePtr.hpp"
#include "Parallelization/Parallelization.hpp"
#include "Parallelization/ThreadTwoSided.hpp"
#include "Parallelization/ThreadSharedMemory.hpp"

// initialization of a random number generator
std::mt19937 mersenne_twister_generator(std::random_device{}());

int main(int argc, char **argv)
{
  LaunchParallelSession(argc, argv);
  {
#if defined(MPI_PARAMETER_SCAN)
    Thread thread(argc, argv);
    SimulationEngine engine(&thread);
    engine.RunSimulation();
#elif defined(MPI_FAST_INTERACTION)
    ThreadTwoSided thread(argc, argv);
    SimulationEngine engine(&thread);
    engine.RunSimulation();
#elif defined(MPI_FAST_INTERACTION_SHARED_MEMORY)
    ThreadSharedMemory thread(argc, argv);
    SimulationEnginePtr engine(&thread);
    engine.RunChimeraSimulation();
//    engine.RunVortexArraySimulation();
#else
    Thread thread(argc, argv);
    SimulationEngine engine(&thread);
    engine.RunSimulation();
#endif
  }
  FinalizeParallelSession();

  return 0;
}