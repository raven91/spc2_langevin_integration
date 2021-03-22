//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_SIMULATIONENGINEPTR_HPP
#define SPC2ODEINTEGRATION_SIMULATIONENGINEPTR_HPP

#include "../Definitions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"

#include <string>

class SimulationEnginePtr
{
 public:

  explicit SimulationEnginePtr(ThreadSharedMemory *thread);
  ~SimulationEnginePtr();

  void RunChimeraSimulation();
  void RunVortexArraySimulation();

 private:

  ThreadSharedMemory *thread_;
  PeriodicBoundaryConditions pbc_config_;
  Real *system_state_;
  long system_state_size_;

  void InitializeRandomSystemState();
  void InitializeSystemStateFromPreviousSolution(const std::string &file_name);

};

#endif //SPC2ODEINTEGRATION_SIMULATIONENGINEPTR_HPP
