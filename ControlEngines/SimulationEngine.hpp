//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_SIMULATIONENGINE_HPP
#define SPC2ODEINTEGRATION_SIMULATIONENGINE_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/Thread.hpp"

#include <string>

class SimulationEngine
{
 public:

  explicit SimulationEngine(Thread *thread);
  ~SimulationEngine();

  void RunSimulation();

 private:

  Thread *thread_;
  PeriodicBoundaryConditions pbc_config_;
  std::vector<Real> system_state_;

  void InitializeRandomSystemState();
  void InitializeSystemStateFromFile(const std::string &file_name);

};

#endif //SPC2ODEINTEGRATION_SIMULATIONENGINE_HPP
