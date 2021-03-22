//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_RUNGEKUTTA4STEPPERPTR_HPP
#define SPC2ODEINTEGRATION_RUNGEKUTTA4STEPPERPTR_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ChimeraSystemPtr.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"
#include "../DynamicalSystems/VortexArraysSystemPtr.hpp"

#include <vector>

class RungeKutta4StepperPtr
{
 public:

  explicit RungeKutta4StepperPtr(ThreadSharedMemory *thread);
  ~RungeKutta4StepperPtr();

  void DoStep(ChimeraSystemPtr &chimera_system, Real *const system_state, Real t, Real dt);
  void DoStep(VortexArraysSystemPtr &vortex_array_system, Real *const system_state, Real t, Real dt);

 private:

  ThreadSharedMemory *thread_;
  std::vector<Real> k_1_;
  std::vector<Real> k_2_;
  std::vector<Real> k_3_;
  std::vector<Real> k_4_;

};

#endif //SPC2ODEINTEGRATION_RUNGEKUTTA4STEPPERPTR_HPP
