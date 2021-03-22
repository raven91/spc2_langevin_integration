//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP
#define SPC2ODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ChimeraSystem.hpp"
#include "../Parallelization/Thread.hpp"

#include <vector>

class RungeKutta4Stepper
{
 public:

  explicit RungeKutta4Stepper(Thread *thread);
  ~RungeKutta4Stepper();

  void DoStep(ChimeraSystem &chimera_system, std::vector<Real> &system_state, Real t, Real dt);

 private:

  Thread *thread_;
  std::vector<Real> k_1_;
  std::vector<Real> k_2_;
  std::vector<Real> k_3_;
  std::vector<Real> k_4_;

};

#endif //SPC2ODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP
