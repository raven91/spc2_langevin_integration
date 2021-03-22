//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_STOCHASTICRUNGEKUTTASTEPPER_HPP
#define SPC2ODEINTEGRATION_STOCHASTICRUNGEKUTTASTEPPER_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ChimeraSystem.hpp"
#include "../Parallelization/Thread.hpp"

#include <cmath>
#include <algorithm> // std::fill

/**
 * Numerical Solution of Stochastic Differential Equations with Jumps in Finance
 * Strong Order 1.5 Taylor Scheme
 */
class StochasticRungeKuttaStepper
{
 public:

  explicit StochasticRungeKuttaStepper(Thread *thread) :
      thread_(thread),
      buffer_(kS * kN, 0.0)
  {

  }
  ~StochasticRungeKuttaStepper()
  {
    buffer_.clear();
  }

  void DoStep(ChimeraSystem &system, std::vector<Real> &system_state, Real t, Real dt)
  {
    static std::vector<Real> derivative(kS * kN, 0.0);
    std::fill(derivative.begin(), derivative.end(), 0.0);
    static std::vector<std::vector<Real>> additional_derivative(kN, std::vector<Real>(2, 0.0));
    std::fill(additional_derivative.begin(), additional_derivative.end(), std::vector<Real>(2, 0.0));
    system.EvaluateRhs(system_state, derivative, additional_derivative, dt);

    std::normal_distribution<Real> norm_dist(0.0, 1.0);
    Real U_1 = 0.0, U_2 = 0.0, delta_W = 0.0, delta_Z = 0.0;

    const std::vector<int> &loop_indices = thread_->GetLoopIndices();
    for (int i : loop_indices)
    {
      U_1 = norm_dist(mersenne_twister_generator);
      U_2 = norm_dist(mersenne_twister_generator);
      delta_W = U_1 * std::sqrt(dt);
      delta_Z = 0.5 * dt * std::sqrt(dt) * (U_1 + U_2 / std::sqrt(3.0));

      int ii = kS * i;
      system_state[ii] += derivative[ii] * dt - std::sqrt(2.0 * system.D_phi()) * derivative[ii + 1] * delta_Z
          - 0.5 * (derivative[ii + 2] * derivative[ii + 1] + system.D_phi() * derivative[ii]) * dt * dt;
      system_state[ii + 1] += derivative[ii + 1] * dt + std::sqrt(2.0 * system.D_phi()) * derivative[ii] * delta_Z
          + 0.5 * (derivative[ii + 2] * derivative[ii] - system.D_phi() * derivative[ii + 1]) * dt * dt;
      system_state[ii + 2] += derivative[ii + 2] * dt + std::sqrt(2.0 * system.D_phi()) * delta_W
          + additional_derivative[i][0] * delta_Z + 0.5 * additional_derivative[i][1] * dt * dt;
    } // i
  }

 private:

  Thread *thread_;
  std::vector<Real> buffer_;

};

#endif //SPC2ODEINTEGRATION_STOCHASTICRUNGEKUTTASTEPPER_HPP
