//
// Created by Nikita Kruk on 2019-03-13.
//

#ifndef SPC2ODEINTEGRATION_STOCHASTICEULERSTEPPERPTR_HPP
#define SPC2ODEINTEGRATION_STOCHASTICEULERSTEPPERPTR_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ChimeraSystemPtr.hpp"
#include "../DynamicalSystems/VortexArraysSystemPtr.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

#include <cmath>
#include <algorithm> // std::fill
#include <fstream>

/**
 * Numerical Solution of Stochastic Differential Equations with Jumps in Finance
 * Strong Order 0.5 Taylor Scheme
 */
class StochasticEulerStepperPtr
{
 public:

  explicit StochasticEulerStepperPtr(ThreadSharedMemory *thread) :
      thread_(thread)
  {

  }
  ~StochasticEulerStepperPtr() = default;

  void DoStep(ChimeraSystemPtr &system, Real *const system_state, Real t, Real dt)
  {
    static std::vector<Real> derivative(kS * kN, 0.0);
    std::fill(derivative.begin(), derivative.end(), 0.0);
    static std::vector<std::vector<Real>> additional_derivative(kN, std::vector<Real>(2, 0.0));
    std::fill(additional_derivative.begin(), additional_derivative.end(), std::vector<Real>(2, 0.0));
    system.EvaluateRhs(system_state, derivative, additional_derivative, dt);

    std::normal_distribution<Real> norm_dist(0.0, 1.0);
    Real delta_W = 0.0;
//    static std::vector<double> angular_velocity(kN, 0.0);

    for (int i : thread_->GetLoopIndices())
    {
      delta_W = std::sqrt(dt) * norm_dist(mersenne_twister_generator);

      int ii = kS * i;
      system_state[ii] += derivative[ii] * dt;
      system_state[ii + 1] += derivative[ii + 1] * dt;
      system_state[ii + 2] += derivative[ii + 2] * dt + std::sqrt(2.0 * system.D_phi()) * delta_W;
//      angular_velocity[i] = derivative[ii + 2] * dt + std::sqrt(2.0 * system.D_phi()) * delta_W;
    } // i
  }

  void DoStep(VortexArraysSystemPtr &system, Real *const system_state, Real t, Real dt)
  {
    static std::vector<Real> derivative(kS * kN, 0.0);
    std::fill(derivative.begin(), derivative.end(), 0.0);
    static std::vector<std::vector<Real>> additional_derivative(kN, std::vector<Real>(2, 0.0));
    std::fill(additional_derivative.begin(), additional_derivative.end(), std::vector<Real>(2, 0.0));
    system.EvaluateRhs(system_state, derivative, additional_derivative, dt);

    std::normal_distribution<Real> norm_dist(0.0, 1.0);
    Real delta_W = 0.0;

    for (int i : thread_->GetLoopIndices())
    {
      delta_W = std::sqrt(dt) * norm_dist(mersenne_twister_generator);

      int ii = kS * i;
      system_state[ii] += derivative[ii] * dt;
      system_state[ii + 1] += derivative[ii + 1] * dt;
      system_state[ii + 2] += derivative[ii + 2] * dt + std::sqrt(2.0 * system.D_phi()) * delta_W;
    } // i
  }

 private:

  ThreadSharedMemory *thread_;

};

#endif //SPC2ODEINTEGRATION_STOCHASTICEULERSTEPPERPTR_HPP
