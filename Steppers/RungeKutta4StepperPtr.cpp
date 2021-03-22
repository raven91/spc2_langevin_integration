//
// Created by Nikita Kruk on 2019-02-01.
//

#include "RungeKutta4StepperPtr.hpp"

RungeKutta4StepperPtr::RungeKutta4StepperPtr(ThreadSharedMemory *thread) :
    thread_(thread),
    k_1_(kS * kN, 0.0),
    k_2_(kS * kN, 0.0),
    k_3_(kS * kN, 0.0),
    k_4_(kS * kN, 0.0)
{

}

RungeKutta4StepperPtr::~RungeKutta4StepperPtr()
{
  k_1_.clear();
  k_2_.clear();
  k_3_.clear();
  k_4_.clear();
}

void RungeKutta4StepperPtr::DoStep(ChimeraSystemPtr &chimera_system, Real *const system_state, Real t, Real dt)
{
  chimera_system.EvaluateRhs(system_state, k_1_, k_1_, 0.0, dt);
  chimera_system.EvaluateRhs(system_state, k_1_, k_2_, 0.5, dt);
  chimera_system.EvaluateRhs(system_state, k_2_, k_3_, 0.5, dt);
  chimera_system.last_coefficient_ = true;
  chimera_system.EvaluateRhs(system_state, k_3_, k_4_, 1.0, dt);
  chimera_system.last_coefficient_ = false;

  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    int ii = kS * i;
    system_state[ii] += (k_1_[ii] + 2.0 * k_2_[ii] + 2.0 * k_3_[ii] + k_4_[ii]) * dt / 6.0;
    system_state[ii + 1] += (k_1_[ii + 1] + 2.0 * k_2_[ii + 1] + 2.0 * k_3_[ii + 1] + k_4_[ii + 1]) * dt / 6.0;
    system_state[ii + 2] += (k_1_[ii + 2] + 2.0 * k_2_[ii + 2] + 2.0 * k_3_[ii + 2] + k_4_[ii + 2]) / 6.0 * dt;
//    system_state[ii + 3] += (k_1_[ii + 3] + 2.0 * k_2_[ii + 3] + 2.0 * k_3_[ii + 3] + k_4_[ii + 3]) / 6.0 * dt;
  } // i
}

void RungeKutta4StepperPtr::DoStep(VortexArraysSystemPtr &vortex_array_system,
                                   Real *const system_state,
                                   Real t,
                                   Real dt)
{
  vortex_array_system.EvaluateRhs(system_state, k_1_, k_1_, 0.0, dt);
  vortex_array_system.EvaluateRhs(system_state, k_1_, k_2_, 0.5, dt);
  vortex_array_system.EvaluateRhs(system_state, k_2_, k_3_, 0.5, dt);
  vortex_array_system.last_coefficient_ = true;
  vortex_array_system.EvaluateRhs(system_state, k_3_, k_4_, 1.0, dt);
  vortex_array_system.last_coefficient_ = false;

  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    int ii = kS * i;
    system_state[ii] += (k_1_[ii] + 2.0 * k_2_[ii] + 2.0 * k_3_[ii] + k_4_[ii]) * dt / 6.0;
    system_state[ii + 1] += (k_1_[ii + 1] + 2.0 * k_2_[ii + 1] + 2.0 * k_3_[ii + 1] + k_4_[ii + 1]) * dt / 6.0;
    system_state[ii + 2] += (k_1_[ii + 2] + 2.0 * k_2_[ii + 2] + 2.0 * k_3_[ii + 2] + k_4_[ii + 2]) / 6.0 * dt;
  } // i
}