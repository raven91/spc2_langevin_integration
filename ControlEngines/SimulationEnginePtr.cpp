//
// Created by Nikita Kruk on 2019-02-01.
//

#include "SimulationEnginePtr.hpp"
#include "../Steppers/RungeKutta4StepperPtr.hpp"
#include "../Steppers/StochasticRungeKuttaStepperPtr.hpp"
#include "../Observers/BinaryObserverPtr.hpp"
#include "../DynamicalSystems/ChimeraSystemPtr.hpp"
#include "../DynamicalSystems/VortexArraysSystemPtr.hpp"

#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm> // std::copy
#include <iostream>

SimulationEnginePtr::SimulationEnginePtr(ThreadSharedMemory *thread) :
    thread_(thread),
    pbc_config_(thread, kL, kL),
    system_state_size_(kS * kN)
{
  system_state_ = new Real[system_state_size_];
}

SimulationEnginePtr::~SimulationEnginePtr()
{
  delete[] system_state_;
}

void SimulationEnginePtr::InitializeRandomSystemState()
{
  if (thread_->IsRoot())
  {
    const Real two_pi = 2.0 * M_PI;

    std::uniform_real_distribution<Real> unif_real_dist_0L(0.0, kL);
    std::uniform_real_distribution<Real> unif_real_dist_02pi(0.0, two_pi);

    for (int i = 0; i < kN; ++i)
    {
      system_state_[kS * i] = unif_real_dist_0L(mersenne_twister_generator);
      system_state_[kS * i + 1] = unif_real_dist_0L(mersenne_twister_generator);
      system_state_[kS * i + 2] = unif_real_dist_02pi(mersenne_twister_generator);
    } // i

    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_, system_state_size_);
}

void SimulationEnginePtr::InitializeSystemStateFromPreviousSolution(const std::string &file_name)
{
  if (thread_->IsRoot())
  {
    std::ifstream prev_solution_file(file_name, std::ios::binary | std::ios::in);
    assert(prev_solution_file.is_open());

//	std::streampos size = prev_solution_file.tellg();
    prev_solution_file.seekg((1 + kS * kN) * 2900 * sizeof(float), std::ios::beg);

    float t = 0.0f;
    std::vector<float> system_state_float(kS * kN, 0.0f);

    prev_solution_file.read((char *) &t, sizeof(float));
    prev_solution_file.read((char *) &system_state_float[0], kS * kN * sizeof(float));
    std::copy(system_state_float.begin(), system_state_float.end(), &system_state_[0]);

    prev_solution_file.close();
    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_, system_state_size_);
}

void SimulationEnginePtr::RunChimeraSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

  Real v_0 = 2.0;
  Real sigma = 1.0;
  Real rho = 1.0;
  Real alpha = 0.0;
  Real D_phi = 0.4;
  Real rho_0 = 2.0;

  Real t_0 = 0.0;
  Real t_1 = 100000.0;
  Real dt = 0.1;

  int trial = 0;
//  for (int trial = 0; trial < 10; ++trial)
  {
    InitializeRandomSystemState();
//    InitializeSystemStateFromPreviousSolution("/Volumes/Kruk/spc2/spc2OdeIntegration/from_uniform_initial_condition/v0_0.01_sigma_1_rho_0.01_alpha_1.3_Dphi_0.02_N_50000_0_0.bin");
//    RungeKutta4StepperPtr stepper(thread_);
    StochasticRungeKuttaStepperPtr stepper(thread_);
    ChimeraSystemPtr chimera_system(thread_, v_0, sigma, rho, alpha, D_phi, rho_0, pbc_config_);
    BinaryObserverPtr binary_observer(thread_, v_0, sigma, rho, alpha, D_phi, rho_0, pbc_config_, dt, trial);
//    IntegrateConst(stepper, chimera_system, t_0, t_1, dt, binary_observer);

    Real t = t_0;
    binary_observer.SaveSystemState(system_state_, system_state_size_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(chimera_system, system_state_, t, dt);
      // with the linked list technique keep all positions under periodic boundaries
      pbc_config_.ApplyPeriodicBoundaryConditions(system_state_, system_state_size_);
//    thread_->SynchronizeVector(system_state_);
      binary_observer.SaveSystemState(system_state_, system_state_size_, t);
      binary_observer.SaveSummaryStatistics(system_state_, system_state_size_, t);
    } // t
  } // trial

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}

void SimulationEnginePtr::RunVortexArraySimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

  Real L = 10.0;
  Real velocity = 1.0 / L;
  Real mu_plus = 100.0;
  Real mu_minus = 100.0;
  Real xi_a = 0.4 / L;
  Real xi_r = 0.2 / L;
  Real D_phi = 0.01;
  Real kappa = 20.0;
  Real rho = 2.0 / L;

  Real t_0 = 0.0;
  Real t_1 = 1000.0;
  Real dt = 0.01;

  int trial = 0;
//  for (int trial = 0; trial < 10; ++trial)
  {
    InitializeRandomSystemState();
//    RungeKutta4StepperPtr stepper(thread_);
//    StochasticEulerStepperPtr stepper(thread_);
    StochasticRungeKuttaStepperPtr stepper(thread_);
    VortexArraysSystemPtr
        vortex_array_system(thread_, velocity, mu_plus, mu_minus, xi_a, xi_r, D_phi, kappa, rho, pbc_config_);
    BinaryObserverPtr
        binary_observer(thread_, velocity, mu_plus, mu_minus, xi_a, xi_r, D_phi, kappa, rho, pbc_config_, dt, trial);
//    IntegrateConst(stepper, vortex_array_system, t_0, t_1, dt, binary_observer);

    Real t = t_0;
    binary_observer.SaveSystemState(system_state_, system_state_size_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(vortex_array_system, system_state_, t, dt);
      // with the linked list technique keep all positions under periodic boundaries
      pbc_config_.ApplyPeriodicBoundaryConditions(system_state_, system_state_size_);
//    thread_->SynchronizeVector(system_state_);
      binary_observer.SaveSystemState(system_state_, system_state_size_, t);
      binary_observer.SaveSummaryStatistics(system_state_, system_state_size_, t);
    } // t
  } // trial

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}