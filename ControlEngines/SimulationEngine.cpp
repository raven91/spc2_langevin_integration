//
// Created by Nikita Kruk on 2019-02-01.
//

#include "SimulationEngine.hpp"
#include "../Steppers/StochasticEulerStepper.hpp"
#include "../Steppers/StochasticRungeKuttaStepper.hpp"
#include "../Steppers/RungeKutta4Stepper.hpp"
#include "../DynamicalSystems/ChimeraSystem.hpp"
#include "../Observers/BinaryObserver.hpp"

#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm> // std::copy
#include <iostream>

SimulationEngine::SimulationEngine(Thread *thread) :
    thread_(thread),
    pbc_config_(thread, 1.0, 1.0),
    system_state_(kS * kN, 0.0)
{

}

SimulationEngine::~SimulationEngine()
{
  system_state_.clear();
}

void SimulationEngine::InitializeRandomSystemState()
{
  if (thread_->IsRoot())
  {
    const Real two_pi = 2.0 * M_PI;

    std::uniform_real_distribution<Real> unif_real_dist_01(0.0, 1.0);
    std::uniform_real_distribution<Real> unif_real_dist_02pi(0.0, two_pi);

    for (int i = 0; i < kN; ++i)
    {
      system_state_[kS * i] = unif_real_dist_01(mersenne_twister_generator);
      system_state_[kS * i + 1] = unif_real_dist_01(mersenne_twister_generator);
      system_state_[kS * i + 2] = unif_real_dist_02pi(mersenne_twister_generator);
    } // i

    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_);
}

void SimulationEngine::InitializeSystemStateFromFile(const std::string &file_name)
{
  if (thread_->IsRoot())
  {
    std::ifstream init_cond_file(file_name, std::ios::binary | std::ios::in);
    assert(init_cond_file.is_open());

    //	std::streampos size = init_cond_file.tellg();
    init_cond_file.seekg(0, std::ios::beg);

//		float t = 0.0f;
    std::vector<float> system_state_float(kS * kN, 0.0f);

//		init_cond_file.read((char *)&t, sizeof(float));
    init_cond_file.read((char *) &system_state_float[0], kS * kN * sizeof(float));
    std::copy(system_state_float.begin(), system_state_float.end(), system_state_.begin());

    init_cond_file.close();
    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_);
}

void SimulationEngine::RunSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

#if defined(MPI_PARAMETER_SCAN)
  Real v_0 = 1.0;
  Real sigma = 1.0;
    Real rho = 0.3;// + (mpi_thread_rank % 10) * 0.05;
    Real alpha = 1.54;// + (mpi_thread_rank / 10) * 0.02;
    Real D_phi = 0.001 + (mpi_thread_rank) * 0.001;//0.0005 + (mpi_thread_rank) * 0.0005;
#else
  Real v_0 = 1.0;
  Real sigma = 4.0;
  Real rho = 0.05;
  Real alpha = 1.54;
  Real D_phi = 0.01;
#endif

  Real t_0 = 0.0;
  Real t_1 = 1000.0;
  Real dt = 0.01;
//	Real abs_err = 1.0e-10;
//	Real rel_err = 1.0e-6;

  int trial = 0;
//  for (int trial = 0; trial < 10; ++trial)
  {
    InitializeRandomSystemState();
//	InitializeSystemStateFromFile(std::string("/Users/nikita/Documents/spc2/spc2ContinuationMethod/startup/localized_initial_condition_sigma_1_rho_0.3_alpha_1.54_Dphi_0_N_1000_0.bin"));
    StochasticRungeKuttaStepper stepper(thread_);
//    RungeKutta4Stepper stepper(thread_);
    ChimeraSystem chimera_system(thread_, v_0, sigma, rho, alpha, D_phi, pbc_config_);
    BinaryObserver binary_observer(thread_, v_0, sigma, rho, alpha, D_phi, pbc_config_, dt, trial);

    Real t = t_0;
    binary_observer.SaveSystemState(system_state_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(chimera_system, system_state_, t, dt);
      // with the linked list technique keep all positions under periodic boundaries
      pbc_config_.ApplyPeriodicBoundaryConditions(system_state_);
      thread_->SynchronizeVector(system_state_);
      if (thread_->IsRoot())
      {
        binary_observer.SaveSystemState(system_state_, t);
        binary_observer.SaveSummaryStatistics(system_state_, t);
      }
    } // t

#if defined(MPI_PARAMETER_SCAN)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}