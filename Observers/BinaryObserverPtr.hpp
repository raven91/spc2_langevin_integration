//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_BINARYOBSERVERPTR_HPP
#define SPC2ODEINTEGRATION_BINARYOBSERVERPTR_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

#include <string>
#include <fstream>
#include <chrono>

class BinaryObserverPtr
{
 public:

  explicit BinaryObserverPtr(ThreadSharedMemory *thread,
                             Real v_0,
                             Real sigma,
                             Real rho,
                             Real alpha,
                             Real D_phi,
                             Real rho_0,
                             PeriodicBoundaryConditions &pbc_config,
                             Real dt,
                             int trial);
  explicit BinaryObserverPtr(ThreadSharedMemory *thread,
                             Real velocity,
                             Real mu_plus,
                             Real mu_minus,
                             Real xi_a,
                             Real xi_r,
                             Real D_phi,
                             Real kappa,
                             Real rho,
                             PeriodicBoundaryConditions &pbc_config,
                             Real dt,
                             int trial);
  ~BinaryObserverPtr();

  void SaveSystemState(Real *const system_state, long size, Real t);
  void SaveSummaryStatistics(const Real *const system_state, long size, Real t);

 private:

  ThreadSharedMemory *thread_;
  PeriodicBoundaryConditions &pbc_config_;
  std::chrono::time_point<std::chrono::system_clock> integration_step_timer_;
  int output_time_counter_[2];
  int output_time_threshold_[2];

  std::string simulation_file_name_;
  std::ofstream simulation_file_;

  std::string summary_statistics_file_name_;
  std::ofstream summary_statistics_file_;

};

#endif //SPC2ODEINTEGRATION_BINARYOBSERVERPTR_HPP
