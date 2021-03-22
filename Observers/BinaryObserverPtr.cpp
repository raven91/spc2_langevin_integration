//
// Created by Nikita Kruk on 2019-02-01.
//

#include "BinaryObserverPtr.hpp"

#include <mpi.h>
#include <mpio.h>
#include <sstream>
#include <cassert>
#include <cmath>
#include <complex>
#include <numeric> // std::accumulate
#include <iostream>
#include <iterator> // std::back_inserter, std::copy

BinaryObserverPtr::BinaryObserverPtr(ThreadSharedMemory *thread,
                                     Real v_0,
                                     Real sigma,
                                     Real rho,
                                     Real alpha,
                                     Real D_phi,
                                     PeriodicBoundaryConditions &pbc_config,
                                     Real dt,
                                     int trial) :
    thread_(thread),
    pbc_config_(pbc_config),
    output_time_counter_{0, 0},
    output_time_threshold_{100, 100} // mod 1 - save at every dt
{
  integration_step_timer_ = std::chrono::system_clock::now();

#if defined(__linux__) && defined(LICHTENBERG)
  std::string file_folder("/work/scratch/nk59zoce/cpp/spc2OdeIntegration/srk/");
#elif defined(__linux__) && (defined(BCS_CLUSTER))
  std::string file_folder("/home/nkruk/cpp/spc2OdeIntegration/output/");
#elif defined(__linux__)
  std::string file_folder("/home/nikita/Documents/spc2OdeIntegration/");
#elif defined(__APPLE__)
  std::string file_folder("/Users/nikita/Documents/Projects/spc2/spc2OdeIntegration/");
//  std::string file_folder("/Volumes/Kruk/spc2/spc2OdeIntegration/continued/");
#endif

  std::ostringstream simulation_file_name_buffer;
  simulation_file_name_buffer << file_folder << "v0_" << v_0 << "_sigma_" << sigma << "_rho_" << rho
                              << "_alpha_" << alpha << "_Dphi_" << D_phi << "_N_" << kN;
  simulation_file_name_buffer << "_" << 0;
  simulation_file_name_buffer << "_" << trial << ".bin";
  simulation_file_name_ = simulation_file_name_buffer.str();
  if (thread_->IsRoot())
  {
    std::remove(simulation_file_name_.c_str());
  }
  simulation_file_.open(simulation_file_name_, std::ios::binary | std::ios::out | std::ios::app);
  assert(simulation_file_.is_open());

  std::ostringstream summary_statistics_file_name_buffer;
  summary_statistics_file_name_buffer << file_folder << "summary_statistics_" << "v0_" << v_0 << "_sigma_" << sigma
                                      << "_rho_" << rho << "_alpha_" << alpha << "_Dphi_" << D_phi << "_N_" << kN;
  summary_statistics_file_name_buffer << "_" << 0;
  summary_statistics_file_name_buffer << "_" << trial << ".txt";
  summary_statistics_file_name_ = summary_statistics_file_name_buffer.str();
  if (thread_->IsRoot())
  {
    std::remove(summary_statistics_file_name_.c_str());
  }
  summary_statistics_file_.open(summary_statistics_file_name_, std::ios::out | std::ios::app);
  assert(summary_statistics_file_.is_open());

  MPI_Barrier(MPI_COMM_WORLD);
}

BinaryObserverPtr::BinaryObserverPtr(ThreadSharedMemory *thread,
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
                                     int trial) :
    thread_(thread),
    pbc_config_(pbc_config),
    output_time_counter_{0, 0},
    output_time_threshold_{1, 1} // mod 1 - save at every dt
{
  integration_step_timer_ = std::chrono::system_clock::now();

#if defined(__linux__) && defined(LICHTENBERG)
  std::string file_folder("/work/scratch/nk59zoce/cpp/spc2OdeIntegration/srk/");
#elif defined(__linux__) && (defined(BCS_CLUSTER))
  std::string file_folder("/home/nkruk/cpp/spc2OdeIntegration/output/");
#elif defined(__linux__)
  std::string file_folder("/home/nikita/Documents/spc2OdeIntegration/");
#elif defined(__APPLE__)
  std::string file_folder("/Users/nikita/Documents/Projects/spc2/spc2OdeIntegration/");
#endif

  std::ostringstream simulation_file_name_buffer;
  simulation_file_name_buffer << file_folder << "v_" << velocity << "_muplus_" << mu_plus << "_muminus_" << mu_minus
                              << "_xia_" << xi_a << "_xir_" << xi_r << "_Dphi_" << D_phi << "_kappa_" << kappa
                              << "_rho_" << rho << "_N_" << kN;
  simulation_file_name_buffer << "_" << 0;
  simulation_file_name_buffer << "_" << trial << ".bin";
  simulation_file_name_ = simulation_file_name_buffer.str();
  if (thread_->IsRoot())
  {
    std::remove(simulation_file_name_.c_str());
  }
  simulation_file_.open(simulation_file_name_, std::ios::binary | std::ios::out | std::ios::app);
  assert(simulation_file_.is_open());

  std::ostringstream summary_statistics_file_name_buffer;
  summary_statistics_file_name_buffer << file_folder << "summary_statistics_" << "v_" << velocity << "_muplus_"
                                      << mu_plus << "_muminus_" << mu_minus
                                      << "_xia_" << xi_a << "_xir_" << xi_r << "_Dphi_" << D_phi << "_kappa_" << kappa
                                      << "_rho_" << rho << "_N_" << kN;
  summary_statistics_file_name_buffer << "_" << 0;
  summary_statistics_file_name_buffer << "_" << trial << ".txt";
  summary_statistics_file_name_ = summary_statistics_file_name_buffer.str();
  if (thread_->IsRoot())
  {
    std::remove(summary_statistics_file_name_.c_str());
  }
  summary_statistics_file_.open(summary_statistics_file_name_, std::ios::out | std::ios::app);
  assert(summary_statistics_file_.is_open());

  MPI_Barrier(MPI_COMM_WORLD);
}

BinaryObserverPtr::~BinaryObserverPtr()
{
  if (thread_->IsRoot())
  {
    if (simulation_file_.is_open())
    {
      simulation_file_.close();
    }
    if (summary_statistics_file_.is_open())
    {
      summary_statistics_file_.close();
    }
  }
}

void BinaryObserverPtr::SaveSystemState(Real *const system_state, long size, Real t)
{
  // Sequential file output
  thread_->SynchronizeVector(system_state, size);
  if (thread_->IsRoot())
  {
    if (!(output_time_counter_[0] % output_time_threshold_[0]))
    {
      std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
      std::cout << "value at t = " << t << " integrated in " << elapsed_seconds.count() << "s" << std::endl;
      integration_step_timer_ = std::chrono::system_clock::now();

      float t_float = t;
      static std::vector<float> system_state_float(size, 0.0);//(&system_state[0], &system_state[size]);
      std::copy(&system_state[0], &system_state[size], &system_state_float[0]);

      simulation_file_.write((char *) (&t_float), sizeof(float));
      simulation_file_.write((char *) (&system_state_float[0]), kS * kN * sizeof(float));
    }
    ++output_time_counter_[0];
  }

  // Parallel file output
  /*if (!(output_time_counter_[0] % output_time_threshold_[0]))
  {
    if (thread_->IsRoot())
    {
      std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
      std::cout << "t:" << t << " | integration time:" << elapsed_seconds.count() << "s" << std::endl;
      integration_step_timer_ = std::chrono::system_clock::now();
    }

    RealOutput t_float = (RealOutput) t;
    if (thread_->IsRoot())
    {
      MPI_File output_file;
      MPI_File_open(MPI_COMM_SELF,
                    simulation_file_name_.c_str(),
                    MPI_MODE_APPEND | MPI_MODE_WRONLY,
                    MPI_INFO_NULL,
                    &output_file);
      MPI_File_write(output_file, &t_float, 1, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_close(&output_file);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_File output_file;
    MPI_File_open(MPI_COMM_WORLD,
                  simulation_file_name_.c_str(),
                  MPI_MODE_APPEND | MPI_MODE_WRONLY,
                  MPI_INFO_NULL,
                  &output_file);
    static std::vector<RealOutput> system_state_float(thread_->GetNumberOfParticlesPerMpichThread() * kS, 0.0);
    std::copy(&system_state[thread_->GetNumberOfParticlesPerMpichThread() * kS * thread_->GetRank()],
              &system_state[thread_->GetNumberOfParticlesPerMpichThread() * kS * (thread_->GetRank() + 1)],
              system_state_float.begin());
//              std::back_inserter(system_state_float));
    MPI_Offset position;
    MPI_File_get_position(output_file, &position);
    MPI_File_set_view(output_file,
                      position + thread_->GetNumberOfParticlesPerMpichThread() * kS * thread_->GetRank()
                          * sizeof(RealOutput),
                      MPI_FLOAT,
                      MPI_FLOAT,
                      "native",
                      MPI_INFO_NULL);
    MPI_Request request;
    MPI_File_iwrite(output_file,
                    &system_state_float[0],
                    thread_->GetNumberOfParticlesPerMpichThread() * kS,
                    MPI_FLOAT,
                    &request);
    MPI_Waitall(1, &request, MPI_STATUS_IGNORE);

    MPI_File_close(&output_file);
  }
  ++output_time_counter_[0];*/
}

void BinaryObserverPtr::SaveSummaryStatistics(const Real *const system_state, long size, Real t)
{
  // Sequential computation
  if (thread_->IsRoot())
  {
    if (!(output_time_counter_[1] % output_time_threshold_[1]))
    {
      float t_float = t;
      std::complex<float> complex_order_parameter(0.0f, 0.0f);
      for (int i = 0; i < kN; ++i)
      {
        complex_order_parameter +=
            std::complex<float>((float) std::cos(system_state[kS * i + 2]), (float) std::sin(system_state[kS * i + 2]));
      }
      complex_order_parameter /= kN;

      summary_statistics_file_ << t_float << '\t' << std::abs(complex_order_parameter) << '\t'
                               << std::arg(complex_order_parameter);
      summary_statistics_file_ << std::endl;
    }
    ++output_time_counter_[1];
  }

  // Parallel computation
  /*if (!(output_time_counter_[1] % output_time_threshold_[1]))
  {
    float t_float = (float) t;
    float polar_order_subparameter_x = 0.0f, polar_order_subparameter_y = 0.0f;
    int i = 0, j = 0, k = 0;
    for (int i : thread_->GetLoopIndices())
    {
      polar_order_subparameter_x += (float) std::cos(system_state[kS * i + 2]);
      polar_order_subparameter_y += (float) std::sin(system_state[kS * i + 2]);
    } // idx

    float polar_order_parameter_x = 0.0f, polar_order_parameter_y = 0.0f;
    MPI_Reduce(&polar_order_subparameter_x, &polar_order_parameter_x, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&polar_order_subparameter_y, &polar_order_parameter_y, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (thread_->IsRoot())
    {
      std::complex<float> polar_order_parameter(polar_order_parameter_x / kN, polar_order_parameter_y / kN);
      summary_statistics_file_ << t_float << '\t' << std::abs(polar_order_parameter) << '\t'
                               << std::arg(polar_order_parameter);
      summary_statistics_file_ << std::endl;
    }
  }
  ++output_time_counter_[1];*/
}