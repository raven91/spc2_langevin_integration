//
// Created by Nikita Kruk on 2019-02-01.
//

#include "BinaryObserver.hpp"

#include <sstream>
#include <cassert>
#include <cmath>
#include <complex>
#include <numeric> // std::accumulate
#include <iostream>
#include <iterator> // std::copy

BinaryObserver::BinaryObserver(Thread *thread,
                               Real v_0,
                               Real sigma,
                               Real rho,
                               Real alpha,
                               Real D_phi,
                               Real rho_0,
                               PeriodicBoundaryConditions &pbc_config,
                               Real dt,
                               int trial) :
    thread_(thread),
    pbc_config_(pbc_config),
    output_time_counter_{0, 0},
    output_time_threshold_{1, 1} // mod 1 - save at every dt
{
//#if defined(MPI_FAST_INTERACTION)
  if (thread_->IsRoot())
  {
//#endif

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
    simulation_file_name_buffer << file_folder << "v0_" << v_0 << "_sigma_" << sigma << "_rho_" << rho
                                << "_alpha_" << alpha << "_Dphi_" << D_phi << "_N_" << kN << "_rho0_" << rho_0;
#if defined(MPI_PARAMETER_SCAN)
    simulation_file_name_buffer << "_" << mpi_thread_rank;
#else
    simulation_file_name_buffer << "_" << 0;
#endif
    simulation_file_name_buffer << "_" << trial << ".bin";
    simulation_file_name_ = simulation_file_name_buffer.str();
    std::remove(simulation_file_name_.c_str());

    simulation_file_.open(simulation_file_name_, std::ios::binary | std::ios::out | std::ios::app);
    assert(simulation_file_.is_open());

    std::ostringstream summary_statistics_file_name_buffer;
    summary_statistics_file_name_buffer << file_folder << "summary_statistics_" << "v0_" << v_0 << "_sigma_" << sigma
                                        << "_rho_" << rho << "_alpha_" << alpha << "_Dphi_" << D_phi << "_N_" << kN
                                        << "_rho0_" << rho_0;
#if defined(MPI_PARAMETER_SCAN)
    summary_statistics_file_name_buffer << "_" << mpi_thread_rank;
#else
    summary_statistics_file_name_buffer << "_" << 0;
#endif
    summary_statistics_file_name_buffer << "_" << trial << ".txt";
    summary_statistics_file_name_ = summary_statistics_file_name_buffer.str();
    std::remove(summary_statistics_file_name_.c_str());

    summary_statistics_file_.open(summary_statistics_file_name_, std::ios::out | std::ios::app);
    assert(summary_statistics_file_.is_open());

//#if defined(MPI_FAST_INTERACTION)
  }
//#endif
}

BinaryObserver::BinaryObserver(Thread *thread,
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
//#if defined(MPI_FAST_INTERACTION)
  if (thread_->IsRoot())
  {
//#endif

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
#if defined(MPI_PARAMETER_SCAN)
    simulation_file_name_buffer << "_" << mpi_thread_rank;
#else
    simulation_file_name_buffer << "_" << 0;
#endif
    simulation_file_name_buffer << "_" << trial << ".bin";
    simulation_file_name_ = simulation_file_name_buffer.str();
    std::remove(simulation_file_name_.c_str());

    simulation_file_.open(simulation_file_name_, std::ios::binary | std::ios::out | std::ios::app);
    assert(simulation_file_.is_open());

    std::ostringstream summary_statistics_file_name_buffer;
    summary_statistics_file_name_buffer << file_folder << "summary_statistics_" << "v_" << velocity << "_muplus_"
                                        << mu_plus << "_muminus_" << mu_minus
                                        << "_xia_" << xi_a << "_xir_" << xi_r << "_Dphi_" << D_phi << "_kappa_"
                                        << kappa << "_rho_" << rho << "_N_" << kN;
#if defined(MPI_PARAMETER_SCAN)
    summary_statistics_file_name_buffer << "_" << mpi_thread_rank;
#else
    summary_statistics_file_name_buffer << "_" << 0;
#endif
    summary_statistics_file_name_buffer << "_" << trial << ".txt";
    summary_statistics_file_name_ = summary_statistics_file_name_buffer.str();
    std::remove(summary_statistics_file_name_.c_str());

    summary_statistics_file_.open(summary_statistics_file_name_, std::ios::out | std::ios::app);
    assert(summary_statistics_file_.is_open());

//#if defined(MPI_FAST_INTERACTION)
  }
//#endif
}

BinaryObserver::~BinaryObserver()
{
//#if defined(MPI_FAST_INTERACTION)
  if (thread_->IsRoot())
  {
//#endif

    if (simulation_file_.is_open())
    {
      simulation_file_.close();
    }

    if (summary_statistics_file_.is_open())
    {
      summary_statistics_file_.close();
    }

//#if defined(MPI_FAST_INTERACTION)
  }
//#endif
}

void BinaryObserver::SaveSystemState(std::vector<Real> &system_state, Real t)
{
//	if (t >= 9900.0)
  if (!(output_time_counter_[0] % output_time_threshold_[0]))
  {
    std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
    std::cout << "value at t = " << t << " integrated in " << elapsed_seconds.count() << "s" << std::endl;
    integration_step_timer_ = std::chrono::system_clock::now();

    float t_float = t;
    static std::vector<float> system_state_float(system_state.size(), 0.0);//(system_state.begin(), system_state.end());
    std::copy(system_state.begin(), system_state.end(), system_state_float.begin());

    simulation_file_.write((char *) (&t_float), sizeof(float));
    simulation_file_.write((char *) (&system_state_float[0]), kS * kN * sizeof(float));
  }
  ++output_time_counter_[0];
}

void BinaryObserver::SaveSummaryStatistics(std::vector<Real> &system_state, Real t)
{
  if (!(output_time_counter_[1] % output_time_threshold_[1]))
  {
    float t_float = t;
    std::complex<float> complex_order_parameter(0.0f, 0.0f);
    for (int i = 0; i < kN; ++i)
    {
      complex_order_parameter +=
          std::complex<float>((float) std::cos(system_state[kS * i + 2]), (float) std::sin(system_state[kS * i + 2]));
    } // i
    complex_order_parameter /= kN;

    summary_statistics_file_ << t_float << '\t' << std::abs(complex_order_parameter) << '\t'
                             << std::arg(complex_order_parameter);

//		float rotational_order_parameter = 0.0f;
//		std::vector<float> r_i(2, 0.0f), r_g(2, 0.0f), v_i(2, 0.0f);
//		for (size_t i = 0; i < kN; ++i)
//		{
//			r_g[0] += system_state[kS * i    ];
//			r_g[1] += system_state[kS * i + 1];
//		}
//		r_g[0] /= kN;
//		r_g[1] /= kN;
//		for (size_t i = 0; i < kN; ++i)
//		{
//			r_i[0] = system_state[kS * i    ];
//			r_i[1] = system_state[kS * i + 1];
//			v_i[0] = std::cos(system_state[kS * i + 2]);
//			v_i[1] = std::sin(system_state[kS * i + 2]);
//
//			rotational_order_parameter += (r_i[0] - r_g[0]) * v_i[1] - (r_i[1] - r_g[1]) * v_i[0];
//		}
//		rotational_order_parameter /= kN;
//
//		summary_statistics_file_ << '\t' << rotational_order_parameter;

    summary_statistics_file_ << std::endl;
  }
  ++output_time_counter_[1];
}
