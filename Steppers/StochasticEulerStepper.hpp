//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_STOCHASTICEULERSTEPPER_HPP
#define SPC2ODEINTEGRATION_STOCHASTICEULERSTEPPER_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ChimeraSystem.hpp"
#include "../Parallelization/Thread.hpp"

#include <cmath>
#include <algorithm> // std::fill

class StochasticEulerStepper
{
 public:

  explicit StochasticEulerStepper(Thread *thread) :
      thread_(thread),
      buffer_(kS * kN, 0.0)
  {

  }

  ~StochasticEulerStepper()
  {
    buffer_.clear();
  }

  void DoStep(ChimeraSystem &chimera_system, std::vector<Real> &x, Real t, Real dt)
  {
    static std::vector<Real> deterministic(x.size(), 0.0), stochastic(x.size(), 0.0);
    std::fill(deterministic.begin(), deterministic.end(), 0.0);
    std::fill(stochastic.begin(), stochastic.end(), 0.0);

    chimera_system.EvaluateRhs(x, deterministic, deterministic, 0.0, dt);
    chimera_system.AddNoise(x, stochastic, t);

    const std::vector<int> &loop_indices = thread_->GetLoopIndices();
    for (int i : loop_indices)
    {
      x[kS * i] += dt * deterministic[kS * i] + std::sqrt(dt) * stochastic[kS * i];
      x[kS * i + 1] += dt * deterministic[kS * i + 1] + std::sqrt(dt) * stochastic[kS * i + 1];
      x[kS * i + 2] += dt * deterministic[kS * i + 2] + std::sqrt(dt) * stochastic[kS * i + 2];
    } // i
  }

//  void SynchronizeVectorThroughBuffer(std::vector<Real> &vec)
//  {
//    SynchronizeVectorThroughBuffer(vec, buffer_);
//  }

 private:

  Thread *thread_;
  std::vector<Real> buffer_;

//  void SynchronizeVectorThroughBuffer(std::vector<Real> &vec, std::vector<Real> &buf)
//  {
//#if defined(MPI_FAST_INTERACTION)
//    MPI_Gather(&vec[mpi_thread_rank * mpi_num_particles_per_thread * kS],
//               mpi_num_particles_per_thread * kS,
//               MPI_DOUBLE,
//               &buf[0],
//               mpi_num_particles_per_thread * kS,
//               MPI_DOUBLE,
//               mpi_root_rank,
//               MPI_COMM_WORLD);
//    if (mpi_thread_rank == mpi_root_rank)
//    {
//      vec.swap(buf);
//    }
//    MPI_Bcast(&vec[0], (int) vec.size(), MPI_DOUBLE, mpi_root_rank, MPI_COMM_WORLD);
//
//    /*	for (int i = 0; i < mpi_num_threads; ++i)
//     {
//     if (i != mpi_thread_rank)
//     {
//     MPI_Send(&vec[mpi_thread_rank * mpi_num_particles_per_thread * kS], mpi_num_particles_per_thread * kS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
//     }
//     }
//     MPI_Status status;
//     for (int i = 0; i < mpi_num_threads; ++i)
//     {
//     if (i != mpi_thread_rank)
//     {
//     MPI_Recv(&vec[i * mpi_num_particles_per_thread * kS], mpi_num_particles_per_thread * kS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
//     }
//     }
//     */
//#endif
//  }
};

#endif //SPC2ODEINTEGRATION_STOCHASTICEULERSTEPPER_HPP
