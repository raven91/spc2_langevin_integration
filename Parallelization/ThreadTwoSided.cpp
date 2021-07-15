//
// Created by Nikita Kruk on 2019-02-04.
//

#include "ThreadTwoSided.hpp"

#include <mpi.h>
#include <vector>
#include <numeric>  // std::iota
#include <cassert>

ThreadTwoSided::ThreadTwoSided(int argc, char **argv) :
    Thread(argc, argv)
{
  root_rank_ = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_mpich_threads_);
  assert(!(kN
      % number_of_mpich_threads_));  // number of particles must be divisible by the number of threads for this kind of parallelization
  number_of_particles_per_mpich_thread_ = kN / number_of_mpich_threads_;

  loop_indices_ = std::vector<int>(number_of_particles_per_mpich_thread_, 0);
  std::iota(loop_indices_.begin(), loop_indices_.end(), rank_ * number_of_particles_per_mpich_thread_);
}

ThreadTwoSided::~ThreadTwoSided()
{

}

bool ThreadTwoSided::IsRoot() const
{
  return (rank_ == root_rank_);
}

void ThreadTwoSided::SynchronizeVector(std::vector<Real> &vec)
{
  std::vector<Real> buf(kN * kS, 0.0);
  SynchronizeVectorThroughBuffer(vec, buf);
}

void ThreadTwoSided::SynchronizeVectorThroughBuffer(std::vector<Real> &vec, std::vector<Real> &buf)
{
  MPI_Gather(&vec[rank_ * number_of_particles_per_mpich_thread_ * kS],
             number_of_particles_per_mpich_thread_ * kS,
             MPI_DOUBLE,
             &buf[0],
             number_of_particles_per_mpich_thread_ * kS,
             MPI_DOUBLE,
             root_rank_,
             MPI_COMM_WORLD);
  if (IsRoot())
  {
    vec.swap(buf);
  }
  MPI_Bcast(&vec[0], (int) vec.size(), MPI_DOUBLE, root_rank_, MPI_COMM_WORLD);

  /*	for (int i = 0; i < number_of_mpich_threads_; ++i)
	{
		if (i != rank_)
		{
			MPI_Send(&vec[rank_ * number_of_particles_per_mpich_thread_ * kS], number_of_particles_per_mpich_thread_ * kS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Status status;
	for (int i = 0; i < number_of_mpich_threads_; ++i)
	{
		if (i != rank_)
		{
			MPI_Recv(&vec[i * number_of_particles_per_mpich_thread_ * kS], number_of_particles_per_mpich_thread_ * kS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
		}
	}
*/
}

void ThreadTwoSided::BroadcastVector(std::vector<Real> &vec)
{
  MPI_Bcast(&vec[0], (int) vec.size(), MPI_DOUBLE, root_rank_, MPI_COMM_WORLD);
}

void ThreadTwoSided::SynchronizePrelinkedList(std::vector<int> &pre_linked_list) const
{
  if (!IsRoot())
  {
    MPI_Gather(&pre_linked_list[rank_ * number_of_particles_per_mpich_thread_],
               number_of_particles_per_mpich_thread_,
               MPI_INT,
               &pre_linked_list[0],
               number_of_particles_per_mpich_thread_,
               MPI_INT,
               root_rank_,
               MPI_COMM_WORLD);
  } else
  {
    MPI_Gather(MPI_IN_PLACE,
               number_of_particles_per_mpich_thread_,
               MPI_INT,
               &pre_linked_list[0],
               number_of_particles_per_mpich_thread_,
               MPI_INT,
               root_rank_,
               MPI_COMM_WORLD);
  }
  MPI_Bcast(&pre_linked_list[0], (int) pre_linked_list.size(), MPI_INT, root_rank_, MPI_COMM_WORLD);
}