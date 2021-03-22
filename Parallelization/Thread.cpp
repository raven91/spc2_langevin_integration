//
// Created by Nikita Kruk on 2019-02-04.
//

#include "Thread.hpp"

#include <cassert>
#include <iostream>
#include <numeric>  // std::iota

Thread::Thread(int argc, char **argv)
{
  root_rank_ = 0;
  rank_ = 0;
  number_of_mpich_threads_ = 1;
  assert(!(kN
      % number_of_mpich_threads_));  // number of particles must be divisible by the number of threads for this kind of parallelization
  number_of_particles_per_mpich_thread_ = kN / number_of_mpich_threads_;

  loop_indices_ = std::vector<int>(number_of_particles_per_mpich_thread_, 0);
  std::iota(loop_indices_.begin(), loop_indices_.end(), rank_ * number_of_particles_per_mpich_thread_);
}

Thread::~Thread()
{
  loop_indices_.clear();
}

int Thread::GetRank()
{
  return rank_;
}

int Thread::GetNumberOfMpichThreads()
{
  return number_of_mpich_threads_;
}

int Thread::GetNumberOfParticlesPerMpichThread()
{
  return number_of_particles_per_mpich_thread_;
}

const std::vector<int> &Thread::GetLoopIndices() const
{
  return loop_indices_;
}

bool Thread::IsRoot() const
{
  return true;
}

void Thread::SynchronizeVector(std::vector<Real> &vec)
{

}

void Thread::SynchronizeVectorThoughBuffer(std::vector<Real> &vec, std::vector<Real> &buf)
{

}

void Thread::BroadcastVector(std::vector<Real> &vec)
{

}

void Thread::SynchronizePrelinkedList(std::vector<int> &pre_linked_list) const
{

}