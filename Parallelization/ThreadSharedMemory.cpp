//
// Created by Nikita Kruk on 2019-02-04.
//

#include "ThreadSharedMemory.hpp"

#include <mpi.h>
#include <iostream>
#include <cassert>
#include <algorithm> // std::copy

ThreadSharedMemory::ThreadSharedMemory(int argc, char **argv) :
    Thread(argc, argv)
{
  root_rank_ = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_mpich_threads_);
  assert(!(kN
      % number_of_mpich_threads_));  // width by height must be divisible by the number of threads for this kind of parallelization
  number_of_particles_per_mpich_thread_ = kN / number_of_mpich_threads_;

  loop_indices_ = std::vector<int>(number_of_particles_per_mpich_thread_, 0);
  std::iota(loop_indices_.begin(), loop_indices_.end(), rank_ * number_of_particles_per_mpich_thread_);

  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shared_communicator_);
  MPI_Comm_rank(shared_communicator_, &shared_rank_);
  MPI_Comm_size(shared_communicator_, &number_of_shared_mpich_threads_);
  shared_root_rank_ = 0;
  shared_ranks_ = std::vector<int>((unsigned long) number_of_mpich_threads_, 0);
  MPI_Gather(&shared_rank_, 1, MPI_INT, &shared_ranks_[0], 1, MPI_INT, root_rank_, MPI_COMM_WORLD);
  MPI_Bcast(&shared_ranks_[0], (int) shared_ranks_.size(), MPI_INT, root_rank_, MPI_COMM_WORLD);
  numbers_of_shared_mpich_threads_ = std::vector<int>((unsigned long) number_of_mpich_threads_, 0);
  MPI_Gather(&number_of_shared_mpich_threads_,
             1,
             MPI_INT,
             &numbers_of_shared_mpich_threads_[0],
             1,
             MPI_INT,
             root_rank_,
             MPI_COMM_WORLD);
  MPI_Bcast(&numbers_of_shared_mpich_threads_[0],
            (int) numbers_of_shared_mpich_threads_.size(),
            MPI_INT,
            root_rank_,
            MPI_COMM_WORLD);
  for (int i = 0; i < shared_ranks_.size(); ++i)
  {
    if (shared_ranks_[i] == shared_root_rank_)
    {
      shared_root_rank_indexes_.push_back(i);
    }
  } // i
}

ThreadSharedMemory::~ThreadSharedMemory()
{
  for (std::pair<const std::string, MPI_Win> &window : windows_)
  {
    MPI_Win_free(&(window.second));
  } // window
  MPI_Comm_free(&shared_communicator_);
}

int ThreadSharedMemory::GetRootRank()
{
  return root_rank_;
}

int ThreadSharedMemory::GetSharedRootRank()
{
  return shared_root_rank_;
}

MPI_Comm &ThreadSharedMemory::GetSharedCommunicator()
{
  return shared_communicator_;
}

MPI_Win &ThreadSharedMemory::GetWindow(const std::string &window_name)
{
  return windows_[window_name];
}

bool ThreadSharedMemory::IsRoot() const
{
  return (rank_ == root_rank_);
}

bool ThreadSharedMemory::IsSharedRoot()
{
  return (shared_rank_ == shared_root_rank_);
}

void ThreadSharedMemory::SynchronizeVector(std::vector<Real> &vec)
{
  Real *const p = &vec[0];
  SynchronizeVector(p, vec.size());
}

void ThreadSharedMemory::SynchronizeVector(Real *const vec, long size)
{
  MPI_Win win;
  if (rank_ == root_rank_)
  {
    MPI_Win_create(&vec[0], size * sizeof(Real), sizeof(Real), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Win_fence(0, win);
    MPI_Win_fence(0, win);
  } else
  {
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Put(&vec[rank_ * number_of_particles_per_mpich_thread_ * kS],
            number_of_particles_per_mpich_thread_ * kS,
            MPI_DOUBLE,
            root_rank_,
            rank_ * number_of_particles_per_mpich_thread_ * kS,
            number_of_particles_per_mpich_thread_ * kS,
            MPI_DOUBLE,
            win);
    MPI_Win_fence(0, win);
    MPI_Get(&vec[0], size, MPI_DOUBLE, root_rank_, 0, size, MPI_DOUBLE, win);
    MPI_Win_fence(0, win);
  }
  MPI_Win_free(&win);
}

void ThreadSharedMemory::SynchronizeVectorThroughBuffer(std::vector<Real> &vec, std::vector<Real> &buf)
{

}

void ThreadSharedMemory::SynchronizeVectorThroughoutClusters(Real *const vec)
{
  // if there are threads on different nodes
  if (shared_root_rank_indexes_.size() > 1)
  {
    if (IsRoot())
    {
      std::vector<MPI_Request> requests(shared_root_rank_indexes_.size() - 1);
      for (int i = 1; i < shared_root_rank_indexes_.size(); ++i)
      {
//        MPI_Recv(&vec[number_of_particles_per_mpich_thread_ * kS * shared_root_rank_indexes_[i]],
//                 number_of_particles_per_mpich_thread_ * kS * numbers_of_shared_mpich_threads_[shared_root_rank_indexes_[i]],
//                 MPI_DOUBLE,
//                 shared_root_rank_indexes_[i],
//                 0,
//                 MPI_COMM_WORLD,
//                 MPI_STATUS_IGNORE);
        MPI_Irecv(&vec[number_of_particles_per_mpich_thread_ * kS * shared_root_rank_indexes_[i]],
                  number_of_particles_per_mpich_thread_ * kS
                      * numbers_of_shared_mpich_threads_[shared_root_rank_indexes_[i]],
                  MPI_DOUBLE,
                  shared_root_rank_indexes_[i],
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else if (IsSharedRoot())
    {
      MPI_Request request;
//      MPI_Send(&vec[number_of_particles_per_mpich_thread_ * kS * rank_],
//               number_of_particles_per_mpich_thread_ * kS * numbers_of_shared_mpich_threads_[rank_],
//               MPI_DOUBLE,
//               root_rank_,
//               0,
//               MPI_COMM_WORLD);
      MPI_Isend(&vec[number_of_particles_per_mpich_thread_ * kS * rank_],
                number_of_particles_per_mpich_thread_ * kS * numbers_of_shared_mpich_threads_[rank_],
                MPI_DOUBLE,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }

    if (IsRoot())
    {
      std::vector<MPI_Request> requests(shared_root_rank_indexes_.size() - 1);
      for (int i = 1; i < shared_root_rank_indexes_.size(); ++i)
      {
//        MPI_Send(&vec[0],
//                 number_of_particles_per_mpich_thread_ * kS * number_of_mpich_threads_,
//                 MPI_DOUBLE,
//                 shared_root_rank_indexes_[i],
//                 0,
//                 MPI_COMM_WORLD);
        MPI_Isend(&vec[0],
                  number_of_particles_per_mpich_thread_ * kS * number_of_mpich_threads_,
                  MPI_DOUBLE,
                  shared_root_rank_indexes_[i],
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else if (IsSharedRoot())
    {
      MPI_Request request;
//      MPI_Recv(&vec[0],
//               number_of_particles_per_mpich_thread_ * kS * number_of_mpich_threads_,
//               MPI_DOUBLE,
//               root_rank_,
//               0,
//               MPI_COMM_WORLD,
//               MPI_STATUS_IGNORE);
      MPI_Irecv(&vec[0],
                number_of_particles_per_mpich_thread_ * kS * number_of_mpich_threads_,
                MPI_DOUBLE,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD); // for all the non-root processes
  }
}

void ThreadSharedMemory::SynchronizeVectorThroughoutClusters(int *const vec)
{
  // if there are threads on different nodes
  if (shared_root_rank_indexes_.size() > 1)
  {
    if (IsRoot())
    {
      std::vector<MPI_Request> requests(shared_root_rank_indexes_.size() - 1);
      for (int i = 1; i < shared_root_rank_indexes_.size(); ++i)
      {
//        MPI_Recv(&vec[number_of_particles_per_mpich_thread_ * shared_root_rank_indexes_[i]],
//                 number_of_particles_per_mpich_thread_ * numbers_of_shared_mpich_threads_[shared_root_rank_indexes_[i]],
//                 MPI_INT,
//                 shared_root_rank_indexes_[i],
//                 0,
//                 MPI_COMM_WORLD,
//                 MPI_STATUS_IGNORE);
        MPI_Irecv(&vec[number_of_particles_per_mpich_thread_ * shared_root_rank_indexes_[i]],
                  number_of_particles_per_mpich_thread_
                      * numbers_of_shared_mpich_threads_[shared_root_rank_indexes_[i]],
                  MPI_INT,
                  shared_root_rank_indexes_[i],
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else if (IsSharedRoot())
    {
      MPI_Request request;
//      MPI_Send(&vec[number_of_particles_per_mpich_thread_ * rank_],
//               number_of_particles_per_mpich_thread_ * numbers_of_shared_mpich_threads_[rank_],
//               MPI_INT,
//               root_rank_,
//               0,
//               MPI_COMM_WORLD);
      MPI_Isend(&vec[number_of_particles_per_mpich_thread_ * rank_],
                number_of_particles_per_mpich_thread_ * numbers_of_shared_mpich_threads_[rank_],
                MPI_INT,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }

    if (IsRoot())
    {
      std::vector<MPI_Request> requests(shared_root_rank_indexes_.size() - 1);
      for (int i = 1; i < shared_root_rank_indexes_.size(); ++i)
      {
//        MPI_Send(&vec[0],
//                 number_of_particles_per_mpich_thread_ * number_of_mpich_threads_,
//                 MPI_INT,
//                 shared_root_rank_indexes_[i],
//                 0,
//                 MPI_COMM_WORLD);
        MPI_Isend(&vec[0],
                  number_of_particles_per_mpich_thread_ * number_of_mpich_threads_,
                  MPI_INT,
                  shared_root_rank_indexes_[i],
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else if (IsSharedRoot())
    {
      MPI_Request request;
//      MPI_Recv(&vec[0],
//               number_of_particles_per_mpich_thread_ * number_of_mpich_threads_,
//               MPI_INT,
//               root_rank_,
//               0,
//               MPI_COMM_WORLD,
//               MPI_STATUS_IGNORE);
      MPI_Irecv(&vec[0],
                number_of_particles_per_mpich_thread_ * number_of_mpich_threads_,
                MPI_INT,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD); // for all the non-root processes
  }
}

void ThreadSharedMemory::BroadcastVector(std::vector<Real> &vec)
{
  MPI_Win win;
  if (rank_ == root_rank_)
  {
    MPI_Win_create(&vec[0], vec.size() * sizeof(Real), sizeof(Real), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Win_fence(0, win);
  } else
  {
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Get(&vec[0], vec.size(), MPI_DOUBLE, root_rank_, 0, vec.size(), MPI_DOUBLE, win);
    MPI_Win_fence(0, win);
  }
  MPI_Win_free(&win);
}

void ThreadSharedMemory::BroadcastVector(Real *const vec, long size)
{
  MPI_Win win;
  if (IsRoot())
  {
    MPI_Win_create(&vec[0], size * sizeof(Real), sizeof(Real), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Win_fence(0, win);
  } else
  {
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Get(&vec[0], size, MPI_DOUBLE, root_rank_, 0, size, MPI_DOUBLE, win);
    MPI_Win_fence(0, win);
  }
  MPI_Win_free(&win);
}

void ThreadSharedMemory::SynchronizePrelinkedList(std::vector<int> &pre_linked_list) const
{

}

void ThreadSharedMemory::AllocateSharedWindow(MPI_Aint size, Real *&vec, const std::string &window_name)
{
  if (windows_.find(window_name) == windows_.end())
  {
    MPI_Win win;
    MPI_Win_allocate_shared(size * sizeof(Real),
                            sizeof(Real),
                            MPI_INFO_NULL,
                            shared_communicator_,
                            &vec,
                            &win);
    windows_[window_name] = win;
  } else
  {
    std::cout << "Window '" << window_name << "' already exists. The previous window is kept." << std::endl;
  }
}

void ThreadSharedMemory::AllocateSharedWindow(MPI_Aint size, int *&vec, const std::string &window_name)
{
  if (windows_.find(window_name) == windows_.end())
  {
    MPI_Win win;
    MPI_Win_allocate_shared(size * sizeof(int),
                            sizeof(int),
                            MPI_INFO_NULL,
                            shared_communicator_,
                            &vec,
                            &win);
    windows_[window_name] = win;
  } else
  {
    std::cout << "Window '" << window_name << "' already exists. The previous window is kept." << std::endl;
  }
}

void ThreadSharedMemory::FreeSharedWindow(const std::string &window_name)
{
  if (windows_.find(window_name) != windows_.end())
  {
    MPI_Win_free(&(windows_[window_name]));
    windows_.erase(window_name);
  } else
  {
    std::cout << "Window '" << window_name << "' does not exist. Nothing is erased." << std::endl;
  }
}