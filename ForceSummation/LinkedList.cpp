//
// Created by Nikita Kruk on 22.01.20.
//

#include "LinkedList.hpp"

LinkedList::LinkedList(const Thread *const thread,
                       Real x_size,
                       Real y_size,
                       int num_cells_x,
                       int num_cells_y,
                       PeriodicBoundaryConditions &pbc_config) :
    thread_(thread),
    x_size_(x_size),
    y_size_(y_size),
    num_cells_x_(num_cells_x),
    num_cells_y_(num_cells_y),
    pbc_config_(pbc_config)
{
  pre_linked_list_ = std::vector<int>(kN, 0);
  linked_list_ = std::vector<std::vector<int>>(num_cells_x_ * num_cells_y_, std::vector<int>());

  //all neighbors
  neighboring_cells_ =
      {
          {-1, -1}, {0, -1}, {1, -1},
          {-1, 0}, {0, 0}, {1, 0},
          {-1, 1}, {0, 1}, {1, 1}
      };
  //half of neighbors
//	neighboring_cells_ =
//	{
//			{0, 0}, {1, 0},
//			{-1, 1}, {0, 1}, {1, 1}
//	};
}

LinkedList::~LinkedList()
{
  pre_linked_list_.clear();
  for (std::vector<int> &linked_list : linked_list_)
  {
    linked_list.clear();
  } // linked_list
  linked_list_.clear();
  neighboring_cells_.clear();
}

void LinkedList::ConstructLinkedList(const std::vector<Real> &system_state)
{
  for (std::vector<int> &ll : linked_list_)
  {
    ll.clear();
  } // ll

  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int ii = 0;
  Real x_i = 0.0, y_i = 0.0;
//	if (first_coefficient_)
//	{
//#if defined(MPI_FAST_INTERACTION)
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    ii = kS * i;
    x_i = system_state[ii];
    y_i = system_state[ii + 1];
    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_cells_x_);
    i_cell_y = int(y_i / y_size_ * num_cells_y_);
    i_cell = i_cell_y * num_cells_x_ + i_cell_x;
    pre_linked_list_[i] = i_cell;
  } // i

//  if (mpi_thread_rank != mpi_root_rank)
//  {
//    MPI_Gather(&pre_linked_list_[mpi_thread_rank * mpi_num_particles_per_thread],
//               mpi_num_particles_per_thread,
//               MPI_INT,
//               &pre_linked_list_[0],
//               mpi_num_particles_per_thread,
//               MPI_INT,
//               mpi_root_rank,
//               MPI_COMM_WORLD);
//  } else
//  {
//    MPI_Gather(MPI_IN_PLACE,
//               mpi_num_particles_per_thread,
//               MPI_INT,
//               &pre_linked_list_[mpi_root_rank * mpi_num_particles_per_thread],
//               mpi_num_particles_per_thread,
//               MPI_INT,
//               mpi_root_rank,
//               MPI_COMM_WORLD);
//  }
//  MPI_Bcast(&pre_linked_list_[0], (int) pre_linked_list_.size(), MPI_INT, mpi_root_rank, MPI_COMM_WORLD);
  thread_->SynchronizePrelinkedList(pre_linked_list_);

  for (int i = 0; i < kN; ++i)
  {
    linked_list_[pre_linked_list_[i]].push_back(i);
  } // i
//#else
//  for (int i = 0; i < kN; ++i)
//    {
//        ii = kS * i;
//
//        x_i   = system_state[ii    ] + k_coef * dt * k_prev[ii    ];
//        y_i   = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
////		phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
//
//        if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
//        {
//            pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
//        }
//
//        i_cell_x = int(x_i / x_size_ * num_subcells_x_);
//        i_cell_y = int(y_i / y_size_ * num_subcells_y_);
//        i_cell = i_cell_y * num_subcells_x_ + i_cell_x;
//
//        linked_list_[i_cell].push_back(i);
//    }
//#endif
//	}
}

const std::vector<std::vector<int> > &LinkedList::GetLinkedList() const
{
  return linked_list_;
}