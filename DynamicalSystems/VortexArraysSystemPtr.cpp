//
// Created by Nikita Kruk on 2019-03-12.
//

#include "VortexArraysSystemPtr.hpp"

#include <iostream>

VortexArraysSystemPtr::VortexArraysSystemPtr(ThreadSharedMemory *thread,
                                             Real velocity,
                                             Real mu_plus,
                                             Real mu_minus,
                                             Real xi_a,
                                             Real xi_r,
                                             Real D_phi,
                                             Real kappa,
                                             Real rho,
                                             PeriodicBoundaryConditions &pbc_config) :
    thread_(thread),
    pbc_config_(pbc_config),
    velocity_(velocity),
    mu_plus_(mu_plus),
    mu_minus_(mu_minus),
    xi_a_(xi_a),
    xi_r_(xi_r),
    D_phi_(D_phi),
    kappa_(kappa),
    rho_(rho),
    alignment_force_(kN, 0.0),
    antialignment_force_(kN, 0.0),
    cohesion_force_(kN, 0.0),
    neighborhood_cardinality_(kN, std::vector<Real>(3, 0.0)),
    x_size_(pbc_config.GetXSize()),
    y_size_(pbc_config.GetYSize()),
    num_subcells_x_(int(pbc_config.GetXSize() / rho)),
    num_subcells_y_(int(pbc_config.GetYSize() / rho))
{
  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(kS * kN, rk_system_state_, std::string("rk_system_state_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread_->AllocateSharedWindow(0, rk_system_state_, std::string("rk_system_state_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("rk_system_state_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &rk_system_state_);
  }

  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(kN, pre_linked_list_, std::string("pre_linked_list_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread_->AllocateSharedWindow(0, pre_linked_list_, std::string("pre_linked_list_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("pre_linked_list_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &pre_linked_list_);
  }
  linked_list_ = std::vector<std::vector<int>>(num_subcells_x_ * num_subcells_y_, std::vector<int>());

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

  last_coefficient_ = true;
}

VortexArraysSystemPtr::~VortexArraysSystemPtr()
{
  alignment_force_.clear();
  cohesion_force_.clear();
  neighborhood_cardinality_.clear();
  neighboring_cells_.clear();
  for (int i = 0; i < linked_list_.size(); ++i)
  {
    linked_list_[i].clear();
  }
  linked_list_.clear();

  MPI_Barrier(thread_->GetSharedCommunicator());
  thread_->FreeSharedWindow(std::string("rk_system_state_window"));
  thread_->FreeSharedWindow(std::string("pre_linked_list_window"));
}

void VortexArraysSystemPtr::EvaluateRhs(const Real *const system_state,
                                        const std::vector<Real> &k_prev,
                                        std::vector<Real> &k_next,
                                        Real k_coef,
                                        Real dt)
{
  if (rho_ <= 0.25 * x_size_)
  {
    EvaluateInteractionsWithLinkedList(system_state, k_prev, k_next, k_coef, dt);
  } else
  {
    EvaluateInteractionsWithAllPairs(system_state, k_prev, k_next, k_coef, dt);
  }
}

void VortexArraysSystemPtr::EvaluateRhs(const Real *const system_state,
                                        std::vector<Real> &derivative,
                                        std::vector<std::vector<Real>> &additional_derivative,
                                        Real dt)
{
  if (rho_ <= 0.25 * x_size_)
  {
    EvaluateInteractionsWithLinkedList(system_state, derivative, additional_derivative, dt);
  } else
  {
    EvaluateInteractionsWithAllPairs(system_state, derivative, additional_derivative, dt);
  }
}

void VortexArraysSystemPtr::EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                                             const std::vector<Real> &k_prev,
                                                             std::vector<Real> &k_next,
                                                             Real k_coef,
                                                             Real dt)
{
  MPI_Win win_rk_sys = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_rk_sys);
  for (int i : thread_->GetLoopIndices())
  {
    // Note system_state is intentionally kept unsynchronized
    rk_system_state_[i * kS] = system_state[i * kS] + k_coef * dt * k_prev[i * kS];
    rk_system_state_[i * kS + 1] = system_state[i * kS + 1] + k_coef * dt * k_prev[i * kS + 1];
    rk_system_state_[i * kS + 2] = system_state[i * kS + 2] + k_coef * dt * k_prev[i * kS + 2];
  } // i
  MPI_Win_sync(win_rk_sys);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(rk_system_state_);

  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(antialignment_force_.begin(), antialignment_force_.end(), 0.0);
  std::fill(cohesion_force_.begin(), cohesion_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), std::vector<Real>(3, 0.0));

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr = 0.0;

  MPI_Win_lock_all(0, win_rk_sys);
  for (int i : thread_->GetLoopIndices())
  {
    int ii = kS * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
    phi_i = rk_system_state_[ii + 2];

    // j ~ i
    for (int j = 0; j < kN; ++j)
    {
      int jj = kS * j;
      if (j != i)
      {
        x_j = rk_system_state_[jj];
        y_j = rk_system_state_[jj + 1];
        phi_j = rk_system_state_[jj + 2];

        pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
        dr = std::sqrt(dx * dx + dy * dy);

        if (dr <= xi_a_)
        {
          alignment_force_[i] += AlignmentStrength(dr) * std::sin(phi_j - phi_i);
          ++neighborhood_cardinality_[i][0];
        } else if (dr <= rho_)
        {
          antialignment_force_[i] += AlignmentStrength(dr) * std::sin(phi_j - phi_i);
          ++neighborhood_cardinality_[i][1];
        }
        if (dr <= xi_r_)
        {
          cohesion_force_[i] += -kappa_ * std::sin(std::atan2(dy, dx) - phi_i);
          ++neighborhood_cardinality_[i][2];
        }
      }
    } // j

    k_next[ii] = velocity_ * std::cos(phi_i);
    k_next[ii + 1] = velocity_ * std::sin(phi_i);
//    k_next[ii + 2] = alignment_force_[i] + cohesion_force_[i];
    if ((neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1]) > 0)
    {
      k_next[ii + 2] = (alignment_force_[i] + antialignment_force_[i])
          / (neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1]);
    }
    if (neighborhood_cardinality_[i][2] > 0)
    {
      k_next[ii + 2] += cohesion_force_[i] / neighborhood_cardinality_[i][2];
    }
  } // i
  MPI_Win_unlock_all(win_rk_sys);
}

void VortexArraysSystemPtr::EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                                             std::vector<Real> &derivative,
                                                             std::vector<std::vector<Real>> &additional_derivative,
                                                             Real dt)
{
  MPI_Win win_rk_sys = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_rk_sys);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    rk_system_state_[i * kS] = system_state[i * kS];
    rk_system_state_[i * kS + 1] = system_state[i * kS + 1];
    rk_system_state_[i * kS + 2] = system_state[i * kS + 2];
  } // i
  MPI_Win_sync(win_rk_sys);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(rk_system_state_);

  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(antialignment_force_.begin(), antialignment_force_.end(), 0.0);
  std::fill(cohesion_force_.begin(), cohesion_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), std::vector<Real>(3, 0.0));

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr = 0.0;

  MPI_Win_lock_all(0, win_rk_sys);
  for (int i : loop_indices)
  {
    int ii = kS * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
    phi_i = rk_system_state_[ii + 2];

    // j ~ i
    for (int j = 0; j < kN; ++j)
    {
      int jj = kS * j;
      if (j != i)
      {
        x_j = rk_system_state_[jj];
        y_j = rk_system_state_[jj + 1];
        phi_j = rk_system_state_[jj + 2];

        pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
        dr = std::sqrt(dx * dx + dy * dy);

        if (dr <= xi_a_)
        {
          alignment_force_[i] += AlignmentStrength(dr) * std::sin(phi_j - phi_i);
          ++neighborhood_cardinality_[i][0];
        } else if (dr <= rho_)
        {
          antialignment_force_[i] += AlignmentStrength(dr) * std::sin(phi_j - phi_i);
          ++neighborhood_cardinality_[i][1];
        }
        if (dr <= xi_r_)
        {
          cohesion_force_[i] += -kappa_ * std::sin(std::atan2(dy, dx) - phi_i);
          ++neighborhood_cardinality_[i][2];
        }
      }
    } // j

    derivative[ii] = velocity_ * std::cos(phi_i);
    derivative[ii + 1] = velocity_ * std::sin(phi_i);
//    derivative[ii + 2] = alignment_force_[i] + cohesion_force_[i];
    if ((neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1]) > 0)
    {
      derivative[ii + 2] = (alignment_force_[i] + antialignment_force_[i])
          / (neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1]);
    }
    if (neighborhood_cardinality_[i][2] > 0)
    {
      derivative[ii + 2] += cohesion_force_[i] / neighborhood_cardinality_[i][2];
    }
  } // i

  Real sin_of_phase_difference = 0.0, cos_of_phase_difference = 0.0;
  Real sin_of_positional_phase_difference = 0.0, cos_of_positional_phase_difference = 0.0;
  Real alignment_strength = 0.0, alignment_strength_derivative = 0.0;
  Real norm_1 = 0.0, norm_2 = 0.0;
  thread_->SynchronizeVector(derivative);
  for (int i : loop_indices)
  {
    int ii = kS * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
    phi_i = rk_system_state_[ii + 2];

    // j ~ i
    for (int j = 0; j < kN; ++j)
    {
      int jj = kS * j;
      if (j != i)
      {
        x_j = rk_system_state_[jj];
        y_j = rk_system_state_[jj + 1];
        phi_j = rk_system_state_[jj + 2];

        pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
        dr = std::sqrt(dx * dx + dy * dy);

        norm_1 = neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1];
        norm_2 = neighborhood_cardinality_[i][2];
        if ((dr <= rho_) && (norm_1 > 0.0))
        {
          cos_of_phase_difference = std::cos(phi_j - phi_i);
          sin_of_phase_difference = std::sin(phi_j - phi_i);
          alignment_strength = AlignmentStrength(dr);
          alignment_strength_derivative = AlignmentStrengthDerivative(dr);
          additional_derivative[i][0] += -alignment_strength * cos_of_phase_difference / norm_1;
          additional_derivative[i][1] += ((derivative[jj] - derivative[ii]) * alignment_strength_derivative * dx / dr
              + (derivative[jj + 1] - derivative[ii + 1]) * alignment_strength_derivative * dy / dr)
              * sin_of_phase_difference / norm_1
              + (derivative[jj + 2] - derivative[ii + 2]) * alignment_strength * cos_of_phase_difference / norm_1
              - 2.0 * D_phi_ * sin_of_phase_difference * alignment_strength / norm_1;
        }
        if ((dr <= xi_r_) && (norm_2 > 0.0))
        {
          cos_of_positional_phase_difference = std::cos(std::atan2(dy, dx) - phi_i);
          sin_of_positional_phase_difference = std::sin(std::atan2(dy, dx) - phi_i);
          additional_derivative[i][0] += kappa_ * cos_of_positional_phase_difference / norm_2;
          additional_derivative[i][1] += -kappa_ * ((derivative[jj] - derivative[ii]) * (-dy) / dr / dr
              + (derivative[jj + 1] - derivative[ii + 1]) * dx / dr / dr - derivative[ii + 2])
              * cos_of_positional_phase_difference / norm_2
              + kappa_ * D_phi_ * sin_of_positional_phase_difference / norm_2;
        }
      }
    } // j
    additional_derivative[i][0] *= std::sqrt(2.0 * D_phi_);
  } // i
  MPI_Win_unlock_all(win_rk_sys);
}

void VortexArraysSystemPtr::EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                                               const std::vector<Real> &k_prev,
                                                               std::vector<Real> &k_next,
                                                               Real k_coef,
                                                               Real dt)
{
  MPI_Win win_rk_sys = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_rk_sys);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    rk_system_state_[i * kS] = system_state[i * kS] + k_coef * dt * k_prev[i * kS];
    rk_system_state_[i * kS + 1] = system_state[i * kS + 1] + k_coef * dt * k_prev[i * kS + 1];
    rk_system_state_[i * kS + 2] = system_state[i * kS + 2] + k_coef * dt * k_prev[i * kS + 2];
  } // i
  MPI_Win_sync(win_rk_sys);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(rk_system_state_);

  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(antialignment_force_.begin(), antialignment_force_.end(), 0.0);
  std::fill(cohesion_force_.begin(), cohesion_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), std::vector<Real>(3, 0.0));

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr = 0.0;

  for (int i = 0; i < linked_list_.size(); ++i)
  {
    linked_list_[i].clear();
  } // i

  // construct a linked list
  MPI_Win win_pre_linked_list = thread_->GetWindow(std::string("pre_linked_list_window"));
  MPI_Win_lock_all(0, win_pre_linked_list);
  MPI_Win_lock_all(0, win_rk_sys);
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int ii = 0;
//	if (first_coefficient_)
//	{
  for (int i : loop_indices)
  {
    ii = kS * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
//		phi_i = rk_system_state_[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    pre_linked_list_[i] = i_cell;
  } // i
  MPI_Win_sync(win_pre_linked_list);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(pre_linked_list_);
  // obtain all the linked_list
  for (int i = 0; i < kN; ++i)
  {
    linked_list_[pre_linked_list_[i]].push_back(i);
  } // i
  MPI_Win_unlock_all(win_pre_linked_list);

  // loop using the linked list
  MPI_Win_lock_all(0, win_rk_sys);
  int j = 0, jj = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  for (int i : loop_indices)
  {
    ii = kS * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
    phi_i = rk_system_state_[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    for (int neighboring_cell = 0; neighboring_cell < neighboring_cells_.size(); ++neighboring_cell)
    {
      j_cell_x = i_cell_x + neighboring_cells_[neighboring_cell][0];
      j_cell_y = i_cell_y + neighboring_cells_[neighboring_cell][1];
      AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
      j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

      for (std::vector<int>::iterator neighbr_iter = linked_list_[j_cell].begin(),
               neighbr_iter_end = linked_list_[j_cell].end(); neighbr_iter != neighbr_iter_end; ++neighbr_iter)
      {
        j = *neighbr_iter;
        jj = kS * j;
        if (j != i)
        {
          x_j = rk_system_state_[jj];
          y_j = rk_system_state_[jj + 1];
          phi_j = rk_system_state_[jj + 2];

          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
          {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
          }

          pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
          dr = std::sqrt(dx * dx + dy * dy);

          if (dr <= xi_a_)
          {
            alignment_force_[i] += AlignmentStrength(dr) * std::sin(phi_j - phi_i);
            ++neighborhood_cardinality_[i][0];
          } else if (dr <= rho_)
          {
            antialignment_force_[i] += AlignmentStrength(dr) * std::sin(phi_j - phi_i);
            ++neighborhood_cardinality_[i][1];
          }
          if (dr <= xi_r_)
          {
            cohesion_force_[i] += -kappa_ * std::sin(std::atan2(dy, dx) - phi_i);
            ++neighborhood_cardinality_[i][2];
          }
        }
      } // neighbr_iter
    } // neighboring_cell

    k_next[ii] = velocity_ * std::cos(phi_i);
    k_next[ii + 1] = velocity_ * std::sin(phi_i);
//    k_next[ii + 2] = alignment_force_[i] + cohesion_force_[i];
    if ((neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1]) > 0)
    {
      k_next[ii + 2] = (alignment_force_[i] + antialignment_force_[i])
          / (neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1]);
    }
    if (neighborhood_cardinality_[i][2] > 0)
    {
      k_next[ii + 2] += cohesion_force_[i] / neighborhood_cardinality_[i][2];
    }
  } // i
  MPI_Win_unlock_all(win_rk_sys);
}

void VortexArraysSystemPtr::EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                                               std::vector<Real> &derivative,
                                                               std::vector<std::vector<Real>> &additional_derivative,
                                                               Real dt)
{
  MPI_Win win_rk_sys = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_rk_sys);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    rk_system_state_[i * kS] = system_state[i * kS];
    rk_system_state_[i * kS + 1] = system_state[i * kS + 1];
    rk_system_state_[i * kS + 2] = system_state[i * kS + 2];
  } // i
  MPI_Win_sync(win_rk_sys);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(rk_system_state_);

  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(antialignment_force_.begin(), antialignment_force_.end(), 0.0);
  std::fill(cohesion_force_.begin(), cohesion_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), std::vector<Real>(3, 0.0));

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr = 0.0;

  for (int i = 0; i < linked_list_.size(); ++i)
  {
    linked_list_[i].clear();
  } // i

  // construct a linked list
  MPI_Win win_pre_linked_list = thread_->GetWindow(std::string("pre_linked_list_window"));
  MPI_Win_lock_all(0, win_pre_linked_list);
  MPI_Win_lock_all(0, win_rk_sys);
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int ii = 0;
//	if (first_coefficient_)
//	{
  for (int i : loop_indices)
  {
    ii = kS * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
//		phi_i = rk_system_state_[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    pre_linked_list_[i] = i_cell;
  } // i
  MPI_Win_sync(win_pre_linked_list);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(pre_linked_list_);
  // obtain all the linked_list
  for (int i = 0; i < kN; ++i)
  {
    linked_list_[pre_linked_list_[i]].push_back(i);
  } // i
  MPI_Win_unlock_all(win_pre_linked_list);

  // loop using the linked list
  MPI_Win_lock_all(0, win_rk_sys);
  int j = 0, jj = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  for (int i : loop_indices)
  {
    ii = kS * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
    phi_i = rk_system_state_[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    for (int neighboring_cell = 0; neighboring_cell < neighboring_cells_.size(); ++neighboring_cell)
    {
      j_cell_x = i_cell_x + neighboring_cells_[neighboring_cell][0];
      j_cell_y = i_cell_y + neighboring_cells_[neighboring_cell][1];
      AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
      j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

      for (std::vector<int>::iterator neighbr_iter = linked_list_[j_cell].begin(),
               neighbr_iter_end = linked_list_[j_cell].end(); neighbr_iter != neighbr_iter_end; ++neighbr_iter)
      {
        j = *neighbr_iter;
        jj = kS * j;
        if (j != i)
        {
          x_j = rk_system_state_[jj];
          y_j = rk_system_state_[jj + 1];
          phi_j = rk_system_state_[jj + 2];

          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
          {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
          }

          pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
          dr = std::sqrt(dx * dx + dy * dy);

          if (dr <= xi_a_)
          {
            alignment_force_[i] += AlignmentStrength(dr) * std::sin(phi_j - phi_i);
            ++neighborhood_cardinality_[i][0];
          } else if (dr <= rho_)
          {
            antialignment_force_[i] += AlignmentStrength(dr) * std::sin(phi_j - phi_i);
            ++neighborhood_cardinality_[i][1];
          }
          if (dr <= xi_r_)
          {
            cohesion_force_[i] += -kappa_ * std::sin(std::atan2(dy, dx) - phi_i);
            ++neighborhood_cardinality_[i][2];
          }
        }
      } // neighbr_iter
    } // neighboring_cell

    derivative[ii] = velocity_ * std::cos(phi_i);
    derivative[ii + 1] = velocity_ * std::sin(phi_i);
    if ((neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1]) > 0)
    {
      derivative[ii + 2] = (alignment_force_[i] + antialignment_force_[i])
          / (neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1]);
    }
    if (neighborhood_cardinality_[i][2] > 0)
    {
      derivative[ii + 2] += cohesion_force_[i] / neighborhood_cardinality_[i][2];
    }
  } // i

  Real sin_of_phase_difference = 0.0, cos_of_phase_difference = 0.0;
  Real sin_of_positional_phase_difference = 0.0, cos_of_positional_phase_difference = 0.0;
  Real alignment_strength = 0.0, alignment_strength_derivative = 0.0;
  Real norm_1 = 0.0, norm_2 = 0.0;
  thread_->SynchronizeVector(derivative);
  for (int i : loop_indices)
  {
    ii = kS * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
    phi_i = rk_system_state_[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    for (int neighboring_cell = 0; neighboring_cell < neighboring_cells_.size(); ++neighboring_cell)
    {
      j_cell_x = i_cell_x + neighboring_cells_[neighboring_cell][0];
      j_cell_y = i_cell_y + neighboring_cells_[neighboring_cell][1];
      AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
      j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

      for (std::vector<int>::iterator neighbr_iter = linked_list_[j_cell].begin(),
               neighbr_iter_end = linked_list_[j_cell].end(); neighbr_iter != neighbr_iter_end; ++neighbr_iter)
      {
        j = *neighbr_iter;
        jj = kS * j;
        if (j != i)
        {
          x_j = rk_system_state_[jj];
          y_j = rk_system_state_[jj + 1];
          phi_j = rk_system_state_[jj + 2];

          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
          {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
          }

          pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
          dr = std::sqrt(dx * dx + dy * dy);

          norm_1 = neighborhood_cardinality_[i][0] + neighborhood_cardinality_[i][1];
          norm_2 = neighborhood_cardinality_[i][2];
          if ((dr <= rho_) && (norm_1 > 0.0))
          {
            cos_of_phase_difference = std::cos(phi_j - phi_i);
            sin_of_phase_difference = std::sin(phi_j - phi_i);
            alignment_strength = AlignmentStrength(dr);
            alignment_strength_derivative = AlignmentStrengthDerivative(dr);
            additional_derivative[i][0] += -alignment_strength * cos_of_phase_difference / norm_1;
            additional_derivative[i][1] += ((derivative[jj] - derivative[ii]) * alignment_strength_derivative * dx / dr
                + (derivative[jj + 1] - derivative[ii + 1]) * alignment_strength_derivative * dy / dr)
                * sin_of_phase_difference / norm_1
                + (derivative[jj + 2] - derivative[ii + 2]) * alignment_strength * cos_of_phase_difference / norm_1
                - 2.0 * D_phi_ * sin_of_phase_difference * alignment_strength / norm_1;
          }
          if ((dr <= xi_r_) && (norm_2 > 0.0))
          {
            cos_of_positional_phase_difference = std::cos(std::atan2(dy, dx) - phi_i);
            sin_of_positional_phase_difference = std::sin(std::atan2(dy, dx) - phi_i);
            additional_derivative[i][0] += kappa_ * cos_of_positional_phase_difference / norm_2;
            additional_derivative[i][1] += -kappa_ * ((derivative[jj] - derivative[ii]) * (-dy) / dr / dr
                + (derivative[jj + 1] - derivative[ii + 1]) * dx / dr / dr - derivative[ii + 2])
                * cos_of_positional_phase_difference / norm_2
                + kappa_ * D_phi_ * sin_of_positional_phase_difference / norm_2;
          }
        }
      } // neighbr_iter
    } // neighboring_cell
    additional_derivative[i][0] *= std::sqrt(2.0 * D_phi_);
  } // i
  MPI_Win_unlock_all(win_rk_sys);
}

Real VortexArraysSystemPtr::AlignmentStrength(Real r)
{
  if (r < xi_a_)
  {
    return mu_plus_ * (1.0 - r * r / (xi_a_ * xi_a_));
  } else if (r < rho_)
  {
    return -mu_minus_ * 4.0 * (r - xi_a_) * (rho_ - r) / (rho_ - xi_a_) / (rho_ - xi_a_);
  } else
  {
    return 0.0;
  }
}

Real VortexArraysSystemPtr::AlignmentStrengthDerivative(Real r)
{
  if (r < xi_a_)
  {
    return mu_plus_ * (-2.0 * r / xi_a_ / xi_a_);
  } else if (r < rho_)
  {
    return -mu_minus_ * 4.0 * (-2.0 * r + rho_ + xi_a_) / (rho_ - xi_a_) / (rho_ - xi_a_);
  } else
  {
    return 0.0;
  }
}

void VortexArraysSystemPtr::AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y)
{
  if (cell_x < 0)
  {
    cell_x = num_subcells_x_ - 1;
  } else if (cell_x >= num_subcells_x_)
  {
    cell_x = 0;
  }

  if (cell_y < 0)
  {
    cell_y = num_subcells_y_ - 1;
  } else if (cell_y >= num_subcells_y_)
  {
    cell_y = 0;
  }
}