//
// Created by Nikita Kruk on 2019-03-12.
//

#ifndef SPC2ODEINTEGRATION_VORTEXARRAYSSYSTEMPTR_HPP
#define SPC2ODEINTEGRATION_VORTEXARRAYSSYSTEMPTR_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

class VortexArraysSystemPtr
{
 public:

  explicit VortexArraysSystemPtr(ThreadSharedMemory *thread,
                                 Real velocity,
                                 Real mu_plus,
                                 Real mu_minus,
                                 Real xi_a,
                                 Real xi_r,
                                 Real D_phi,
                                 Real kappa,
                                 Real rho,
                                 PeriodicBoundaryConditions &pbc_config);
  ~VortexArraysSystemPtr();

  void EvaluateRhs(const Real *const system_state,
                   const std::vector<Real> &k_prev,
                   std::vector<Real> &k_next,
                   Real k_coef,
                   Real dt);
  void EvaluateRhs(const Real *const system_state,
                   std::vector<Real> &derivative,
                   std::vector<std::vector<Real>> &additional_derivative,
                   Real dt);

  void EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                        const std::vector<Real> &k_prev,
                                        std::vector<Real> &k_next,
                                        Real k_coef,
                                        Real dt);
  void EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                        std::vector<Real> &derivative,
                                        std::vector<std::vector<Real>> &additional_derivative,
                                        Real dt);
  void EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                          const std::vector<Real> &k_prev,
                                          std::vector<Real> &k_next,
                                          Real k_coef,
                                          Real dt);
  void EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                          std::vector<Real> &derivative,
                                          std::vector<std::vector<Real>> &additional_derivative,
                                          Real dt);

  Real D_phi() const
  { return D_phi_; }

  bool last_coefficient_;

 private:

  ThreadSharedMemory *thread_;
  PeriodicBoundaryConditions &pbc_config_;
  Real velocity_;
  Real mu_plus_;
  Real mu_minus_;
  Real xi_a_;
  Real xi_r_;
  Real D_phi_;
  Real kappa_;
  Real rho_;

  Real *rk_system_state_;
  std::vector<Real> alignment_force_;
  std::vector<Real> antialignment_force_;
  std::vector<Real> cohesion_force_;
  std::vector<std::vector<Real>> neighborhood_cardinality_;

  // Cell Subdivision Routine
  Real x_size_;
  Real y_size_;
  int num_subcells_x_;
  int num_subcells_y_;
  int *pre_linked_list_;
  std::vector<std::vector<int>> linked_list_;
  std::vector<std::vector<int>> neighboring_cells_;

  Real AlignmentStrength(Real r);
  Real AlignmentStrengthDerivative(Real r);
  void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y);

};

#endif //SPC2ODEINTEGRATION_VORTEXARRAYSSYSTEMPTR_HPP
