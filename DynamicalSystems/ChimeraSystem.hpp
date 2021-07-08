//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_CHIMERASYSTEM_HPP
#define SPC2ODEINTEGRATION_CHIMERASYSTEM_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/Thread.hpp"

class ChimeraSystem
{
 public:

  explicit ChimeraSystem(Thread *thread,
                         Real v_0,
                         Real sigma,
                         Real rho,
                         Real alpha,
                         Real D_phi,
                         Real rho_0_,
                         PeriodicBoundaryConditions &pbc_config);
  ~ChimeraSystem();

  void EvaluateRhs(std::vector<Real> &system_state,
                   const std::vector<Real> &k_prev,
                   std::vector<Real> &k_next,
                   Real k_coef,
                   Real dt);
  void EvaluateRhs(std::vector<Real> &system_state,
                   std::vector<Real> &derivative,
                   std::vector<std::vector<Real>> &additional_derivative,
                   Real dt);

  void EvaluateInteractionsWithAllPairs(std::vector<Real> &system_state,
                                        const std::vector<Real> &k_prev,
                                        std::vector<Real> &k_next,
                                        Real k_coef,
                                        Real dt);
  void EvaluateInteractionsWithAllPairs(std::vector<Real> &system_state,
                                        std::vector<Real> &derivative,
                                        std::vector<std::vector<Real>> &additional_derivative,
                                        Real dt);
  void EvaluateInteractionsWithLinkedList(std::vector<Real> &system_state,
                                          const std::vector<Real> &k_prev,
                                          std::vector<Real> &k_next,
                                          Real k_coef,
                                          Real dt);
  void EvaluateInteractionsWithLinkedList(std::vector<Real> &system_state,
                                          std::vector<Real> &derivative,
                                          std::vector<std::vector<Real>> &additional_derivative,
                                          Real dt);
  void EvaluateInteractionsWithVerletNeighborList(std::vector<Real> &system_state,
                                                  const std::vector<Real> &k_prev,
                                                  std::vector<Real> &k_next,
                                                  Real k_coef,
                                                  Real dt);

  // stochastic part of the Stochastic Euler method
  void AddNoise(const std::vector<Real> &system_state, std::vector<Real> &derivative, Real t) const;

  Real rho() const
  { return rho_; }
  Real D_phi() const
  { return D_phi_; }
  Real alpha() const
  { return alpha_; }

  bool last_coefficient_;

 private:

  Thread *thread_;
  PeriodicBoundaryConditions &pbc_config_;
  Real v_0_;
  Real sigma_;
  Real rho_;
  Real rho_squared_;
  Real alpha_;
  Real D_phi_;
  Real rho_0_;

  std::vector<Real> alignment_force_;
  std::vector<Real> neighborhood_cardinality_;

  // Cell Subdivision Routine
  Real x_size_;
  Real y_size_;
  int num_subcells_x_;
  int num_subcells_y_;
  std::vector<int> pre_linked_list_;
  std::vector<std::vector<int>> linked_list_;
  std::vector<std::vector<int>> neighboring_cells_;

  Real r_max_;
  Real r_min_;
  std::vector<std::vector<int>> verlet_list_;
  Real accumulated_displacement_;
  bool should_update_lists_;

  void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y) const;

};

#endif //SPC2ODEINTEGRATION_CHIMERASYSTEM_HPP
