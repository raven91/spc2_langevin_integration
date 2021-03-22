//
// Created by Nikita Kruk on 22.01.20.
//

#ifndef SPC2ODEINTEGRATION_LINKEDLIST_HPP
#define SPC2ODEINTEGRATION_LINKEDLIST_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/Thread.hpp"

#include <vector>

class LinkedList
{
 public:

  explicit LinkedList(const Thread *const thread,
                      Real x_size,
                      Real y_size,
                      int num_cells_x,
                      int num_cells_y,
                      PeriodicBoundaryConditions &pbc_config);
  ~LinkedList();

  void ConstructLinkedList(const std::vector<Real> &system_state);
  [[nodiscard]] const std::vector<std::vector<int>> &GetLinkedList() const;

 private:

  const Thread *const thread_;
  PeriodicBoundaryConditions &pbc_config_;
  Real x_size_;
  Real y_size_;
  int num_cells_x_;
  int num_cells_y_;
  std::vector<int> pre_linked_list_; // for parallelization of list computation
  std::vector<std::vector<int>> linked_list_;
  std::vector<std::vector<int>> neighboring_cells_;

};

#endif //SPC2ODEINTEGRATION_LINKEDLIST_HPP
