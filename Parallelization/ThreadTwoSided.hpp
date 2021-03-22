//
// Created by Nikita Kruk on 2019-02-04.
//

#ifndef SPC2ODEINTEGRATION_THREADTWOSIDED_HPP
#define SPC2ODEINTEGRATION_THREADTWOSIDED_HPP

#include "Thread.hpp"

class ThreadTwoSided : public Thread
{
 public:

  ThreadTwoSided(int argc, char **argv);
  ~ThreadTwoSided();

  bool IsRoot() const override;

  void SynchronizeVector(std::vector<Real> &vec) override;
  void SynchronizeVectorThoughBuffer(std::vector<Real> &vec, std::vector<Real> &buf) override;
  void BroadcastVector(std::vector<Real> &vec) override;
  void SynchronizePrelinkedList(std::vector<int> &pre_linked_list) const override;

};

#endif //SPC2ODEINTEGRATION_THREADTWOSIDED_HPP
