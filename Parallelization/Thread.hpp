//
// Created by Nikita Kruk on 2019-02-04.
//

#ifndef SPC2ODEINTEGRATION_THREAD_HPP
#define SPC2ODEINTEGRATION_THREAD_HPP

#include "../Definitions.hpp"

#include <vector>

class Thread
{
 public:

  Thread(int argc, char **argv);
  virtual ~Thread();

  int GetRank();
  int GetNumberOfMpichThreads();
  int GetNumberOfParticlesPerMpichThread();
  [[nodiscard]] const std::vector<int> &GetLoopIndices() const;

  virtual bool IsRoot() const;

  virtual void SynchronizeVector(std::vector<Real> &vec);
  virtual void SynchronizeVectorThroughBuffer(std::vector<Real> &vec, std::vector<Real> &buf);
  virtual void BroadcastVector(std::vector<Real> &vec);
  virtual void SynchronizePrelinkedList(std::vector<int> &pre_linked_list) const;

 protected:

  int root_rank_;
  int rank_;
  int number_of_mpich_threads_;
  int number_of_particles_per_mpich_thread_;
  std::vector<int> loop_indices_;

};

#endif //SPC2ODEINTEGRATION_THREAD_HPP
