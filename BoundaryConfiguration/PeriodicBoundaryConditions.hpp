//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP
#define SPC2ODEINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP

#include "../Definitions.hpp"
#include "../Parallelization/Thread.hpp"

class PeriodicBoundaryConditions
{
 public:

  explicit PeriodicBoundaryConditions(Thread *thread, Real x_size, Real y_size);
  ~PeriodicBoundaryConditions();

  void ClassAEffectiveParticleDistance(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);

  void ClassBEffectiveParticleDistanceSigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  void ClassBEffectiveParticleDistanceUnsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);

  void ClassCEffectiveParticleDistanceSigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  void ClassCEffectiveParticleDistanceUnsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);

  void ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc);
  void ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state);
  void ApplyPeriodicBoundaryConditions(Real *const system_state, long size);

  Real GetXSize();
  Real GetYSize();

 private:

  Thread *thread_;
  Real x_size_;
  Real y_size_;
  Real x_rsize_;
  Real y_rsize_;

};

#endif //SPC2ODEINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP
