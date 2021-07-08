//
// Created by Nikita Kruk on 2019-02-01.
//

#include "PeriodicBoundaryConditions.hpp"

PeriodicBoundaryConditions::PeriodicBoundaryConditions(Thread *thread,
                                                       Real x_size,
                                                       Real y_size) :
    thread_(thread),
    x_size_(x_size),
    y_size_(y_size),
    x_rsize_(1.0 / x_size),
    y_rsize_(1.0 / y_size)
{

}

PeriodicBoundaryConditions::~PeriodicBoundaryConditions() = default;

// periodic signed distance if all interaction sites are in the same simulation box
void PeriodicBoundaryConditions::ClassAEffectiveParticleDistance(Real x_i,
                                                                 Real y_i,
                                                                 Real x_j,
                                                                 Real y_j,
                                                                 Real &dx,
                                                                 Real &dy)
{
  dx = x_j - x_i;
  dx -= static_cast<int>(dx * 2.0 * x_rsize_) * x_size_;

  dy = y_j - y_i;
  dy -= static_cast<int>(dy * 2.0 * y_rsize_) * y_size_;
}

// if the centers or primary sites of all molecules are in the same box, but the other sites can be away from the centers
void PeriodicBoundaryConditions::ClassBEffectiveParticleDistanceSigned(Real x_i,
                                                                       Real y_i,
                                                                       Real x_j,
                                                                       Real y_j,
                                                                       Real &dx,
                                                                       Real &dy)
{
  dx = x_j - x_i;
  if (dx > 0.5 * x_size_)
  {
    dx -= x_size_;
  } else if (dx < -0.5 * x_size_)
  {
    dx += x_size_;
  }

  dy = y_j - y_i;
  if (dy > 0.5 * y_size_)
  {
    dy -= y_size_;
  } else if (dy < -0.5 * y_size_)
  {
    dy += y_size_;
  }
}

// if the centers or primary sites of all molecules are in the same box, but the other sites can be away from the centers
void PeriodicBoundaryConditions::ClassBEffectiveParticleDistanceUnsigned(Real x_i,
                                                                         Real y_i,
                                                                         Real x_j,
                                                                         Real y_j,
                                                                         Real &dx,
                                                                         Real &dy)
{
  dx = std::fabs(x_j - x_i);
  if (dx > 0.5 * x_size_)
  {
    dx -= x_size_;
  }

  dy = std::fabs(y_j - y_i);
  if (dy > 0.5 * y_size_)
  {
    dy -= y_size_;
  }
}

// calculate remaindars at any distance
//if the sign of the distance is relevant
void PeriodicBoundaryConditions::ClassCEffectiveParticleDistanceSigned(Real x_i,
                                                                       Real y_i,
                                                                       Real x_j,
                                                                       Real y_j,
                                                                       Real &dx,
                                                                       Real &dy)
{
  dx = x_j - x_i;
  dx -= x_size_ * std::nearbyint(dx * x_rsize_);

  dy = y_j - y_i;
  dy -= y_size_ * std::nearbyint(dy * y_rsize_);
}

// calculate remaindars at any distance
//if the sign of the distance is not relevant
//#pragma acc routine seq
void PeriodicBoundaryConditions::ClassCEffectiveParticleDistanceUnsigned(Real x_i,
                                                                         Real y_i,
                                                                         Real x_j,
                                                                         Real y_j,
                                                                         Real &dx,
                                                                         Real &dy)
{
  dx = std::fabs(x_j - x_i);
  dx -= static_cast<int>(dx * x_rsize_ + 0.5) * x_size_;

  dy = std::fabs(y_j - y_i);
  dy -= static_cast<int>(dy * y_rsize_ + 0.5) * y_size_;
}

void PeriodicBoundaryConditions::ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc)
{
  x_pbc = x - std::floor(x * x_rsize_) * x_size_;
  y_pbc = y - std::floor(y * y_rsize_) * y_size_;
}

void PeriodicBoundaryConditions::ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state)
{
  for (int i : thread_->GetLoopIndices())
  {
    int ii = i * kS;
    system_state[ii] -= std::floor(system_state[ii] * x_rsize_) * x_size_;
    system_state[ii + 1] -= std::floor(system_state[ii + 1] * y_rsize_) * y_size_;
  } // i
}

void PeriodicBoundaryConditions::ApplyPeriodicBoundaryConditions(Real *const system_state, long size)
{
  for (int i : thread_->GetLoopIndices())
  {
    int ii = i * kS;
    system_state[ii] -= std::floor(system_state[ii] * x_rsize_) * x_size_;
    system_state[ii + 1] -= std::floor(system_state[ii + 1] * y_rsize_) * y_size_;
  } // i
}

Real PeriodicBoundaryConditions::GetXSize() const
{
  return x_size_;
}

Real PeriodicBoundaryConditions::GetYSize() const
{
  return y_size_;
}