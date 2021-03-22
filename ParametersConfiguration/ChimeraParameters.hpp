//
// Created by Nikita Kruk on 2019-03-12.
//

#ifndef SPC2ODEINTEGRATION_CHIMERAPARAMETERS_HPP
#define SPC2ODEINTEGRATION_CHIMERAPARAMETERS_HPP

#include "../Definitions.hpp"

class ChimeraParameters
{
 public:

  ChimeraParameters();
  ~ChimeraParameters();

 private:

  Real sigma_;
  Real rho_;
  Real alpha_;
  Real D_phi_;

};

#endif //SPC2ODEINTEGRATION_CHIMERAPARAMETERS_HPP
