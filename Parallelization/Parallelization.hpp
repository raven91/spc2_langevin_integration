//
// Created by Nikita Kruk on 2019-02-04.
//

#ifndef SPC2ODEINTEGRATION_PARALLELIZATION_HPP
#define SPC2ODEINTEGRATION_PARALLELIZATION_HPP

#include "../Definitions.hpp"

void LaunchParallelSession(int argc, char **argv);
void FinalizeParallelSession();

#endif //SPC2ODEINTEGRATION_PARALLELIZATION_HPP
