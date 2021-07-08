//
// Created by Nikita Kruk on 2019-02-01.
//

#ifndef SPC2ODEINTEGRATION_DEFINITIONS_HPP
#define SPC2ODEINTEGRATION_DEFINITIONS_HPP

//#define LICHTENBERG
//#define BCS_CLUSTER

//#define MPI_PARAMETER_SCAN
//#define MPI_FAST_INTERACTION
#define MPI_FAST_INTERACTION_SHARED_MEMORY

#include <random>

typedef double Real;
typedef float RealOutput;

const int kN = 8192;
const int kS = 3;
const Real kL = 64.0;

extern std::mt19937 mersenne_twister_generator;

#endif //SPC2ODEINTEGRATION_DEFINITIONS_HPP
