//
// Created by Nikita Kruk on 2019-02-04.
//

#include "../Definitions.hpp"

#if defined(MPI_PARAMETER_SCAN) \
 || defined(MPI_FAST_INTERACTION) || defined(MPI_FAST_INTERACTION_SHARED_MEMORY)
#include <mpi.h>
#endif

#if defined(OPENMP_FAST_INTERACTION)
#include <omp.h>
#endif

#if defined(MPI_PARAMETER_SCAN) \
 || defined(MPI_FAST_INTERACTION) || defined(MPI_FAST_INTERACTION_SHARED_MEMORY)
void LaunchParallelSession(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
}
#else
void LaunchParallelSession(int argc, char **argv)
{

}
#endif

#if defined(MPI_PARAMETER_SCAN) \
 || defined(MPI_FAST_INTERACTION) || defined(MPI_FAST_INTERACTION_SHARED_MEMORY)
void FinalizeParallelSession()
{
  MPI_Finalize();
}
#else
void FinalizeParallelSession()
{

}
#endif