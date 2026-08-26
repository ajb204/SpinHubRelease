/*
 * Abort.cpp
 *
 *  Created on: Jun 19, 2009
 *      Author: guillaume
 */

#include <standard.h>
#include <Abort.h>

#ifdef DFH_MPI
#include <mpi.h>
#endif

void Abort(int errn) {
  #ifdef _MPI_H
    MPI::Finalize();
  #endif
  std::cerr<<" Abort has been called with error status: "<<errn<<std::endl;
  exit(errn);
  return;
};
