/*Legal_notice:

Copyright 1996-2001 by Soren M. Kristensen, Department of Chemistry,
University of Copenhagen, DK-2100 Copenhagen O, Denmark.  All rights
reserved.

This software and its related software and documentation is provided 'as
is' without express or implied warrenty. The Departmet of Chemistry, The
University of Copenhagen and the author make no warrenties as to any matter
whatsoever with respect to the program and the related software and
documentation. In particular, any and all warranties of merchantability and
fitness for any particular purpose are expressly excluded.

Legal_notice*/
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdio>
#include "fuda_classes.H"
#include "fuda_utils.H"


///////////////////////////////////////////////////////////////
// FuncCall member functions.
///////////////////////////////////////////////////////////////

// Constructor for FuncCall data structure.
FuncCall::FuncCall()
{
  // Initialization.
  defined = 0; 
  p = NULL;
  dp = NULL;
  dp_flg = NULL;
  dp_flg_false = NULL;
  psize = 0;
}


// Destructor.
FuncCall::~FuncCall()
{
  if (defined and psize < 0)
    {
      delete[] p;
      delete[] dp;
      delete[] dp_flg;
      delete[] dp_flg_false;      
    }
}



// Set new allocation size.
void FuncCall::set_size(unsigned int size)
{
  if (size<1) throw Uferr::FuncCallInvalSize(size);

  // We only resize if size>psize.
  if (size>psize)
    {
      // If already allocated, we deallocate old data. 
      if (defined)
	{
	  delete[] p;
	  delete[] dp;
	  delete[] dp_flg;
	  delete[] dp_flg_false;
	}
      
      // Allocate new space.
      p = new double[size];
      dp = new double[size];
      dp_flg = new int[size];
      dp_flg_false = new int[size];
      
      // Initialize dp_flg_false with zeros.
      for(unsigned int i=0; i<size; i++) dp_flg_false[i] = 0;

      // Set new size.
      psize = size;

      // Set defined flag.
      defined = 1;
    }
}


// Return the allocation size.
unsigned int FuncCall::get_size()
{
  return (psize);
}



