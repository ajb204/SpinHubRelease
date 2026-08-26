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
//#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdio>
#include "fuda_classes.H"
#include "fuda_utils.H"



/* Data class member functions */

Data::Data(Func *a_func, unsigned int size) : 
  func(a_func), x(size, 0.0), init(0.0), value(0.0), u(1.0) { }

Data::~Data()
{
  // Deallocate x.
  x.clear();
}


Func *Data::get_func() 
{
  return (func);
}


Dtype *Data::get_dtype() 
{
  return (func->get_dtype());
}


unsigned int Data::get_dim()
{
  return (x.size());
}

double Data::get_ix(unsigned int i) 
{
  if (i<0 || i>=x.size()) throw Uferr::DataIndexInvalid(i);
  return (x[i]);
}

void Data::set_ix(unsigned int i, double val)
{
  if (i<0 || i>=x.size()) throw Uferr::DataIndexInvalid(i);
  x[i]=val;
  func->get_fuda()->data_clear_sync();
}

double Data::get_init()
{
  return (init); 
}

void Data::set_init(double a_value)
{ 
  // We both set init and value.
  init=value=a_value;
  func->get_fuda()->data_clear_sync();
}

void Data::init_value()
{ 
  value=init;
  func->get_fuda()->data_clear_sync();
}

void Data::montecarlo_value(double random_factor)
{ 
  value=init+u*random_factor;
  func->get_fuda()->data_clear_sync();
}

double Data::get_value()
{
  return (value); 
}

void Data::set_value(double a_value)
{ 
  value=a_value;
  func->get_fuda()->data_clear_sync();
}

double Data::get_u()
{
  return (u); 
}

void Data::set_u(double a_u)
{
  // The uncertainty must be positive.
  if (a_u<=0.0) throw Uferr::DataUInvalid(a_u);
  u = a_u;
  func->get_fuda()->data_clear_sync();
}

void Data::print()
{
  /* Print data record to standard output */
  for(unsigned int i=0; i<x.size(); i++)
    std::cout << "x[" << i << "] = " << x[i] << "\n";

  // Print value and uncertainty.
  std::cout << "init : " << init << "\n";
  std::cout << "value: " << value << "\n";
  std::cout << "u    : " << u << "\n";
  std::cout << "\n";
}


bool DataPtr_less::operator()(Data *a, Data *b)
{
  int adim = a->x.size();
  int bdim = b->x.size();
  if (adim!=bdim) throw Uferr::DataDimMismatch();

  bool rtn = 0;
  
  // Compare records x[adim-1] first, x[0] last.
  for (int i=adim-1; i>=0; i--)
    if ((a->x[i])<(b->x[i]))
      {
	rtn = 1;
	break;
      }
    else if ((b->x[i])<(a->x[i]))
      {
	rtn = 0;
	break;
      }
  
  return (rtn);
}


