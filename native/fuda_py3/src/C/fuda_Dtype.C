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


/* Dtype class member functions */

Dtype::Dtype(Fuda *a_fuda, std::string& a_name, unsigned int dim) :
  fuda(a_fuda), name(a_name), pvec(dim), purge_flg(0),
  purge_radius(dim) {}

// Destructor.
Dtype::~Dtype()
{
  // Deallocate vectors.
  pvec.clear();
  purge_radius.clear();
  eval.funcs.clear();
  eval.data.clear();
}

Fuda *Dtype::get_fuda() 
{
  return (fuda);
}

void Dtype::get_name(std::string& a_name)
{
  a_name = name;
}

unsigned int Dtype::get_dim()
{
  return (pvec.size());
}

void Dtype::set_ip(unsigned int i, Param *p)
{
  if (i<0 || i>=get_dim()) throw Uferr::DtypeExplIndexInvalid(i);
  pvec[i]=p; 
  fuda->dtype_clear_sync();
}

Param *Dtype::get_ip(unsigned int i)
{
  if (i<0 || i>=get_dim()) throw Uferr::DtypeExplIndexInvalid(i);
  return (pvec[i]);
}

/* Print object contents */
void Dtype::print()
{
  std::cout << "name  : " << name << "\n";
  std::cout << "dim   : " << pvec.size() << "\n";
  for(unsigned int i=0; i<pvec.size(); i++)
    {
      std::string str;
      pvec[i]->get_name(str);
      std::cout << "pvec[" << i << "]  : " << str << "\n";
    }

  unsigned int i=0;
  if (fuda->eval_is_rsync())
    {
      for(std::list<Func*>::iterator fi=eval.funcs.begin();
	  fi!=eval.funcs.end(); fi++)
	{
	  Func& fn = **fi;
	  std::string str;
	  fn.get_name(str);
	  std::cout << "eval.funcs[" << i++ << "] : " << str << "\n";
	}
      std::cout << "eval.data.size()  : " << eval.data.size() << "\n";
    }
  std::cout << "\n";
}


void Dtype::set_purge(bool flg)
{
  purge_flg = flg;
}

void Dtype::set_purge_radius(unsigned int ip, double radius)
{
  if (ip<0 || ip>=get_dim()) throw Uferr::DtypeExplIndexInvalid(ip);
  purge_radius[ip] = radius;
}

bool Dtype::is_purge()
{
  return (purge_flg);
}

double Dtype::get_purge_radius(unsigned int ip)
{
  if (ip<0 || ip>=get_dim()) throw Uferr::DtypeExplIndexInvalid(ip);
  return (purge_radius[ip]);
}


void Dtype::purge()
{
  // Do we purge?
  if (purge_flg)
    { 
      // Proceede if we got more than one point.
      if (eval.data.size()>1) 
	{
	  // Sort the eval.data list.
	  eval.data.sort(DataPtr_less());

	  // Move first data record from eval.data to dl list.
	  std::list<Data*> dl;
	  dl.push_back(eval.data.front());
	  eval.data.pop_front();

	  // Loop over rest of elements and purge while moving.
	  while (eval.data.size()>0)
	    {
	      // Get last data pointer on dl.
	      Data *pre_data = dl.back();
	      
	      // Pop first data pointer from eval.data.
	      Data *this_data = eval.data.front();
	      eval.data.pop_front();
	      
	      // Compare.
	      bool remove_flg = 1;
	      for (unsigned int i=0; i<get_dim(); i++)
		if (fabs(pre_data->x[i] - this_data->x[i]) > purge_radius[i])
		  {
		    remove_flg = 0;
		    break;
		  }
	      
	      // Evt. move record pointer.
	      if (!remove_flg) dl.push_back(this_data);	  
	    }

	  // Finally, move the data back.
	  eval.data.splice(eval.data.end(), dl, dl.begin(), dl.end());
	}      
    }
}


bool Dtype::is_referenced()
{
  // Return true if the Dtype is referenced by at least one Func
  // object.

  bool referenced = 0;
  
  // Loop over all functions.
  for(Func_iterator fi = fuda->func_begin(); fi != fuda->func_end(); fi++)
    {
      // Compare reference to self refererence.
      if (this == (*fi)->get_dtype())
	{
	  referenced = 1;
	  break;
	}
    }
  return (referenced);
}


