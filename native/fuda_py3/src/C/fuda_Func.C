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


/* Func class member functions */

/* Constructor */
Func::Func(Fuda *a_fuda, std::string& a_name,
	   Ftype *a_ftype, Dtype *a_dtype)
  : fuda(a_fuda), name(a_name),
    ftype(a_ftype), dtype(a_dtype)
{
  // Set use flag to true.
  use = 1;

  /* Initialize paramter list with NULLs */
  params.resize(ftype->get_nparam());
  for(unsigned int i=0; i<ftype->get_nparam(); i++)
    params[i] = NULL;
}

// Destructor.
Func::~Func()
{
  // Deallocate vectors.
  params.clear();
  data.clear();
  eval.x.clear();
  eval.p.clear();
  eval.c.clear();
}

Fuda *Func::get_fuda() 
{
  return (fuda);
}

void Func::get_name(std::string& a_name)
{ 
  a_name = name;
}

bool Func::get_use()
{
  return (use);
}

void Func::set_use(bool a_use)
{
  // Set the use flag.
  use=a_use;

  // Clear the func sync flag.
  fuda->func_clear_sync();
}

Param *Func::get_param(unsigned int i)
{
  if (i<0 || i>=params.size()) throw Uferr::ParamIndexInvalid(i);
  return (params[i]);
}

void Func::set_param(unsigned int i, Param *a_param) 
{
  // Check bounds.
  if (i<0 || i>=params.size()) throw Uferr::ParamIndexInvalid(i);

  // If the ftype parameter is not a variable, check that we got a
  // constant parameter.
  if (!ftype->is_p_var(i) && a_param->get_kind()!=FUDA::PKIND_CONST)
    throw Uferr::ParamNotConst();

  // Set param.
  params[i] = a_param;

  // Clear sync flag.
  fuda->func_clear_sync();
}


Param *Func::get_var(unsigned int i)
{
  if (i<0 || i>=ftype->get_nvar()) throw Uferr::VarIndexInvalid(i);
  return (params[ftype->get_var_index(i)]);
}


void Func::set_var(unsigned int i, Param *a_param) 
{
  // Check bounds.
  if (i<0 || i>=ftype->get_nvar()) throw Uferr::VarIndexInvalid(i);

  // Set param.
  params[ftype->get_var_index(i)] = a_param;

  // Clear sync flag.
  fuda->func_clear_sync();
}


Ftype *Func::get_ftype() 
{
  return (ftype);
}


Dtype *Func::get_dtype()
{
  return (dtype); 
}


unsigned int Func::get_nexpl()
{
  return (eval.x.size());
}


unsigned int Func::get_ndata()
{
  return (data.size());
}


Data_iterator Func::data_begin()
{
  return (data.begin());
}


Data_iterator Func::data_end()
{
  return (data.end());
}


// Add data record to function data list.
Data *Func::add_data()
{
  // Create, add and return ref to new data record.
  Data new_data(this, dtype->get_dim());
  data.push_back(new_data);
  Data& da = data.back();
  fuda->func_clear_sync();
  return (&da);
}


// Delete all data records in func record data list.
void Func::delete_data()
{
  // delete all data.
  data.clear();

  // Clear the data sync and rsync flags.
  fuda->data_clear_sync();
  fuda->data_clear_rsync();

  // Clear sync flg.
  fuda->func_clear_sync();
}

// Initialize the function eval.x vector.
void Func::eval_x_init()
{
  // Clear eval vector.
  eval.x.clear();
  
  // Count up number of explanatory parameters and add entries to
  // eval.x vector.
  for (unsigned int i=0; i<ftype->get_nparam(); i++)
    {
      Param *pm = get_param(i);
      FUDA::Param_kind kind = pm->get_kind();
      IntIndex ii;
      if (kind==FUDA::PKIND_EXPL)
	{
          // Loop over dtype explanatory parmeters to find
          // the corresponding parameter.
	  bool found = 0;
          for (unsigned int ep=0; ep<dtype->get_dim(); ep++)
            {
              if (pm==dtype->get_ip(ep))
                {
		  found = 1;		  
		  ii.iv = i;
                  ii.val = ep;
                  break;
                }
            }

	  // Did we find it? - if not we have an internal
	  // inconsistency in uf, because Fuda::func_add and
	  // Fuda::func_set_param are responsible for checking the
	  // validity of the parameters.
	  if (!found)
	    {
	      std::string pname;
	      pm->get_name(pname);
	      throw Uferr::ExplVarInvalid(pname);
	    }
	  
	  // Add to eval.x
	  eval.x.push_back(ii);
	  
	}
      else if (kind==FUDA::PKIND_PARAM)
	{
	}
      else if (kind==FUDA::PKIND_CONST)
	{
	}
      else throw Uferr::ParamKindInvalid();
    }  
}



// Initialize the function eval.p and eval.c vectors.
void Func::eval_pc_init()
{
  // Clear eval vectors.
  eval.p.clear();
  eval.c.clear();
  
  // free parameters and constant-value parameters (either fixed or
  // constant) and add entries to eval vectors.
  for (unsigned int i=0; i<ftype->get_nparam(); i++)
    {
      Param *pm = get_param(i);
      FUDA::Param_kind kind = pm->get_kind();
      IntIndex ii;
      if (kind==FUDA::PKIND_EXPL)
	{
	}
      else if (kind==FUDA::PKIND_PARAM)
	{
	  if (pm->is_free())
	    {
	      
	      // Add to eval.p
	      ii.iv = i;
	      eval.p.push_back(ii);
	      
		    }
	  else
	    {
	      // Add to eval.c
	      ii.iv = i;
	      eval.c.push_back(ii);
	      
	    }
	}      
      else if (kind==FUDA::PKIND_CONST)
	{
	  // Add to eval.c
	  ii.iv = i;
	  eval.c.push_back(ii);
	  
	}
      else throw Uferr::ParamKindInvalid();
    }  
}



// Calculate function value as a function of explanatory parameters in dtype.
double Func::call(std::vector<double>& expl_vec)
{
  // Vector size must match number of explanatory parameters in dtype.
  if (expl_vec.size()!=dtype->get_dim())
    throw Uferr::NumExplInvalid(expl_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(fuda->fcall_get_ref());

  // Setup func_call structure.

  // Set all parameter values including explanatory parameters which
  // have no valid value. The explanatory parameter values are
  // subsequently overwritten with the correct ones.
  for (unsigned int i=0; i<params.size(); i++)
    func_call.p[i] = params[i]->get_value();
  
  // Set (overwrite) all explanatory parameter values with values from
  // argument vector.
  for (unsigned int i=0; i<eval.x.size(); i++)
    func_call.p[eval.x[i].iv] = expl_vec[eval.x[i].val];
  
  // Call the function.
  double value;
  ftype->call(func_call.p,
	      func_call.dp_flg_false,
	      func_call.dp, &value);

  return (value);
}


// Calculate function value as a function of the explanatory
// parameters of the function. Note that this call differs from
// Func::call because a functions explanatory parameters need not
// include all the explanatory parameters of the dtype but may be a
// subset of these and differ in their order.
double Func::call_by_expl(std::vector<double>& expl_vec)
{
  // Vector size must match number of explanatory parameters in function.
  if (expl_vec.size()!=eval.x.size())
    throw Uferr::NumExplInvalid(expl_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(fuda->fcall_get_ref());

  // Setup func_call structure.

  // Set all parameter values including explanatory parameters which
  // have no valid value. The explanatory parameter values are
  // subsequently overwritten with the correct ones.
  for (unsigned int i=0; i<params.size(); i++)
    func_call.p[i] = params[i]->get_value();
  
  // Set (overwrite) all explanatory parameter values with values from
  // argument vector.
  for (unsigned int i=0; i<eval.x.size(); i++)
    func_call.p[eval.x[i].iv] = expl_vec[i];
  
  // Call the function.
  double value;
  ftype->call(func_call.p,
	      func_call.dp_flg_false,
	      func_call.dp, &value);

  return (value);
}


// Calculate function as a function of variable parameters.
double Func::call_by_var(std::vector<double>& v_vec)
{
  // Vector size must match number of variable parameters in function.
  if (v_vec.size()!=ftype->get_nvar())
    throw Uferr::NumVarInvalid(v_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(fuda->fcall_get_ref());

  // Setup func_call structure.

  // Set all constant values.
  for (unsigned int i=0; i<ftype->get_nparam(); i++)
    if (!ftype->is_p_var(i)) 
      func_call.p[i] = params[i]->get_value();

  // Set variable values.
  for (unsigned int i=0; i<ftype->get_nvar(); i++)
    func_call.p[ftype->get_var_index(i)] = v_vec[i];

  // Call the function.
  double value;
  ftype->call(func_call.p,
	      func_call.dp_flg_false,
	      func_call.dp, &value);
  
  return (value);
}


// Calculate function as a function of variable parameters and return
// derivatives in vector argument.
double Func::call_by_var(std::vector<double>& v_vec,
			 std::vector<double>& dv_vec)
{
  // Vector size must match number of variable parameters in function.
  if (v_vec.size()!=ftype->get_nvar())
    throw Uferr::NumVarInvalid(v_vec.size());

  // dv_vec vector size must match number of variable parameters in function.
  if (dv_vec.size()!=ftype->get_nvar())
    throw Uferr::NumVarInvalid(dv_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(fuda->fcall_get_ref());

  // Setup func_call structure.

  // Set all constant values.
  for (unsigned int i=0; i<ftype->get_nparam(); i++)
    if (!ftype->is_p_var(i)) 
      func_call.p[i] = params[i]->get_value();

  // Set variable values.
  for (unsigned int i=0; i<ftype->get_nvar(); i++)
    func_call.p[ftype->get_var_index(i)] = v_vec[i];

  // Set dp_flg flags to true for all vars and false to all consts.
  for (unsigned int i=0; i<ftype->get_nparam(); i++)
    func_call.dp_flg[i] = ftype->is_p_var(i);

  // Call the function.
  double value;
  ftype->call(func_call.p,
	      func_call.dp_flg,
	      func_call.dp, &value);
  
  // Return derivatives.
  for (unsigned int i=0; i<ftype->get_nvar(); i++)
    dv_vec[i] = func_call.dp[ftype->get_var_index(i)];

  return (value);
}


// Calculate function as a function of all parameters.
double Func::call_by_param(std::vector<double>& p_vec)
{
  // Vector size must match number of parameters in function.
  if (p_vec.size()!=ftype->get_nparam())
    throw Uferr::NumParamInvalid(p_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(fuda->fcall_get_ref());

  // Setup func_call structure.

  // Set  values.
  for (unsigned int i=0; i<ftype->get_nparam(); i++)
    func_call.p[i] = p_vec[i];

  // Call the function.
  double value;
  ftype->call(func_call.p,
	      func_call.dp_flg_false,
	      func_call.dp, &value);
  
  return (value);
}


/* Print object contents */
void Func::print()
{
  std::cout << "name  : " << name << "\n";
  std::string str;
  ftype->get_name(str);
  std::cout << "ftype : " << str << "\n";
  dtype->get_name(str);
  std::cout << "dtype : " << str << "\n";
  std::cout << "use   : " << use << "\n";
  std::cout << "params :";
  for(unsigned int i=0; i<params.size(); i++)
    {
      std::string str;
      params[i]->get_name(str);
      std::cout << " " << str;
    }
  std::cout << "\n";

  // Print eval structure.
  std::cout << "eval.x :";
  for (unsigned int i=0; i<eval.x.size(); i++)
    std::cout << " (" << eval.x[i].iv << "," << eval.x[i].val << ")";
  std::cout << "\n";

  std::cout << "eval.p :";
  for (unsigned int i=0; i<eval.p.size(); i++)
    std::cout << " (" << eval.p[i].iv << "," << eval.p[i].val << ")";
  std::cout << "\n";

  std::cout << "eval.c :";
  for (unsigned int i=0; i<eval.c.size(); i++)
    std::cout << " (" << eval.c[i].iv << "," << eval.c[i].val << ")";
  std::cout << "\n";

  std::cout << "ndata : " << get_ndata() << "\n";

  std::cout << "\n";
}


/* Print function data contents */
void Func::print_data()
{
  unsigned int count = 0;
  for(Data_iterator di=data_begin(); di!=data_end(); di++)
    {
      std::cout << "[" << count++ << "]\n";
      di->print();
    }
}

