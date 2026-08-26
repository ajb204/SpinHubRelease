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

#include <cstring>
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
// Fuda member functions.
///////////////////////////////////////////////////////////////

// Constructor for Fuda data structure.
Fuda::Fuda()
{
  // Workspace initialization.
  wspace = NULL;
  wsize = 0;
  wlocked = 0;

  // List sync flag initialization.
  param_sync = 0;
  param_rsync = 0;
  dtype_sync = 0;
  dtype_rsync = 0;
  func_sync = 0;
  func_rsync = 0;
  data_sync = 0;
  data_rsync = 0;
  
  // Eval initialization.
  eval_set_mod_count(0); 
  eval_clear();

  // fcall initialization.
  fcall = new FuncCall();
}


// Return pointer to func call data structure.
FuncCall *Fuda::fcall_get_ref()
{
  return (fcall);
}


// Set eval modification counter.
void Fuda::eval_set_mod_count(unsigned int value)
{
  eval.mod_count = value;
}


// Increment eval modification counter.
void Fuda::eval_incr_mod_count()
{
  eval.mod_count++;
}


// Return eval modification counter.
unsigned int Fuda::eval_get_mod_count()
{
  return (eval.mod_count);
}


// Deallocate wspace.
void Fuda::wspace_deallocate()
{
  if (wspace!=0) delete[] wspace;
  wsize = 0;
  wlocked = 0;
}


// Get wspace size.
unsigned int Fuda::wspace_size()
{
  return (wsize);
}


// Lock wspace, and return true if successful.
bool Fuda::wspace_lock()
{
  if (wlocked) return (0);
  else
    {
      wlocked = 1;
      return (1);
    }
}


// Unlock wspace.
void Fuda::wspace_unlock()
{
  wlocked = 0;
}


// Get wspace ref to workspace at least of size minsize.
double *Fuda::wspace_ref(unsigned int minsize)
{
  wspace_allocate(minsize);
  return (wspace);
}


// Allocate wspace.
void Fuda::wspace_allocate(unsigned int minsize)
{
  if (minsize>wsize)
    {
      // Eventually deallocate previous wspace.
      if (wspace!=NULL) delete[] wspace;

      // Allocate new wspace.
      wsize = minsize;
      wspace = new double[minsize];
    }
}


// Ptype iterator begin and end functions.
Ptype_iterator Fuda::ptype_begin() { return (ptypes.begin()); }
Ptype_iterator Fuda::ptype_end() { return (ptypes.end()); }
unsigned int Fuda::get_nptype() { return (ptypes.size()); }


// Add a ptype record.
Ptype *Fuda::ptype_add(std::string name, FUDA::Param_kind a_kind)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::PtypeNameInvalid(name);

  /* Is the name free? */
  if (ptype_find(name) != NULL) throw Uferr::PtypeNameAlreadyUsed(name);

  /* Add link */
  Ptype *ptype = new Ptype(this, name, a_kind);
  ptypes.push_back(ptype);

  /* Return ref to added link */
  return (ptype);
}


// Delete a Ptype object in the ptypes list.
void Fuda::ptype_del(std::string& name)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::PtypeNameInvalid(name);

  // Find the ptype by looping over ptypes list.
  Ptype_iterator pti;
  bool ptiempty = true;
  std::string iname;
  for(Ptype_iterator i = ptype_begin(); i != ptype_end(); i++)
    {
      // Compare names and if they match, we set pti and break out.
      (*i)->get_name(iname);
      if (iname == name) 
	{
	  pti = i;
    ptiempty = false;
	  break;
	}
    }

  //  Did we find the ptype?
  if (ptiempty == true) throw Uferr::PtypeNameInvalid(name);

  // Is the ptype referenced by any parameter?
  if ((*pti)->is_referenced()) throw Uferr::PtypeReferencedByParam(name);

  // First destroy the ptype object.
  delete *pti;

  // Then remove the empty pointer from the params list.
  ptypes.erase(pti);
}


// Delete a Ptype object in the ptypes list.
void Fuda::ptype_del(const char name[])
{
  std::string nm = name;
  ptype_del(nm);
}


// Delete all declared ptypes.
void Fuda::ptype_del_all()
{
  /* No ptype may be referenced by a parameter. This strange
     construction of looking up the first ptype referenced by the
     first parameter is just to get the name of a ptype to return
     when throwing the exception. */
  if (get_nparam()>0)
    {
      // Loop over parameters.
      for(Param_iterator pi = param_begin(); pi != param_end(); pi++)
	if ((*pi)->get_ptype() != NULL)
	  {
	    // Get name of first ptype referenced by param obj.
	    std::string ptname;
	    (*pi)->get_ptype()->get_name(ptname);

	    // Finally throw exception with this ptype name.
	    throw Uferr::PtypeReferencedByParam(ptname);
	  }
    }
  
  /* Loop over ptypes list entries */
  for(Ptype_iterator i = ptype_begin(); i != ptype_end(); i++)
    {
      // Destroy the Ptype object.
      delete *i;
    }
  
  // Then remove all the empty pointers from the ptypes list.
  ptypes.clear();
}


// Find a Ptype with a specified name.
Ptype *Fuda::ptype_find(std::string& name)
{
  Ptype *ptype = NULL;

  /* Loop over list entries */
  for(Ptype_iterator i = ptypes.begin(); i != ptypes.end(); i++)
    {
      Ptype& pt = **i;
      if (pt.name == name) 
	{
	  ptype = &pt;
	  break;
	}
    }
  return (ptype);
}


// Find a ptype with a specified name.
Ptype *Fuda::ptype_find(const char name[])
{
  std::string nm = name;
  return (ptype_find(nm));
}


// Param iterator begin and end functions.
Param_iterator Fuda::param_begin() { return (params.begin()); }
Param_iterator Fuda::param_end() { return (params.end()); }
unsigned int Fuda::get_nparam() { return (params.size()); }
void Fuda::param_set_sync() { param_sync = 1; }
void Fuda::param_set_rsync() { param_rsync = 1; }
void Fuda::param_clear_sync() { param_sync = 0; }
void Fuda::param_clear_rsync() { param_rsync = 0; }
bool Fuda::param_is_sync() { return (param_sync); }
bool Fuda::param_is_rsync() { return (param_rsync); }


// Add a parameter with name and ptype.
Param *Fuda::param_add(std::string name,
			std::string ptype_name)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::ParamNameInvalid(name);

  /* Is the name free? */
  if (param_find(name) != NULL) throw Uferr::ParamNameAlreadyUsed(name);

  /* Is the ptype ok? */
  Ptype *ptype = ptype_find(ptype_name);
  if (ptype == NULL) throw Uferr::PtypeNameInvalid(name);

  /* Add link */
  Param *param = new Param(this, name, ptype);
  params.push_back(param);

  /* Return iterator to added link */
  return (param);
}


// Add a parameter with name and parameter kind.
Param *Fuda::param_add(std::string name, FUDA::Param_kind a_kind)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::ParamNameInvalid(name);

  /* Is the name free? */
  if (param_find(name) != NULL) throw Uferr::ParamNameAlreadyUsed(name);

  /* Add link */
  Param *param = new Param(this, name, a_kind);
  params.push_back(param);

  /* Return iterator to added link */
  return (param);
}


// Delete a Param object in the params list.
void Fuda::param_del(std::string& name)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::ParamNameInvalid(name);

  // Find the parameter by looping over params list.
  Param_iterator pi;
  bool piempty = true;
  std::string iname;
  for(Param_iterator i = param_begin(); i != param_end(); i++)
    {
      // Compare names and if they match, we set pi and break out.
      (*i)->get_name(iname);
      if (iname == name) 
	{
	  pi = i;
    piempty=false;
	  break;
	}
    }

  //  Did we find the parameter?
  if (piempty == true) throw Uferr::ParamNameInvalid(name);

  // Is the parameter referenced by any function?
  if ((*pi)->is_referenced()) throw Uferr::ParamReferencedByFunc(name);

  // First destroy the param object.
  delete *pi;

  // Then remove the empty pointer from the params list.
  params.erase(pi);

  // Clear the sync flags.
  param_clear_sync();
  param_clear_rsync();  
}


// Delete a Param object in the params list.
void Fuda::param_del(const char name[])
{
  std::string nm = name;
  param_del(nm);
}


// Delete all declared parameters.
void Fuda::param_del_all()
{
  /* No parameters may be referenced by a function. This strange
     construction of looking up the first parameter referenced by the
     first function is just to get the name of a parameter to return
     when throwing the exception. */
  if (get_nfunc()>0)
    {
      // Loop over functions.
      for(Func_iterator fi = func_begin(); fi != func_end(); fi++)
	if ((*fi)->get_ftype()->get_nparam()>0)
	  {
	    // Get parameter name of first paramter referenced by func obj.
	    std::string pname;
	    (*fi)->get_param(0)->get_name(pname);

	    // Finally throw exception with this param name.
	    throw Uferr::ParamReferencedByFunc(pname);
	  }
    }
  
  /* Loop over params list entries */
  for(Param_iterator i = param_begin(); i != param_end(); i++)
    {
      // Destroy the Param object.
      delete *i;
    }
  
  // Then remove all the empty pointers from the params list.
  params.clear();

  // Clear sync flags.
  param_clear_sync();
  param_clear_rsync();
}


// Find a parameter with a name.
Param *Fuda::param_find(std::string& name)
{
  Param *param = NULL;

  /* Loop over list entries */
  for(Param_iterator i = param_begin(); i != param_end(); i++)
    {
      Param& pm = **i;
      if (pm.name == name) 
	{
	  param = &pm;
	  break;
	}
    }
  return (param);
}


// find parameter with a name.
Param *Fuda::param_find(const char name[])
{
  std::string nm = name;
  return (param_find(nm));
}


// Dtype iterator begin and end functions.
Dtype_iterator Fuda::dtype_begin() { return (dtypes.begin()); }
Dtype_iterator Fuda::dtype_end() { return (dtypes.end()); }
unsigned int Fuda::get_ndtype() { return (dtypes.size()); }
void Fuda::dtype_set_sync() { dtype_sync = 1; }
void Fuda::dtype_set_rsync() { dtype_rsync = 1; }
void Fuda::dtype_clear_sync() { dtype_sync = 0; }
void Fuda::dtype_clear_rsync() { dtype_rsync = 0; }
bool Fuda::dtype_is_sync() { return (dtype_sync); }
bool Fuda::dtype_is_rsync() { return (dtype_rsync); }


// Add a Dtype record.
Dtype *Fuda::dtype_add(std::string name, std::list<std::string>& p_names)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::DtypeNameInvalid(name);

  /* Is the name free? */
  if (dtype_find(name) != NULL) throw Uferr::DtypeNameAlreadyUsed(name);

  // Get dim.
  unsigned int dim = p_names.size();

  /* Make a dtype record */
  Dtype *dtype = new Dtype(this, name, dim);

  /* Check paramters and set parameter references */
  unsigned int p_count = 0;
  for(std::list<std::string>::iterator pi=p_names.begin(); pi!=p_names.end(); pi++)
    {
      Param *p_ref = param_find(*pi);
      if (p_ref != NULL)
	{
	  // This must be an explanatory variable.
	  if (p_ref->get_kind()!=FUDA::PKIND_EXPL)
	    {
	      // Deallocate dtype and throw exception.
	      delete dtype;
	      throw Uferr::ParamNotExplanatory(*pi);
	    }
	  
	  dtype->set_ip((p_count++), p_ref);
	}
      else
	{
	  // Deallocate dtype and throw exception.
	  delete dtype;
	  throw Uferr::ParamNameInvalid(*pi);
	}
      
    }

  /* Add record to list */
  dtypes.push_back(dtype);

  /* Return ref to added link */
  return (dtype);
}


// Delete a Dtype object in the dtypes list.
void Fuda::dtype_del(std::string& name)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::DtypeNameInvalid(name);

  // Find the dtype by looping over dtypes list.
  Dtype_iterator dti;
  bool dtiempty = true;
  std::string iname;
  for(Dtype_iterator i = dtype_begin(); i != dtype_end(); i++)
    {
      // Compare names and if they match, we set dti and break out.
      (*i)->get_name(iname);
      if (iname == name) 
	{
	  dti = i;
    dtiempty = false;
	  break;
	}
    }

  //  Did we find the dtype?
  if (dtiempty == true) throw Uferr::DtypeNameInvalid(name);

  // Is the dtype referenced by any Func object?
  if ((*dti)->is_referenced()) throw Uferr::DtypeReferencedByFunc(name);

  // First destroy the dtype object.
  delete *dti;

  // Then remove the empty pointer from the dtypes list.
  dtypes.erase(dti);

  // Clear the sync flags.
  dtype_clear_sync();
  dtype_clear_rsync();
}


// Delete a Dtype object in the ptypes list.
void Fuda::dtype_del(const char name[])
{
  std::string nm = name;
  dtype_del(nm);
}


// Delete all declared dtypes.
void Fuda::dtype_del_all()
{
  /* No dtype may be referenced by a Func object. This strange
     construction of looking up the first dtype referenced by the
     first function is just to get the name of a dtype to return
     when throwing the exception. */
  if (get_nfunc()>0)
    {
      // Get dtype name of first function.
      std::string dtname;
      Func* fn = *func_begin();
      fn->get_dtype()->get_name(dtname);

      // Finally throw exception with this ptype name.
      throw Uferr::DtypeReferencedByFunc(dtname);
    }
  
  /* Loop over dtypes list entries */
  for(Dtype_iterator i = dtype_begin(); i != dtype_end(); i++)
    {
      // Destroy the Dtype object.
      delete *i;
    }
  
  // Then remove all the empty pointers from the dtypes list.
  dtypes.clear();
}


// Find dtype by name. If not found return NULL.
Dtype *Fuda::dtype_find(std::string& name)
{
  Dtype *dtype = NULL;

  /* Loop over list entries */
  for(Dtype_iterator i = dtypes.begin(); i != dtypes.end(); i++)
    {
      Dtype& dt = **i;
      if (dt.name == name) 
	{
	  dtype = &dt;
	  break;
	}
    }
  return (dtype);
}


// Find Dtype by name.
Dtype *Fuda::dtype_find(const char name[])
{
  std::string nm = name;
  return (dtype_find(nm));
}


// Data sync and rsysnc set/get/clear procedures.
void Fuda::data_set_sync() { data_sync = 1; }
void Fuda::data_set_rsync() { data_rsync = 1; }
void Fuda::data_clear_sync() { data_sync = 0; }
void Fuda::data_clear_rsync() { data_rsync = 0; }
bool Fuda::data_is_sync() { return (data_sync); }
bool Fuda::data_is_rsync() { return (data_rsync); }


// Ftype iterator begin and end functions.
Ftype_iterator Fuda::ftype_begin() { return (ftypes.begin()); }
Ftype_iterator Fuda::ftype_end() { return (ftypes.end()); }
unsigned int Fuda::get_nftype() { return (ftypes.size()); }


// Add preconfigured ftype record.
Ftype *Fuda::ftype_add(Ftype *ftype_ptr)
{
  // Check that the fuda referference in the ftype matches.
  if (ftype_ptr->get_fuda() != this) throw Uferr::FtypeFudaRefMismatch();

  // Get the ftype name.
  std::string name;
  ftype_ptr->get_name(name);
  
  // Is the name valid?
  if (name.size()<1) throw Uferr::FtypeNameInvalid(name);

  // Is the name free or already in use?
  if (ftype_find(name) != NULL) throw Uferr::FtypeNameAlreadyUsed(name);

  // Get the number of parameters and variables in ftype.
  unsigned int nparam = ftype_ptr->get_nparam();
  unsigned int nvar = ftype_ptr->get_nvar();

  // Check param. Must be at least one.
  if (nparam < 1) throw Uferr::NumParamInvalid(nparam);

  // Check nvar. Must be at least one and not larger than nparam.
  if (nvar < 1 || nvar > nparam) throw Uferr::NumVarInvalid(nvar);

  /* Update fcall record which must be sized according to the number
     of parameters of the ftype with the largest number of parameters */
  fcall->set_size(nparam);

  /* Add link reference */
  ftypes.push_back(ftype_ptr);

  /* Return ref to added link */
  return (ftype_ptr);
}


// Add ftype record with auto generation of p_name, p_descr, p_deriv
// and p_var.
Ftype *Fuda::ftype_add_cfunc(std::string name,
			      unsigned int nparam,
			      FUDA::CFUNC func,
			      void *fs,
			      unsigned int fs_size)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::FtypeNameInvalid(name);

  /* Is the name free? */
  if (ftype_find(name) != NULL) throw Uferr::FtypeNameAlreadyUsed(name);

  // We set number of variable parameters to the number of parameters.
  unsigned int nvar = nparam;

  // Check nvar, must be at least one.
  if (nvar < 1) throw Uferr::NumVarInvalid(nvar);

  /* Update fcall record which must be sized according to the number
     of parameters of the ftype with the largest number of parameters */
  fcall->set_size(nparam);

  // Create ftype record and fill it.
  FtypeCFUNC *ftype = new FtypeCFUNC(this, name, nparam, func, fs, fs_size);

  /* Add link reference */
  ftypes.push_back(ftype);

  /* Return ref to added link */
  return (ftype);
}


// Add product ftype.
Ftype *Fuda::ftype_add_prod(std::string name,
			     bool a_scale_flg,
			     std::vector<Ftype*>& ftvec)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::FtypeNameInvalid(name);

  /* Is the name free? */
  if (ftype_find(name) != NULL) throw Uferr::FtypeNameAlreadyUsed(name);

  // We need at least one ftype in ftvec to make a product.
  if (ftvec.size()<1) throw Uferr::FtypeEmptyFtypeVec();

  // Loop over Ftypes and check validity.
  for (unsigned int i=0; i<ftvec.size(); i++)
    if (ftvec[i]==0 || ftvec[i]->get_fuda()!=this)
      throw Uferr::FtypeRefInvalid();
    
  // Create the product ftype structure.
  FtypePROD *ftype = new FtypePROD(this, name, a_scale_flg, ftvec);

  /* Add link reference */
  ftypes.push_back(ftype);

  // Here goes some post processing.

  // todo: to avoid memory lekage in case of an exception in the next
  // section, all the rest should be in a try{} with deallocation of
  // ftype in case of exception.
  

  /* Update fcall record which must be sized according to the number
     of parameters of the ftype with the largest number of parameters */
  fcall->set_size(ftype->get_nparam());

  return (ftype);
}


// Add sum ftype.
Ftype *Fuda::ftype_add_sum(std::string name,
			   bool a_scale_flg,
			   std::vector<Ftype*>& ftvec)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::FtypeNameInvalid(name);

  /* Is the name free? */
  if (ftype_find(name) != NULL) throw Uferr::FtypeNameAlreadyUsed(name);

  // We need at least one ftype in ftvec to make a product.
  if (ftvec.size()<1) throw Uferr::FtypeEmptyFtypeVec();

  // Loop over Ftypes and check validity.
  for (unsigned int i=0; i<ftvec.size(); i++)
    if (ftvec[i]==0 || ftvec[i]->get_fuda()!=this)
      throw Uferr::FtypeRefInvalid();
    
  // Create the sum ftype structure.
  FtypeSUM *ftype = new FtypeSUM(this, name, a_scale_flg, ftvec);

  /* Add link reference */
  ftypes.push_back(ftype);

  // Here goes some post processing.

  // todo: to avoid memory lekage in case of an exception in the next
  // section, all the rest should be in a try{} with deallocation of
  // ftype in case of exception.
  

  /* Update fcall record which must be sized according to the number
     of parameters of the ftype with the largest number of parameters */
  fcall->set_size(ftype->get_nparam());

  return (ftype);
}


// Add composite ftype.
Ftype *Fuda::ftype_add_comp(std::string name,
			    bool a_scale_flg,
			    unsigned int p_index,
			    Ftype *f,
			    Ftype *g)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::FtypeNameInvalid(name);

  /* Is the name free? */
  if (ftype_find(name) != NULL) throw Uferr::FtypeNameAlreadyUsed(name);

  // Check validity of ftypes.
  if (f==0 || f->get_fuda()!=this)
    throw Uferr::FtypeRefInvalid();
  if (g==0 || g->get_fuda()!=this)
    throw Uferr::FtypeRefInvalid();

  // Create the composite ftype structure.
  FtypeCOMP *ftype = new FtypeCOMP(this, name, a_scale_flg, p_index, f, g);

  /* Add link reference */
  ftypes.push_back(ftype);

  /* Update fcall record which must be sized according to the number
     of parameters of the ftype with the largest number of parameters */
  fcall->set_size(ftype->get_nparam());

  return (ftype);
}


// Find ftype by name. Return NULL if not found.
Ftype *Fuda::ftype_find(std::string& name)
{
  Ftype *ftype = NULL;

  /* Loop over list entries */
  for(Ftype_iterator i = ftypes.begin(); i != ftypes.end(); i++)
    {
      Ftype& ft = **i;
      std::string ftname;
      ft.get_name(ftname);
      if (ftname == name) 
	{
	  ftype = &ft;
	  break;
	}
    }
  return (ftype);
}


// Find ftype by name.
Ftype *Fuda::ftype_find(const char name[])
{
  std::string nm = name;
  return (ftype_find(nm));
}


// Func iterator begin and end functions.
Func_iterator Fuda::func_begin() { return (funcs.begin()); }
Func_iterator Fuda::func_end() { return (funcs.end()); }
unsigned int Fuda::get_nfunc() { return (funcs.size()); }
void Fuda::func_set_sync() { func_sync = 1; }
void Fuda::func_set_rsync() { func_rsync = 1; }
void Fuda::func_clear_sync() { func_sync = 0; }
void Fuda::func_clear_rsync() { func_rsync = 0; }
bool Fuda::func_is_sync() { return (func_sync); }
bool Fuda::func_is_rsync() { return (func_rsync); }


// Add Func record.
Func *Fuda::func_add(std::string name,
		      std::string ftype_name,
		      std::string dtype_name,
		      std::list<std::string>& p_names)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::FuncNameInvalid(name);

  /* Is the name free? */
  if (func_find(name) != NULL) throw Uferr::FuncNameAlreadyUsed(name);

  /* Does the Ftype exist? */
  Ftype *ftype = ftype_find(ftype_name);
  if (ftype == NULL) throw Uferr::FtypeNameInvalid(ftype_name);

  /* Does the Dtype exist? */
  Dtype *dtype = dtype_find(dtype_name);
  if (dtype == NULL) throw Uferr::DtypeNameInvalid(dtype_name);

  /* Create a func record */
  Func *func = new Func(this, name, ftype, dtype);

  /* Do we have the right number of variable parameter names */
  if (p_names.size() != ftype->get_nparam())
    {
      delete func;
      throw Uferr::NparamInvalid(p_names.size());
    }
  
  /* Check variable paramters and set parameter references */
  unsigned int p_count = 0;
  for(std::list<std::string>::iterator pi=p_names.begin(); pi!=p_names.end(); pi++)
    {
      Param *p_ref = param_find(*pi);
      if (p_ref != NULL)
	{
	  // If corresponding ftype param is not variable, we check.
	  if (!(func->get_ftype()->is_p_var(p_count)) &&
	      p_ref->get_kind()!=FUDA::PKIND_CONST)
	    {
	      delete func;
	      throw Uferr::ParamNotConst(p_count);
	    }
	  func->set_param(p_count++, p_ref);
	}
      else
	{
	  delete func;
	  throw Uferr::ParamNameInvalid(*pi);
	}
    }
  
  // Initialize the eval.x vector of the function.
  func->eval_x_init();

  /* Add link */
  funcs.push_back(func);

  /* Return ref to added link */
  return (func);
}


// Delete a Func object in the func list.
void Fuda::func_del(std::string& name)
{
  // Is the name valid?
  if (name.size()<1) throw Uferr::FuncNameInvalid(name);

  /* Loop over list entries and search for function */
  Func_iterator fi;
  bool fiempty = true;
  std::string iname;
  for(Func_iterator i = func_begin(); i != func_end(); i++)
    {
      // Compare names and if they match, we set pi and break out.
      (*i)->get_name(iname);
      if ((*i)->name == name) 
	{
	  fi = i;
    fiempty = false;
	  break;
	}
    }

  /* Does the func exist? */
  if (fiempty == true) throw Uferr::FuncNameInvalid(name);

  // First destroy the function object.
  delete *fi;

  // Then remove the empty pointer from the funcs list.
  funcs.erase(fi);

  // Unsync functions and data.
  func_clear_sync();
  func_clear_rsync();
  data_clear_sync();
  data_clear_rsync();
}


// Delete a Func object in the funcs list.
void Fuda::func_del(const char name[])
{
  std::string nm = name;
  func_del(nm);
}


// Delete all declared functions.
void Fuda::func_del_all()
{
  /* Loop over list entries */
  for(Func_iterator i = func_begin(); i != func_end(); i++)
    {
      // Destroy the function object.
      delete *i;
    }
  
  // Then remove all the empty pointers from the funcs list.
  funcs.clear();

  // Clear sync flags of functions and data.
  func_clear_sync();
  func_clear_rsync();
  data_clear_sync();
  data_clear_rsync();
}


// Set a parameter associated with a Func record.
void Fuda::func_set_param(std::string name, unsigned int pindex,
			   std::string pname)
{
  /* The func must exist. */
  Func *fn = func_find(name);
  if (fn==NULL) throw Uferr::FuncNameInvalid(name);

  /* Get the ftype */
  Ftype *ftype = fn->get_ftype();

  /* Get the dtype */
  Dtype *dtype = fn->get_dtype();

  /* Is the parameter index within bounds */
  if (pindex<0 || pindex >= ftype->get_nparam())
    throw Uferr::ParamIndexInvalid(pindex);

  // Does the parameter name exist.
  Param *pm = param_find(pname);
  if (pm==NULL) throw Uferr::ParamNameInvalid(pname);
  
  // Get new parameter kind.
  FUDA::Param_kind kind = pm->get_kind();

  // We cannot change a parameter to an explanatory paramter which is not
  // in the dtype. 
  if (kind==FUDA::PKIND_EXPL)
    {
      // Search dtype for the explanatory variable.
      bool found = 0;
      for (unsigned int i=0; i<dtype->get_dim(); i++)
	if (pm==dtype->get_ip(i)) found = 1;
      if (!found) throw Uferr::ExplVarInvalid(pname);
    }
  
  // We cannot mutate a non-variable parameter to PARAM-kind.
  if (kind==FUDA::PKIND_PARAM && !ftype->is_p_var(pindex))
    throw Uferr::ParamKindInvalid();

  // Make the change.
  fn->set_param(pindex,pm);

  // Update the eval structure of the function.
  fn->eval_x_init();

  // Clear func sync flag.
  func_clear_sync();

  // Clear eval stucture which may refer to old parameter.
  eval_clear();
}


// Find func record by name. Return NULL if not found.
Func *Fuda::func_find(std::string& name)
{
  Func *func = NULL;

  /* Loop over list entries */
  for(Func_iterator i = func_begin(); i != func_end(); i++)
    {
      Func& fn = **i;
      if (fn.name == name) 
	{
	  func = &fn;
	  break;
	}
    }
  return (func);
}


// Find func record by name.
Func *Fuda::func_find(const char name[])
{
  std::string nm = name;
  return func_find(nm);
}


// Return true if eval structure and relevant lists are in sync.
bool Fuda::eval_is_sync()
{
  // We return true if eval.sync and all lists are in sync.
  return (eval.sync && param_is_sync() && dtype_is_sync() && 
	  data_is_sync() && func_is_sync());
}


// Return true if eval structure and relevant lists are in rsync.
bool Fuda::eval_is_rsync()
{
  // We return true if eval.sync and all lists are in rsync.
  return (eval.sync && param_is_rsync() && dtype_is_rsync() && 
	  data_is_rsync() && func_is_rsync());
}


// Set eval and list sync flags.
void Fuda::eval_set_sync()
{
  /* Set all sync and rsync flags. */
  eval.sync = 1;
  param_set_sync();
  param_set_rsync();
  dtype_set_sync(); 
  dtype_set_rsync(); 
  data_set_sync();
  data_set_rsync();
  func_set_sync();
  func_set_rsync();
}


// Clear the eval structure and clear eval sync flag.
void Fuda::eval_clear()
{
  // Clear eval.sync and eval data structures plus dependent structures.
  eval.sync = 0;
  eval.x.clear();
  eval.p.clear();
  eval.c.clear();
  eval.d.clear();

  // Increment modification counter.
  eval_incr_mod_count();

  // Clear dtype eval structures.
  for (Dtype_iterator dti=dtype_begin(); dti!=dtype_end(); dti++)
    {
      Dtype& dt = **dti;
      dt.eval.funcs.clear();
      dt.eval.data.clear();
    }

  // Clear param eval structures.
  for (Param_iterator pi=param_begin(); pi!=param_end(); pi++)
    {
      Param& pm = **pi;
      pm.eval.active = 0;
      pm.eval.index = 0;
    }
}


// Initialise/setup eval structure.
void Fuda::eval_init()
{
  // If we are in sync, we do nothing.
  if (eval_is_sync() && eval_is_rsync()) return;

  // Clear fuda eval data structures.
  eval_clear();

  // Declare and zero data counter.
  unsigned int d_count = 0;

  /* Loop over functions, distribute them on dtype eval.funcs list and
     count data points. */
  for (Func_iterator fi=func_begin(); fi!=func_end(); fi++)
    {
      Func& fn = **fi;

      // We only take functions in use.
      if (fn.use)
	{
	  // Add to dtype eval.funcs list.
	  fn.dtype->eval.funcs.push_back(&fn);

	  // Add data to dtype eval.data list.
	  for(Data_iterator di=fn.data_begin(); di!=fn.data_end(); di++)
	    {
	      Data& data = *di;
	      fn.dtype->eval.data.push_back(&data);
	      d_count++;
	    }	  
	}
    }

  // Declare expl. variable, indep. parameter, and constant counters.
  unsigned int x_count = 0, p_count = 0, c_count = 0;

  /* Loop over dtypes and functions and count number of explanatory,
     free and constant parameters */
  for (Dtype_iterator dti=dtype_begin(); dti!=dtype_end(); dti++)
    {
      Dtype& dt = **dti;

      // Setup functions with dtype containing data.
      if (dt.eval.data.size()>0)
	{
	  // Loop over functions.
	  for (std::list<Func*>::iterator fi=dt.eval.funcs.begin();
	       fi!=dt.eval.funcs.end(); fi++)
	    {
	      Func& fn = **fi;

	      // Loop over function parameters.
	      for (unsigned int i=0; i<fn.ftype->get_nparam(); i++)
		{
		  // If not set, set the parameter eval.active flag and
		  // step expl. var, independ. param or const. counter.
		  Param *pm = fn.get_param(i);
		  if (!(pm->eval.active))
		    {
		      /* set active flag and increment expl, free
			 param. or c counter. */
		      pm->eval.active = 1;
		      if (pm->kind==FUDA::PKIND_EXPL)
			x_count++;
		      else if (pm->kind==FUDA::PKIND_PARAM)
			if (pm->free_flg) p_count++;
			else c_count++;
		      else if (pm->kind==FUDA::PKIND_CONST)
			c_count++;
		      else throw Uferr::ParamKindInvalid();
		    }
		}
	    }
	}
      else 
	{
	  // Clear eval.funcs lists when no data are present.
	  dt.eval.funcs.clear();
	}
    }

  // Loop over Dtypes, purge data and recalc d_count.
  d_count = 0;
  for (Dtype_iterator dti=dtype_begin(); dti!=dtype_end(); dti++)
    {
      Dtype& dt = **dti;
      dt.purge();
      d_count += dt.eval.data.size();
    }

  // We need at least one data point.
  if (d_count<1) throw Uferr::NoData();

  // Resize eval.d vector to number of data points.
  eval.d.resize(d_count);

  // Loop over data of all dtypes and setup eval.d Data* vector.
  unsigned int idata = 0;
  for (Dtype_iterator dti=dtype_begin(); dti!=dtype_end(); dti++)
    {
      Dtype& dt = **dti;

      // Loop over data.
      for (std::list<Data*>::iterator di=dt.eval.data.begin();
	   di!=dt.eval.data.end(); di++)
	{
	  // Set eval.d vector element.
	  eval.d[idata] = *di;
          idata++;
	}
    }

  // Allocate eval.x, eval.p and eval.c vectors.
  eval.x.resize(x_count);
  eval.p.resize(p_count);
  eval.c.resize(c_count);

  // Declare x, p and c indexes (for x, p and c vectors).
  unsigned int x_index = 0, p_index = 0, c_index = 0;

  // Loop over params and fill fuda.eval x, p and c vectors.
  for (Param_iterator pi=param_begin(); pi!=param_end(); pi++)
    {
      Param& pm = **pi;
      if (pm.eval.active)
	{
	  // Branch out according to expl, free param or constant.
	  if (pm.kind==FUDA::PKIND_EXPL)
	    {
	      // Save explanatory variable index in parameter.
	      pm.eval.index = x_index;
	      
	      // Save pointer to parameter in eval.x vector.
	      eval.x[x_index] = &pm;
	      
	      // Increment index.
	      x_index++;
	    }
	  else if (pm.kind==FUDA::PKIND_PARAM)
	    {
	      if (pm.free_flg)
		{
		  // Save dependent parameter index in parameter.
		  pm.eval.index = p_index;
		  
		  // Save pointer to parameter in eval.p vector.
		  eval.p[p_index] = &pm;
		  
		  // Increment index.
		  p_index++;
		}
	      else
		{
		  // Save constant parameter index in parameter.
		  pm.eval.index = c_index;
		  
		  // Save pointer to parameter in eval.c vector.
		  eval.c[c_index] = &pm;
		  
		  // Increment index.
		  c_index++;
		}
	    }
	  else if (pm.kind==FUDA::PKIND_CONST)
	    {
	      // Save constant parameter index in parameter.
	      pm.eval.index = c_index;
	      
	      // Save pointer to parameter in eval.c vector.
	      eval.c[c_index] = &pm;
	      
	      // Increment index.
	      c_index++;
	    }
	  else throw Uferr::ParamKindInvalid();
	}
    }


  /* Loop over dtypes and eval.funcs to setup eval.p and eval.c for
     each function. The eval.x vector is fully setup when the function
     is created or modified. */
  for (Dtype_iterator dti=dtype_begin(); dti!=dtype_end(); dti++)
    {
      Dtype& dt = **dti;

      // Loop over functions.
      for (std::list<Func*>::iterator fi=dt.eval.funcs.begin();
	   fi!=dt.eval.funcs.end(); fi++)
	{
	  Func& fn = **fi;

	  // Initialize the eval.p and eval.c structure of the function.
	  fn.eval_pc_init();

	  // eval.p setup.
	  for (unsigned int i=0; i<fn.eval.p.size(); i++)
	    {
	      fn.eval.p[i].val = 
		fn.params[fn.eval.p[i].iv]->get_eval_index();
	    }
	  
	  // eval.c setup.
	  for (unsigned int i=0; i<fn.eval.c.size(); i++)
	    {	      
	      fn.eval.c[i].val = 
		fn.params[fn.eval.c[i].iv]->get_eval_index();
	    }	  
	}
    }

  // Set the sync flags.
  eval_set_sync();
}


// Get number of free parameters in fit.
unsigned int Fuda::eval_get_nfree()
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  return (eval.p.size());
}


// Return Pointer to free param.
Param *Fuda::eval_get_free(unsigned int i)
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  if (i<0 || i>=eval_get_nfree()) throw Uferr::EvalFreeIndexInvalid(i);
  return (eval.p[i]);
}


// Get number of explanatory parameters in fit.
unsigned int Fuda::eval_get_nexpl() 
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  return (eval.x.size());
}


// Return Pointer to explanatory param.
Param *Fuda::eval_get_expl(unsigned int i)
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  if (i<0 || i>=eval_get_nexpl()) throw Uferr::EvalExplIndexInvalid(i);
  return (eval.x[i]);
}


// Get number of constant (fixed) parameters in fit.
unsigned int Fuda::eval_get_nconst() 
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  return (eval.c.size());
}


// Return Pointer to constant param.
Param *Fuda::eval_get_const(unsigned int i)
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  if (i<0 || i>=eval_get_nconst()) throw Uferr::EvalConstIndexInvalid(i);
  return (eval.c[i]);
}


// Get number of data points.
unsigned int Fuda::eval_get_ndata() 
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  return (eval.d.size());
}


// Return Pointer to Data record.
Data *Fuda::eval_get_data(unsigned int i)
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  if (i<0 || i>=eval_get_ndata()) throw Uferr::EvalDataIndexInvalid(i);
  return (eval.d[i]);
}


// Return Pointer to Dtype record.
Dtype *Fuda::eval_get_dtype(unsigned int i)
{
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  if (i<0 || i>=eval_get_ndata()) throw Uferr::EvalDataIndexInvalid(i);
  return (eval.d[i]->get_func()->get_dtype());
}


/* Evaluate overall function value and partial derrivatives with
   respect to the parameters for a specified data type and a vector of
   explanatory variables. */
double Fuda::eval_call(Dtype& dtype, std::vector<double>& xvec,
		       double pvec[], double dpvec[], bool dp_flg)
{
  // Eval structure must be in rsync.
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  
  // Get reference to fcall record.
  FuncCall& func_call = *(fcall_get_ref());

  // If dp_flg set, we clear dpvec.
  if (dp_flg)
    for(unsigned int i=0; i<eval.p.size(); i++)
      dpvec[i] = 0.0;

  // Accumulated (calculated) function value.
  double calc = 0.0;
  
  // Loop over functions for dtype.
  for (std::list<Func*>::iterator fi=dtype.eval.funcs.begin();
       fi!=dtype.eval.funcs.end(); fi++)
    {
      // Get ref. to func record.
      Func& fn = **fi;
      
      // Setup arguments for function call (in f_call).
      
      // Loop over explanatory variables.
      for(unsigned int i=0; i<fn.eval.x.size(); i++)
        func_call.p[fn.eval.x[i].iv] = xvec[fn.eval.x[i].val];
      
      // Loop over dependent parameters.
      for(unsigned int i=0; i<fn.eval.p.size(); i++)
        func_call.p[fn.eval.p[i].iv] = pvec[fn.eval.p[i].val];
      
      // Loop over constant parameters.
      for(unsigned int i=0; i<fn.eval.c.size(); i++)
        func_call.p[fn.eval.c[i].iv] =
          eval.c[fn.eval.c[i].val]->get_value();
      
      // Setup function call.
      double value;

      /* Branch out according to whether derrivatives are to be calculated
         or not and call function */
      if (dp_flg)
        {
          // We setup func_call.dp_flg array.
          
          // Zero func_call.dp_flg.
          for(unsigned int i=0; i<func_call.psize; i++)
            func_call.dp_flg[i] = 0;
      
          // Set dp_flg for dependent parameters.
          for(unsigned int i=0; i<fn.eval.p.size(); i++)
            func_call.dp_flg[fn.eval.p[i].iv] = 1;

          // Function call.
          fn.ftype->call(func_call.p, func_call.dp_flg,
			 func_call.dp, &value);

          // Accumulate calculated dp values.
          for(unsigned int i=0; i<fn.eval.p.size(); i++)
            dpvec[fn.eval.p[i].val] += func_call.dp[fn.eval.p[i].iv];
        }
      else
        {
          fn.ftype->call(func_call.p, func_call.dp_flg_false,
			 func_call.dp, &value);
        }

      // Accumulate value.
      calc += value;
    }

  return (calc);
}


/* Evaluate overall function value and partial derrivatives with
   respect to the parameters for a data point referenced by the data
   index in the eval structure. */
double Fuda::eval_call(unsigned int d_index, double pvec[],
                        double dpvec[], bool dp_flg)
{
  // Eval structure must be in rsync.
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  
  // Get references to data and dtype records.
  Data& data = *(eval.d[d_index]);
  Dtype& dtype = *(eval.d[d_index]->get_func()->get_dtype());

  // Call general eval_call taking dtype and xvec arguments.
  double calc = eval_call(dtype, data.x, pvec, dpvec, dp_flg);

  return (calc);
}


/* Evaluate overall (obs-calc)/u and partial derrivatives with
   respect to the parameters for a data point referenced by the data
   index in the eval structure. */
double Fuda::eval_obs_calc(unsigned int d_index, double pvec[],
			    double dpvec[], bool dp_flg)
{
  // Get references to data record.
  Data& data = *(eval.d[d_index]);

  // We let eval_call handle all exceptions, so no checking.

  // First evaluate calculated value and possibly derrivatives.
  double calc = eval_call(d_index, pvec, dpvec, dp_flg);
  
  /* If dp_flg set, we scale dpvec with -1/uncertainty for the data
     point. The minus is because we need to return the derrivatives of
     (obs-calc)/u. */
  if (dp_flg)
    for(unsigned int i=0; i<eval.p.size(); i++)
      dpvec[i] /= -data.u;

  // Return (obs-calc)/uncertainty.
  double obs_calc = (data.value-calc)/data.u;

  return (obs_calc);
}


/* Evaluate euclidian norm for the free parameters given in pvec. */
double Fuda::eval_enorm(double pvec[])
{
  // Get number of data points.
  unsigned int ndata = eval_get_ndata();

  // Get reference to wspace.
  double *dvec;
  if (wspace_lock())
    dvec = wspace_ref(ndata);
  else
    exit(1);

  // Loop over data records and setup obs-calc vector.
  for (unsigned int d_index=0; d_index<ndata; d_index++)
    {
      // Calculate (obs-calc)/u).
      dvec[d_index] = eval_obs_calc(d_index, pvec, 0, 0);
    }

  // Unlock workspace.
  wspace_unlock();

  // Return euclidian norm.
  return (FUDA::calc_enorm(ndata,pvec));
}


/* Evaluate euclidian norm for the free parameters. */
double Fuda::eval_enorm()
{
  // Get number of free paramters.
  unsigned int nfree = eval_get_nfree();

  // Get number of data points.
  unsigned int ndata = eval_get_ndata();

  // Get reference to wspace.
  double *dvec, *pvec;
  if (wspace_lock())
    {
      pvec = wspace_ref(nfree+ndata);
      dvec = &pvec[nfree];
    }
  else
    exit(1);

  // Setup pvec.
  for(unsigned int p_index=0; p_index<nfree; p_index++)
    pvec[p_index] = eval_get_free(p_index)->get_value();
  
  // Loop over data records and setup obs-calc vector.
  for (unsigned int d_index=0; d_index<ndata; d_index++)
    {
      // Calculate (obs-calc)/u).
      dvec[d_index] = eval_obs_calc(d_index, pvec, 0, 0);
    }

  // Calculate Euclidian norm.
  double enorm = FUDA::calc_enorm(ndata,dvec);
  
  // Unlock workspace.
  wspace_unlock();

  // Return enorm.
  return (enorm);
}


/* Recalculate the value of all data points for the current parameters
   and initialize the data records with the calculated values. Be ware
   that the previous initial value and value is overwritten in this
   process. */
void Fuda::eval_data_recalc()
{
  // Eval structure must be in rsync.
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  
  // Get number of free paramters.
  unsigned int nfree = eval_get_nfree();

  // Get number of data points.
  unsigned int ndata = eval_get_ndata();

  // Get reference to wspace.
  double *pvec;
  if (wspace_lock())
    {
      pvec = wspace_ref(nfree);
    }
  else
    exit(1);

  // Setup pvec with current parameter values.
  for(unsigned int p_index=0; p_index<nfree; p_index++)
    pvec[p_index] = eval_get_free(p_index)->get_value();
  
  // Loop over data records and set evaluated value in data record.
  for (unsigned int d_index=0; d_index<ndata; d_index++)
    {
      // Calculate function value and set it in the data record.
      eval.d[d_index]->set_init(eval_call(d_index, pvec, 0, 0));
    }

  // Unlock workspace.
  wspace_unlock();
}


/* Generate montecarlo data with a gaussian distribution of noise from
   initial data value and uncertainty for each data point. */
void Fuda::eval_data_random()
{
  // Eval structure must be in rsync.
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  
  // Loop over data records and set evaluated value in data record.
  for (unsigned int d_index=0; d_index<eval_get_ndata(); d_index++)
    {
      // Calculate function value and set it in the data record.
      eval.d[d_index]->montecarlo_value(FUDA::rand_gauss(0.0,1.0));
    }
}


/* Set value to initial value for all data in eval structure. */
void Fuda::eval_data_init()
{
  // Eval structure must be in rsync.
  if (!eval_is_rsync()) throw Uferr::EvalNotSync();
  
  // Loop over data records and set initial value.
  for (unsigned int d_index=0; d_index<eval_get_ndata(); d_index++)
    {
      // Calculate function value and set it in the data record.
      eval.d[d_index]->init_value();
    }
}


void Fuda::get_version(std::string& version_string)
{
  char str[32];
  snprintf(str, sizeof(str),"%d.%d",get_version_major(),get_version_minor());
  version_string = str;
}


unsigned int Fuda::get_version_major()
{
  return (FUDA_VERSION_MAJOR);
}


unsigned int Fuda::get_version_minor()
{
  return (FUDA_VERSION_MINOR);
}


// Print contents of all (most) fuda data structures.
void Fuda::print()
{
  std::cout << "***************\n";
  std::cout << "*** ptypes  ***\n";
  std::cout << "***************\n";
  for(Ptype_iterator pt=ptype_begin(); pt!=ptype_end(); pt++)
    {
      (*pt)->print();
    }
  std::cout << "\n";

  std::cout << "***************\n";
  std::cout << "*** params  ***\n";
  std::cout << "***************\n";
  std::cout << "param_sync  : " << param_is_sync() << "\n";
  std::cout << "param_rsync : " << param_is_rsync() << "\n\n";
  for(Param_iterator pm=param_begin(); pm!=param_end(); pm++)
    {
      (*pm)->print();
    }
  std::cout << "\n";

  std::cout << "***************\n";
  std::cout << "*** dtypes  ***\n";
  std::cout << "***************\n";
  std::cout << "dtype_sync  : " << dtype_is_sync() << "\n";
  std::cout << "dtype_rsync : " << dtype_is_rsync() << "\n\n";
  for(Dtype_iterator dt=dtype_begin(); dt!=dtype_end(); dt++)
    {
      (*dt)->print();
    }
  std::cout << "\n";

  std::cout << "***************\n";
  std::cout << "*** ftypes  ***\n";
  std::cout << "***************\n";
  std::cout << "fcall->defined : " << fcall->defined << "\n";
  std::cout << "fcall->psize   : " << fcall->psize << "\n";
  std::cout << "\n";

  for(Ftype_iterator ft=ftype_begin(); ft!=ftype_end(); ft++)
    {
      (*ft)->print();
    }
  std::cout << "\n";

  std::cout << "***************\n";
  std::cout << "*** funcs   ***\n";
  std::cout << "***************\n";
  std::cout << "func_sync  : " << func_is_sync() << "\n";
  std::cout << "func_rsync : " << func_is_rsync() << "\n\n";
  for(Func_iterator fn=func_begin(); fn!=func_end(); fn++)
    {
      (*fn)->print();
      // fn->print_data();
    }
  std::cout << "\n";

  std::cout << "***************\n";
  std::cout << "***  data   ***\n";
  std::cout << "***************\n";
  std::cout << "data_sync  : " << data_is_sync() << "\n";
  std::cout << "data_rsync : " << data_is_rsync() << "\n\n";
  std::cout << "\n";

  std::cout << "***************\n";
  std::cout << "*** eval    ***\n";
  std::cout << "***************\n";
  std::cout << "sync       : " << eval_is_sync() << "\n";
  std::cout << "mod_count  : " << eval_get_mod_count() << "\n";

  std::cout << "eval.x :";
  for(unsigned int i=0; i<eval.x.size(); i++)
    {
      Param& pm = *(eval.x[i]);
      std::string str;
      pm.get_name(str);
      std::cout << " [" << i << "]:" << str;
    }
  std::cout << "\n";

  std::cout << "eval.p :";
  for(unsigned int i=0; i<eval.p.size(); i++)
    {
      Param& pm = *(eval.p[i]);
      std::string str;
      pm.get_name(str);
      std::cout << " [" << i << "]:" << str;
    }
  std::cout << "\n";

  std::cout << "eval.c :";
  for(unsigned int i=0; i<eval.c.size(); i++)
    {
      Param& pm = *(eval.c[i]);
      std::string str;
      pm.get_name(str);
      std::cout << " [" << i << "]:" << str;
    }
  std::cout << "\n";

  std::cout << "eval.d.size() : " << eval.d.size() << "\n";

  std::cout << "\n";

      
  std::cout << "\n";
}








