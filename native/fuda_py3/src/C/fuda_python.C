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

#include <Python.h>
#include <numpy/arrayobject.h>
#include <cstring>
#include "fuda_classes.H"
#include "fuda_utils.H"
#include "fuda_MinpackLM.H"
#include "fuda_cftypes.H"


// Some global variables.
static Fuda *fuda;
static MinpackLM *lm;
static unsigned int seed;
static char fstr[512];
static std::string current_func_name = "";
static Func *current_func = NULL;

// Variables used by ftype_find and func_find.
static std::string last_func_name = "";
static Func *last_func = NULL;
static std::string last_ftype_name = "";
static Ftype *last_ftype = NULL;


bool opt_arg(PyObject *args, PyObject *kw, char *key , int iarg)
{
  /* This function checks for the presence of an optional argument
     number iarg with the keyword key by checking the args and kw
     objects. */
  return (kw!=NULL &&
	  (PyTuple_Size(args)>=iarg+1 ||
	   PyDict_GetItemString (kw, key)!=NULL));
}


Ftype *ftype_find(std::string &name)
{
  /* This procedure checks the last found ftype to speed up the search
     for the peviously found ftype. */

  // Check last found ftype.
  if (last_ftype != NULL && name==last_ftype_name)
    return (last_ftype);
  else
    {
      Ftype *ft = fuda->ftype_find(name);
      if (ft != NULL)
	{
	  last_ftype = ft;
	  last_ftype_name = name;
	}
      return (ft);
    }
}


Ftype *ftype_find(const char name[])
{
  std::string nm = name;
  return (ftype_find(nm));
}


Func *func_find(std::string &name)
{
  /* This procedure checks the last found func to speed up the search
     for the peviously found func. */

  // Check last found func.
  if (last_func != NULL && name==last_func_name)
    return (last_func);
  else
    {
      Func *fn = fuda->func_find(name);
      if (fn != NULL)
	{
	  last_func = fn;
	  last_func_name = name;
	}
      return (fn);
    }
}


Func *func_find(const char name[])
{
  std::string nm = name;
  return (func_find(nm));
}

//------------------------------------------------------------
//------------------------------------------------------------
// fudaPYTHON begin

class FtypePYTHON : public Ftype
{
  PyObject *pyfunc;       /* Pointer to callable python object*/
  PyObject *pyfs;         /* Pointer to python object which is given
			     to the pyfunc and so provides a way to
			     supply extra information to func */
  // Python arguments for calling python function.
  PyArrayObject *arg_p;     
  PyArrayObject *arg_dp_flg;
  PyArrayObject *arg_dp;
  PyArrayObject *arg_v;
  PyObject      *arg_list;
  
 public:
  friend class Fuda;

  // Constructor and destructor.
  FtypePYTHON (Fuda *a_fuda, std::string a_name, unsigned int nparam,
	       PyObject *a_pyfunc, PyObject *a_pyfs);
  ~FtypePYTHON();

  PyObject *get_pyfunc();
  PyObject *get_pyfs();
  double call(std::vector<double>& p_vec);
  void call(double p[], int dp_flg[], double dp[], double *value);
  FUDA::Ftype_type get_type();
  void print();
};

// Constructor.
FtypePYTHON::FtypePYTHON(Fuda *a_fuda, std::string a_name,
			 unsigned int a_nparam,
			 PyObject *a_pyfunc, PyObject *a_pyfs)
  : Ftype(a_fuda, a_name)
{
  // Save pyfunction reference and reference to pyfs object and
  // increment reference counts.
  pyfunc = a_pyfunc;
  Py_XINCREF(pyfunc);
  pyfs = a_pyfs;
  Py_XINCREF(pyfs);

  // Dimension the ftype.
  set_nparam(a_nparam);

  // Check that python function object is callable.
  if (!PyCallable_Check(pyfunc))
    {
      std::string msg = "FtypePython: Python function object is not callable\n";
      throw Uferr::FtypeCallError(msg);
    }

  // Check that python function object is callable.
  if (PyCallable_Check(pyfs))
    {
      std::string msg = "FtypePYTHON: Python function argument object is callable\n";
      throw Uferr::FtypeCallError(msg);
    }

  // Construct the python function argument object.

  /* The python function must takes five arguments:
     p (mutable list), dp_flg (mutable list), dp (mutable list),
     v (mutable list with one element which is the return value),
     pyfs (PyObject reference to object passed to the function).
  */

  // Create Numeric python array objects for calling python function.
  unsigned int dim = 1;
  //int *dimsize_ptr = new int[dim];
  npy_intp *dimsize_ptr = new npy_intp[dim];
  dimsize_ptr[0] = a_nparam;
  //arg_p = 
  //(PyArrayObject*) PyArray_FromDims(dim,dimsize_ptr,PyArray_DOUBLE);
  //arg_dp_flg = 
  //(PyArrayObject*) PyArray_FromDims(dim,dimsize_ptr,PyArray_INT);
  //arg_dp = 
  //(PyArrayObject*) PyArray_FromDims(dim,dimsize_ptr,PyArray_DOUBLE);

  arg_p = 
    (PyArrayObject*) PyArray_SimpleNew(dim,dimsize_ptr,NPY_DOUBLE);
  arg_dp_flg = 
    (PyArrayObject*) PyArray_SimpleNew(dim,dimsize_ptr,NPY_INT32);
  arg_dp = 
  (PyArrayObject*) PyArray_SimpleNew(dim,dimsize_ptr,NPY_DOUBLE);
  dimsize_ptr[0] = 1;
  //arg_v = 
  //(PyArrayObject*) PyArray_FromDims(dim,dimsize_ptr,PyArray_DOUBLE);

  arg_v = 
    (PyArrayObject*) PyArray_SimpleNew(dim,dimsize_ptr,NPY_DOUBLE);
  
  delete[] dimsize_ptr;

  // Create argument list for calling python function.
  arg_list = PyTuple_New (5);
  
  PyTuple_SetItem(arg_list, 0, (PyObject*) arg_p);
  Py_XINCREF(arg_p);
  PyTuple_SetItem(arg_list, 1, (PyObject*) arg_dp_flg);
  Py_XINCREF(arg_dp_flg);
  PyTuple_SetItem(arg_list, 2, (PyObject*) arg_dp);
  Py_XINCREF(arg_dp);
  PyTuple_SetItem(arg_list, 3, (PyObject*) arg_v);
  Py_XINCREF(arg_v);
  PyTuple_SetItem(arg_list, 4, (PyObject*) pyfs);
  Py_XINCREF(pyfs);
}


FtypePYTHON::~FtypePYTHON()
{
  std::cout << "FtypePYTHON::~Ftype() not implemented\n";
  exit(0);  
}


PyObject *FtypePYTHON::get_pyfunc() 
{
  return (pyfunc);
}


PyObject *FtypePYTHON::get_pyfs()
{
  return (pyfs);
}


FUDA::Ftype_type FtypePYTHON::get_type()
{
  return (FUDA::FTYPE_PYTHON);
}


// Make direct call to underlying Ftype function.
void FtypePYTHON::call(double p[], int dp_flg[], double dp[], double *value)
{
  unsigned int np = get_nparam();

  // Copy values to python argument objects.
  double *arg_p_ref = (double*) arg_p->data;
  int *arg_dp_flg_ref = (int*) arg_dp_flg->data;
  for (unsigned int i=0; i<np; i++)
    {
      arg_p_ref[i] = p[i];
      arg_dp_flg_ref[i] = dp_flg[i];
    }

  // Call the python function.
  //PyObject *return_obj = PyEval_CallObject(pyfunc, arg_list);
  PyObject *return_obj = PyObject_Call(pyfunc, arg_list,NULL);

  // Check return value. NULL means that the call failed.
  if (return_obj==NULL)
    {
      // Clear python error state.
      PyErr_Clear();      
      
      std::string msg;
      msg="FtypePYTHON: Python function call failed\n";
      throw Uferr::FtypeCallError(msg);
    }

    // Copy return value.
  double *arg_v_ref = (double*) arg_v->data;
  *value = arg_v_ref[0];

  // Copy derivatives.
  double *arg_dp_ref = (double*) arg_dp->data;
  for (unsigned int i=0; i<np; i++)
      dp[i] = arg_dp_ref[i];
  
  // Check that return object is an integer object.
  if (!PyLong_Check(return_obj))  // CB: changed from PyInt_Check to PyLong_Check in line with updated API
    {
      // Clear python error state.
      PyErr_Clear();      
      
      std::string msg;
      msg="FtypePYTHON: Python function did not return an integer value\n";
      throw Uferr::FtypeCallError(msg);
    }

  // Get the integer return state.
  int return_val = PyLong_AsLong(return_obj); // CB: changed from PyInt_AsLong to PyLong_AsLong in line with updated API
  Py_XDECREF(return_obj);

  // Check that return value is zero.
  if (return_val != 0)
    {
      std::string msg;
      msg="FtypePYTHON: Python function returned non-zero value for ftype ";
      std::string name;
      get_name(name);
      msg+=name+"\n";
      throw Uferr::FtypeCallError(msg);
    }
}


// Calculate ftype function as a function of all parameters.
double FtypePYTHON::call(std::vector<double>& p_vec)
{
  // Get number of parameters.
  unsigned int np = get_nparam();

  // Vector size must match number of parameters in function.
  if (p_vec.size()!=np)
    throw Uferr::NumParamInvalid(p_vec.size());

  // Copy values to python argument objects.
  double *arg_p_ref = (double*) arg_p->data;
  int *arg_dp_flg_ref = (int*) arg_dp_flg->data;
  for (unsigned int i=0; i<np; i++)
    {
      arg_p_ref[i] = p_vec[i];
      arg_dp_flg_ref[i] = 0;
    }

  // Call the python function.
  //PyObject *return_obj = PyEval_CallObject(pyfunc, arg_list);
  PyObject *return_obj = PyObject_Call(pyfunc, arg_list,NULL);

  // Check return value.
  if (return_obj==NULL)
    {
      PyErr_Clear();      
      std::string msg = "FtypePYTHON: Python function call failed\n";
      throw Uferr::FtypeCallError(msg);
    }

  // Decrease reference count for return object, as we don't use it.
  Py_XDECREF(return_obj);

  // Copy function return value.
  double *arg_v_ref = (double*) arg_v->data;
  double value = arg_v_ref[0];

  // Return the function value.
  return (value);
}


void FtypePYTHON::print()
{
  Ftype::print();
}


extern "C" PyObject *
fuda_ftype_python(PyObject *self, PyObject *args)
{
  /* fuda_ftype_python declares a new python function type. The
     required arguments are: ftype_name, function reference and python
     data object. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=4)
    {
      std::string err =  "fuda_ftype_python: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 4)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_ftype_python: 1. arg (name) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj); // CB: upadted PyUnicode_AsUTF8() with PyUnicode_asUTF8 which should encode the same. this function returns a const char so this is updated too

  // 2nd argument must be the parameter count of the function.
  PyObject *nparam_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(nparam_obj)) // CB: PyInt_Check --> PyLong_Check update
    {
      std::string err =  "fuda_ftype_python: 2. arg (nparam) must be an integer";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get nparam and check bounds.
  int nparam = PyLong_AsLong(nparam_obj); // CB: PyInt_AsLong --> PyLong_AsLong update
  if (nparam<1 || nparam>50)
    {
      std::string err =  "fuda_ftype_python: 2. arg (nparam) invalid";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 3rd argument must be the python function to declare.
  PyObject *pyfunc_obj = PyTuple_GetItem(args, 2);
  if (!PyCallable_Check(pyfunc_obj))
    {
      std::string err =  "fuda_ftype_python: 3. arg (pyfunc) must be a callable object";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 4th argument must be some non-callable python object.
  PyObject *pyfs_obj = PyTuple_GetItem(args, 3);
  if (PyCallable_Check(pyfs_obj))
    {
      std::string err =  "fuda_ftype_python: 4. arg (pyfs) must be a non-callable python object";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }


  // Add/register the ftype.
  FtypePYTHON *ftype_ptr;  
  try {
    ftype_ptr = new FtypePYTHON(fuda, name, nparam, pyfunc_obj, pyfs_obj);
    fuda->ftype_add(ftype_ptr);
  }
  catch (Uferr::FtypeNameInvalid& e) {
    std::string err =  "fuda_ftype_python: name invalid: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FtypeNameAlreadyUsed& e) {
    std::string err =  "fuda_ftype_python: name already used: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_ftype_python: failed with unexpected exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}





// fudaPYTHON end
//------------------------------------------------------------
//------------------------------------------------------------


extern "C" PyObject *
fuda_fuda_print(PyObject *self, PyObject *args)
{
  /* fuda_fuda_print simply calls the fuda print method which in a
     rather brute force manner lists the fuda data structure. Meant
     for debugging. */
  fuda->print();
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_ptype(PyObject *self, PyObject *args, PyObject *kw)
{
  /* fuda_ptype declares a new ptype or assigns the value or other
     attribute for an already assigned ptype. It always takes the
     ptype name as the first argument followed by a key=keyval
     list. When attributes are set for already declared ptypes, the
     same attributes are set for all parameters of the same ptype. */

  // Argument keyword list.
  static char kw_name[] = "name", kw_value[] = "value", kw_free[] = "free";
  static char kw_bounds[] = "bounds", kw_lower[] = "lower", kw_upper[] = "upper";
  static char kw_delta[] = "delta", kw_kind[] = "kind", kw_esd[] = "esd";
  static char *kwlist[] = {kw_name, kw_value, kw_free, kw_bounds, kw_lower,
                           kw_upper, kw_delta, kw_kind, kw_esd, NULL};
  char *name, *kind_name;
  double value, lower, upper, delta, esd;
  int free, bounds;

  // Check for optional keyword arguments.
  int iarg = 1;
  bool value_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool free_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool bounds_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool lower_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool upper_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool delta_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool kind_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool esd_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;

  // Pass argument list.
  if (!PyArg_ParseTupleAndKeywords(args, kw, "s|diidddsd", kwlist, 
				   &name, &value, &free,
				   &bounds, &lower, &upper, &delta,
				   &kind_name, &esd))
    return NULL; 

  // Find the ptype - if not found pm will be NULL.
  Ptype *pt = fuda->ptype_find(name);
  
  // The kind cannot be set for an existing ptype.
  if (pt!=0 and kind_flg)
    {
      std::string err = "fuda_ptype: the kind cannot be changed: ";
      err += kind_name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Set the parameter kind for the ptype.
  FUDA::Param_kind kind;
  if (kind_flg)
    if (strcmp(kind_name,"CONST")==0) 
      {
	kind = FUDA::PKIND_CONST;
      }
    else if (strcmp(kind_name,"PARAM")==0) 
      {
	kind = FUDA::PKIND_PARAM;
      }
    else if (strcmp(kind_name,"EXPL")==0) 
      {
	kind = FUDA::PKIND_EXPL;
      }
    else
      {
	std::string err = "fuda_param: invalid parameter kind: ";
	err += kind_name;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
  else
    {
      // Default parameter kind.
      kind = FUDA::PKIND_PARAM;
    }
  
  // If the ptype does not exist, we create it.
  if (pt==NULL)
    {
      // Declare ptype?
      try {
	pt = fuda->ptype_add(name,kind);
      }
      catch (Uferr::PtypeNameInvalid) {
	std::string err = "fuda_ptype: ptype name (type) invalid";
	PyErr_SetString(PyExc_Exception, err.c_str());
	std::cout << "fuda_param: param_add failed\n";
	return NULL;
      }	
      catch (...) {
	std::string err = "fuda_ptype: unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	std::cout << "fuda_param: param_add failed\n";
	return NULL;
      }	
    }
  else
    {
      if (kind_flg)
	{
	  std::string err = "fuda_ptype: cannot change parameter kind ";
	  err += "of existing ptype";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;	  
	}
    }
  
  // Set optional arguments.
  if (value_flg)
    {
      // Set value for ptype.
      pt->set_value(value);

      // Set value for matching parameters.
      for (Param_iterator pi=fuda->param_begin(); pi!=fuda->param_end(); pi++)
	if ((*pi)->get_ptype()==pt) (*pi)->set_value(value);
    }
  
  if (free_flg) {
    if (free==0)
      {
	// Set ptype to fixed.
	pt->fix();

	// Fix matching parameters.
	for (Param_iterator pi=fuda->param_begin();
	     pi!=fuda->param_end(); pi++)
	  if ((*pi)->get_ptype()==pt) (*pi)->fix();
      }
    else if (free==1)
      {
	// Set ptype to free
	pt->free();

	// Free matching parameters.
	for (Param_iterator pi=fuda->param_begin();
	     pi!=fuda->param_end(); pi++)
	  if ((*pi)->get_ptype()==pt) (*pi)->free();
      }
    else
      {
	std::string err = "fuda_ptype: free value invalid: ";
	snprintf(fstr, sizeof(fstr),"%d",free);
	err += fstr;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
  }

  if (bounds_flg) {
    if (bounds>=0 && bounds<=3)
      {
	// Set bounds for ptype.
	pt->set_bounds(bounds);

	// Set bounds for matching parameters.
	for (Param_iterator pi=fuda->param_begin();
	     pi!=fuda->param_end(); pi++)
	  if ((*pi)->get_ptype()==pt) (*pi)->set_bounds(bounds);
      }
    else
      {
	std::string err = "fuda_ptype: bounds value is invalid: ";
	snprintf(fstr, sizeof(fstr),"%d",bounds);
	err += fstr;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
  }

  if (lower_flg)
    {
      // Set lower for ptype.
      pt->set_lower(lower);

      // Set lower for matching parameters.
      for (Param_iterator pi=fuda->param_begin(); pi!=fuda->param_end(); pi++)
	if ((*pi)->get_ptype()==pt) (*pi)->set_lower(lower);
      
    }
  
  if (upper_flg)
    {
      // Set upper for ptype.
      pt->set_upper(upper);

      // Set upper for matching parameters.
      for (Param_iterator pi=fuda->param_begin(); pi!=fuda->param_end(); pi++)
	if ((*pi)->get_ptype()==pt) (*pi)->set_upper(upper);
    }
  
  if (delta_flg) {
    if (delta > 0.0)
      {
	// Set delta for ptype.
	pt->set_delta(delta);

	// Set delta for matching parameters.
	for (Param_iterator pi=fuda->param_begin(); 
	     pi!=fuda->param_end(); pi++)
	  if ((*pi)->get_ptype()==pt) (*pi)->set_delta(delta);
      }
    else
      {
	std::string err = "fuda_ptype: delta value invalid: ";
	snprintf(fstr, sizeof(fstr),"%f",delta);
	err += fstr;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
  }

  if (esd_flg)
    {
      // Set esd for ptype.
      pt->set_esd(esd);

      // Set esd for matching parameters.
      for (Param_iterator pi=fuda->param_begin(); pi!=fuda->param_end(); pi++)
	if ((*pi)->get_ptype()==pt) (*pi)->set_esd(esd);
    }
  

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_param(PyObject *self, PyObject *args, PyObject *kw)
{
  /* fuda_param declares a new parameter or assigns the value or other
     characteristics for an already assigned parameter. It always
     takes the parameter name as the first argument followed by a
     key=keyval list */

  // Argument keyword list.
  static char kw_name[] = "name", kw_value[] = "value", kw_free[] = "free";
  static char kw_bounds[] = "bounds", kw_lower[] = "lower", kw_upper[] = "upper";
  static char kw_delta[] = "delta", kw_kind[] = "kind", kw_type[] = "type", kw_esd[] = "esd";
  static char *kwlist[] = {kw_name, kw_value, kw_free, kw_bounds, kw_lower,
                           kw_upper, kw_delta, kw_kind, kw_type, kw_esd, NULL};
  char *name, *type_name, *kind_name;
  double value, lower, upper, delta, esd;
  int free, bounds;

  // Check for optional keyword arguments.
  int iarg = 1;
  bool value_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool free_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool bounds_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool lower_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool upper_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool delta_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool kind_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool type_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool esd_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;

  // Pass argument list.
  if (!PyArg_ParseTupleAndKeywords(args, kw, "s|diidddssd", kwlist, 
				   &name, &value, &free,
				   &bounds, &lower, &upper, &delta,
				   &kind_name, &type_name, &esd))
    return NULL; 

  // Find the parameter - if not found pm will be NULL.
  Param *pm = fuda->param_find(name);

  // The kind cannot be set for an existing parameter.
  if (pm!=0 and kind_flg)
    {
      std::string err = "fuda_param: the kind cannot be changed: ";
      err += kind_name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // The type cannot be set for an existing parameter.
  if (pm!=0 and type_flg)
    {
      std::string err = "fuda_param: the type cannot be changed: ";
      err += kind_name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // The type determins the kind, so type and kind are exclusive.
  if (type_flg and kind_flg)
    {
      std::string err = "fuda_param: type and kind are exclusive keywords: ";
      err += kind_name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Set the parameter kind.
  FUDA::Param_kind kind;
  if (kind_flg)
    if (strcmp(kind_name,"CONST")==0) 
      {
	kind = FUDA::PKIND_CONST;
      }
    else if (strcmp(kind_name,"PARAM")==0) 
      {
	kind = FUDA::PKIND_PARAM;
      }
    else if (strcmp(kind_name,"EXPL")==0) 
      {
	kind = FUDA::PKIND_EXPL;
      }
    else
      {
	std::string err = "fuda_param: invalid parameter kind: ";
	err += kind_name;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
  else
    {
      // Default parameter kind.
      kind = FUDA::PKIND_PARAM;
    }
  
  // If the parameter does not exist, we create it.
  if (pm==NULL)
    {
      // Declare by kind or ptype?
      try {
	if (type_flg)
	  {
	    pm = fuda->param_add(name,type_name);
	  }
	else
	  {
	    pm = fuda->param_add(name,kind);
	  }
      }
      catch (Uferr::ParamNameInvalid) {
	std::string err = "fuda_param: param name invalid";
	PyErr_SetString(PyExc_Exception, err.c_str());
	std::cout << "fuda_param: param_add failed\n";
	return NULL;
      }	
      catch (Uferr::PtypeNameInvalid) {
	std::string err = "fuda_param: ptype name (type) invalid";
	PyErr_SetString(PyExc_Exception, err.c_str());
	std::cout << "fuda_param: param_add failed\n";
	return NULL;
      }	
      catch (...) {
	std::string err = "fuda_param: unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	std::cout << "fuda_param: param_add failed\n";
	return NULL;
      }	
    }
  else
    {
      if (kind_flg)
	{
	  std::string err = "fuda_param: cannot change parameter kind ";
	  err += "of existing parameter";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;	  
	}
    }
  
  // Set optional arguments.
  if (value_flg)
    pm->set_value(value);

  if (free_flg) {
    if (free==0) pm->fix();
    else if (free==1) pm->free();
    else
      {
	std::string err = "fuda_param: free value invalid: ";
	snprintf(fstr, sizeof(fstr),"%d",free);
	err += fstr;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
  }

  if (bounds_flg) {
    if (bounds>=0 && bounds<=3) pm->set_bounds(bounds);
    else
      {
	std::string err = "fuda_param: bounds value is invalid: ";
	snprintf(fstr, sizeof(fstr),"%d",bounds);
	err += fstr;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
  }

  if (lower_flg)
    pm->set_lower(lower);

  if (upper_flg)
    pm->set_upper(upper);

  if (delta_flg) {
    if (delta > 0.0) pm->set_delta(delta);
    else
      {
	std::string err = "fuda_param: delta value invalid: ";
	snprintf(fstr, sizeof(fstr),"%f",delta);
	err += fstr;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
  }

  if (esd_flg)
    pm->set_esd(esd);

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_param_del(PyObject *self, PyObject *args)
{
  /* param_del deletes a parameter. */

  char *name;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  try {
    // Delete the param.
    fuda->param_del(name);
  }
  catch (Uferr::ParamNameInvalid& e) {
    std::string err =  "fuda_param_del: param name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamReferencedByFunc& e) {
    std::string err =  "fuda_param_del: param referenced by func: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_param_del: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_param_del_all(PyObject *self, PyObject *args)
{
  /* param_del_all deletes all parameters. */
  
  if (PyTuple_Size(args)>0)
    {
      std::string err =  "fuda_param_del_all: No arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  try {
    fuda->param_del_all();
  }
  catch (Uferr::ParamReferencedByFunc& e) {
    std::string err = 
      "fuda_param_del_all: param referenced by func: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_param_del_all: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_param_get(PyObject *self, PyObject *args)
{
  /* This routine returns either a single object or, if more entities
     are asked for, a tuple of objects. The first argument must be the
     parameter name and all remaining arguments must be strings
     specifying a keyword for an entity to return in the tuple. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg<2)
    {
      std::string err =  "fuda_param_get: at least 2 arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj)) // CB: PyUnicode_Check --> PyLong_Check update
    {
      std::string err =  "fuda_param_get: Second argument must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Find the parameter - if not found pm will be NULL.
  Param *pm = fuda->param_find(name);
  if (pm==NULL)
    {
      // Parameter not found.
      std::string err =  "fuda_param_get: parameter not found: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Create the output tuple.
  int nkey = narg-1;
  PyObject* ot_obj = PyTuple_New (nkey);

  // Loop over remaining arguments which are keywords.
  for (int iarg=1; iarg<narg; iarg++)
    {
      int ikey = iarg-1;

      // All keys must be of type string.
      PyObject *key_obj = PyTuple_GetItem(args, iarg);
      if (!PyUnicode_Check(key_obj))
	{
	  std::string err =  "fuda_param_get: ";
	  snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
          err += ". argument must be a keyword string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
      
      // Get the key.
      const char *key = PyUnicode_AsUTF8(key_obj);

      // Lookup key.
      if (strcmp(key,"name")==0)
	{
	  // Create and set a string opject.
	  std::string str;
	  pm->get_name(str);
	  PyObject* key_obj = PyUnicode_FromString(str.c_str());
 
	  // Insert in tup_obj.
	  if (PyTuple_SetItem(ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"value")==0)
	{
	  // Create and set a float opject.
	  double d = pm->get_value();
	  PyObject* key_obj = PyFloat_FromDouble(d);
 
	  // Insert in ot_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"init")==0)
	{
	  // Create and set a float opject.
	  double d = pm->get_init();
	  PyObject* key_obj = PyFloat_FromDouble(d);
 
	  // Insert in ot_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"esd")==0)
	{
	  // Create and set a float opject.
	  double d = pm->get_esd();
	  PyObject* key_obj = PyFloat_FromDouble(d);
 
	  // Insert in ot_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"lower")==0)
	{
	  // Create and set a float opject.
	  double d = pm->get_lower();
	  PyObject* key_obj = PyFloat_FromDouble(d);
 
	  // Insert in ot_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"upper")==0)
	{
	  // Create and set a float opject.
	  double d = pm->get_upper();
	  PyObject* key_obj = PyFloat_FromDouble(d);
 
	  // Insert in ot_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"delta")==0)
	{
	  // Create and set a float opject.
	  double d = pm->get_delta();
	  PyObject* key_obj = PyFloat_FromDouble(d);
 
	  // Insert in ot_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"free")==0)
	{
	  // Create and set an int opject.
	  int i = pm->is_free();
	  PyObject* key_obj = PyLong_FromLong (i);
 
	  // Insert in tup_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"bounds")==0)
	{
	  // Create and set an int opject.
	  int i = pm->get_bounds();
	  PyObject* key_obj = PyLong_FromLong (i);
 
	  // Insert in tup_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"eval_active")==0)
	{
	  // Create and set an int opject.
	  int i = pm->is_eval_active();
	  PyObject* key_obj = PyLong_FromLong (i);
 
	  // Insert in tup_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"eval_index")==0)
	{
	  // Create and set an int opject.
	  int i = pm->get_eval_index();
	  PyObject* key_obj = PyLong_FromLong (i);
 
	  // Insert in tup_obj.
	  if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"ptype")==0)
	{
	  // Create and set a string opject.
	  std::string str;
	  pm->get_ptype_name(str);
	  PyObject *key_obj;
	  if (str=="")
	    {
	      key_obj = Py_None;
	      Py_INCREF(Py_None);
	    }
	  else
	    key_obj = PyUnicode_FromString(str.c_str());
 
	  // Insert in tup_obj.
	  if (PyTuple_SetItem(ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else if (strcmp(key,"kind")==0)
	{
	  // Create and set a string opject.
	  std::string str;
	  FUDA::Param_kind kind = pm->get_kind();
	  if (kind==FUDA::PKIND_PARAM) str = "PARAM";
	  else if (kind==FUDA::PKIND_CONST) str = "CONST";
	  else if (kind==FUDA::PKIND_EXPL) str = "EXPL";
	  else str = "UNKNOWN";

	  PyObject* key_obj = PyUnicode_FromString(str.c_str());
 
	  // Insert in tup_obj.
	  if (PyTuple_SetItem(ot_obj, ikey, key_obj)!=0)
	    {
	      std::string err =  "fuda_param_get: PyTuple_SetItem failed on ";
	      snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	      err += ". argument";
	      PyErr_SetString(PyExc_Exception, err.c_str());
	      return NULL;
	    }
	}
      else
	{
	  // Invalid keyword.

	  // Deallocate the output tuple.
	  Py_DECREF(ot_obj);

	  // Report error and return.
	  std::string err =  "fuda_param_get: ";
	  snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	  err += ". argument is an invalid keyword: ";
	  err += key;
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
    }

  // If nkey==1, we don't wrap in a tuple.
  if (nkey==1)
    {
      // Extract the only element there is.
      PyObject *out_obj = PyTuple_GetItem(ot_obj,0);
      Py_INCREF(out_obj);

      // Deallocate tuple object which we don't return.
      Py_DECREF(ot_obj);

      // Return single object.
      return out_obj;
    }
  else
    {
      // Return output tuple.
      return ot_obj;
    }
}


extern "C" PyObject *
fuda_param_init_value(PyObject *self, PyObject *args)
{
  /* This routine initializes the value of a parameter to it initial
     value by a call to the Param init_value() method. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=1)
    {
      std::string err =  "fuda_param_init_value: exactly 1 argument required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_param_init_value: Second argument must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Find the parameter - if not found pm will be NULL.
  Param *pm = fuda->param_find(name);
  if (pm==NULL)
    {
      // Parameter not found.
      std::string err =  "fuda_param_init_value: parameter not found: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Initialize the value.
  pm->init_value();

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_param_get_all(PyObject *self, PyObject *args)
{
  /* This routine returns a tuple with all parameters defined in
     fuda. The routine takes no arguments */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=0)
    {
      std::string err =  "fuda_param_get_all: no arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Create the output tuple.
  int size = fuda->get_nparam();
  PyObject* tuple_obj = PyTuple_New (size);

  // Loop over parameters.
  int ituple = 0;
  for (Param_iterator pi=fuda->param_begin();
       pi!=fuda->param_end(); pi++)
    {
      // Get parameter name.
      std::string name;
      (*pi)->get_name(name);

      // Make a string object of the name.
      PyObject* obj = PyUnicode_FromString(name.c_str());
      
      // Insert in tuple.
      if (PyTuple_SetItem(tuple_obj, ituple, obj)!=0)
	{
	  // Deallocate tuple object and abort.
	  Py_DECREF(tuple_obj);
	  return NULL;
	}

      ituple++;
    }

  // Return output tuple.
  return tuple_obj;
}


extern "C" PyObject *
fuda_ftype_product(PyObject *self, PyObject *args)
{
  /* fuda_ftype_product declares a new product function type. The
     required arguments are: ftype_name, scale_flag and ftype_seq. The
     ftype_seq is a sequence with the ftypes used to generate the new
     ftype. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=3)
    {
      std::string err =  "fuda_ftype_product: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 3)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_ftype_product: 1. arg (name) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // 2nd argument must be the scale flag.
  PyObject *scale_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(scale_obj))
    {
      std::string err =  "fuda_ftype_product: 2. arg (scale_flg) must be an integer";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the scale flag.
  bool scale_flg = PyLong_AsLong(scale_obj);

  /* 3rd argument must be a sequence with the names of the
     ftypes in the product ftype. */
  PyObject *ftypes_obj = PyTuple_GetItem(args,2);
  if (!PySequence_Check(ftypes_obj))
    {
      std::string err =  "fuda_ftype_product: 3. arg (ftypes) must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be at least one item in the ftype sequence.
  if (PySequence_Length(ftypes_obj)<1)
    {
      std::string err = "fuda_ftype_product: 3. arg must be a sequence ";
      err += "with one or more ftype names";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Declare and dimension ftypes vector with size of ftype_obj.
  std::vector<Ftype*> ftvec;
  ftvec.resize(PySequence_Length(ftypes_obj));

  // Loop over the items of the ftype sequence and get ftypes.
  for(int i=0; i<PyTuple_Size(ftypes_obj); i++)
    {
      // Get item from sequence.
      PyObject *ftype_name_obj = PySequence_GetItem(ftypes_obj,i);

      // Check it is a string.
      if (!PyUnicode_Check(ftype_name_obj))
	{
	  std::string err =  "fuda_ftype_product: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 3. arg must be of type string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the ftype name.
      const char *ftype_name = PyUnicode_AsUTF8(ftype_name_obj);

      // Look for an ftype with that name.
      Ftype *ft = ftype_find(ftype_name);

      // If we did not find any ftype by that name, we abort.
      if (ft==0)
	{
	  std::string err =  "fuda_ftype_product: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 3rd arg. is not a valid ftype: ";
	  err += ftype_name;
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Add the ftype reference to the ftype vector.
      ftvec[i] = ft;
    }

  // Create the ftype.
  Ftype *ftype;
  try {
    ftype = fuda->ftype_add_prod(name, scale_flg, ftvec);
  }
  catch (Uferr::FtypeNameInvalid& e) {
    std::string err =  "fuda_ftype_product: name invalid: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FtypeNameAlreadyUsed& e) {
    std::string err =  "fuda_ftype_product: name already used: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FtypeEmptyFtypeVec) {
    std::string err =  "fuda_ftype_product: empty ftype vec";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_ftype_product: failed";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_ftype_sum(PyObject *self, PyObject *args)
{
  /* fuda_ftype_sum declares a new sum function type. The
     required arguments are: ftype_name, scale_flag and ftype_seq. The
     ftype_seq is a sequence with the ftypes used to generate the new
     ftype. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=3)
    {
      std::string err =  "fuda_ftype_sum: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 3)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_ftype_sum: 1. arg (name) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // 2nd argument must be the scale flag.
  PyObject *scale_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(scale_obj))
    {
      std::string err =  "fuda_ftype_sum: 2. arg (scale_flg) must be an integer";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the scale flag.
  bool scale_flg = PyLong_AsLong(scale_obj);

  /* 3rd argument must be a sequence with the names of the
     ftypes in the product ftype. */
  PyObject *ftypes_obj = PyTuple_GetItem(args,2);
  if (!PySequence_Check(ftypes_obj))
    {
      std::string err =  "fuda_ftype_sum: 3. arg (ftypes) must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be at least one item in the ftype sequence.
  if (PySequence_Length(ftypes_obj)<1)
    {
      std::string err = "fuda_ftype_sum: 3. arg must be a sequence ";
      err += "with one or more ftype names";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Declare and dimension ftypes vector with size of ftype_obj.
  std::vector<Ftype*> ftvec;
  ftvec.resize(PySequence_Length(ftypes_obj));

  // Loop over the items of the ftype sequence and get ftypes.
  for(int i=0; i<PyTuple_Size(ftypes_obj); i++)
    {
      // Get item from sequence.
      PyObject *ftype_name_obj = PySequence_GetItem(ftypes_obj,i);

      // Check it is a string.
      if (!PyUnicode_Check(ftype_name_obj))
	{
	  std::string err =  "fuda_ftype_sum: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 3. arg must be of type string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the ftype name.
      const char *ftype_name = PyUnicode_AsUTF8(ftype_name_obj);

      // Look for an ftype with that name.
      Ftype *ft = ftype_find(ftype_name);

      // If we did not find any ftype by that name, we abort.
      if (ft==0)
	{
	  std::string err =  "fuda_ftype_sum: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 3rd arg. is not a valid ftype: ";
	  err += ftype_name;
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Add the ftype reference to the ftype vector.
      ftvec[i] = ft;
    }

  // Create the ftype.
  Ftype *ftype;
  try {
    ftype = fuda->ftype_add_sum(name, scale_flg, ftvec);
  }
  catch (Uferr::FtypeNameInvalid& e) {
    std::string err =  "fuda_ftype_sum: name invalid: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FtypeNameAlreadyUsed& e) {
    std::string err =  "fuda_ftype_sum: name already used: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FtypeEmptyFtypeVec) {
    std::string err =  "fuda_ftype_sum: empty ftype vec";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_ftype_sum: failed";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_ftype_composite(PyObject *self, PyObject *args)
{
  /* fuda_ftype_composite declares a new composite function type. The
     required arguments are: ftype_name, scale_flag, p_index, ftype_name_f,
     and ftype_name_g */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=5)
    {
      std::string err =  "fuda_ftype_composite: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 5)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_ftype_composite: 1. arg (name) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // 2nd argument must be the scale flag.
  PyObject *scale_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(scale_obj))
    {
      std::string err =  "fuda_ftype_composite: 2. arg (scale_flg) must be an integer";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the scale flag.
  bool scale_flg = PyLong_AsLong(scale_obj);

  // 3rd argument must be the p_index.
  PyObject *p_index_obj = PyTuple_GetItem(args, 2);
  if (!PyLong_Check(p_index_obj))
    {
      std::string err =  "fuda_ftype_composite: 3. arg (p_index) must be an integer";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the p_index.
  int p_index = PyLong_AsLong(p_index_obj);

  // 4th argument must be the name of f function ftype.
  PyObject *fname_obj = PyTuple_GetItem(args, 3);
  if (!PyUnicode_Check(fname_obj))
    {
      std::string err =  "fuda_ftype_composite: 4. arg (fname) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the name.
  const char *fname = PyUnicode_AsUTF8(fname_obj);

  // Look for an ftype with that name.
  Ftype *f = ftype_find(fname);

  // If we did not find any ftype by that name, we abort.
  if (f==0)
    {
      std::string err =  "fuda_ftype_composite: ";
      err += fname;
      err += " is not a valid ftype";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 5th argument must be the name of g function ftype.
  PyObject *gname_obj = PyTuple_GetItem(args, 4);
  if (!PyUnicode_Check(gname_obj))
    {
      std::string err =  "fuda_ftype_composite: 5. arg (gname) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the name.
  const char *gname = PyUnicode_AsUTF8(gname_obj);

  // Look for an ftype with that name.
  Ftype *g = ftype_find(gname);

  // If we did not find any ftype by that name, we abort.
  if (g==0)
    {
      std::string err =  "fuda_ftype_composite: ";
      err += gname;
      err += " is not a valid ftype";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // p_index must not exeede the number of parameters in f.
  if (p_index<0 || p_index >= (int) f->get_nparam())
    {
      std::string err =  "fuda_ftype_composite: 3. arg (p_index) is invalid";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  Ftype *ftype;
  try {
    ftype = fuda->ftype_add_comp(name, scale_flg, p_index, f, g);
  }
  catch (Uferr::FtypeNameInvalid& e) {
    std::string err =  "fuda_ftype_composite: name invalid: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FtypeNameAlreadyUsed& e) {
    std::string err =  "fuda_ftype_composite: name already used: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FtypeEmptyFtypeVec) {
    std::string err =  "fuda_ftype_composite: empty ftype vec";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamIndexInvalid& e) {
    std::string err =  "fuda_ftype_composite: parameter index invalid:";
    err += e.i;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_ftype_composite: failed";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_ftype_get_all(PyObject *self, PyObject *args)
{
  /* This routine returns a tuple with all defined ftypes declared in
     uf. The routine takes no arguments */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=0)
    {
      std::string err =  "fuda_ftype_get_all: no arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Create the output tuple.
  int size = fuda->get_nftype();
  PyObject* tuple_obj = PyTuple_New (size);

  // Loop over ftypes.
  int ituple = 0;
  for (Ftype_iterator fi=fuda->ftype_begin();
       fi!=fuda->ftype_end(); fi++)
    {
      // Get ftype name.
      std::string name;
      (*fi)->get_name(name);

      // Make a string object of the name.
      PyObject* obj = PyUnicode_FromString(name.c_str());
      
      // Insert in tuple.
      if (PyTuple_SetItem(tuple_obj, ituple, obj)!=0)
	{
	  // Deallocate tuple object and abort.
	  Py_DECREF(tuple_obj);
	  return NULL;
	}
      ituple++;
    }

  // Return output tuple.
  return tuple_obj;
}


extern "C" PyObject *
fuda_ftype_exists(PyObject *self, PyObject *args)
{
  /* ftype_exists returns true (1) if the specified ftype exists. */

  char *name;
  int exists_flg;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Get var flag for specified parameter.
  if (ftype_find(name)==NULL)
    exists_flg = 0;
  else
    exists_flg = 1;

  // Return value.
  return Py_BuildValue("i", exists_flg);
}


extern "C" PyObject *
fuda_func_exists(PyObject *self, PyObject *args)
{
  /* func_exists returns true (1) if the specified func exists. */

  char *name;
  int exists_flg;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Get var flag for specified parameter.
  if (func_find(name)==NULL)
    exists_flg = 0;
  else
    exists_flg = 1;

  // Return value.
  return Py_BuildValue("i", exists_flg);
}


extern "C" PyObject *
fuda_func_del(PyObject *self, PyObject *args)
{
  /* func_del deletes a function. */

  char *name;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  try {
    // Delete the func.
    fuda->func_del(name);

    // If we deleted the current function, we zero the current function.
    if (name==current_func_name)
      {
	current_func = NULL;
	current_func_name = "";
      }    

    // If we deleted the last found function, we zero it.
    if (name==last_func_name)
      {
	last_func = NULL;
	last_func_name = "";
      }    
  }
  catch (Uferr::FuncNameInvalid& e) {
    std::string err =  "fuda_func_del: func name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_func_del: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_func_del_all(PyObject *self, PyObject *args)
{
  /* func_del deletes all functions. */
  
  if (PyTuple_Size(args)>0)
    {
      std::string err =  "fuda_func_del_all: No arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  try {
    fuda->func_del_all();
  }
  catch (...) {
    std::string err = 
      "fuda_func_del_all: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Forget about current function and last found function.
  current_func = NULL;
  current_func_name = "";
  last_func = NULL;
  last_func_name = "";

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_param_exists(PyObject *self, PyObject *args)
{
  /* param_exists returns true (1) if the specified param exists. */

  char *name;
  int exists_flg;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Get var flag for specified parameter.
  if (fuda->param_find(name)==NULL)
    exists_flg = 0;
  else
    exists_flg = 1;

  // Return value.
  return Py_BuildValue("i", exists_flg);
}


extern "C" PyObject *
fuda_param_is_referenced(PyObject *self, PyObject *args)
{
  /* param_is_referenced returns true (1) if the specified param is
     referenced by any function. */

  char *name;
  int referenced_flg;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Get referenced flag for specified parameter.
  try {
    Param* pm = fuda->param_find(name);
    if (pm==NULL)
      {
	std::string err = 
	  "fuda_param_is_referenced: param name invalid: ";
	err += name;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
    referenced_flg = pm->is_referenced();
  }
  catch (...) {
    std::string err = 
      "fuda_param_is_referenced: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  // Return value.
  return Py_BuildValue("i", referenced_flg);
}


extern "C" PyObject *
fuda_dtype_exists(PyObject *self, PyObject *args)
{
  /* dtype_exists returns true (1) if the specified dtype exists. */

  char *name;
  int exists_flg;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Get var flag for specified parameter.
  if (fuda->dtype_find(name)==NULL)
    exists_flg = 0;
  else
    exists_flg = 1;

  // Return value.
  return Py_BuildValue("i", exists_flg);
}


extern "C" PyObject *
fuda_dtype_is_referenced(PyObject *self, PyObject *args)
{
  /* dtype_is_referenced returns true (1) if the specified dtype is
     referenced by any function. */

  char *name;
  int referenced_flg;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Get referenced flag for specified dtype.
  try {
    Dtype* dt = fuda->dtype_find(name);
    if (dt == NULL)
      {
	std::string err = 
	  "fuda_dtype_is_referenced: dtype name invalid: ";
	err += name;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }    
    referenced_flg = dt->is_referenced();
  }
  catch (...) {
    std::string err = 
      "fuda_dtype_is_referenced: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  // Return value.
  return Py_BuildValue("i", referenced_flg);
}


extern "C" PyObject *
fuda_ptype_is_referenced(PyObject *self, PyObject *args)
{
  /* ptype_is_referenced returns true (1) if the specified ptype is
     referenced by any parameter. */

  char *name;
  int referenced_flg;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Get referenced flag for specified ptype.
  try {
    Ptype* pt = fuda->ptype_find(name);
    if (pt == NULL)
      {
	std::string err = 
	  "fuda_ptype_is_referenced: ptype name invalid: ";
	err += name;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }    
    referenced_flg = pt->is_referenced();
  }
  catch (...) {
    std::string err = 
      "fuda_ptype_is_referenced: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  // Return value.
  return Py_BuildValue("i", referenced_flg);
}


extern "C" PyObject *
fuda_ptype_exists(PyObject *self, PyObject *args)
{
  /* ptype_exists returns true (1) if the specified ptype exists. */

  char *name;
  int exists_flg;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Get var flag for specified parameter.
  if (fuda->ptype_find(name)==NULL)
    exists_flg = 0;
  else
    exists_flg = 1;

  // Return value.
  return Py_BuildValue("i", exists_flg);
}


extern "C" PyObject *
fuda_ptype_del(PyObject *self, PyObject *args)
{
  /* ptype_del deletes a ptype. */

  char *name;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  try {
    // Delete the ptype.
    fuda->ptype_del(name);
  }
  catch (Uferr::PtypeNameInvalid& e) {
    std::string err =  "fuda_ptype_del: ptype name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::PtypeReferencedByParam& e) {
    std::string err =  "fuda_ptype_del: ptype referenced by param: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_ptype_del: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_ptype_del_all(PyObject *self, PyObject *args)
{
  /* ptype_del_all deletes all parameters. */
  
  if (PyTuple_Size(args)>0)
    {
      std::string err =  "fuda_ptype_del_all: No arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  try {
    fuda->ptype_del_all();
  }
  catch (Uferr::PtypeReferencedByParam& e) {
    std::string err = 
      "fuda_ptype_del_all: ptype referenced by param: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_ptype_del_all: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_ftype_call(PyObject *self, PyObject *args)
{
  /* fuda_ftype_call call the specified function with the
     specified parameter values. It takes two arguments: ftype_name and
     param_value_seq, where param_value_seq is a sequence of double
     precision numbers. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=2)
    {
      std::string err =  "fuda_ftype_call: invalid number of args: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 2)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the ftype name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_ftype_call: 2. arg. must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Lookup ftype.
  Ftype *ft = ftype_find(name);
  if (ft==NULL)
    {
      std::string err =  "fuda_ftype_call: Invalid ftype name: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get number of parameters in ftype.
  int nparam = ft->get_nparam();

  // Second argument must be a sequence with the values of the parameters.
  PyObject *varseq_obj = PyTuple_GetItem(args,1);
  if (!PySequence_Check(varseq_obj))
    {
      std::string err =  "fuda_ftype_call: 2. arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be nparam items in the varseq.
  if (PySequence_Length(varseq_obj)!=nparam)
    {
      std::string err = "fuda_ftype_call: 2. arg. sequence size invalid: ";
      snprintf(fstr, sizeof(fstr), "%zd", PySequence_Length(varseq_obj));
      err += fstr;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the varseq and setup vector with values.
  std::vector<double> value_vector(nparam);
  for(int i=0; i<nparam; i++)
    {
      // Get item from tuple.
      PyObject *varval_obj = PySequence_GetItem(varseq_obj,i);

      // Check it is a float.
      if (!PyFloat_Check(varval_obj))
	{
	  std::string err =  "fuda_ftype_call: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg. must be of type float";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      value_vector[i] = PyFloat_AsDouble(varval_obj);
    }

  // Evaluate function.
  double fval;
  try {
    fval = ft->call(value_vector);
  }
  catch (Uferr::NumParamInvalid& e) {
    std::string err =
      std::string("fuda_ftype_call: Invalid number of parameters: ") + std::to_string(e.n);
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_ftype_call: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("d", fval);
}



extern "C" PyObject *
fuda_ptype_get_all(PyObject *self, PyObject *args)
{
  /* This routine returns a tuple with all defined ptypes in uf. The
     routine takes no arguments. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=0)
    {
      std::string err =  "fuda_ptype_get_all: no arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Create the output tuple.
  int size = fuda->get_nptype();
  PyObject* tuple_obj = PyTuple_New (size);

  // Loop over ptypes.
  int ituple = 0;
  for (Ptype_iterator pti=fuda->ptype_begin();
       pti!=fuda->ptype_end(); pti++)
    {
      // Get ptype name.
      std::string name;
      (*pti)->get_name(name);

      // Make a string object of the name.
      PyObject* obj = PyUnicode_FromString(name.c_str());
      
      // Insert in tuple.
      if (PyTuple_SetItem(tuple_obj, ituple, obj)!=0)
	{
	  // Deallocate tuple object and abort.
	  Py_DECREF(tuple_obj);
	  return NULL;
	}
      ituple++;
    }

  // Return output tuple.
  return tuple_obj;
}


extern "C" PyObject *
fuda_ftype_get_param(PyObject *self, PyObject *args)
{
  /* ftype_get_param returns the i'th parameter name for a specified
     ftype. It takes two arguments: the ftype name and the index (i)
     of the parameter for which to return the name of. */

  char *name;
  int index;
  
  if (!PyArg_ParseTuple(args, "si", &name, &index))
    return NULL;

  // Get param name.
  std::string pname;
  try {
    // find ftype.
    Ftype *ft = ftype_find(name);
    if (ft==NULL)
      {
	std::string err =  "fuda_ftype_get_param: Invalid ftype";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Get parameter name.
    ft->get_p_name(index,pname);
  }
  catch (Uferr::ParamIndexInvalid& e) {
    std::string err =  "fuda_ftype_get_param: Invalid param index";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_param: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", pname.c_str());
}


extern "C" PyObject *
fuda_ftype_get_param_descr(PyObject *self, PyObject *args)
{
  /* ftype_get_param_descr returns the i'th parameter descriptor
     string for a specified ftype. It takes two arguments: the ftype
     name and the index (i) of the parameter for which to return the
     descriptor string of. */

  char *name;
  int index;
  
  if (!PyArg_ParseTuple(args, "si", &name, &index))
    return NULL;

  // Get param descriptor string.
  std::string descr;
  try {
    // find ftype.
    Ftype *ft = ftype_find(name);
    if (ft==NULL)
      {
	std::string err =  "fuda_ftype_get_param_descr: Invalid ftype";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Get descriptor.
    ft->get_p_descr(index,descr);
  }
  catch (Uferr::ParamIndexInvalid& e) {
    std::string err =  "fuda_ftype_get_param_descr: Invalid param index";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_param_descr: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", descr.c_str());
}


extern "C" PyObject *
fuda_ftype_get_var(PyObject *self, PyObject *args)
{
  /* ftype_get_var returns true (1) if the i'th parameter for a
     specified ftype is a variable. It takes two arguments: the ftype
     name and the index (i) of the parameter for which to return the
     variable status of. */

  char *name;
  int index;
  
  if (!PyArg_ParseTuple(args, "si", &name, &index))
    return NULL;

  // Get var flag for specified parameter.
  int var_flg;
  try {
    // find ftype.
    Ftype *ft = ftype_find(name);
    if (ft==NULL)
      {
	std::string err =  "fuda_ftype_get_var: Invalid ftype";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Get var flag.
    var_flg = ft->is_p_var(index);
  }
  catch (Uferr::ParamIndexInvalid& e) {
    std::string err =  "fuda_ftype_get_var: Invalid param index";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_var: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("i", var_flg);
}


extern "C" PyObject *
fuda_ftype_get_deriv(PyObject *self, PyObject *args)
{
  /* ftype_get_deriv returns true (1) if the i'th parameter for a
     specified ftype support derivative calculation. It takes two
     arguments: the ftype name and the index (i) of the parameter for
     which to return the derivative status of. */

  char *name;
  int index;
  
  if (!PyArg_ParseTuple(args, "si", &name, &index))
    return NULL;

  // Get var flag for specified parameter.
  int deriv_flg;
  try {
    // find ftype.
    Ftype *ft = ftype_find(name);
    if (ft==NULL)
      {
	std::string err =  "fuda_ftype_get_deriv: Invalid ftype";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Get var flag.
    deriv_flg = ft->is_p_deriv(index);
  }
  catch (Uferr::ParamIndexInvalid& e) {
    std::string err =  "fuda_ftype_get_var: Invalid param index";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_var: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("i", deriv_flg);
}


extern "C" PyObject *
fuda_ftype_get_var_index(PyObject *self, PyObject *args)
{
  /* ftype_get_var_index returns the parameter index for the i'th
     variable for the specified ftype. It takes two arguments: the
     ftype name and the variable index for which to return the
     corresponding parameter index. */

  char *name;
  int index;
  
  if (!PyArg_ParseTuple(args, "si", &name, &index))
    return NULL;

  // Get parameter index for specified variable.
  int var_index;
  try {
    // find ftype.
    Ftype *ft = ftype_find(name);
    if (ft==NULL)
      {
	std::string err =  "fuda_ftype_get_var_index: Invalid ftype";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Get var flag.
    var_index = ft->get_var_index(index);
  }
  catch (Uferr::VarIndexInvalid& e) {
    std::string err =  "fuda_ftype_get_var_index: Invalid param index";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_var_index: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("i", var_index);
}


extern "C" PyObject *
fuda_dtype_get_all(PyObject *self, PyObject *args)
{
  /* This routine returns a tuple with all defined dtypes. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=0)
    {
      std::string err =  "fuda_dtype_get_all: no arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Create the output tuple.
  int size = fuda->get_ndtype();
  PyObject* tuple_obj = PyTuple_New (size);

  // Loop over parameters.
  int ituple = 0;
  for (Dtype_iterator dti=fuda->dtype_begin();
       dti!=fuda->dtype_end(); dti++)
    {
      // Get parameter name.
      std::string name;
      (*dti)->get_name(name);

      // Make a string object of the name.
      PyObject* obj = PyUnicode_FromString(name.c_str());
      
      // Insert in tuple.
      if (PyTuple_SetItem(tuple_obj, ituple, obj)!=0)
	{
	  // Deallocate tuple object and abort.
	  Py_DECREF(tuple_obj);
	  return NULL;
	}

      ituple++;
    }

  // Return output tuple.
  return tuple_obj;
}


extern "C" PyObject *
fuda_func_get_all(PyObject *self, PyObject *args)
{
  /* This routine returns a tuple with all defined functions. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=0)
    {
      std::string err =  "fuda_func_get_all: no arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Create the output tuple.
  int size = fuda->get_nfunc();
  PyObject* tuple_obj = PyTuple_New (size);

  // Loop over functions.
  int ituple = 0;
  for (Func_iterator fi=fuda->func_begin();
       fi!=fuda->func_end(); fi++)
    {
      // Get parameter name.
      std::string name;
      (*fi)->get_name(name);

      // Make a string object of the name.
      PyObject* obj = PyUnicode_FromString(name.c_str());
      
      // Insert in tuple.
      if (PyTuple_SetItem(tuple_obj, ituple, obj)!=0)
	{
	  // Deallocate tuple object and abort.
	  Py_DECREF(tuple_obj);
	  return NULL;
	}

      ituple++;
    }

  // Return output tuple.
  return tuple_obj;
}


extern "C" PyObject *
fuda_func_use(PyObject *self, PyObject *args)
{
  /* fuda_func_use sets or clears the use status of a function. It takes
     two arguments, the name of the function and the use status which
     can be 0 or 1. */

  char *name;
  int use;
  
  if (!PyArg_ParseTuple(args, "si", &name, &use))
    return NULL;

  // Set use.
  std::string pname;
  try {
    // find func.
    Func *fn = func_find(name);
    if (fn==NULL)
      {
	std::string err =  "fuda_func_use: Invalid func";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Set use.
    fn->set_use(use);
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_param: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_func_use_all(PyObject *self, PyObject *args)
{
  /* fuda_func_use_all sets or clears the use status of a function. It
     a single argument, the use status, which can be 0 or 1. */

  int use;
  
  if (!PyArg_ParseTuple(args, "i", &use))
    return NULL;

  // Set use.
  try {
    // Loop over functions.
    for (Func_iterator fi=fuda->func_begin();
	 fi!=fuda->func_end(); fi++)
      {
	// Set use flag.
	std::string name;
	(*fi)->set_use(use);
      }
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_param: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_func_set_param(PyObject *self, PyObject *args)
{
  /* fuda_func_set_param sets i'th parameter in the specified function
     to a partucular parameter and so changes the parameters
     associated with a given function. It takes three arguemts:
     func_name, param_index and param_name. */

  char *name, *pname;
  int pindex;
  
  if (!PyArg_ParseTuple(args, "sis", &name, &pindex, &pname))
    return NULL;

  try {
    fuda->func_set_param(name, pindex, pname);
  }
  catch (Uferr::FuncNameInvalid& e) {
    std::string err =  "fuda_func_set_param: func name invalid: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamNameInvalid& e) {
    std::string err =  "fuda_func_set_param: param name invalid: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ExplVarInvalid& e) {
    std::string err =  "fuda_func_set_param: explanatory variable ";
    err += "not member of dtype: ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamKindInvalid) {
    std::string err =  "fuda_func_set_param: param kind invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamIndexInvalid& e) {
    std::string err =  "fuda_ftype_get_param: Invalid param index : ";
    snprintf(fstr, sizeof(fstr), "%d", e.i); err += fstr;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_func_set_param: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_func_get_param(PyObject *self, PyObject *args)
{
  /* fuda_func_get_param returns the name of the i'th parameter for a
     given function. It takes two arguemnts: func_name and
     param_index. */

  char *name;
  int index;
  
  if (!PyArg_ParseTuple(args, "si", &name, &index))
    return NULL;

  // Get param name.
  std::string pname;
  try {
    // find func.
    Func *fn = func_find(name);
    if (fn==NULL)
      {
	std::string err =  "fuda_func_get_param: Invalid func";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Get parameter name.
    Param *pm = fn->get_param(index);
    pm->get_name(pname);
  }
  catch (Uferr::ParamIndexInvalid& e) {
    std::string err =  "fuda_ftype_get_param: Invalid param index";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_param: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", pname.c_str());
}


extern "C" PyObject *
fuda_func_get_var(PyObject *self, PyObject *args)
{
  /* fuda_func_get_var returns the name of the i'th variable for a
     given function. It takes two arguemnts: func_name and
     var_index. */

  char *name;
  int index;
  
  if (!PyArg_ParseTuple(args, "si", &name, &index))
    return NULL;

  // Get param name.
  std::string pname;
  try {
    // find func.
    Func *fn = func_find(name);
    if (fn==NULL)
      {
	std::string err =  "fuda_func_get_param: Invalid func";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Get parameter name.
    Param *pm = fn->get_var(index);
    pm->get_name(pname);
  }
  catch (Uferr::VarIndexInvalid& e) {
    std::string err =  "fuda_ftype_get_param: Invalid param index";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_ftype_get_param: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", pname.c_str());
}



extern "C" PyObject *
fuda_dtype(PyObject *self, PyObject *args)
{
  /* fuda_dtype declare a new dtype. It takes two arguments: dtype_name
     and expl_seq, where expl_seq is a tuple with the explanatory
     variables for this dtype. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=2)
    {
      std::string err =  "fuda_dtype: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 2)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_dtype: Second argument must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  /* Second argument must be a sequence with the names of the
     explanatory variables. */
  PyObject *xvar_obj = PyTuple_GetItem(args,1);
  if (!PySequence_Check(xvar_obj))
    {
      std::string err =  "fuda_dtype: Second argument must be a seqence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the xvar sequence and get parameter names.
  std::list<std::string> p_names;
  for(int i=0; i<PySequence_Length(xvar_obj); i++)
    {
      // Get item from tuple.
      PyObject *xvar_name_obj = PySequence_GetItem(xvar_obj,i);

      // Check it is a string.
      if (!PyUnicode_Check(xvar_name_obj))
	{
	  std::string err =  "fuda_dtype: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of second argument must be of type string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the parameter name.
      const char *xvar_name = PyUnicode_AsUTF8(xvar_name_obj);
      p_names.push_back(xvar_name);
    }

  // Create the dtype.
  try {
    fuda->dtype_add(name, p_names);
  }
  catch (Uferr::DtypeNameAlreadyUsed& e) {
    std::string err = 
      "fuda_dtype: dtype name already used: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::DtypeNameInvalid& e) {
    std::string err = 
      "fuda_dtype: dtype name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamNameInvalid& e) {
    std::string err = 
      "fuda_dtype: parameter name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamNotExplanatory& e) {
    std::string err = 
      "fuda_dtype: parameter is not explanatory: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_dtype: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}



extern "C" PyObject *
fuda_dtype_del(PyObject *self, PyObject *args)
{
  /* dtype_del deletes a parameter. */

  char *name;
  
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  try {
    // Delete the dtype.
    fuda->dtype_del(name);
  }
  catch (Uferr::DtypeNameInvalid& e) {
    std::string err =  "fuda_dtype_del: dtype name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::DtypeReferencedByFunc& e) {
    std::string err =  "fuda_dtype_del: dtype referenced by func: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_dtype_del: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_dtype_del_all(PyObject *self, PyObject *args)
{
  /* dtype_del_all deletes all parameters. */
  
  if (PyTuple_Size(args)>0)
    {
      std::string err =  "fuda_dtype_del_all: No arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  try {
    fuda->dtype_del_all();
  }
  catch (Uferr::DtypeReferencedByFunc& e) {
    std::string err = 
      "fuda_dtype_del_all: dtype referenced by func: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_dtype_del_all: failed - unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_dtype_get(PyObject *self, PyObject *args)
{
  /* fuda_dtype_get returns either a single object or, if more entities
     are asked for, a tuple of objects. Arguments must be strings
     specifying an entity to return in the tuple. The first argument
     must be the name of the dtype to inquire and the remainder of the
     arguments must be keyvals fo entities to return. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg<2)
    {
      std::string err =  "fuda_dtype_get: at least 2 argument required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_dtype_get: First argument must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Find the dtype - if not found ft will be NULL.
  Dtype *dt = fuda->dtype_find(name);
  if (dt==NULL)
    {
      // Ftype not found.
      std::string err =  "fuda_dtype_get: dtype not found: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Create the output tuple.
  int nkey = narg-1;
  PyObject* ot_obj = PyTuple_New (nkey);

  // Loop over arguments which are keywords.
  for (int iarg=1; iarg<narg; iarg++)
    {
      int ikey = iarg-1;

      // All keys must be of type string.
      PyObject *key_obj = PyTuple_GetItem(args, iarg);
      if (!PyUnicode_Check(key_obj))
	{
	  std::string err =  "fuda_dtype_get: ";
	  snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
          err += ". argument must be a keyword string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
      
      // Get the key.
      const char *key = PyUnicode_AsUTF8(key_obj);

      // Lookup key.
      try {
	if (strcmp(key,"dim")==0)
	  {
	    // Create and set an int opject.
	    int i = dt->get_dim();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_dtype_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"purge")==0)
	  {
	    // Create and set an int opject.
	    int i = dt->is_purge();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_dtype_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else
	  {
	    // Invalid keyword.
	    
	    // Deallocate the output tuple.
	    Py_DECREF(ot_obj);
	    
	    // Report error and return.
	    std::string err =  "fuda_dtype_get: ";
	    snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	    err += ". argument is an invalid keyword: ";
	    err += key;
	    PyErr_SetString(PyExc_Exception, err.c_str());
	    return NULL;
	  }
      }
      catch (...) {
	std::string err =  "fuda_dtype_get: unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }      
    }

  // If nkey==1, we don't wrap in a tuple.
  if (nkey==1)
    {
      // Extract the only element there is.
      PyObject *out_obj = PyTuple_GetItem(ot_obj,0);
      Py_INCREF(out_obj);

      // Deallocate tuple object which we don't return.
      Py_DECREF(ot_obj);

      // Return single object.
      return out_obj;
    }
  else
    {
      // Return output tuple.
      return ot_obj;
    }
}


extern "C" PyObject *
fuda_dtype_set_purge(PyObject *self, PyObject *args)
{
  /* fuda_dtype_set_purge sets or clears the purge flag of a dtype. It
     takes two arguments, the name of the dtype and the purge flag
     which can be 0 or 1. */

  char *name;
  int purge_flg;
  
  if (!PyArg_ParseTuple(args, "si", &name, &purge_flg))
    return NULL;

  // Set purge_flg.
  try {
    // find func.
    Dtype *dt = fuda->dtype_find(name);
    if (dt==NULL)
      {
	std::string err =  "fuda_dtype_set_purge: Invalid dtype name";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }

    // Set purge flag.
    dt->set_purge(purge_flg);
  }
  catch (...) {
      std::string err =  "fuda_dtype_set_purge: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_dtype_set_purge_radius(PyObject *self, PyObject *args)
{
  /* fuda_dtype_set_purge_radius sets i'th purge radius in the
     specified dtype. It takes three arguemts: dtype_name, param_index
     and purge_radius. */

  char *name;
  int pindex;
  double purge_radius;
  
  if (!PyArg_ParseTuple(args, "sid", &name, &pindex, &purge_radius))
    return NULL;

  // find dtype.
  Dtype *dt = fuda->dtype_find(name);
  if (dt==NULL)
    {
      std::string err =  "fuda_dtype_set_purge_radius: Invalid dtype name";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  try {
    dt->set_purge_radius(pindex, purge_radius);
  }
  catch (Uferr::DtypeExplIndexInvalid& e) {
    std::string err =  "fuda_dtype_set_purge_radius: Invalid param index : ";
    snprintf(fstr, sizeof(fstr), "%d", e.i); err += fstr;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_dtype_set_purge_radius: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}



extern "C" PyObject *
fuda_dtype_get_purge_radius(PyObject *self, PyObject *args)
{
  /* fuda_dtype_get_param gets the i'th purge radius in the specified
     dtype. It takes 2 arguemts: dtype_name and  param_index and it returns
     the purge radius for the specified parameter index. */

  char *name;
  int pindex;
  double purge_radius;
  
  if (!PyArg_ParseTuple(args, "si", &name, &pindex))
    return NULL;

  // find dtype.
  Dtype *dt = fuda->dtype_find(name);
  if (dt==NULL)
    {
      std::string err =  "fuda_dtype_get_purge_radius: Invalid dtype name";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  try {
    purge_radius = dt->get_purge_radius(pindex);
  }
  catch (Uferr::DtypeExplIndexInvalid& e) {
    std::string err =  "fuda_dtype_get_purge_radius: Invalid param index : ";
    snprintf(fstr, sizeof(fstr), "%d", e.i); err += fstr;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "fuda_dtype_get_purge_radius: Unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return purge radius.
  return Py_BuildValue("d", purge_radius);
}



extern "C" PyObject *
fuda_ftype_get(PyObject *self, PyObject *args)
{
  /* fuda_ftype_get returns either a single object or, if more entities
     are asked for, a tuple of objects. Arguments must be strings
     specifying an entity to return in the tuple. The first argument
     must be the name of the ftype to inquire and the remainder of the
     arguments must be keyvals fo entities to return. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg<2)
    {
      std::string err =  "fuda_ftype_get: at least 2 argument required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_ftype_get: First argument must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Find the ftype - if not found ft will be NULL.
  Ftype *ft = ftype_find(name);
  if (ft==NULL)
    {
      // Ftype not found.
      std::string err =  "fuda_ftype_get: ftype not found: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Create the output tuple.
  int nkey = narg-1;
  PyObject* ot_obj = PyTuple_New (nkey);

  // Loop over arguments which are keywords.
  for (int iarg=1; iarg<narg; iarg++)
    {
      int ikey = iarg-1;

      // All keys must be of type string.
      PyObject *key_obj = PyTuple_GetItem(args, iarg);
      if (!PyUnicode_Check(key_obj))
	{
	  std::string err =  "fuda_ftype_get: ";
	  snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
          err += ". argument must be a keyword string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
      
      // Get the key.
      const char *key = PyUnicode_AsUTF8(key_obj);

      // Lookup key.
      try {
	if (strcmp(key,"nparam")==0)
	  {
	    // Create and set an int opject.
	    int i = ft->get_nparam();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_ftype_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"nvar")==0)
	  {
	    // Create and set an int opject.
	    int i = ft->get_nvar();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_ftype_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"descr")==0)
	  {
	    // Create and set a string opject.
	    std::string str;
	    ft->get_descr(str);
	    PyObject* key_obj = PyUnicode_FromString (str.c_str());
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_ftype_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else
	  {
	    // Invalid keyword.
	    
	    // Deallocate the output tuple.
	    Py_DECREF(ot_obj);
	    
	    // Report error and return.
	    std::string err =  "fuda_ftype_get: ";
	    snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	    err += ". argument is an invalid keyword: ";
	    err += key;
	    PyErr_SetString(PyExc_Exception, err.c_str());
	    return NULL;
	  }
      }
      catch (...) {
	std::string err =  "fuda_ftype_get: unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }      
    }

  // If nkey==1, we don't wrap in a tuple.
  if (nkey==1)
    {
      // Extract the only element there is.
      PyObject *out_obj = PyTuple_GetItem(ot_obj,0);
      Py_INCREF(out_obj);

      // Deallocate tuple object which we don't return.
      Py_DECREF(ot_obj);

      // Return single object.
      return out_obj;
    }
  else
    {
      // Return output tuple.
      return ot_obj;
    }
}




extern "C" PyObject *
fuda_func(PyObject *self, PyObject *args)
{
  /* fuda_func declares a new function. The required arguemnts are:
     func_name, ftype_name, dtype_name and param_seq. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=4)
    {
      std::string err =  "fuda_func: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 4)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_func: 1. arg (name) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Second argument must be the ftype_name.
  PyObject *ftype_name_obj = PyTuple_GetItem(args, 1);
  if (!PyUnicode_Check(ftype_name_obj))
    {
      std::string err =  "fuda_func: 2. arg (ftype) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the ftype name.
  const char *ftype_name = PyUnicode_AsUTF8(ftype_name_obj);

  // 3. argument must be the dtype_name.
  PyObject *dtype_name_obj = PyTuple_GetItem(args, 2);
  if (!PyUnicode_Check(dtype_name_obj))
    {
      std::string err =  "fuda_func: 3. arg (dtype) must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  // Get the dtype name.
  const char *dtype_name = PyUnicode_AsUTF8(dtype_name_obj);

  /* 4. argument must be a sequence with the names of the
     parameters. */
  PyObject *params_obj = PyTuple_GetItem(args,3);
  if (!PySequence_Check(params_obj))
    {
      std::string err =  "fuda_func: 4. arg (parameters) must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be at least one item in the par sequence.
  if (PySequence_Length(params_obj)<1)
    {
      std::string err = "fuda_func: 4. arg must be a sequence ";
      err += "with at one or more parameter names";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the params sequence and get parameter names.
  std::list<std::string> p_names;
  for(int i=0; i<PySequence_Length(params_obj); i++)
    {
      // Get item from tuple.
      PyObject *param_name_obj = PySequence_GetItem(params_obj,i);

      // Check it is a string.
      if (!PyUnicode_Check(param_name_obj))
	{
	  std::string err =  "fuda_func: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 4. arg must be of type string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the parameter name.
      const char *param_name = PyUnicode_AsUTF8(param_name_obj);
      p_names.push_back(param_name);
    }

  // Create the function.
  Func *fn;
  try {
    fn = fuda->func_add(name, ftype_name, dtype_name, p_names);
  }
  catch (Uferr::FuncNameInvalid& e) {
    std::string err = 
      "fuda_func: func name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FuncNameAlreadyUsed& e) {
    std::string err = 
      "fuda_func: func name already used: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::FtypeNameInvalid& e) {
    std::string err = 
      "fuda_func: ftype name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  } 
  catch (Uferr::DtypeNameInvalid& e) {
    std::string err = 
      "fuda_func: dtype name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::NparamInvalid) {
    std::string err = 
      "fuda_func: Invalid number of parameters";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamNameInvalid& e) {
    std::string err = 
      "fuda_func: parameter name invalid: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamNotConst& e) {
    std::string err = 
      "fuda_func: parameter number ";
    snprintf(fstr, sizeof(fstr),"%d",e.i);
    err = err + fstr;
    err = err + " must be of the CONST kind";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::DtypeExplVarMissing& e) {
    std::string err = 
      "fuda_func: dtype explanatory variable missing: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_func: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  // Set current function.
  current_func = fn;
  current_func_name = name;
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_func_get(PyObject *self, PyObject *args)
{
  /* This routine returns either a single object or, if more entities
     are asked for, a tuple of objects. Arguments must be strings
     specifying an entity to return in the tuple. The first argument
     must be th function name and all subsequent arguemnts are
     interpreted as keyvals for which to return entities. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg<2)
    {
      std::string err =  "fuda_func_get: at least 2 arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_func_get: First argument must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Find the func - if not found fn will be NULL.
  Func *fn = func_find(name);
  if (fn==NULL)
    {
      // Func not found.
      std::string err =  "fuda_func_get: func not found: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Create the output tuple.
  int nkey = narg-1;
  PyObject* ot_obj = PyTuple_New (nkey);

  // Loop over arguments which are keywords.
  for (int iarg=1; iarg<narg; iarg++)
    {
      int ikey = iarg-1;

      // All keys must be of type string.
      PyObject *key_obj = PyTuple_GetItem(args, iarg);
      if (!PyUnicode_Check(key_obj))
	{
	  std::string err =  "fuda_func_get: ";
	  snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
          err += ". argument must be a keyword string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
      
      // Get the key.
      const char *key = PyUnicode_AsUTF8(key_obj);

      // Lookup key.
      try {
	if (strcmp(key,"ftype")==0)
	  {
	    // Create and set a string opject.
	    std::string str;
	    fn->get_ftype()->get_name(str);
	    PyObject* key_obj = PyUnicode_FromString (str.c_str());
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_func_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"dtype")==0)
	  {
	    // Create and set a string opject.
	    std::string str;
	    fn->get_dtype()->get_name(str);
	    PyObject* key_obj = PyUnicode_FromString (str.c_str());
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_func_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"nparam")==0)
	  {
	    // Create and set an int opject.
	    int i = fn->get_ftype()->get_nparam();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_func_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"nvar")==0)
	  {
	    // Create and set an int opject.
	    int i = fn->get_ftype()->get_nvar();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_func_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"use")==0)
	  {
	    // Create and set an int opject.
	    int i = fn->get_use();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_func_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"ndata")==0)
	  {
	    // Create and set an int opject.
	    int i = fn->get_ndata();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_func_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else
	  {
	    // Invalid keyword.
	    
	    // Deallocate the output tuple.
	    Py_DECREF(ot_obj);
	    
	    // Report error and return.
	    std::string err =  "fuda_func_get: ";
	    snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	    err += ". argument is an invalid keyword: ";
	    err += key;
	    PyErr_SetString(PyExc_Exception, err.c_str());
	    return NULL;
	  }
      }
      catch (...) {
	std::string err =  "fuda_func_get: unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }      
    }

  // If nkey==1, we don't wrap in a tuple.
  if (nkey==1)
    {
      // Extract the only element there is.
      PyObject *out_obj = PyTuple_GetItem(ot_obj,0);
      Py_INCREF(out_obj);

      // Deallocate tuple object which we don't return.
      Py_DECREF(ot_obj);

      // Return single object.
      return out_obj;
    }
  else
    {
      // Return output tuple.
      return ot_obj;
    }
}


extern "C" PyObject *
fuda_func_call(PyObject *self, PyObject *args)
{
  /* fuda_func_call calls the specified function with the specified
     explanatory variable values which must correspond to the
     explanatory variables of the functions dtype. It takes two
     arguments: func_name and var_value_seq, where var_value_seq
     is a sequence of double precision numbers. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=2)
    {
      std::string err =  "fuda_func_call: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 2)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the func name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_func_call: 2. arg. must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Lookup func.
  Func *fn = func_find(name);
  if (fn==NULL)
    {
      std::string err =  "fuda_func_call: Invalid function name: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get number of explanatory parameters in functions data type.
  int nexpl = fn->get_dtype()->get_dim();

  // Second argument must be a sequence with the values of the parameters.
  PyObject *varseq_obj = PyTuple_GetItem(args,1);
  if (!PySequence_Check(varseq_obj))
    {
      std::string err =  "fuda_func_call: 2. arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be nexpl items in the varseq.
  if (PySequence_Length(varseq_obj)!=nexpl)
    {
      std::string err = "fuda_func_call: 2. arg. sequence size invalid: ";
      snprintf(fstr, sizeof(fstr), "%zd", PySequence_Length(varseq_obj));
      err += fstr;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the varseq and setup vector with values.
  std::vector<double> value_vector(nexpl);
  for(int i=0; i<nexpl; i++)
    {
      // Get item from sequence.
      PyObject *varval_obj = PySequence_GetItem(varseq_obj,i);

      // Check it is a float.
      if (!PyFloat_Check(varval_obj))
	{
	  std::string err =  "fuda_func_call: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg. must be of type float";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      value_vector[i] = PyFloat_AsDouble(varval_obj);
    }

  // Evaluate function.
  double fval;
  try {
    fval = fn->call_by_expl(value_vector);
  }
  catch (Uferr::NumVarInvalid& e) {
    std::string err = 
      std::string("fuda_func_call: Invalid number of variables: ") + std::to_string(e.n);
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalCallFuncError& e) {
    std::string err = "fuda_func_call: eval_call catched exception ";
    err += "from compiled function associated with function : ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_func_call: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("d", fval);
}



extern "C" PyObject *
fuda_func_call_by_expl(PyObject *self, PyObject *args)
{
  /* fuda_func_call_by_expl calls the specified function with the specified
     explanatory variable values which must correspond to the
     explanatory variables of the function. It takes two
     arguments: func_name and var_value_seq, where var_value_seq
     is a sequence of double precision numbers. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=2)
    {
      std::string err =  "fuda_func_call_by_expl: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 2)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the func name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_func_call_by_expl: 2. arg. must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Lookup func.
  Func *fn = func_find(name);
  if (fn==NULL)
    {
      std::string err =  "fuda_func_call_by_expl: Invalid function name: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get number of explanatory parameters in function.
  int nexpl = fn->get_nexpl();

  // Second argument must be a sequence with the values of the parameters.
  PyObject *varseq_obj = PyTuple_GetItem(args,1);
  if (!PySequence_Check(varseq_obj))
    {
      std::string err =  "fuda_func_call_by_expl: 2. arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be nexpl items in the varseq.
  if (PySequence_Length(varseq_obj)!=nexpl)
    {
      std::string err = "fuda_func_call_by_expl: 2. arg. sequence size invalid: ";
      snprintf(fstr, sizeof(fstr), "%zd", PySequence_Length(varseq_obj));
      err += fstr;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the varseq and setup vector with values.
  std::vector<double> value_vector(nexpl);
  for(int i=0; i<nexpl; i++)
    {
      // Get item from tuple.
      PyObject *varval_obj = PySequence_GetItem(varseq_obj,i);

      // Check it is a float.
      if (!PyFloat_Check(varval_obj))
	{
	  std::string err =  "fuda_func_call_by_expl: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg. must be of type float";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      value_vector[i] = PyFloat_AsDouble(varval_obj);
    }

  // Evaluate function.
  double fval;
  try {
    fval = fn->call_by_expl(value_vector);
  }
  catch (Uferr::NumVarInvalid& e) {
    std::string err = 
      std::string("fuda_func_call: Invalid number of variables: ") + std::to_string(e.n);
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalCallFuncError& e) {
    std::string err = "fuda_func_call_by_expl: eval_call catched exception ";
    err += "from compiled function associated with function : ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_func_call_by_expl: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("d", fval);
}



extern "C" PyObject *
fuda_func_call_by_var(PyObject *self, PyObject *args)
{
  /* fuda_func_call_by_var call the specified function with the
     specified variable values. It takes two arguments: func_name and
     var_value_seq, where var_value_seq is a sequence of double
     precision numbers. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=2)
    {
      std::string err =  "fuda_func_call_by_var: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 2)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the func name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_func_call_by_var: 2. arg. must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Lookup func.
  Func *fn = func_find(name);
  if (fn==NULL)
    {
      std::string err =  "fuda_func_call_by_var: Invalid function name: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get number of variables in function.
  int nvar = fn->get_ftype()->get_nvar();

  // Second argument must be a sequence with the values of the variables.
  PyObject *varseq_obj = PyTuple_GetItem(args,1);
  if (!PySequence_Check(varseq_obj))
    {
      std::string err =  "fuda_func_call_by_var: 2. arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be nvar items in the varseq.
  if (PySequence_Length(varseq_obj)!=nvar)
    {
      std::string err = "fuda_func_call_by_var: 2. arg. tuple size invalid: ";
      snprintf(fstr, sizeof(fstr), "%zd", PyTuple_Size(varseq_obj));
      err += fstr;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the varseq and setup vector with values.
  std::vector<double> value_vector(nvar);
  for(int i=0; i<nvar; i++)
    {
      // Get item from sequence.
      PyObject *varval_obj = PySequence_GetItem(varseq_obj,i);

      // Check it is a float.
      if (!PyFloat_Check(varval_obj))
	{
	  std::string err =  "fuda_func_call_by_var: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg. must be of type float";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      value_vector[i] = PyFloat_AsDouble(varval_obj);
    }

  // Evaluate function.
  double fval;
  try {
    fval = fn->call_by_var(value_vector);
  }
  catch (Uferr::NumVarInvalid& e) {
    std::string err = 
      std::string("fuda_func_call_by_var: Invalid number of variables: ") + std::to_string(e.n);
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalCallFuncError& e) {
    std::string err = "fuda_func_call_by_var: eval_call catched exception from ";
    err += "compiled function associated with function : ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_func_call_by_var: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("d", fval);
}



extern "C" PyObject *
fuda_func_deriv_by_var(PyObject *self, PyObject *args)
{
  /* fuda_func_deriv_by_var evaluate the derrivatives of the specified
     function with with respect to the variables at the specified
     variable values. It takes two arguments: func_name and
     var_value_seq, where var_value_seq is a sequence of double
     precision numbers. It returns a tuple with the derivatves of the
     function with respect to the variables. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=2)
    {
      std::string err =  "fuda_func_deriv_by_var: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 2)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the func name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_func_deriv_by_var: 2. arg. must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Lookup func.
  Func *fn = func_find(name);
  if (fn==NULL)
    {
      std::string err =  "fuda_func_deriv_by_var: Invalid function name: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get number of variables in function.
  int nvar = fn->get_ftype()->get_nvar();

  // Second argument must be a sequence with the values of the variables.
  PyObject *varseq_obj = PyTuple_GetItem(args,1);
  if (!PySequence_Check(varseq_obj))
    {
      std::string err =  "fuda_func_deriv_by_var: 2. arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be nvar items in the varseq.
  if (PySequence_Length(varseq_obj)!=nvar)
    {
      std::string err = "fuda_func_deriv_by_var: 2. arg. sequence size invalid: ";
      snprintf(fstr, sizeof(fstr), "%zd", PySequence_Length(varseq_obj));
      err += fstr;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the varseq and setup vector with values.
  std::vector<double> value_vector(nvar);
  for(int i=0; i<nvar; i++)
    {
      // Get item from sequence.
      PyObject *varval_obj = PySequence_GetItem(varseq_obj,i);

      // Check it is a float.
      if (!PyFloat_Check(varval_obj))
	{
	  std::string err =  "fuda_deriv_call_by_var: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg. must be of type float";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      value_vector[i] = PyFloat_AsDouble(varval_obj);
    }

  // Evaluate function.
  double fval;
  std::vector<double> deriv_vector(nvar);
  try {
    fval = fn->call_by_var(value_vector, deriv_vector);
  }
  catch (Uferr::NumVarInvalid& e) {
    std::string err = 
      std::string("fuda_func_deriv_by_var: Invalid number of variables: ") + std::to_string(e.n);
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalCallFuncError& e) {
    std::string err = "fuda_func_deriv_by_var: eval_call catched exception from ";
    err += "compiled function associated with function : ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_func_deriv_by_var: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  // We return the derivative vector as a tuple.

  // Create the output tuple.
  PyObject* ot_obj = PyTuple_New (nvar);

  // Loop over derivatives and fill tuple.
  for (int ivar=0; ivar<nvar; ivar++)
    {
      // Create and set a float opject.
      PyObject* deriv_obj = PyFloat_FromDouble(deriv_vector[ivar]);
 
      // Insert in ot_obj.
      if (PyTuple_SetItem (ot_obj, ivar, deriv_obj)!=0)
	return (NULL);
    }

  // Return output tuple.
  return ot_obj;
}


extern "C" PyObject *
fuda_func_call_by_param(PyObject *self, PyObject *args)
{
  /* fuda_func_call_by_param call the specified function with the
     specified parameter values. It takes two arguments: func_name and
     param_value_seq, where param_value_seq is a sequence of double
     precision numbers. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=2)
    {
      std::string err =  "fuda_func_call_by_param: invalid number of args: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 2)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the func name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_func_call_by_param: 2. arg. must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Lookup func.
  Func *fn = func_find(name);
  if (fn==NULL)
    {
      std::string err =  "fuda_func_call_by_param: Invalid function name: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get number of parameters in function.
  int nparam = fn->get_ftype()->get_nparam();

  // Second argument must be a sequence with the values of the parameters.
  PyObject *varseq_obj = PyTuple_GetItem(args,1);
  if (!PySequence_Check(varseq_obj))
    {
      std::string err =  "fuda_func_call_by_param: 2. arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // There must be nparam items in the varseq.
  if (PySequence_Length(varseq_obj)!=nparam)
    {
      std::string err = "fuda_func_call_by_param: 2. arg. sequence size invalid: ";
      snprintf(fstr, sizeof(fstr), "%zd", PySequence_Length(varseq_obj));
      err += fstr;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the varseq and setup vector with values.
  std::vector<double> value_vector(nparam);
  for(int i=0; i<nparam; i++)
    {
      // Get item from tuple.
      PyObject *varval_obj = PySequence_GetItem(varseq_obj,i);

      // Check it is a float.
      if (!PyFloat_Check(varval_obj))
	{
	  std::string err =  "fuda_func_call_by_param: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg. must be of type float";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      value_vector[i] = PyFloat_AsDouble(varval_obj);
    }

  // Evaluate function.
  double fval;
  try {
    fval = fn->call_by_param(value_vector);
  }
  catch (Uferr::NumParamInvalid& e) {
    std::string err =
      std::string("fuda_func_call_by_param: Invalid number of parameters: ") + std::to_string(e.n);
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_func_call_by_param: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("d", fval);
}



extern "C" PyObject *
fuda_data(PyObject *self, PyObject *args)
{
  /* fuda_data declares a data point. The data point will be associated
     with the current function. The dtype is derived from the current
     function. It takes nexpl+2 arguments exactly corresponding to the
     explanatory variable values, the point value and the
     corresponding uncertainty. */

  // We need a current function to be set.
  if (current_func==NULL)
    {
      std::string err =  "fuda_data: no current function";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get dimension of dtype (number of explanatory variables).
  int nexpl = current_func->get_dtype()->get_dim();

  // We expect nexpl+2 args (explanatory vars, value and uncertainty).
  int narg = PyTuple_Size(args);
  if (narg!=nexpl+2)
    {
      std::string err =  "fuda_data: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d (expected %d)",narg, nexpl+2);
      err += fstr;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over explanatory variable values to check types.
  for(int i=0; i<nexpl+2; i++)
    {
      // Get the object.
      PyObject *expl_obj = PyTuple_GetItem(args, i);
      if (!PyFloat_Check(expl_obj))
	{
	  std::string err =  "fuda_data: arguments must be of type double";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
    }

  // Add a data record to current function.
  Data *data = current_func->add_data();

  // Loop over explanatory variable values and set them.
  for(int i=0; i<nexpl; i++)
    {
      // Get the object.
      PyObject *expl_obj = PyTuple_GetItem(args, i);

      // Get the explanatory variable value.
      double expl_val = PyFloat_AsDouble(expl_obj);

      // Set the value.
      data->set_ix(i,expl_val);
    }

  // Set data init and value.
  PyObject *val_obj = PyTuple_GetItem(args, nexpl);
  double val = PyFloat_AsDouble(val_obj);
  data->set_init(val);

  // Set data uncertainty.
  PyObject *u_obj = PyTuple_GetItem(args, nexpl+1);
  double u = PyFloat_AsDouble(u_obj);
  data->set_u(u);
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_data_del(PyObject *self, PyObject *args)
{
  /* fuda_data_del will delete all data for the current function. */

  // We expect zero args.
  if (PyTuple_Size(args)!=0)
    {
      std::string err =  "fuda_data_del: requires no arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Do we have a current function set.
  if (current_func==NULL)
    {
      // No current function.
      std::string err =  "fuda_data_del: no current function";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  else
    {
      try {  
        // Delete the data.
        current_func->delete_data();
      }
      catch (...) {
        std::string err =  "fuda_data_del: unspecified exception";
        PyErr_SetString(PyExc_Exception, err.c_str());
        return NULL;
      }
    }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_data_del_all(PyObject *self, PyObject *args)
{
  /* fuda_data_del_all will delete all data for all functions. */

  // We expect zero args.
  if (PyTuple_Size(args)!=0)
    {
      std::string err =  "fuda_data_del_all: requires no arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  try {  
    // Loop over func records.
    for (Func_iterator fi=fuda->func_begin();
         fi!=fuda->func_end(); fi++)
      {
        // Delete the data.
        (*fi)->delete_data();
      }
  }
  catch (...) {
    std::string err =  "fuda_data_del_all: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_func_get_current(PyObject *self, PyObject *args)
{
  /* fuda_func_get_current returns the name of the current function. */

  // We expect zero args.
  if (PyTuple_Size(args)!=0)
    {
      std::string err =  "fuda_func_get_current: requires no arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Return the current function name.
  if (current_func==NULL)
    {
      // No current function.
      Py_INCREF(Py_None);
      return Py_None;
    }
  else
    {
      return Py_BuildValue("s", current_func_name.c_str());
    }
  
}


extern "C" PyObject *
fuda_func_current(PyObject *self, PyObject *args)
{
  /* fuda_func_current explicitly sets the current function. */

  char *name;
  
  // Get name.
  if (!PyArg_ParseTuple(args, "s", &name))
    return NULL;

  // Look for function name.
  Func *fn = func_find(name);
  
  // Did we find one?
  if (fn==NULL)
    {
      std::string err =  "fuda_func_current: function not found: ";
      err += name;
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  else
    {
      // Set current function.
      current_func = fn;
      current_func_name = name;
    }
  
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_eval_init(PyObject *self, PyObject *args)
{
  /* fuda_eval_init call the fuda eval_init method which initializes
     the fuda eval structure. An explicit call to this function is
     rarely needed as the fuda eval_init method is automatically
     called when running a minimization. */

  // We expect zero args.
  if (PyTuple_Size(args)!=0)
    {
      std::string err =  "fuda_eval_init: requires no arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Setup eval structure.
  try {
    fuda->eval_init();
  }
  catch (Uferr::ParamKindInvalid) {
    std::string err = 
      "fuda_eval_init: A parameter with an invalid kind specifier\n";
    err += "was addressed. This is an internal inconsistency uf error.\n";
    err += "Please report to the author";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::NoData) {
    std::string err = "fuda_eval_init: there are no data to fit";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_eval_init: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_eval_data_recalc(PyObject *self, PyObject *args)
{
  /* fuda_eval_data_recalc call the fuda eval_data_recalc method which
     initializes the value of all data points to the value calculated
     for the data point for the current parameters. */

  // We expect zero args.
  if (PyTuple_Size(args)!=0)
    {
      std::string err =  "fuda_eval_data_recalc: requires no arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Setup eval structure.
  try {
    fuda->eval_data_recalc();
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_data_recalc: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_eval_data_recalc: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_eval_data_random(PyObject *self, PyObject *args)
{
  /* fuda_eval_data_random call the fuda eval_data_random
     method which initializes the value of all data points to their
     initial value plus gaussian noise with a standard deviation taken
     from the uncertainty of the data points. */

  // We expect zero args.
  if (PyTuple_Size(args)!=0)
    {
      std::string err =  "fuda_eval_data_random: requires no arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Setup eval structure.
  try {
    fuda->eval_data_random();
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_data_random: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_eval_data_random: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_eval_get(PyObject *self, PyObject *args)
{
  /* This routine returns either a single object or, if more entities
     are asked for, a tuple of objects. Arguments must be strings
     specifying an entity to return in the tuple. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg<1)
    {
      std::string err =  "fuda_eval_get: at least 1 argument required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Create the output tuple.
  int nkey = narg;
  PyObject* ot_obj = PyTuple_New (nkey);

  // Loop over arguments which are keywords.
  for (int iarg=0; iarg<narg; iarg++)
    {
      int ikey = iarg;

      // All keys must be of type string.
      PyObject *key_obj = PyTuple_GetItem(args, iarg);
      if (!PyUnicode_Check(key_obj))
	{
	  std::string err =  "fuda_eval_get: ";
	  snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
          err += ". argument must be a keyword string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
      
      // Get the key.
      const char *key = PyUnicode_AsUTF8(key_obj);

      // Lookup key.
      try {
	if (strcmp(key,"sync")==0)
	  {
	    // Create and set an int opject.
	    int i = fuda->eval_is_sync();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_eval_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"nfree")==0)
	  {
	    // Create and set an int opject.
	    int i = fuda->eval_get_nfree();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_eval_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"nexpl")==0)
	  {
	    // Create and set an int opject.
	    int i = fuda->eval_get_nexpl();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_eval_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"nconst")==0)
	  {
	    // Create and set an int opject.
	    int i = fuda->eval_get_nconst();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_eval_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"ndata")==0)
	  {
	    // Create and set an int opject.
	    int i = fuda->eval_get_ndata();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_eval_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else
	  {
	    // Invalid keyword.
	    
	    // Deallocate the output tuple.
	    Py_DECREF(ot_obj);
	    
	    // Report error and return.
	    std::string err =  "fuda_eval_get: ";
	    snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	    err += ". argument is an invalid keyword: ";
	    err += key;
	    PyErr_SetString(PyExc_Exception, err.c_str());
	    return NULL;
	  }
      }
      catch (Uferr::EvalNotSync) {
	std::string err =  "fuda_eval_get: eval structure not in sync";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (...) {
	std::string err =  "fuda_eval_get: unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }      
    }

  // If nkey==1, we don't wrap in a tuple.
  if (nkey==1)
    {
      // Extract the only element there is.
      PyObject *out_obj = PyTuple_GetItem(ot_obj,0);
      Py_INCREF(out_obj);

      // Deallocate tuple object which we don't return.
      Py_DECREF(ot_obj);

      // Return single object.
      return out_obj;
    }
  else
    {
      // Return output tuple.
      return ot_obj;
    }
}



extern "C" PyObject *
fuda_eval_get_free(PyObject *self, PyObject *args)
{
  /* Return the name of the i'th free variable in a fitting run. */

  int index=0;
  
  if (!PyArg_ParseTuple(args, "i", &index))
    return NULL;

  // Get value.
  std::string name;
  try {
    Param *pm = fuda->eval_get_free(index);
    pm->get_name(name);
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_get_free: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalFreeIndexInvalid& e) {
      std::string err =  "fuda_eval_get_free: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", name.c_str());
}


extern "C" PyObject *
fuda_eval_get_expl(PyObject *self, PyObject *args)
{
  /* Return the name of the i'th explanatory variable in a fitting
     run. */

  int index=0;
  
  if (!PyArg_ParseTuple(args, "i", &index))
    return NULL;

  // Get value.
  std::string name;
  try {
    Param *pm = fuda->eval_get_expl(index);
    pm->get_name(name);
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_get_expl: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalExplIndexInvalid& e) {
      std::string err =  "fuda_eval_get_expl: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", name.c_str());
}


extern "C" PyObject *
fuda_eval_get_const(PyObject *self, PyObject *args)
{
  /* Return the name of the i'th constant parameter in a fitting run. */

  int index=0;
  
  if (!PyArg_ParseTuple(args, "i", &index))
    return NULL;

  // Get value.
  std::string name;
  try {
    Param *pm = fuda->eval_get_const(index);
    pm->get_name(name);
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_get_const: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalConstIndexInvalid& e) {
      std::string err =  "fuda_eval_get_const: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", name.c_str());
}


extern "C" PyObject *
fuda_eval_get_data(PyObject *self, PyObject *args)
{
  /* Return the i'th data point as a tuple. */

  int index=0;
  
  if (!PyArg_ParseTuple(args, "i", &index))
    return NULL;

  // Declare the output tuple pointer.
  PyObject* tup_obj = NULL;

  // Get the data point.
  try {
    // Get reference to Data object.
    Data *d = fuda->eval_get_data(index);

    // Get number of explanatory variables.
    int xdim = d->get_dim();
  
    // Create the output tuple.
    tup_obj = PyTuple_New (xdim+2);

    // Loop over explanatory variables.
    for (int ix=0; ix<xdim; ix++)
      {
	// Create and set a float opject.
	double xval = d->get_ix(ix);
	PyObject* x_obj = PyFloat_FromDouble (xval);

	// Insert explanatory variable in tup_obj.
	if (PyTuple_SetItem (tup_obj, ix, x_obj)!=0)
	  {
	    std::string err =  "fuda_eval_get_data: PyTuple_SetItem failed on ";
	    snprintf(fstr, sizeof(fstr),"%d",ix+1); err += fstr;
	    err += ". argument";
	    PyErr_SetString(PyExc_Exception, err.c_str());

	    // Deallocate tuple object which we don't return.
	    Py_DECREF(tup_obj);

	    return NULL;
	  }
      }

    // Create and set value.
    double val = d->get_value();
    PyObject* val_obj = PyFloat_FromDouble (val);

    // Insert value in tup_obj.
    if (PyTuple_SetItem (tup_obj, xdim, val_obj)!=0)
      {
	std::string err =  "fuda_eval_get_data: PyTuple_SetItem failed on ";
	snprintf(fstr, sizeof(fstr),"%d",xdim+1); err += fstr;
	err += ". argument";
	PyErr_SetString(PyExc_Exception, err.c_str());

	// Deallocate tuple object which we don't return.
	Py_DECREF(tup_obj);

	return NULL;
      }

    // Create and set uncertainty.
    double u = d->get_u();
    PyObject* u_obj = PyFloat_FromDouble (u);

    // Insert uncertainty in tup_obj.
    if (PyTuple_SetItem (tup_obj, xdim+1, u_obj)!=0)
      {
	std::string err =  "fuda_eval_get_data: PyTuple_SetItem failed on ";
	snprintf(fstr, sizeof(fstr),"%d",xdim+2); err += fstr;
	err += ". argument";
	PyErr_SetString(PyExc_Exception, err.c_str());

	// Deallocate tuple object which we don't return.
	Py_DECREF(tup_obj);

	return NULL;
      }
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_get_const: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());

    // Deallocate tuple object which we don't return.
    Py_DECREF(tup_obj);

    return NULL;
  }
  catch (Uferr::EvalDataIndexInvalid& e) {
      std::string err =  "fuda_eval_get_data: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());

      // Deallocate tuple object which we don't return.
      Py_DECREF(tup_obj);

      return NULL;
  }

  // Return data tuple.
  return tup_obj;
}


extern "C" PyObject *
fuda_eval_get_func(PyObject *self, PyObject *args)
{
  /* Return the name of the function the i'th data point belongs to */

  int index=0;
  
  if (!PyArg_ParseTuple(args, "i", &index))
    return NULL;

  // Get name.
  std::string name;
  try {
    Data *d = fuda->eval_get_data(index);
    d->get_func()->get_name(name);
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_get_const: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalDataIndexInvalid& e) {
      std::string err =  "fuda_eval_get_const: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", name.c_str());
}


extern "C" PyObject *
fuda_eval_get_dtype(PyObject *self, PyObject *args)
{
  /* Return the name of the dtype of the i'th data point. */

  int index=0;
  
  if (!PyArg_ParseTuple(args, "i", &index))
    return NULL;

  // Get name.
  std::string name;
  try {
    Data *d = fuda->eval_get_data(index);
    d->get_dtype()->get_name(name);
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_get_const: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalDataIndexInvalid& e) {
      std::string err =  "fuda_eval_get_const: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("s", name.c_str());
}


extern "C" PyObject *
fuda_eval_enorm(PyObject *self, PyObject *args)
{
  /* Get enorm for eval structure. */

  // Check # args.
  int narg = PyTuple_Size(args);
  if (narg!=0)
    {
      std::string err =  "fuda_eval_enorm: no arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get enorm.
  double value;
  try {
    value = fuda->eval_enorm();
  }
  catch (...) {
      std::string err =  "fuda_eval_enorm: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("d", value);
}


extern "C" PyObject *
fuda_eval_call(PyObject *self, PyObject *args)
{
  /* fuda_eval_call calculates the value of a data point when a data
     type, a vector with the explanatory variable values and a vector
     with the free parameters for the eval structure is given as
     arguments.  declare a new dtype. It takes two arguments:
     dtype_name and expl_seq, where expl_seq is a tuple with the
     explanatory variables for this dtype. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=3)
    {
      std::string err =  "fuda_eval_call: invalid number of arguments: ";
      snprintf(fstr, sizeof(fstr),"%d",narg);
      err += fstr;
      err += " (expected 3)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the data type name.
  PyObject *name_obj = PyTuple_GetItem(args, 0);
  if (!PyUnicode_Check(name_obj))
    {
      std::string err =  "fuda_dtype: Second argument must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the name.
  const char *name = PyUnicode_AsUTF8(name_obj);

  // Get pointer to dtype if exists or abort.
  Dtype* dtype = fuda->dtype_find(name);
  if (dtype==NULL)
    {
      std::string err =  "fuda_eval_call: data type does not exist";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Get number of explanatory variables for specified data type. This
  // is the number of elements expected in the xvec sequence argument.
  unsigned int xdim = dtype->get_dim();
  
  // Get number of free parameters in the eval structure. This is the
  // number of elements expected in the pvec sequence argument.
  unsigned int pdim = fuda->eval_get_nfree();
  
  /* Second argument must be a sequence with xdim float values */
  PyObject *xvec_obj = PyTuple_GetItem(args,1);
  if (!PySequence_Check(xvec_obj))
    {
      std::string err =  "fuda_eval_call: Second argument must be a seqence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // The length of the xvec sequence must be exactly xdim.
  if ((int) xdim!=PySequence_Length(xvec_obj))
    {
      std::string err =  "fuda_eval_call: Second argument sequence has invalid number of elements";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the xvec sequence and fill xvec.
  std::vector<double> xvec(xdim, 0.0);
  for(unsigned int i=0; i<xdim; i++)
    {
      // Get item from sequence.
      PyObject *xval_obj = PySequence_GetItem(xvec_obj, (int)i);

      // Check it is a double.
      if (!PyFloat_Check(xval_obj))
	{
	  std::string err =  "fuda_eval_call: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of second argument must be of type float";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      xvec[i] = PyFloat_AsDouble(xval_obj);
    }

  /* 3rd argument must be a sequence with pdim float values */
  PyObject *pvec_obj = PyTuple_GetItem(args,2);
  if (!PySequence_Check(pvec_obj))
    {
      std::string err =  "fuda_eval_call: 3rd argument must be a seqence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // The length of the pvec sequence must be exactly pdim.
  if ((int)pdim!=PySequence_Length(pvec_obj))
    {
      std::string err =  "fuda_eval_call: 3rd argument sequence has invalid number of elements";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Loop over the items of the pvec sequence and fill pvec.
  double pvec[pdim];
  for(unsigned int i=0; i<pdim; i++)
    {
      // Get item from sequence.
      PyObject *pval_obj = PySequence_GetItem(pvec_obj,(int)i);

      // Check that it is a double.
      if (!PyFloat_Check(pval_obj))
	{
	  std::string err =  "fuda_eval_call: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 3rd argument must be of type float";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      pvec[i] = PyFloat_AsDouble(pval_obj);
    }

  // Make the function call.
  try {
    fuda->eval_call(*dtype, xvec, pvec, (double*) NULL, (bool) 0);
  }
  catch (Uferr::EvalNotSync) {
    std::string err =  "fuda_eval_call: eval structure not in sync";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_eval_call: unexpected exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}





extern "C" PyObject *
fuda_lm_minimize(PyObject *self, PyObject *args)
{
  /* Execute a Levenberg-Marquart (minpack) minimization. */

  // We expect zero args.
  if (PyTuple_Size(args)!=0)
    {
      std::string err =  "fuda_lm_minimize: requires no arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Minimize.
  try {
    lm->minimize();
  }
  catch (Uferr::NoFreeParam) {
    std::string err = 
      "fuda_lm_minimize: No free parameters";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::NparamGtNdata) {
    std::string err = 
      "fuda_lm_minimize: number of data < number of parameters";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::LmTolInvalid& e) {
    std::string err = 
      "fuda_lm_minimize: Invalid tolerance: " + e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::LmTolsInvalid) {
    std::string err = 
      "fuda_lm_minimize: Invalid tolerances";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::LmMaxfevInvalid) {
    std::string err = 
      "fuda_lm_minimize: maxfev invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::LmFactorInvalid) {
    std::string err = 
      "fuda_lm_minimize: factor invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::ParamKindInvalid) {
    std::string err = 
      "fuda_lm_minimize: A parameter with an invalid kind specifier\n";
    err += "was addressed. This is an internal inconsistency uf error.\n";
    err += "Please report to the author";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::NoData) {
    std::string err = "fuda_lm_minimize: there are no data to fit";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (Uferr::EvalCallFuncError& e) {
    std::string err = "fuda_lm_minimize: eval_call catched exception from ";
    err += "compiled function associated with function : ";
    err += e.name;
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
    std::string err = 
      "fuda_lm_minimize: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_lm(PyObject *self, PyObject *args, PyObject *kw)
{
  /* Set control variables for Levenberg-Marquart minimization. */

  // Argument keyword list.
  static char kw_tol[] = "tol", kw_ftol[] = "ftol", kw_xtol[] = "xtol";
  static char kw_gtol[] = "gtol", kw_ctol[] = "ctol", kw_factor[] = "factor";
  static char kw_numderiv_eps[] = "numderiv_eps", kw_numderiv[] = "numderiv";
  static char kw_maxfev[] = "maxfev", kw_nprint[] = "nprint", kw_scale_covar[] = "scale_covar";
  static char *kwlist[] = {kw_tol, kw_ftol, kw_xtol, kw_gtol, kw_ctol, kw_factor,
                           kw_numderiv_eps, kw_numderiv, kw_maxfev, kw_nprint,
                           kw_scale_covar, NULL};
  double tol, ftol, xtol, gtol, ctol, factor, numderiv_eps;
  int numderiv, maxfev, nprint, scale_covar;

  // Check for optional keyword arguments.
  int iarg = 0;
  bool tol_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool ftol_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool xtol_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool gtol_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool ctol_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool factor_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool numderiv_eps_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool numderiv_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool maxfev_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool nprint_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool scale_covar_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;

  // Pass argument list.
  if (!PyArg_ParseTupleAndKeywords(args, kw, "|dddddddiiii", kwlist, 
				   &tol, &ftol, &xtol, &gtol, &ctol,
				   &factor, &numderiv_eps, &numderiv,
				   &maxfev, &nprint, &scale_covar))
    return NULL; 

  // Set arguments.
  if (tol_flg)
    lm->set_tol(tol);
  if (ftol_flg)
    lm->set_ftol(ftol);
  if (xtol_flg)
    lm->set_xtol(xtol);
  if (gtol_flg)
    lm->set_gtol(gtol);
  if (ctol_flg)
    lm->set_ctol(ctol);
  if (factor_flg)
    lm->set_factor(factor);
  if (numderiv_eps_flg)
    lm->set_numderiv_eps(numderiv_eps);
  if (numderiv_flg)
    lm->set_numderiv(numderiv);
  if (maxfev_flg)
    lm->set_maxfev(maxfev);
  if (nprint_flg)
    lm->set_nprint(nprint);
  if (scale_covar_flg)
    lm->set_scale_covar(scale_covar);

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_lm_get(PyObject *self, PyObject *args)
{
  /* This routine returns either a single object or, if more entities
     are asked for, a tuple of objects. Arguments must be strings
     specifying an entity to return in the tuple. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg<1)
    {
      std::string err =  "fuda_lm_get: at least 1 argument required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Create the output tuple.
  int nkey = narg;
  PyObject* ot_obj = PyTuple_New (nkey);

  // Loop over arguments which are keywords.
  for (int iarg=0; iarg<narg; iarg++)
    {
      int ikey = iarg;

      // All keys must be of type string.
      PyObject *key_obj = PyTuple_GetItem(args, iarg);
      if (!PyUnicode_Check(key_obj))
	{
	  std::string err =  "fuda_lm_get: ";
	  snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
          err += ". argument must be a keyword string";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
      
      // Get the key.
      const char *key = PyUnicode_AsUTF8(key_obj);

      // Lookup key.
      try {
	if (strcmp(key,"tol")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_tol();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"ftol")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_ftol();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"xtol")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_xtol();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"gtol")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_gtol();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"ctol")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_ctol();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"factor")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_factor();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"numderiv_eps")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_numderiv_eps();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"enorm")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_enorm();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"sd")==0)
	  {
	    // Create and set a float opject.
	    double d = lm->get_sd();
	    PyObject* key_obj = PyFloat_FromDouble(d);
	    
	    // Insert in ot_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"sync")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->is_sync();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"minimized")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->is_minimized();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"maxfev")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->get_maxfev();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"numderiv")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->get_numderiv();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"nprint")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->get_nprint();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"scale_covar")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->get_scale_covar();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"info")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->get_info();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"nfev")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->get_nfev();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"njev")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->get_njev();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"fit_ok")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->fit_is_ok();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else if (strcmp(key,"fit_converged")==0)
	  {
	    // Create and set an int opject.
	    int i = lm->is_converged();
	    PyObject* key_obj = PyLong_FromLong (i);
	    
	    // Insert in tup_obj.
	    if (PyTuple_SetItem (ot_obj, ikey, key_obj)!=0)
	      {
		std::string err =  "fuda_lm_get: PyTuple_SetItem failed on ";
		snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
		err += ". argument";
		PyErr_SetString(PyExc_Exception, err.c_str());
		return NULL;
	      }
	  }
	else
	  {
	    // Invalid keyword.
	    
	    // Deallocate the output tuple.
	    Py_DECREF(ot_obj);
	    
	    // Report error and return.
	    std::string err =  "fuda_lm_get: ";
	    snprintf(fstr, sizeof(fstr),"%d",iarg+1); err += fstr;
	    err += ". argument is an invalid keyword: ";
	    err += key;
	    PyErr_SetString(PyExc_Exception, err.c_str());
	    return NULL;
	  }
      }
      catch (Uferr::LmNoFit) {
	std::string err =  "fuda_lm_get: no fit";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (...) {
	std::string err =  "fuda_lm_get: unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }      
    }

  // If nkey==1, we don't wrap in a tuple.
  if (nkey==1)
    {
      // Extract the only element there is.
      PyObject *out_obj = PyTuple_GetItem(ot_obj,0);
      Py_INCREF(out_obj);

      // Deallocate tuple object which we don't return.
      Py_DECREF(ot_obj);

      // Return single object.
      return out_obj;
    }
  else
    {
      // Return output tuple.
      return ot_obj;
    }
}




extern "C" PyObject *
fuda_lm_update_param(PyObject *self, PyObject *args)
{
  /* Update value and esd of all free parameters to values estimated
     by last lm_minimize call. */

  // Check number of arguments.
  if (PyTuple_Size(args)>0)
    {
      std::string err =  "fuda_lm_update_param: no arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Update value and esd of all free parameters to values estimated
  // by last lm_minimize call.
  try {
    // Loop over free paramters.
    for (unsigned int i=0; i<fuda->eval_get_nfree(); i++)
      {
	// Get the parameter.
	Param *param = fuda->eval_get_free(i);
	
	// Update value and esd.
	param->set_value(lm->get_value(i));
	param->set_esd(lm->get_esd(i));
      }    
  }
  catch (Uferr::LmParamIndexInvalid& e) {
      std::string err =  "fuda_lm_get_value: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (Uferr::LmNoFit) {
      std::string err =  "fuda_lm_get_value: no fit";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "fuda_lm_get_value: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_lm_get_value(PyObject *self, PyObject *args)
{
  /* Return value of the i'th variable after minimization. */

  int index=0;
  
  if (!PyArg_ParseTuple(args, "i", &index))
    return NULL;

  // Get value.
  double value;
  try {
    value = lm->get_value(index);
  }
  catch (Uferr::LmParamIndexInvalid& e) {
      std::string err =  "fuda_lm_get_value: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (Uferr::LmNoFit) {
      std::string err =  "fuda_lm_get_value: no fit";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "fuda_lm_get_value: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("d", value);
}


extern "C" PyObject *
fuda_lm_get_esd(PyObject *self, PyObject *args)
{
  /* Get estimated uncertainty of the i'th free variable after minimization. */

  int index=0;
  
  if (!PyArg_ParseTuple(args, "i", &index))
    return NULL;

  // Get value.
  double value;
  try {
    value = lm->get_esd(index);
  }
  catch (Uferr::LmParamIndexInvalid& e) {
      std::string err =  "fuda_lm_get_esd: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (Uferr::LmNoFit) {
      std::string err =  "fuda_lm_get_esd: no fit";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (Uferr::LmFitNotConverged) {
      std::string err =  "fuda_lm_get_esd: fit not converged";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "fuda_lm_get_esd: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("d", value);
}


extern "C" PyObject *
fuda_lm_get_covar(PyObject *self, PyObject *args)
{
  /* Return a value from the covariance matrix. */

  int i=0, j=0;
  
  if (!PyArg_ParseTuple(args, "ii", &i, &j))
    return NULL;

  // Get value.
  double value;
  try {
    value = lm->get_covar(i,j);
  }
  catch (Uferr::LmParamIndexInvalid& e) {
      std::string err =  "fuda_lm_get_covar: index out of bounds";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (Uferr::LmNoFit) {
      std::string err =  "fuda_lm_get_covar: no fit";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "fuda_lm_get_covar: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Return value.
  return Py_BuildValue("d", value);
}


extern "C" PyObject *
fuda_rand_set_seed(PyObject *self, PyObject *args)
{
  /* Set rand seed. */

  int new_seed = seed;
  
  if (!PyArg_ParseTuple(args, "i", &new_seed))
    return NULL;

  seed = new_seed;
  FUDA::srand(seed);

  // Return value.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
fuda_rand_get_seed(PyObject *self, PyObject *args)
{
  /* Return the seed. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg<0)
    {
      std::string err =  "fuda_rand_get_seed: no arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Return value.
  return Py_BuildValue("i", seed);
}


extern "C" PyObject *
fuda_rand_uniform(PyObject *self, PyObject *args)
{
  /* Return uniformly [0;1[ distributed random number. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg<0)
    {
      std::string err =  "fuda_rand_uniform: no arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // Return value.
  return Py_BuildValue("d", FUDA::rand_uniform());
}


extern "C" PyObject *
fuda_rand_gauss(PyObject *self, PyObject *args)
{
  /* Return gaussian distributed random number. */

  double mean, sd;
  
  if (!PyArg_ParseTuple(args, "dd", &mean, &sd))
    return NULL;

  // Return value.
  return Py_BuildValue("d", FUDA::rand_gauss(mean,sd));
}


extern "C" PyObject *
fuda_z_distrib_p(PyObject *self, PyObject *args)
{
  /* Return accumulated probability from -oo to z for normal distribution. */

  double z;
  
  if (!PyArg_ParseTuple(args, "d", &z))
    return NULL;

  // Return value.
  return Py_BuildValue("d", FUDA::z_distrib_p(z));
}


extern "C" PyObject *
fuda_z_distrib_crit(PyObject *self, PyObject *args)
{
  /* Return critical value for normal distribution. */

  double p;
  
  if (!PyArg_ParseTuple(args, "d", &p))
    return NULL;

  // Return value.
  return Py_BuildValue("d", FUDA::z_distrib_crit(p));
}


extern "C" PyObject *
fuda_f_distrib_p(PyObject *self, PyObject *args)
{
  /* Return accumulated probability from -oo to f for F distribution. */

  double f;
  int f1, f2;
  
  if (!PyArg_ParseTuple(args, "dii", &f, &f1, &f2))
    return NULL;

  // Return value.
  return Py_BuildValue("d", FUDA::f_distrib_p(f,f1,f2));
}


extern "C" PyObject *
fuda_f_distrib_crit(PyObject *self, PyObject *args)
{
  /* Return critical value for F distribution. */

  double p;
  int f1, f2;
  
  if (!PyArg_ParseTuple(args, "dii", &p, &f1, &f2))
    return NULL;

  // Return value.
  return Py_BuildValue("d", FUDA::f_distrib_crit(p,f1,f2));
}


extern "C" PyObject *
fuda_chi2_distrib_p(PyObject *self, PyObject *args)
{
  /* Return accumulated probability from 0 to x for chi2 distribution. */

  double x;
  int nf;
  
  if (!PyArg_ParseTuple(args, "di", &x, &nf))
    return NULL;

  // Return value.
  return Py_BuildValue("d", FUDA::chi2_distrib_p(x,nf));
}


extern "C" PyObject *
fuda_chi2_distrib_crit(PyObject *self, PyObject *args)
{
  /* Return critical value for chi2 distribution. */

  double p;
  int nf;
  
  if (!PyArg_ParseTuple(args, "di", &p, &nf))
    return NULL;

  // Return value.
  return Py_BuildValue("d", FUDA::chi2_distrib_crit(p,nf));
}



/* Here goes the mapping array for all the python uf functions. */

static PyMethodDef fudaMethods[] = {
  // Ptype.
  {"ptype",  (PyCFunction)fuda_ptype, METH_VARARGS|METH_KEYWORDS},
  {"ptype_del",  fuda_ptype_del, METH_VARARGS},
  {"ptype_del_all",  fuda_ptype_del_all, METH_VARARGS},
  {"ptype_get_all",  fuda_ptype_get_all, METH_VARARGS},
  {"ptype_exists",  fuda_ptype_exists, METH_VARARGS},
  {"ptype_is_referenced",  fuda_ptype_is_referenced, METH_VARARGS},
  //{"ptype",  (PyCFunction)fuda_ptype, METH_VARARGS|METH_KEYWORDS},
  //{"ptype_get",  fuda_ptype_get, METH_VARARGS},

  // Param.
  {"param",  (PyCFunction)fuda_param, METH_VARARGS|METH_KEYWORDS},
  {"param_del",  fuda_param_del, METH_VARARGS},
  {"param_del_all",  fuda_param_del_all, METH_VARARGS},
  {"param_get",  fuda_param_get, METH_VARARGS},
  {"param_get_all",  fuda_param_get_all, METH_VARARGS},
  {"param_init_value",  fuda_param_init_value, METH_VARARGS},
  {"param_exists",  fuda_param_exists, METH_VARARGS},
  {"param_is_referenced",  fuda_param_is_referenced, METH_VARARGS},
  //{"param_get_sync",  fuda_param_get_sync, METH_VARARGS},
  //{"param_get_rsync",  fuda_param_get_rsync, METH_VARARGS},

  // Dtype.
  {"dtype",  fuda_dtype, METH_VARARGS},
  {"dtype_del",  fuda_dtype_del, METH_VARARGS},
  {"dtype_del_all",  fuda_dtype_del_all, METH_VARARGS},
  {"dtype_get",  fuda_dtype_get, METH_VARARGS},
  {"dtype_get_all",  fuda_dtype_get_all, METH_VARARGS},
  {"dtype_exists",  fuda_dtype_exists, METH_VARARGS},
  {"dtype_is_referenced",  fuda_dtype_is_referenced, METH_VARARGS},
  //{"dtype_get_expl",  fuda_dtype_get_expl, METH_VARARGS},
  {"dtype_get_purge_radius",  fuda_dtype_get_purge_radius, METH_VARARGS},
  {"dtype_set_purge_radius",  fuda_dtype_set_purge_radius, METH_VARARGS},
  {"dtype_set_purge",  fuda_dtype_set_purge, METH_VARARGS},

  // Data.
  {"data",  fuda_data, METH_VARARGS},
  {"data_del",  fuda_data_del, METH_VARARGS},
  {"data_del_all",  fuda_data_del_all, METH_VARARGS},

  // Ftype.
  {"ftype_product",  fuda_ftype_product, METH_VARARGS},
  {"ftype_sum",  fuda_ftype_sum, METH_VARARGS},
  {"ftype_composite",  fuda_ftype_composite, METH_VARARGS},
  {"ftype_python",  fuda_ftype_python, METH_VARARGS},
  {"ftype_get",  fuda_ftype_get, METH_VARARGS},
  {"ftype_exists",  fuda_ftype_exists, METH_VARARGS},
  {"ftype_get_all",  fuda_ftype_get_all, METH_VARARGS},
  {"ftype_get_param",  fuda_ftype_get_param, METH_VARARGS},
  {"ftype_get_param_descr",  fuda_ftype_get_param_descr, METH_VARARGS},
  {"ftype_get_var",  fuda_ftype_get_var, METH_VARARGS},
  {"ftype_get_var_index",  fuda_ftype_get_var_index, METH_VARARGS},
  {"ftype_call",  fuda_ftype_call, METH_VARARGS},

  // Func.
  {"func",  fuda_func, METH_VARARGS},
  {"func_exists",  fuda_func_exists, METH_VARARGS},
  {"func_del",  fuda_func_del, METH_VARARGS},
  {"func_del_all",  fuda_func_del_all, METH_VARARGS},
  {"func_use",  fuda_func_use, METH_VARARGS},
  {"func_use_all",  fuda_func_use_all, METH_VARARGS},
  {"func_get_all",  fuda_func_get_all, METH_VARARGS},
  {"func_get_current",  fuda_func_get_current, METH_VARARGS},
  {"func_current",  fuda_func_current, METH_VARARGS},
  {"func_get",  fuda_func_get, METH_VARARGS},
  {"func_get_param",  fuda_func_get_param, METH_VARARGS},
  {"func_get_var",  fuda_func_get_var, METH_VARARGS},
  {"func_set_param",  fuda_func_set_param, METH_VARARGS},
  {"func_call",  fuda_func_call, METH_VARARGS},
  {"func_call_by_expl",  fuda_func_call_by_expl, METH_VARARGS},
  {"func_call_by_var",  fuda_func_call_by_var, METH_VARARGS},
  {"func_call_by_param",  fuda_func_call_by_param, METH_VARARGS},
  {"func_deriv_by_var",  fuda_func_deriv_by_var, METH_VARARGS},

  // Eval.
  {"eval_get",  fuda_eval_get, METH_VARARGS},
  {"eval_init",  fuda_eval_init, METH_VARARGS},
  {"eval_get_free",  fuda_eval_get_free, METH_VARARGS},
  {"eval_get_expl",  fuda_eval_get_expl, METH_VARARGS},
  {"eval_get_const",  fuda_eval_get_const, METH_VARARGS},
  {"eval_get_data",  fuda_eval_get_data, METH_VARARGS},
  {"eval_get_dtype",  fuda_eval_get_dtype, METH_VARARGS},
  {"eval_get_func",  fuda_eval_get_func, METH_VARARGS},
  {"eval_enorm",  fuda_eval_enorm, METH_VARARGS},
  {"eval_data_recalc",  fuda_eval_data_recalc, METH_VARARGS},
  {"eval_data_random",  fuda_eval_data_random, METH_VARARGS},
  //  {"eval_data_get",  fuda_eval_data_get, METH_VARARGS},

  // Lm.
  {"lm",  (PyCFunction)fuda_lm, METH_VARARGS|METH_KEYWORDS},
  {"lm_minimize",  fuda_lm_minimize, METH_VARARGS},
  {"lm_update_param",  fuda_lm_update_param, METH_VARARGS},
  {"lm_get",  fuda_lm_get, METH_VARARGS},
  {"lm_get_value",  fuda_lm_get_value, METH_VARARGS},
  {"lm_get_esd",  fuda_lm_get_esd, METH_VARARGS},
  {"lm_get_covar",  fuda_lm_get_covar, METH_VARARGS},

  // Misc. functions.
  {"rand_set_seed",  fuda_rand_set_seed, METH_VARARGS},
  {"rand_get_seed",  fuda_rand_get_seed, METH_VARARGS},
  {"rand_gauss",  fuda_rand_gauss, METH_VARARGS},
  {"rand_uniform",  fuda_rand_uniform, METH_VARARGS},
  {"z_distrib_p",  fuda_z_distrib_p, METH_VARARGS},
  {"z_distrib_crit",  fuda_z_distrib_crit, METH_VARARGS},
  {"chi2_distrib_p",  fuda_chi2_distrib_p, METH_VARARGS},
  {"chi2_distrib_crit",  fuda_chi2_distrib_crit, METH_VARARGS},
  {"f_distrib_p",  fuda_f_distrib_p, METH_VARARGS},
  {"f_distrib_crit",  fuda_f_distrib_crit, METH_VARARGS},
  {"fuda_print",  fuda_fuda_print, METH_VARARGS},

  // Sentinel.
  {NULL,      NULL}
};

static struct PyModuleDef fudamodule =
{
    PyModuleDef_HEAD_INIT,
    "fudalib", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    fudaMethods
};

PyMODINIT_FUNC PyInit_fudalib(void)
{
   if(PyArray_ImportNumPyAPI() <0 ) {
    return NULL;
  }

  // Create a fuda instance.
   fuda = new Fuda();

  // Create an lm instance.
  lm = new MinpackLM(fuda);
  
  // Initialize random number generator with random seed.
  seed = FUDA::rand_seed();
  FUDA::srand(seed);  

  // Declare ftypes.
  FUDA::declare_cftypes(fuda,"all");

  // Register module functions.
  // (void) Py_InitModule("fudalib", fudaMethods);
  return PyModule_Create(&fudamodule);


  // Import numeric array functions.
//   if(PyArray_API == NULL)
// {
//     import_array(); 
// }
}

