#include <Python.h>
#include <numpy/arrayobject.h>
#include <cstring>
#include "dataIO.H"
#include "dataIO_nmrpipe.H"
#include "dataIO_sparky.H"

// Some global variables.
static std::vector<DataIO*> io;
static char fstr[512];

bool opt_arg(PyObject *args, PyObject *kw, char *key , int iarg)
{
  /* This function checks for the presence of an optional argument
     number iarg with the keyword key by checking the args and kw
     objects. */
  return (kw!=NULL &&
	  (PyTuple_Size(args)>=iarg+1 ||
	   PyDict_GetItemString (kw, key)!=NULL));
}


extern "C" PyObject *
dataIO_open(PyObject *self, PyObject *args, PyObject *kw)
{
  /* dataIO_open opens a new data file */

  // Argument keyword list.
  static char kw_name[] = "name", kw_type[] = "type", kw_mode[] = "mode";
  static char *kwlist[] = {kw_name, kw_type, kw_mode, NULL};
  char *name, *type, *mode;

  // Check for optional keyword arguments.
  int iarg = 1;
  bool type_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;
  bool mode_flg = opt_arg(args, kw, kwlist[iarg],iarg); iarg++;

  // Pass argument list.
  if (!PyArg_ParseTupleAndKeywords(args, kw, "s|ss", kwlist, 
				   &name, &type, &mode)) return NULL; 

  // Get file name.
  std::string fname(name);
  std::cout<< name << std::endl;
  // Get file type.
  std::string ftype;
  if (type_flg)
    ftype = type;
  else
    ftype = "nmrpipe";
  
  std::cout<< ftype << std::endl;

  // Get file access mode.
  std::string fmode;
  if (mode_flg)
    fmode = mode;
  else
    fmode = "r";
  
  // Branch out according to ftype and open file.
  DataIO *io_new = 0;
  if (ftype=="nmrpipe")  
    {      
      try {
	// Mode is not passed as it is not supported - always "r".
	io_new = new DataIO_nmrpipe(fname);
      }
      catch (DataIO::CannotOpenFile& e) {
	std::string err = 
	  "dataIO_open: Cannot open file: " + e.name;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::CannotReadHeader) {
	std::string err = 
	  "dataIO_open: Cannot read file header";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::HeaderInvalid) {
	std::string err = 
	  "dataIO_open: Invalid file header or file format";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::DimIsComplex& e) {
	std::string err = 
	  "dataIO_open: dimension ";
	snprintf(fstr, sizeof(fstr), "%d", e.dim); err += fstr;
	err += " is complex";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::DimIndexInvalid& e) {
	std::string err = 
	  "dataIO_open: DimIndexInvalid: ";
	snprintf(fstr, sizeof(fstr), "%d", e.index); err += fstr;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::InvalidXmapArgs) {
	std::string err = 
	  "dataIO_open: InvalidXmapArgs";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (...) {
	std::string err =  "dataIO_open: Unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
    }
  else if (ftype=="sparky")  
    {      
      try {
	// Mode is not passed as it is not supported - always "r".
	io_new = new DataIO_sparky(fname);
      }
      catch (DataIO::FreadError) {
	std::string err = 
	  "dataIO_open: fread error";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::FseekError) {
	std::string err = 
	  "dataIO_open: fseek error";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::CannotOpenFile& e) {
	std::string err = 
	  "dataIO_open: Cannot open file: " + e.name;
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::CannotReadHeader) {
	std::string err = 
	  "dataIO_open: Cannot read file header";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::HeaderInvalid) {
	std::string err = 
	  "dataIO_open: Invalid file header or file format";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (DataIO::DimInvalid) {
	std::string err = 
	  "dataIO_open: Invalid dimension of data file";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
      catch (...) {
	std::string err =  "dataIO_open: Unspecified exception";
	PyErr_SetString(PyExc_Exception, err.c_str());
	return NULL;
      }
    }
  else
    {
      // Invalid file type.
      std::string err =  "dataIO_open: invalid type: "+ftype+"\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Look for unused io entry and add dataIO reference.
  unsigned int io_index;
  bool found = 0;
  for (unsigned int i=0; i<io.size(); i++)
    if (io[i]==0)
      {
	found = 1;
	io_index = i;	
	io[io_index] = io_new;
	break;
      }

  // If no free entries are found, we extend io by one.
  if (!found) 
    {
      io_index = io.size();
      io.push_back(io_new);
    }
  
  // Return io_index as python integer.
  return Py_BuildValue("i", io_index);
}


extern "C" PyObject *
dataIO_close(PyObject *self, PyObject *args)
{
  /* dataIO_close closes the file. */

  int io_index;

  // Get argument.
  if (!PyArg_ParseTuple(args, "i", &io_index))
    return NULL;

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_close: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Close the file by destructing the io structure.
  try {
    delete io[io_index];
    io[io_index]=0;
  }
  catch (...) {
    std::string err =  "dataIO_close: Unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
dataIO_get_pt_offset(PyObject *self, PyObject *args)
{
  /* dataIO_get_pt_offset gets the point offset for the specified
     io_index file. */

  int io_index, offset;

  // Get argument.
  if (!PyArg_ParseTuple(args, "i", &io_index))
    return NULL;

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_get_pt_offset: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the offset.
  try {
    offset = io[io_index]->get_pt_offset();
  }
  catch (...) {
    std::string err =  "dataIO_get_pt_offset: Unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("i", offset);
}


extern "C" PyObject *
dataIO_get_dim_offset(PyObject *self, PyObject *args)
{
  /* dataIO_get_dim_offset gets the dimension offset for the specified
     io_index file. */

  int io_index, offset;

  // Get argument.
  if (!PyArg_ParseTuple(args, "i", &io_index))
    return NULL;

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_get_pt_offset: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the offset.
  try {
    offset = io[io_index]->get_dim_offset();
  }
  catch (...) {
    std::string err =  "dataIO_get_dim_offset: Unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("i", offset);
}


extern "C" PyObject *
dataIO_get_dim(PyObject *self, PyObject *args)
{
  /* dataIO_get_dim returns the dimension of the specified
     io_index file. */

  int io_index, dim;

  // Get argument.
  if (!PyArg_ParseTuple(args, "i", &io_index))
    return NULL;

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_get_dim: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the dimension.
  try {
    dim = io[io_index]->get_dim();
  }
  catch (...) {
    std::string err =  "dataIO_get_dim: Unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("i", dim);
}


extern "C" PyObject *
dataIO_get_size(PyObject *self, PyObject *args)
{
  /* dataIO_get_size returns the size of for the specified
     io_index file in the specified dimension. */

  int io_index, idim, isize;

  // Get arguments.
  if (!PyArg_ParseTuple(args, "ii", &io_index, &idim))
    return NULL;

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_set_pt_offset: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Check the validity of the idim.
  if (idim<0 || idim>=io[io_index]->get_dim())
    {
      std::string err =  "dataIO_get_size: Invalid dimension\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the size.
  try {
    isize = io[io_index]->get_size(idim);
  }
  catch (...) {
    std::string err =  "dataIO_get_size: Unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return value.
  return Py_BuildValue("i", isize);
}


extern "C" PyObject *
dataIO_set_pt_offset(PyObject *self, PyObject *args)
{
  /* dataIO_set_pt_offset sets the point offset for the specified
     io_index file. */

  int io_index, offset;

  // Get arguments.
  if (!PyArg_ParseTuple(args, "ii", &io_index, &offset))
    return NULL;

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_set_pt_offset: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Set the offset.
  try {
    io[io_index]->set_pt_offset(offset);
  }
  catch (...) {
    std::string err =  "dataIO_set_pt_offset: Unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
dataIO_set_dim_offset(PyObject *self, PyObject *args)
{
  /* dataIO_set_dim_offset sets the dimension offset for the specified
     io_index file. */

  int io_index, offset;

  // Get arguments.
  if (!PyArg_ParseTuple(args, "ii", &io_index, &offset))
    return NULL;

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_set_dim_offset: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Set the offset.
  try {
    io[io_index]->set_dim_offset(offset);
  }
  catch (...) {
    std::string err =  "dataIO_set_dim_offset: Unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  
  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
dataIO_print(PyObject *self, PyObject *args)
{
  /* Call print() for all open dataIO objects */
  for (unsigned int i=0; i<io.size(); i++)
    if (io[i]!=0)
      {
	std::cout << "Python/dataIO file index: " << i << ":\n";
	io[i]->print();
	std::cout << "\n";
      }
  

  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
dataIO_read_pt(PyObject *self, PyObject *args)
{
  /* read_pt reads a data point in the file and returns the value. The
     routine takes two arguments, the io index and a sequence containing the
     coordinate of the point. */

  // First check of number of arguments. We need at least one for the io
  // index and one for the coordinate tuple.
  int narg = PyTuple_Size(args);
  if (narg<2)
    {
      std::string err =  "dataIO_read_pt: at least 2 arguments required: ";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }

  // First argument must be the integer io_index.
  PyObject *io_index_obj = PyTuple_GetItem(args, 0);
  if (!PyLong_Check(io_index_obj))
    {
      std::string err =  "dataIO_read_pt: First argument must "
	"be integer io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Second argument must be the coordinate sequence.
  PyObject *coord_seq_obj = PyTuple_GetItem(args, 1);
  if (!PySequence_Check(coord_seq_obj))
    {
      std::string err =  "dataIO_read_pt: Second arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the io_index.
  int io_index = PyLong_AsLong(io_index_obj);

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_read_pt: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Get the dimension of the data.
  int dim = io[io_index]->get_dim();
  
  // Check size of coordinate sequence.
  int ncoord = PySequence_Length(coord_seq_obj);
  if (ncoord!=dim)
    {
      std::string err =  "dataIO_read_pt: Invalid number of items "
	"in coordinate sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Extract the coordinates from the coordinate sequence.
  // Loop over the items of the params tuple and get parameter names.
  std::vector<unsigned int> coord;
  coord.resize(dim);
  for(int i=0; i<dim; i++)
    {
      // Get item from tuple.
      PyObject *coord_obj = PySequence_GetItem(coord_seq_obj,i);

      // Check it is an integer.
      if (!PyLong_Check(coord_obj))
	{
	  std::string err =  "dataIO_read_pt: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg must be of type integer";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      coord[i] = PyLong_AsLong(coord_obj);
    }

  // Get the point.
  double fval;
  try {
    fval = io[io_index]->read_pt(coord);
  }
  catch (DataIO::PointDimInvalid) {
    std::string err =  "dataIO_read_pt: Point dimension invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (DataIO::PointIndexInvalid) {
      std::string err =  "dataIO_read_pt: Point index invalid";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (DataIO::FseekError) {
      std::string err =  "dataIO_read_pt: Fseek error";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (DataIO::FreadError) {
      std::string err =  "dataIO_read_pt: Fread error";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "dataIO_read_pt: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  

  // Return value.
  return Py_BuildValue("d", fval);
}


extern "C" PyObject *
dataIO_map_default(PyObject *self, PyObject *args)
{
  /* map_default sets the default linear point to user unit mapping
     for a given dimension.  It takes 3 arguments: (int) io_index,
     (int) idim, (string) unit. */

  // First check of number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=3)
    {
      std::string err =  "dataIO_map_default: exactly 3 arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the integer io_index.
  PyObject *io_index_obj = PyTuple_GetItem(args, 0);
  if (!PyLong_Check(io_index_obj))
    {
      std::string err =  "dataIO_map_default: First argument must "
	"be integer io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the io_index.
  int io_index = PyLong_AsLong(io_index_obj);

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_map_default: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Get the dimension of the data and offsets.
  int dim = io[io_index]->get_dim();
  int dim_offset = io[io_index]->get_dim_offset();
  int pt_offset = io[io_index]->get_pt_offset();

  // 2nd argument must be the dimension in which to set the mapping.
  PyObject *idim_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(idim_obj))
    {
      std::string err =  "dataIO_map_default: 2nd argument must "
	"be integer idim";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the idim.
  int idim = PyLong_AsLong(idim_obj);

  // Check the validity of the index.
  if (idim-dim_offset<0 ||
      idim-dim_offset>=dim)
    {
      std::string err =  "dataIO_map_default: Invalid idim\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 3rd argument is the unit.
  PyObject *unit_obj = PyTuple_GetItem(args, 2);
  if (!PyUnicode_Check(unit_obj))
    {
      std::string err =  "dataIO_map_default: 3rd arg. must be a string";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the unit.
  const char *unit_cstr = PyUnicode_AsUTF8(unit_obj);
  std::string unit_str = unit_cstr;

  // Set the linear mapping.
  try {
    io[io_index]->map_default(idim, unit_str);
  }
  catch (DataIO::DimIndexInvalid) {
    std::string err =  "dataIO_map_default: Dimension index invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (DataIO::InvalidUnitString) {
    std::string err =  "dataIO_map_default: Invalid unit string";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (...) {
      std::string err =  "dataIO_map_linear: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  
  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
dataIO_map_linear(PyObject *self, PyObject *args)
{
  /* map_linear sets a linear point to user unit mapping for a given
     function.  It takes 6 arguments: (int) io_index, (int) idim,
     (int) pt0, (double) x0, (int) pt1, (double) x1. */

  // First check of number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=6)
    {
      std::string err =  "dataIO_map_linear: exactly 6 arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the integer io_index.
  PyObject *io_index_obj = PyTuple_GetItem(args, 0);
  if (!PyLong_Check(io_index_obj))
    {
      std::string err =  "dataIO_map_linear: First argument must "
	"be integer io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the io_index.
  int io_index = PyLong_AsLong(io_index_obj);

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_map_linear: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Get the dimension of the data and offsets.
  int dim = io[io_index]->get_dim();
  int dim_offset = io[io_index]->get_dim_offset();
  int pt_offset = io[io_index]->get_pt_offset();

  // 2nd argument must be the dimension in which to set the mapping.
  PyObject *idim_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(idim_obj))
    {
      std::string err =  "dataIO_map_linear: 2nd argument must "
	"be integer idim";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the idim.
  int idim = PyLong_AsLong(idim_obj);

  // Check the validity of the index.
  if (idim-dim_offset<0 ||
      idim-dim_offset>=dim)
    {
      std::string err =  "dataIO_map_linear: Invalid idim\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 3rd argument is pt0.
  PyObject *pt0_obj = PyTuple_GetItem(args, 2);
  if (!PyLong_Check(pt0_obj))
    {
      std::string err =  "dataIO_map_linear: 3nd argument must "
	"be integer pt0";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the pt0.
  int pt0 = PyLong_AsLong(pt0_obj);

  // 4th argument is x0.
  PyObject *x0_obj = PyTuple_GetItem(args, 3);
  if (!PyFloat_Check(x0_obj))
    {
      std::string err =  "dataIO_map_linear: 4th argument must "
	"be double x0";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the x0.
  double x0 = PyFloat_AsDouble(x0_obj);

  // 5th argument is pt1.
  PyObject *pt1_obj = PyTuple_GetItem(args, 4);
  if (!PyLong_Check(pt1_obj))
    {
      std::string err =  "dataIO_map_linear: 5nd argument must "
	"be integer pt1";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the pt1.
  int pt1 = PyLong_AsLong(pt1_obj);

  // 6th argument is x1.
  PyObject *x1_obj = PyTuple_GetItem(args, 5);
  if (!PyFloat_Check(x1_obj))
    {
      std::string err =  "dataIO_map_linear: 6th argument must "
	"be double x1";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the x1.
  double x1 = PyFloat_AsDouble(x1_obj);

  // Set the linear mapping.
  try {
    io[io_index]->map_linear(idim, pt0, x0, pt1, x1);
  }
  catch (DataIO::DimIndexInvalid) {
    std::string err =  "dataIO_map_linear: Dimension index invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (DataIO::InvalidXmapArgs) {
      std::string err =  "dataIO_map_linear: Invalid map arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "dataIO_map_linear: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  
  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
dataIO_map_point(PyObject *self, PyObject *args)
{
  /* map_point sets one entry in the point to user mapping for a given
     dimension.  It takes 4 arguments: (int) io_index, (int) idim,
     (int) pt, (double) x. */

  // First check of number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=4)
    {
      std::string err =  "dataIO_map_point: exactly 4 arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the integer io_index.
  PyObject *io_index_obj = PyTuple_GetItem(args, 0);
  if (!PyLong_Check(io_index_obj))
    {
      std::string err =  "dataIO_xpoint: First argument must "
	"be integer io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the io_index.
  int io_index = PyLong_AsLong(io_index_obj);

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_map_point: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Get the dimension of the data and offsets.
  int dim = io[io_index]->get_dim();
  int dim_offset = io[io_index]->get_dim_offset();
  int pt_offset = io[io_index]->get_pt_offset();

  // 2nd argument must be the dimension in which to set the mapping.
  PyObject *idim_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(idim_obj))
    {
      std::string err =  "dataIO_map_point: 2nd argument must "
	"be integer idim";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the idim.
  int idim = PyLong_AsLong(idim_obj);

  // Check the validity of the index.
  if (idim-dim_offset<0 ||
      idim-dim_offset>=dim)
    {
      std::string err =  "dataIO_map_point: Invalid idim\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 3rd argument is pt.
  PyObject *pt_obj = PyTuple_GetItem(args, 2);
  if (!PyLong_Check(pt_obj))
    {
      std::string err =  "dataIO_map_point: 3nd argument must "
	"be integer pt";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the pt.
  int pt = PyLong_AsLong(pt_obj);

  // 4th argument is x.
  PyObject *x_obj = PyTuple_GetItem(args, 3);
  if (!PyFloat_Check(x_obj))
    {
      std::string err =  "dataIO_map_point: 4th argument must "
	"be double x";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the x.
  double x = PyFloat_AsDouble(x_obj);

  // Set the point.
  try {
    io[io_index]->map_point(idim, pt, x);
  }
  catch (DataIO::DimIndexInvalid) {
    std::string err =  "dataIO_map_point: Dimension index invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (DataIO::InvalidXmapArgs) {
      std::string err =  "dataIO_map_point: Invalid map arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "dataIO_map_point: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  
  // Return none.
  Py_INCREF(Py_None);
  return Py_None;
}


extern "C" PyObject *
dataIO_map_x2p(PyObject *self, PyObject *args)
{
  /* map_x2p returns the closest point index for a given
     x-value. It takes 3 arguments: (int) io_index, (int) idim,
     (double) x. */

  // First check of number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=3)
    {
      std::string err =  "dataIO_map_x2p: exactly 3 arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the integer io_index.
  PyObject *io_index_obj = PyTuple_GetItem(args, 0);
  if (!PyLong_Check(io_index_obj))
    {
      std::string err =  "dataIO_map_x2p: First argument must "
	"be integer io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the io_index.
  int io_index = PyLong_AsLong(io_index_obj);

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_map_x2p: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Get the dimension of the data and offsets.
  int dim = io[io_index]->get_dim();
  int dim_offset = io[io_index]->get_dim_offset();
  int pt_offset = io[io_index]->get_pt_offset();

  // 2nd argument must be the dimension in which to set the mapping.
  PyObject *idim_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(idim_obj))
    {
      std::string err =  "dataIO_map_x2p: 2nd argument must "
	"be integer idim";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the idim.
  int idim = PyLong_AsLong(idim_obj);

  // Check the validity of the index.
  if (idim-dim_offset<0 ||
      idim-dim_offset>=dim)
    {
      std::string err =  "dataIO_map_x2p: Invalid idim\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 3th argument is x.
  PyObject *x_obj = PyTuple_GetItem(args, 2);
  if (!PyFloat_Check(x_obj))
    {
      std::string err =  "dataIO_map_x2p: 3th argument must "
	"be double x";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the x.
  double x = PyFloat_AsDouble(x_obj);
  int pt;

  // Convert to point index.
  try {
    pt = io[io_index]->map_x2p(idim, x);
  }
  catch (DataIO::DimIndexInvalid) {
    std::string err =  "dataIO_map_x2p: Dimension index invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (DataIO::InvalidXmapArgs) {
      std::string err =  "dataIO_map_x2p: Invalid map arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "dataIO_map_x2p: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  
  // Return pt as python integer.
  return Py_BuildValue("i", pt);
}


extern "C" PyObject *
dataIO_map_p2x(PyObject *self, PyObject *args)
{
  /* map_p2x returns the x value of the given point index. It
     takes 3 arguments: (int) io_index, (int) idim, (int) pt. */

  // First check of number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=3)
    {
      std::string err =  "dataIO_map_p2x: exactly 3 arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the integer io_index.
  PyObject *io_index_obj = PyTuple_GetItem(args, 0);
  if (!PyLong_Check(io_index_obj))
    {
      std::string err =  "dataIO_map_p2x: First argument must "
	"be integer io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the io_index.
  int io_index = PyLong_AsLong(io_index_obj);

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_map_p2x: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Get the dimension of the data and offsets.
  int dim = io[io_index]->get_dim();
  int dim_offset = io[io_index]->get_dim_offset();
  int pt_offset = io[io_index]->get_pt_offset();

  // 2nd argument must be the dimension in which to set the mapping.
  PyObject *idim_obj = PyTuple_GetItem(args, 1);
  if (!PyLong_Check(idim_obj))
    {
      std::string err =  "dataIO_map_p2x: 2nd argument must "
	"be integer idim";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the idim.
  int idim = PyLong_AsLong(idim_obj);

  // Check the validity of the index.
  if (idim-dim_offset<0 ||
      idim-dim_offset>=dim)
    {
      std::string err =  "dataIO_map_p2x: Invalid idim\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 3th argument is pt.
  PyObject *pt_obj = PyTuple_GetItem(args, 2);
  if (!PyLong_Check(pt_obj))
    {
      std::string err =  "dataIO_map_p2x: 3th argument must "
	"be integer pt";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the pt.
  int pt = PyLong_AsLong(pt_obj);
  double x;

  // Convert to point index.
  try {
    x = io[io_index]->map_p2x(idim, pt);
  }
  catch (DataIO::DimIndexInvalid) {
    std::string err =  "dataIO_map_p2x: Dimension index invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (DataIO::InvalidXmapArgs) {
      std::string err =  "dataIO_map_p2x: Invalid map arguments";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "dataIO_map_p2x: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  
  // Return x as python float.
  return Py_BuildValue("d", x);
}


extern "C" PyObject *
dataIO_read_mx(PyObject *self, PyObject *args)
{
  /* read_mx reads a submatrix of data points in the file and returns
     the matrix as a PyArrayObject of floats.

     Arguments:

     io_index : integer dataIO io index.

     dim_seq : sequence of integers with the dimensions to read from in
     the file. E.g., if we have a 3D file, then the dim_seq "(1,2)"
     will give us a 2D matrix with data from dimension 1 in the file
     in the first dimension of the matrix and data from dimension 2 in
     the file in the second dimension of the matrix. If "(2,1)" is
     specified the data are transposed in the matrix, so data in
     dimension 2 in the file will be the first dimension in the
     matrix. If "(1,3)" is specified a subplane along dimension 1 and
     3 in the data file is returned. The dim_tuple can have between 1
     and dim elements, where dim is the dimension of the file
     data. E.g., "(1)" returns a 1D slice (along the first dimension
     of the file) from the file. "(1,2,3)" returns a sub cube of data
     from the file.

     dim_size_seq: sequence of integers with the sizes of the matrix to
     return. dim_size_seq must have the same number of elements as
     dim_seq.

     coord_seq: tuple of integers specifying the coordinate in the
     data file to start reading from. */

  // Check number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=4)
    {
      std::string err =  "dataIO_read_mx: exactly 4 arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
      
    }
  
  // First argument must be the integer io_index.
  PyObject *io_index_obj = PyTuple_GetItem(args, 0);
  if (!PyLong_Check(io_index_obj))
    {
      std::string err =  "dataIO_read_mx: First argument must "
	"be integer io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the io_index.
  int io_index = PyLong_AsLong(io_index_obj);

  // Check the validity of the io_index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_read_mx: Invalid io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get reference to the relevant dataIO object.
  DataIO& data_io = *(io[io_index]);
  
  // Get the dimension of the data.
  int dim = data_io.get_dim();

  // Get dim and pt offsets.
  unsigned int dim_offset = data_io.get_dim_offset();
  unsigned int pt_offset = data_io.get_pt_offset();

  // Second argument must be the dim_seq.
  PyObject *dim_seq = PyTuple_GetItem(args, 1);
  if (!PySequence_Check(dim_seq))
    {
      std::string err =  "dataIO_read_mx: Second arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Check size of dim_seq which is the dimension of the matrix to
  // return.
  int mxdim = PySequence_Length(dim_seq);
  if (mxdim<1 || mxdim>dim)
    {
      std::string err =  "dataIO_read_mx: Invalid dim_seq size";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Extract and check the dim_seq values and setup the dimvec vector.
  std::vector<unsigned int> dimvec;
  dimvec.resize(dim);
  for(int i=0; i<mxdim; i++)
    {
      // Get item from tuple.
      PyObject *dim_obj = PySequence_GetItem(dim_seq,i);

      // Check it is an integer.
      if (!PyLong_Check(dim_obj))
	{
	  std::string err =  "dataIO_read_mx: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of dim_seq must be of type integer";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the dim value.
      dimvec[i] = PyLong_AsLong(dim_obj);

      // Check the dim value.
      if (dimvec[i]-dim_offset<0 || dimvec[i]-dim_offset>=dim)
	{
	  std::string err =  "dataIO_read_mx: value of item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " in dim_tuple is invalid";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Check that it has not been used before.
      for (int previous=0;previous<i;previous++)
	if (dimvec[i]==dimvec[previous])
	  {
	    std::string err =  "dataIO_read_mx: value of item ";
	    snprintf(fstr, sizeof(fstr),"%d",i);
	    err += fstr;
	    err += " in dim_tuple is already used in dim_tuple";
	    PyErr_SetString(PyExc_Exception, err.c_str());
	    return NULL;
	  }
    }

  // 3rd argument must be the dim_size_seq.
  PyObject *dimsize_seq = PyTuple_GetItem(args, 2);
  if (!PySequence_Check(dimsize_seq))
    {
      std::string err =  "dataIO_read_mx: 3rd arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Check size of dim_size_seq which must be mxdim.
  if (mxdim!=PySequence_Length(dimsize_seq))
    {
      std::string err =  "dataIO_read_mx: Invalid dimsize_seq size";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Extract and check the dim_size_seq values and setup the
  // dimsize vector.
  std::vector<unsigned int> dimsize;
  dimsize.resize(mxdim);
  for(int i=0; i<mxdim; i++)
    {
      // Get item from sequence.
      PyObject *dimsize_obj = PySequence_GetItem(dimsize_seq,i);

      // Check it is an integer.
      if (!PyLong_Check(dimsize_obj))
	{
	  std::string err =  "dataIO_read_mx: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of dimsize_seq must be of type integer";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the dimsize value.
      dimsize[i] = PyLong_AsLong(dimsize_obj);

      // Check the dimsize value (must be positive).
      if (dimsize[i]<1)
	{
	  std::string err =  "dataIO_read_mx: value of item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " in dimsize_seq is invalid";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
    }

  // 4th argument must be the coord_seq.
  PyObject *coord_seq = PyTuple_GetItem(args, 3);
  if (!PySequence_Check(coord_seq))
    {
      std::string err =  "dataIO_read_mx: 4th arg. must be a sequence";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Check size of coord_seq.
  if (dim!=PySequence_Length(coord_seq))
    {
      std::string err =  "dataIO_read_mx: In valid number of items "
	"in coord_seq";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Extract the coordinates from the coord_seq.
  std::vector<unsigned int> coord;
  coord.resize(dim);
  for(int i=0; i<dim; i++)
    {
      // Get item from tuple.
      PyObject *coord_obj = PySequence_GetItem(coord_seq,i);

      // Check it is an integer.
      if (!PyLong_Check(coord_obj))
	{
	  std::string err =  "dataIO_read_mx: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg must be of type integer";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the coord value.
      coord[i] = PyLong_AsLong(coord_obj);

      // Check that the coord is within bounds.
      if (coord[i]<pt_offset || 
	  coord[i]>=data_io.get_size(i+dim_offset)+pt_offset)
	{
	  std::string err =  "dataIO_read_mx: value of item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " in coord_seq is out of bounds";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
    }

  // Check the upper bounds of the matrix to read in.
  for(int i=0; i<mxdim; i++)
    {
      unsigned int upper = coord[dimvec[i]-dim_offset]+dimsize[i]-1;
      if(upper<1 || upper>data_io.get_size(dimvec[i]))
	{
	  std::string err =  "dataIO_read_mx: matrix is out of bounds in ";
	  snprintf(fstr, sizeof(fstr),"matrix dimension %d (file dimension %d)",
		  i,dimvec[i]);
	  err += fstr;
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
    }

  // Create python float array object.
  npy_intp *dimsize_ptr = new npy_intp[mxdim];
  for (unsigned int i=0; i<mxdim; i++)
    dimsize_ptr[i] = dimsize[i];
  //PyArrayObject *mx = 
  //(PyArrayObject*) PyArray_FromDims(mxdim,dimsize_ptr,PyArray_FLOAT);

  PyArrayObject *mx =
    (PyArrayObject*) PyArray_SimpleNew(mxdim,dimsize_ptr,NPY_FLOAT);
  
  delete[] dimsize_ptr;

  // Setup the matrix strides vector from mx->strides[].
  std::vector<int> mx_strides;
  mx_strides.resize(mxdim);
  int char_per_float = sizeof(float)/sizeof(char);
  for (int i=0; i<mxdim; i++)
    mx_strides[i] = mx->strides[i]/char_per_float;

  try {
    data_io.read_float_mx(dimvec, dimsize, coord,
			  mx_strides, (float*)(mx->data));
  }
  catch (...) {
    std::string err =  "dataIO_read_mx: unspecified exception";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }

  // Return the matrix.
  return (PyArray_Return(mx));  
}



/* Here are some routines that have been renamed or have been
   superseeded by others and will be discontinued in the future. */

extern "C" PyObject *
dataIO_set_xmap(PyObject *self, PyObject *args)
{
  return dataIO_map_linear(self,args);
}


extern "C" PyObject *
dataIO_set_xval(PyObject *self, PyObject *args)
{
  return dataIO_map_point(self,args);
}


extern "C" PyObject *
dataIO_map_x_to_pt(PyObject *self, PyObject *args)
{
  return dataIO_map_x2p(self,args);
}


extern "C" PyObject *
dataIO_map_pt_to_x(PyObject *self, PyObject *args)
{
  return dataIO_map_p2x(self,args);
}


extern "C" PyObject *
dataIO_read_fuda_pt(PyObject *self, PyObject *args)
{
  // read_fuda_pt has been replaced by data_read_pt in the fudaIO module.

  /* read_fuda_pt reads a data point in the file and returns the value
     as part of a tuple which can be used directly as a data tuple in
     fuda. The routine takes three arguments, the io index,a tuple
     containing the integer coordinates of the point and the
     uncertainty of the data point. The routine returns a tuple
     containing the xmapped coordinates, the point value and the
     uncertainty */

  // First check of number of arguments.
  int narg = PyTuple_Size(args);
  if (narg!=3)
    {
      std::string err =  "dataIO_read_pt: exactly 3 arguments required";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // First argument must be the integer io_index.
  PyObject *io_index_obj = PyTuple_GetItem(args, 0);
  if (!PyLong_Check(io_index_obj))
    {
      std::string err =  "dataIO_read_pt: First argument must "
	"be integer io index";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Second argument must be the coordinate tuple.
  PyObject *coord_tup_obj = PyTuple_GetItem(args, 1);
  if (!PyTuple_Check(coord_tup_obj))
    {
      std::string err =  "dataIO_read_pt: Second arg. must be a tuple";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // 3rd argument must be the data uncertainty.
  PyObject *uncertainty_obj = PyTuple_GetItem(args, 2);
  if (!PyFloat_Check(uncertainty_obj))
    {
      std::string err =  "dataIO_read_pt: 3rd argument must "
	"be double (data uncertainty)";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Get the io_index.
  int io_index = PyLong_AsLong(io_index_obj);

  // Check the validity of the index.
  if (io_index<0 || io_index>=io.size() || io[io_index]==0)
    {
      std::string err =  "dataIO_read_pt: Invalid io index\n";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Get the dimension of the data.
  int dim = io[io_index]->get_dim();
  
  // Check size of coordinate tuple.
  int ncoord = PyTuple_Size(coord_tup_obj);
  if (ncoord!=dim)
    {
      std::string err =  "dataIO_read_pt: Invalid number of items "
	"in coordinate tuple";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }

  // Extract the coordinates from the coordinate tuple.
  // Loop over the items of the params tuple and get parameter names.
  std::vector<unsigned int> coord;
  coord.resize(dim);
  for(int i=0; i<PyTuple_Size(coord_tup_obj); i++)
    {
      // Get item from tuple.
      PyObject *coord_obj = PyTuple_GetItem(coord_tup_obj,i);

      // Check it is an integer.
      if (!PyLong_Check(coord_obj))
	{
	  std::string err =  "dataIO_read_pt: item ";
	  snprintf(fstr, sizeof(fstr),"%d",i);
	  err += fstr;
	  err += " of 2. arg must be of type integer";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}

      // Get the value.
      coord[i] = PyLong_AsLong(coord_obj);
    }

  // Get the point.
  double fval;
  try {
    fval = io[io_index]->read_pt(coord);
  }
  catch (DataIO::PointDimInvalid) {
    std::string err =  "dataIO_read_pt: Point dimension invalid";
    PyErr_SetString(PyExc_Exception, err.c_str());
    return NULL;
  }
  catch (DataIO::PointIndexInvalid) {
      std::string err =  "dataIO_read_pt: Point index invalid";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (DataIO::FseekError) {
      std::string err =  "dataIO_read_pt: Fseek error";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (DataIO::FreadError) {
      std::string err =  "dataIO_read_pt: Fread error";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }
  catch (...) {
      std::string err =  "dataIO_read_pt: unspecified exception";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
  }

  // Get the uncertainty.
  double uncertainty = PyFloat_AsDouble(uncertainty_obj);
  
  // Map point coordinates to user values.
  std::vector<double> xcoord;
  xcoord.resize(dim);
  for (int idim=0; idim<dim; idim++)
    xcoord[idim] = io[io_index]->map_p2x(idim,coord[idim]);
  
  // Create the return tuple.
  int nelement = dim+2; // the coordinates + value + uncertainty.
  PyObject* ot_obj = PyTuple_New (nelement);

  // Add coordinates to the return tuple.
  for (int idim=0; idim<dim; idim++)
    {
      PyObject* coord_obj = PyFloat_FromDouble(xcoord[idim]);

      // Insert in ot_obj.
      if (PyTuple_SetItem (ot_obj, idim, coord_obj)!=0)
	{
	  std::string err =  "dataIO_read_fuda_pt : panic - ";
	  err += "PyTuple_SetItem failed. Internal inconsistency";
	  PyErr_SetString(PyExc_Exception, err.c_str());
	  return NULL;
	}
    }

  // Add data value to return tuple.
  PyObject* fval_obj = PyFloat_FromDouble(fval);

  // Insert in ot_obj.
  if (PyTuple_SetItem (ot_obj, dim, fval_obj)!=0)
    {
      std::string err =  "dataIO_read_fuda_pt : panic - ";
      err += "PyTuple_SetItem failed. Internal inconsistency";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Add data uncertainty to return tuple.
  PyObject* u_obj = PyFloat_FromDouble(uncertainty);

  // Insert in ot_obj.
  if (PyTuple_SetItem (ot_obj, dim+1, u_obj)!=0)
    {
      std::string err =  "dataIO_read_fuda_pt : panic - ";
      err += "PyTuple_SetItem failed. Internal inconsistency";
      PyErr_SetString(PyExc_Exception, err.c_str());
      return NULL;
    }
  
  // Return tuple.
  return (ot_obj);
}






/* Here goes the mapping array for all the dataIO functions. */

static PyMethodDef dataIO_Methods[] = {
  {"open",  (PyCFunction)dataIO_open, METH_VARARGS|METH_KEYWORDS},
  {"close",  dataIO_close, METH_VARARGS},
  {"get_dim",  dataIO_get_dim, METH_VARARGS},
  {"get_pt_offset",  dataIO_get_pt_offset, METH_VARARGS},
  {"get_dim_offset",  dataIO_get_dim_offset, METH_VARARGS},
  {"get_size",  dataIO_get_size, METH_VARARGS},
  {"set_pt_offset",  dataIO_set_pt_offset, METH_VARARGS},
  {"set_dim_offset",  dataIO_set_dim_offset, METH_VARARGS},
  {"read_pt",  dataIO_read_pt, METH_VARARGS},
  {"read_mx",  dataIO_read_mx, METH_VARARGS},
  {"map_default",  dataIO_map_default, METH_VARARGS},
  {"map_linear",  dataIO_map_linear, METH_VARARGS},
  {"map_point",  dataIO_map_point, METH_VARARGS},
  {"map_p2x",  dataIO_map_p2x, METH_VARARGS},
  {"map_x2p",  dataIO_map_x2p, METH_VARARGS},

  // Routines to be phased out.
  {"set_xmap",  dataIO_set_xmap, METH_VARARGS},
  {"set_xval",  dataIO_set_xval, METH_VARARGS},
  {"map_pt_to_x",  dataIO_map_pt_to_x, METH_VARARGS},
  {"map_x_to_pt",  dataIO_map_x_to_pt, METH_VARARGS},
  {"read_fuda_pt",  dataIO_read_fuda_pt, METH_VARARGS},
  {"info",  dataIO_print, METH_VARARGS},

  // Sentinel.
  {NULL,      NULL}
};


static struct PyModuleDef dataIOmodule =
{
    PyModuleDef_HEAD_INIT,
    "dataIO", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    dataIO_Methods
};

PyMODINIT_FUNC PyInit_dataIO()
{

  if(PyArray_ImportNumPyAPI() <0 ) {
    return NULL;
  }

  io.resize(1);
  // (void) Py_InitModule("dataIO", dataIO_Methods);
  return PyModule_Create(&dataIOmodule);

}

