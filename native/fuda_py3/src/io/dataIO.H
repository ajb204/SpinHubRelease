#ifndef _DATAIO_H
#define _DATAIO_H

// The largest number of dimensions a file can have. This value can be
// changed to higher values if wanted without any memory penalty.
#define DATAIO_MAXDIM 20

#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <cstdio>
#include <cmath>

/************************************************************************/
/* dataIO is the base (and abstract) class for all classes for reading
   and writing data files. Currently, only real matrixes are
   supported and only read operations are implemented. */
class DataIO
{
  std::string name; // Name of file.
  unsigned int dim; // Number of dimensions.
  std::vector<unsigned int> size; // Dimensions of data file. 
  unsigned int dim_offset; // if =1, dimensions are numbered 1, 2, 3, ..
  unsigned int pt_offset; // if =1, points are numbered 1, 2, 3, 4, ...
  std::vector<std::vector<double> > xmap;

public:

  // Creator and destructor.
  DataIO(std::string a_name);
  virtual ~DataIO();

  // Object methods.
  void set_dim(unsigned int a_dim);
  unsigned int get_dim();

  void set_dim_offset(unsigned int offset);
  unsigned int get_dim_offset();

  void set_pt_offset(unsigned int offset);
  unsigned int get_pt_offset();

  void set_size(unsigned int idim, unsigned int isize);
  unsigned int get_size(unsigned int idim);

  void print_base();
  void map_linear(unsigned int idim, int pt0, double val0,
			int pt1, double val1);
  void map_point(unsigned int idim, unsigned int pt,
			  double xval);
  double map_p2x(unsigned int idim, unsigned int pt_val);
  unsigned int map_x2p(unsigned int idim, double x_val);
  double read_pt(std::vector<unsigned int> coord);
  void read_float_mx(std::vector<unsigned int> dimvec,
		     std::vector<unsigned int> dimsize,
		     std::vector<unsigned int> coord,
		     std::vector<int> mx_strides,
		     float *mx);  

  virtual std::string get_type() = 0;
  virtual float read_pt_nooffset(std::vector<unsigned int> coord) = 0;
  virtual void print() = 0;
  virtual void set_debug(int dbg) = 0;
  virtual void map_default(unsigned int idim, std::string& unit) = 0;

  // Exceptions.
  struct CannotOpenFile {
    std::string name;
    CannotOpenFile(std::string n) : name(n) {}
  };
  struct InvalidUnitString {
    std::string str;
    InvalidUnitString(std::string s) : str(s) {}
  };
  struct CannotReadHeader{};
  struct HeaderInvalid{};
  struct DimIsComplex {
    int dim;
    DimIsComplex(int d) : dim(d) {}
  };
  struct DimIndexInvalid {
    int index;
    DimIndexInvalid(int i) : index(i) {}
  };
  struct DimInvalid {
    int dim;
    DimInvalid(int i) : dim(i) {}
  };
  struct DimSizeInvalid {
    int size;
    DimSizeInvalid(int i) : size(i) {}
  };
  struct CoordInvalid {
    int dim;
    int coord;
    CoordInvalid(int i, int c) : dim(i), coord(c) {}
  };
  struct UpperBoundsInvalid {
    int dim;
    int coord;
    UpperBoundsInvalid(int i, int c) : dim(i), coord(c) {}
  };
  struct MatrixDimInvalid {
    int dim;
    MatrixDimInvalid(int i) : dim(i) {}
  };
  struct VectorDimInvalid {
    int dim;
    VectorDimInvalid(int i) : dim(i) {}
  };
  struct CoordDimInvalid {
    int dim;
    CoordDimInvalid(int i) : dim(i) {}
  };
  struct PointDimInvalid {
    int dim;
    PointDimInvalid(int i) : dim(i) {}
  };
  struct PointIndexInvalid {
    int dim, index;
    PointIndexInvalid(int d, int i) : dim(d), index(i) {}
  };
  struct FseekError {
    int offset;
    FseekError(int off) : offset(off) {}
  };
  struct FreadError {
    int offset, size, rsize;
    FreadError(int off, int siz, int rsiz) : 
      offset(off), size(siz), rsize(rsiz) {}
  };
  struct InvalidXmapArgs{};
    
};

#endif
