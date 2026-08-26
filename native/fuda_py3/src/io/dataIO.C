
#include <string>
#include "dataIO.H"


DataIO::DataIO(std::string a_name) : name(a_name),
				    dim(0),
				    dim_offset(0),
				    pt_offset(0) 
{
}


DataIO::~DataIO()
{
}


void DataIO::set_dim(unsigned int a_dim)
{
  // Check dim.
  if (a_dim<1 || a_dim>DATAIO_MAXDIM)
    throw DataIO::DimInvalid(a_dim);

  // Set dimension of data file.
  dim = a_dim;

  // Resize size and xmap vectors according to dim.
  size.resize(dim);
  xmap.resize(dim);
}


unsigned int DataIO::get_dim()
{
  return (dim);
}


void DataIO::set_dim_offset(unsigned int offset)
{
  dim_offset = offset;
}


unsigned int DataIO::get_dim_offset()
{
  return (dim_offset);
}


void DataIO::set_pt_offset(unsigned int offset)
{
  pt_offset = offset;
}


unsigned int DataIO::get_pt_offset()
{
  return (pt_offset);
}


void DataIO::set_size(unsigned int idim, unsigned int a_size)
{
  if (idim-dim_offset<0 || idim-dim_offset>=dim)
    throw DataIO::DimIndexInvalid(idim);
  size[idim-dim_offset] = a_size;

  // Resize xmap.
  xmap[idim-dim_offset].resize(a_size);
}


unsigned int DataIO::get_size(unsigned int idim)
{
  if (idim-dim_offset<0 || idim-dim_offset>=dim)
    throw DataIO::DimIndexInvalid(idim);
  return (size[idim-dim_offset]);
}


void DataIO::print_base()
{
  std::cout << "name:      " << name << "\n";
  std::cout << "dim:       " << dim << "\n";
  std::cout << "dim_offset:" << dim_offset << "\n";
  std::cout << "pt_offset: " << pt_offset << "\n";
  for (unsigned int i=0; i<dim; i++)
    std::cout << "size[" << i << "]:  " << size[i] << "\n";
}


void DataIO::map_linear(unsigned int idim, 
			int pt0, double val0,
			int pt1, double val1)
{
  // Check idim.
  if (idim-dim_offset<0 || idim-dim_offset>=dim)
    throw DataIO::DimIndexInvalid(idim);
 
  // Check that pt0 and pt1 are not identical.
  if (pt0==pt1) 
    throw DataIO::InvalidXmapArgs();

  // Offset correct all values.
  idim -= dim_offset;
  pt0 -= pt_offset;
  pt1 -= pt_offset;
  
  // Setup xmap[idim].
  for (int i=0; i<xmap[idim].size(); i++)
      xmap[idim][i] = 
	(i-pt0) * (val1-val0)/(pt1-pt0) + val0;
}


void DataIO::map_point(unsigned idim, 
		       unsigned int pt,
		       double xval)
{
  // Check idim.
  if (idim-dim_offset<0 || idim-dim_offset>=dim)
    throw DataIO::DimIndexInvalid(idim);
 
  // Check that pt is within bounds.
  if (pt-pt_offset<0 || pt-pt_offset>=size[idim-dim_offset]) 
    throw DataIO::PointIndexInvalid(idim,pt);

  // Offset correct values.
  idim -= dim_offset;
  pt -= pt_offset;

  // Set the value.
  xmap[idim][pt] = xval;
}


double DataIO::map_p2x(unsigned idim, 
		       unsigned int pt)
{
  // Check idim.
  if (idim-dim_offset<0 || idim-dim_offset>=dim)
    throw DataIO::DimIndexInvalid(idim);
 
  // Check that pt is within bounds.
  if (pt-pt_offset<0 || pt-pt_offset>=size[idim-dim_offset]) 
    throw DataIO::PointIndexInvalid(idim,pt);

  // Offset correct values.
  idim -= dim_offset;
  pt -= pt_offset;

  // Return the value.
  return (xmap[idim][pt]);
}


unsigned int DataIO::map_x2p(unsigned idim, 
			     double x)
{
  // Check idim.
  if (idim-dim_offset<0 || idim-dim_offset>=dim)
    throw DataIO::DimIndexInvalid(idim);
 
  // Search for closest matching point.
  unsigned int pt = 0;
  double min_adif = fabs(xmap[idim][pt]-x);
  for (int i=1; i<size[idim]; i++)
    {
      double adif = fabs(xmap[idim][i]-x);
      if (adif<min_adif)
	{
	  min_adif = adif;
	  pt = i;
	}
    }

  // Return offset corrected point index.
  return (pt+pt_offset);
}


void DataIO::read_float_mx(std::vector<unsigned int> dimvec,
			   std::vector<unsigned int> dimsize,
			   std::vector<unsigned int> coord,
			   std::vector<int> mx_strides,
			   float *mx)  
{
  /* read_float_mx reads a submatrix of data points into dbuf[] which
     is interpreted as a dimvec.size() dimensional data
     matrix. read_mx is a generic routine which calls read_pt many
     times and inherited classes implementing particular matrix file
     formats may supply tailored versions which are more efficient.

     Arguments:

     dimvec : vector of integers with the dimensions to read from in
     the file. E.g., if we have a 3D file, then the dimvec "(1,2)"
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
     from the file. dimvec takes into account dim_offset.

     dimsize: vector of integers with the sizes of the matrix to
     return. dimsize must have the same number of elements as
     dimvec.

     coord: vector of integers specifying the coordinate in the
     data file to start reading from.

     mx_strides: vector of matrix dimension increments. I.e. the number of
     positions (doubles) to move in mx to get to the next address in a
     particular dimension.  */

  // Get dimension of data matrix.
  int mxdim = dimvec.size();
  if (mxdim<1 || mxdim>dim) throw DataIO::MatrixDimInvalid(mxdim);

  // mx_strides must have mxdim elements.
  if (mx_strides.size()!=mxdim) 
    throw DataIO::VectorDimInvalid(mx_strides.size());

  // Check the dimvec values.
  for(int i=0; i<mxdim; i++)
    {
      if (dimvec[i]-dim_offset<0 || dimvec[i]-dim_offset>=dim)
	throw DataIO::DimIndexInvalid(dimvec[i]);

      // Check that it has not been used before.
      for (int previous=0;previous<i;previous++)
	if (dimvec[i]==dimvec[previous])
	  throw DataIO::DimIndexInvalid(dimvec[i]);
    }

  // Correct dimvec for offset.
  for (int i=0; i<dim; i++)
    dimvec[i] -=dim_offset;  

  // dimsize vector must have same dimensions.
  if (dimsize.size()!=mxdim)
    throw DataIO::VectorDimInvalid(dimsize.size());
  
  // Check the dimsize values.
  for(int i=0; i<mxdim; i++)
    if (dimsize[i]<1) throw DataIO::DimSizeInvalid(dimsize[i]);

  // coord vectors must have as many dimensions at the data file.
  if (coord.size()!=dim)
    throw DataIO::CoordDimInvalid(coord.size());
  
  // Check coordinate.
  for(int i=0; i<dim; i++)
    if (coord[i]-pt_offset<0 || 
	coord[i]-pt_offset>=size[i])
      throw DataIO::CoordInvalid(i,coord[i]-pt_offset);

  // Correct coord for offset.
  for (int i=0; i<dim; i++)
    coord[i] -=pt_offset;  

  // Check the upper bounds of the matrix to read in.
  for(int i=0; i<mxdim; i++)
    {
      unsigned int upper = coord[dimvec[i]]+dimsize[i]-1;
      if(upper<0 || upper>=size[dimvec[i]])
	throw DataIO::UpperBoundsInvalid(i, upper);
    }

  // Setup some helper vectors for reading in the matrix.

  // is_mxdim is a vector of same dimension as the file and with
  // elements true if the dimension is a dimension in mx. mx_dimvec is
  // a vector of same dimension as the file and containing the the
  // indices of the corresponding dimensions in mx. The values of
  // mxdimvec are only defined if the given dimension of the file is
  // also a dimension in mx, i.e. if the corresponding element in
  // is_mxvec is true.
  std::vector<bool> is_mxdim(dim,0);
  std::vector<unsigned int> mxdimvec(dim,0);
  for (unsigned int i=0; i<mxdim; i++)
    {
      is_mxdim[dimvec[i]] = 1;
      mxdimvec[dimvec[i]] = i;      
    }

  // last_coord is the coordinate of the last point to read in.
  std::vector<unsigned int> last_coord = coord;
  for (unsigned int i=0; i<dim; i++)
    if (is_mxdim[i]) last_coord[i] += dimsize[mxdimvec[i]]-1;

  // current is the coordinate of the point in the file currently to
  // be read. It changes as the matrix is read. mx_current is the
  // coordinate of the current data point but in the matrix. It
  // changes as the marix is read.
  std::vector<unsigned int> current = coord;
  std::vector<unsigned int> mx_current(mxdim,0);

  // Last dimension in the file which is also a dimension in mx.
  unsigned int lastdim = 0;
  for (unsigned int i=0; i<dim; i++)
    if (is_mxdim[i]) lastdim = i;
  
  // Read the matrix.
  bool done_flg = 0;
  while (!done_flg)
    {
      // Calculate the address of the receiving location in the matrix.
      float *mxadr = mx;
      for (unsigned int i=0; i<dim; i++)
	if (is_mxdim[i])
	  mxadr += mx_current[mxdimvec[i]]*mx_strides[mxdimvec[i]];

      // Read the data point.
      *(float*)mxadr = (float) read_pt_nooffset(current);
      
      // Step to next data point.
      bool reset_flg = 1;
      for (unsigned i=0; i<dim; i++)
	{
	  if (is_mxdim[i])
	    {	      
	      if (reset_flg)
		{
		  // Increment this dimensions coordinate.
		  current[i]++;
		  mx_current[mxdimvec[i]]++;
		  
		  // Have we gone beyond the end in this dimension.
		  if (current[i]>last_coord[i])
		    {
		      // Reset this dimensions coordinate.
		      current[i] = coord[i];
		      mx_current[mxdimvec[i]] = 0;
		      
		      // If this is lastdim, we are done.
		      if (i==lastdim) 
			{
			  done_flg = 1;
			  break;
			}		      
		    }
		  else break;
		}
	    }
	}
    }  
}


double DataIO::read_pt(std::vector<unsigned int> coord)
{
  for (int i=0; i<dim; i++)
    coord[i] -=pt_offset;
  return (read_pt_nooffset(coord));
}
