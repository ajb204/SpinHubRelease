
#include "dataIO.H"
#include "dataIO_nmrpipe.H"

// NMRPipe definitions taken from fdatap.h, dimloc.h, specunit.h.
#define FDATASIZE          512   /* Length of header in 4-byte float values. */
#define MAXDIM      4

#define HDR_OK      0
#define HDR_SWAPPED 1

#define FDDIMCOUNT  9
#define FD_COMPLEX  0

#define NDPARM        1000
#define NDSIZE        (1+NDPARM)  /* Number of points in dimension.          */
#define NDQUADFLAG    (7+NDPARM)  /* Data Type Code (See Below).             */

#define LAB_PTS     "Pts"
#define LAB_HZ      "Hz"
#define LAB_PPM     "ppm"


// Get a parameter from the header.
extern "C" float getParm(const float header[],
		       int paramCode, int origDimCode);

// Test validity of header.
extern "C" int testHdr(const float header[]);

// Byteswap the header.
extern "C" int swapHdr(const float header[]);

// Convert from points to user units (e.g. ppm).
extern "C" float iPnt2spec(const float header[], int dimCode, 
			   int pntVal, const char *specLabel);


/* ************************************************************** */
/* Here goes the class and object methods.                        */


// Creator. Opens a file and reads the parameter header.
DataIO_nmrpipe::DataIO_nmrpipe(std::string fname) : 
  DataIO(fname),
  open_flg(0),
  swap_flg(0),
  debug(0),
  header(0),
  wsize(0),
  wspace(0)
{
  // Open the file.
  stream = fopen(fname.c_str(),"r");
  if (stream==0) throw DataIO::CannotOpenFile(fname);
  open_flg = 1;
  
  // Allocate and read the header.
  header = new float[FDATASIZE];
  size_t count = fread(header, sizeof(float), FDATASIZE, stream);
  std::cout << count << ' ' << FDATASIZE << std::endl;
  if (count!=FDATASIZE) throw DataIO::CannotReadHeader();
  
  // Test the data header and eventually byte swap the header */
  int return_stat = testHdr(header);
  if (return_stat==HDR_OK)
    {
      // The file is not byte swapped.
      swap_flg = 0;
    }
  else if (return_stat==HDR_SWAPPED)
    {
      // The header is swapped.
      swap_flg = 1;
    }
  else
    throw DataIO::HeaderInvalid();
  
  // Get and set number of dimensions.
  int dim = (int) getParm(header, FDDIMCOUNT, 0);
  set_dim(dim);

  // Check that we have a real matrix.
  for (unsigned int idim=1; idim<=dim; idim++)
    if (FD_COMPLEX== (unsigned int) getParm(header, NDQUADFLAG, idim))
      throw DataIO::DimIsComplex(idim);
  
  // Get matrix size.
  for (unsigned int idim=1; idim<=dim; idim++)
    set_size(idim-1, (unsigned int) getParm(header,NDSIZE,idim));

  /* Setup bsize.
     bsize[0] holds the size of one float in bytes.
     bsize[1] holds the size of one row in bytes.
     bsize[2] holds the size of one plane in bytes.
     bsize[4] holds the size of one cube in bytes. */
  bsize.resize(MAXDIM);
  bsize[0] = sizeof(float);
  for(int i=1; i<dim; i++) bsize[i] = bsize[i-1]*get_size(i-1);

  // These are actually not used and therefore quite silly.
  for(int i=dim; i<MAXDIM; i++) bsize[i] = bsize[dim-1]*get_size(dim-1);

  // Default workspace allocation.
  wsize = get_size(0); // Size of one row.
  wspace = new float[wsize];

  // Allocate and initialize wspace contents record.
  wcoord.resize(MAXDIM);
  for(int i=0; i<MAXDIM; i++) wcoord[i] = 1;

  // Make sure we don't think we got a vald row in wspace already.
  wdefined = 0;

  // Initialize map with values from 0.0 to 1.0.
  for (int idim=0; idim<dim; idim++)
    map_linear(idim, 0, 0.0, get_size(idim), 1.0);
}


// Destructor. Closes a file and deallocates workspace.
DataIO_nmrpipe::~DataIO_nmrpipe()
{
  // Close the file.
  if (open_flg) fclose(stream);

  // Deallocate work space.
  if (wspace!=0) delete[] wspace;
}


// Set debug flag.
void DataIO_nmrpipe::set_debug(int dbg)
{
  debug = dbg;
}




void swap_float(unsigned int n, float f[])
{
  /* Byteswap n 4 byte floats */

  int i;
  union swap4 { float f; char s[4]; } in, out;
  
  for(i=0;i<n;i++)
    {
      in.f = f[i];
      out.s[0] = in.s[3];
      out.s[1] = in.s[2];
      out.s[2] = in.s[1];
      out.s[3] = in.s[0];
      f[i] = out.f;
    }
}



float DataIO_nmrpipe::read_pt_nooffset(std::vector<unsigned int> coord)
{
  // Is the vector valid.
  if (get_dim()!=coord.size())
    throw DataIO::PointDimInvalid(coord.size());

  // Is the requested point within bounds.
  for (unsigned int i=0; i<get_dim(); i++)
    if (coord[i]<0 || coord[i]>=get_size(i))
      throw DataIO::PointIndexInvalid(i,coord[i]);

  // Do we have the data row in the buffer?
  bool inside_flg = 1;
  if (wdefined)
    {
      for (unsigned int i=1; i<get_dim(); i++)
	if (coord[i]!=wcoord[i]) 
	  {
	    inside_flg = 0;
	    break;
	  }
    }
  else inside_flg = 0;
  
  // If data is not in the buffer, we read in a row.
  if (!inside_flg)
    {
      /* Determine the absolute byte offset in the file for
	 reading a full row */
      unsigned int offset = FDATASIZE*sizeof(float);
      for(unsigned int i=1; i<get_dim(); i++)
	offset += coord[i]*bsize[i];
      
      /* Position the file */
      if (fseek(stream, offset, SEEK_SET)!=0)
	throw DataIO::FseekError(offset);
      
      /* Read the row */
      unsigned int rtnval = fread(wspace, sizeof(float), get_size(0), stream);
      if (rtnval!=get_size(0)) 
	throw DataIO::FreadError(offset,get_size(0),rtnval);
      
      /* Eventually byteswap data */
      if (swap_flg) swap_float(get_size(0), wspace);      

      /* Set wcoord */
      wcoord[0] = 0;
      for (int i=1; i<coord.size(); i++)
	wcoord[i] = coord[i];

      // Set buffer contents to be defined.
      wdefined = 1;
    }

  // Get the point.
  float fval = wspace[coord[0]];

  if (debug>0)
    {
      /* Determine the absolute byte offset in the file for
	 reading a full row */
      unsigned int offset = FDATASIZE*sizeof(float);
      for(unsigned int i=1; i<get_dim(); i++)
	offset += coord[i]*bsize[i];
      
      std::cout 
	<< "("
	<< coord[0];

      for(int i=1;i<get_dim();i++)
	std::cout
	  << ","
	  << coord[i];
      
      std::cout 
	<< ")   wcoord("
	<< wcoord[0];

      for(int i=1;i<get_dim();i++)
	std::cout
	  << ","
	  << wcoord[i];
      
      std::cout
	<< ")  offset: "
	<< offset
	<< "  fval: "
	<< fval
	<< "\n";
    }

  // Return the wanted point.
  return fval;
}



void DataIO_nmrpipe::print()
{
  print_base();
  std::cout << "type:      " << get_type() << "\n";
  std::cout << "swap_flg:  " << is_swap() << "\n";
  for (unsigned int i=0; i<get_dim(); i++)
    std::cout << "bsize[" << i << "]: " << bsize[i] << "\n";

  std::cout << "wsize:    " << wsize << "\n";
  for (unsigned int i=0; i<get_dim(); i++)
    std::cout << "wcoord[" << i << "]: " << wcoord[i] << "\n";
}



std::string DataIO_nmrpipe::get_type()
{
  return ("nmrpipe");
}


bool DataIO_nmrpipe::is_swap()
{
  return (swap_flg);
}


bool DataIO_nmrpipe::is_open()
{
  return (open_flg);
}


void DataIO_nmrpipe::map_default(unsigned int idim, std::string& unit)  
{
  // Check idim.
  if (idim-get_dim_offset()<0 || idim-get_dim_offset()>=get_dim())
    throw DataIO::DimIndexInvalid(idim);
 
  // Dimensions in nmrpipe are numbered 1,2,3 and 4.
  int idim_nmrpipe = idim-get_dim_offset()+1;

  if ((unit=="ppm") or (unit=="PPM"))
    {
      // Points in nmrpipe are numbered 1,2,3,4,5....,N.

      // First point.
      int pt0_nmrpipe = 1;
      int pt0 = get_pt_offset();
      double val0 = (double) iPnt2spec( header, idim_nmrpipe, 
					pt0_nmrpipe, LAB_PPM );

      // Last point.
      int pt1_nmrpipe = get_size(idim);
      int pt1 = get_pt_offset()+get_size(idim)-1;
      double val1 = (double) iPnt2spec( header, idim_nmrpipe,
					pt1_nmrpipe, LAB_PPM );

      // Set up a linear mapping.
      map_linear(idim, pt0, val0, pt1, val1);
    }
  else 
    throw DataIO::InvalidUnitString(unit);
}
