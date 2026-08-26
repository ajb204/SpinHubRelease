
#include "dataIO.H"
#include "dataIO_sparky.H"

// Return true if the machine is little endian.
static bool little_endian(void)
{
  int i = 1;

  return *(char *)&i;
}


static void swap_float_array(unsigned int n, float f[])
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


static float char2float(char* cbuf, unsigned int pos, bool swap_flg)
{
  /* Retrieve a float from pos in a char array and eventually byteswap. */

  union swap4 { float f; char s[4]; } out;
  
  if (swap_flg)
    {      
      out.s[0] = cbuf[pos+3];
      out.s[1] = cbuf[pos+2];
      out.s[2] = cbuf[pos+1];
      out.s[3] = cbuf[pos+0];
      return out.f;
    }
  else
    {
      out.s[0] = cbuf[pos+0];
      out.s[1] = cbuf[pos+1];
      out.s[2] = cbuf[pos+2];
      out.s[3] = cbuf[pos+3];
      return out.f;
    }  
}


static int char2int(char* cbuf, unsigned int pos, bool swap_flg)
{
  /* Retrieve an int from pos in a char array and eventually byteswap. */

  union swap4 { int i; char s[4]; } out;
  
  if (swap_flg)
    {      
      out.s[0] = cbuf[pos+3];
      out.s[1] = cbuf[pos+2];
      out.s[2] = cbuf[pos+1];
      out.s[3] = cbuf[pos+0];
      return out.i;
    }
  else
    {
      out.s[0] = cbuf[pos+0];
      out.s[1] = cbuf[pos+1];
      out.s[2] = cbuf[pos+2];
      out.s[3] = cbuf[pos+3];
      return out.i;
    }  
}


static unsigned int char2uint(char* cbuf, unsigned int pos, bool swap_flg)
{
  /* Retrieve an unsigned int from pos in a char array and eventually
     byteswap. */

  union swap4 { unsigned int u; char s[4]; } out;
  
  if (swap_flg)
    {      
      out.s[0] = cbuf[pos+3];
      out.s[1] = cbuf[pos+2];
      out.s[2] = cbuf[pos+1];
      out.s[3] = cbuf[pos+0];
      return out.u;
    }
  else
    {
      out.s[0] = cbuf[pos+0];
      out.s[1] = cbuf[pos+1];
      out.s[2] = cbuf[pos+2];
      out.s[3] = cbuf[pos+3];
      return out.u;
    }  
}


/* Some tiny helper classes. */

DataIO_sparky_tile::DataIO_sparky_tile(FILE *stream,
				       unsigned int tile_number,
				       unsigned int bsize,
				       unsigned int bskip) :
  tnumber(tile_number),
  buf(0)
{
  // Allocate the buffer.
  buf = new float[bsize/4];

  // Position the file.
  if (fseek(stream, bskip, SEEK_SET)!=0)
    throw DataIO::FseekError(bskip);
      
  // Read the tile.
  unsigned int rtnval = fread(buf, 1, bsize, stream);
  if (rtnval!=bsize) 
    throw DataIO::FreadError(bskip,bsize,rtnval);
}


DataIO_sparky_tile::~DataIO_sparky_tile()
{
  if (buf != 0) delete[] buf;
}


unsigned int DataIO_sparky_tile::get_tnumber()
{
  return tnumber;
}


float *DataIO_sparky_tile::get_buf()
{
  return buf;
}



/* ************************************************************** */
/* Here goes the class and object methods.                        */


// Creator. Opens a file and reads the parameter header.
DataIO_sparky::DataIO_sparky(std::string fname) : 
  DataIO(fname),
  open_flg(0),
  debug(0),
  swap_flg(0)
{
  // Open the file.
  stream = fopen(fname.c_str(),"r");
  if (stream==0) throw DataIO::CannotOpenFile(fname);
  open_flg = 1;
  
  /* Set swap flag if the machine is little endian as sparky files are
     always saved as big endian */
  swap_flg = little_endian();

  // Read the first header.
  
  // Position the file.
  if (fseek(stream, 0, SEEK_SET)!=0)
    throw DataIO::FseekError(0);
      
  // Read the header.
  unsigned int rtnval = fread(head, 1, SPARKY_HEAD_SIZE, stream);
  if (rtnval!=SPARKY_HEAD_SIZE) 
    throw DataIO::FreadError(0,SPARKY_HEAD_SIZE,rtnval);

  // Get number of dimensions.
  int dim = head[10];

  // Check that we have a valid sparky file.
  if (head[11]!=1 ||head[13]!=2 || dim<1 || dim>4)
    throw DataIO::HeaderInvalid();

  // Set number of dimensions.
  set_dim(dim);

  // Read dimheaders one at a time.
  dimhead.resize(dim);
  for(int i=0;i<dim;i++)
    {
      dimhead[i] = new char[SPARKY_DIMHEAD_SIZE];
      unsigned int rtnval = fread(dimhead[i], 1, SPARKY_DIMHEAD_SIZE, stream);
      if (rtnval!=SPARKY_DIMHEAD_SIZE) 
	throw DataIO::FreadError(0,SPARKY_HEAD_SIZE,rtnval);
    }
  
  
  // Get matrix size.
  for(int i=0;i<dim;i++)
    set_size(i,char2int(dimhead[dim-i-1],8,swap_flg));
  
  // Initialize map with values from 0.0 to 1.0.
  for (int idim=0; idim<dim; idim++)
    map_linear(idim, 0, 0.0, get_size(idim), 1.0);

  // Setup tile data structures.

  // The tile dimensions.
  tilesize.resize(dim);
  for(int i=0;i<dim;i++)
    tilesize[i] = char2int(dimhead[dim-i-1],16,swap_flg);
  
  // Number of tiles in each dimension.
  tilecount.resize(dim);
  for (int idim=0; idim<dim; idim++)
    tilecount[idim] = get_size(idim)/tilesize[idim];
  
  // Number of tiles to skip in each dimension.
  tileskip.resize(dim);
  tileskip[0] = 1;
  for (int idim=1; idim<dim; idim++)
    tileskip[idim] = tileskip[idim-1]*tilecount[idim-1];
  
  // Number of positions to skip in a tile in each dimension.
  tilebufskip.resize(dim);
  tilebufskip[0] = 1;
  for (int idim=1; idim<dim; idim++)
    tilebufskip[idim] = tilebufskip[idim-1]*tilesize[idim-1];
  
  // The size of a tile in bytes.
  tile_bsize = 4;
  for (int idim=0; idim<dim; idim++)
    tile_bsize *= tilesize[idim];

  // This is the number of bytes to skip in the beginning of the file.
  tile_initskip_bsize = SPARKY_HEAD_SIZE + SPARKY_DIMHEAD_SIZE*dim;
  
  // We don't have a current tile.
  current_tile = 0;

  // Adjust the tile list to contain enough tiles to span any dimension.
  max_tilelist_size = 0;
  for(int i=0;i<dim;i++) 
    if (tilecount[i]>max_tilelist_size)
      max_tilelist_size = tilecount[i];

  // Get sfreq, swidth, and ppmcenter (ppm value of point size/2).
  sfreq.resize(dim);
  swidth.resize(dim);
  ppmcenter.resize(dim);
  for(int i=0;i<dim;i++)
    {    
      sfreq[i] = char2float(dimhead[dim-i-1],20,swap_flg);
      swidth[i] = char2float(dimhead[dim-i-1],24,swap_flg);
      ppmcenter[i] = char2float(dimhead[dim-i-1],28,swap_flg);
    }
  
}


// Destructor. Closes a file and deallocates workspace.
DataIO_sparky::~DataIO_sparky()
{
  // Close the file.
  if (open_flg) fclose(stream);

  // Deallocate the current tile.
  delete current_tile;

  // Deallocate the tiles in the tilelist.
  while (tilelist.size()>0)
    {
      DataIO_sparky_tile *tile = tilelist.back();
      tilelist.pop_back();
      delete tile;
    }
  
}


// Set debug flag. 
void DataIO_sparky::set_debug(int dbg)
{
  debug = dbg;
}



unsigned int DataIO_sparky::pt2tnumber(std::vector<unsigned int>& pt)
{
  unsigned int tnumber = 0;
  for(int i=0; i<get_dim(); i++)
    {
      unsigned int tindex = pt[i]/tilesize[i];
      tnumber += tindex*tileskip[i];
    }
  return tnumber;
}

unsigned int DataIO_sparky::pt2tpos(std::vector<unsigned int>& pt)
{
  unsigned int tpos = 0;
  for(int i=0; i<get_dim(); i++)
    {
      unsigned int tindex = pt[i]/tilesize[i];
      unsigned int tcoord = pt[i]-tindex*tilesize[i];
      tpos += tcoord*tilebufskip[i];
    }
  return tpos;
}


unsigned int DataIO_sparky::tnumber2byteskip(unsigned int tnumber)
{
  return (tile_initskip_bsize + tnumber*tile_bsize);
}


float DataIO_sparky::read_pt_nooffset(std::vector<unsigned int> coord)
{
  // Get dim.
  int dim = get_dim();

  // Is the vector valid.
  if (dim!=coord.size())
    throw DataIO::PointDimInvalid(coord.size());

  // Is the requested point within bounds.
  for (unsigned int i=0; i<get_dim(); i++)
    if (coord[i]<0 || coord[i]>=get_size(i))
      throw DataIO::PointIndexInvalid(i,coord[i]);

  // Calculate the tile number and the position in the tile.
  unsigned int tnumber = pt2tnumber(coord);
  unsigned int tpos = pt2tpos(coord);
    
  // First, look in the tile list.
  if (current_tile!=0 && current_tile->get_tnumber()!=tnumber)
    {
      // First we push the current tile on the list.
      tilelist.push_front(current_tile);
      
      // Zero the current tile.
      current_tile = 0;

      // Look in the list for tnumber and if found, set current tile.
      for(std::list<DataIO_sparky_tile*>::iterator ti=tilelist.begin();
	  ti!=tilelist.end();ti++)
	if ((*ti)->get_tnumber()==tnumber)
	  {
	    // Set the current tile and remove it from the tile list.
	    current_tile = *ti;
	    tilelist.erase(ti);
	    break;
	  }

      // If capacity is exceeded, we remove tiles from the end (the oldest).
      while (tilelist.size()>max_tilelist_size)
	{
	  DataIO_sparky_tile *tile = tilelist.back();
	  tilelist.pop_back();
	  delete tile;
	}
      
    }
      
  // If no current tile, we read the tile.
  if (current_tile==0)
    {
      // Get number of bytes to skip before reading tile.
      unsigned int bskip = tnumber2byteskip(tnumber);

      // Allocate and read new tile.
      current_tile = new DataIO_sparky_tile(stream, tnumber,
					    tile_bsize, bskip);

      /* Eventually byteswap data */
      if (swap_flg) swap_float_array(tile_bsize/4, current_tile->get_buf());
    }
  
  // Get the value.
  float fval = (current_tile->get_buf())[tpos];

  if (debug>0)
    {    
      std::cout 
	<< "("
	<< coord[0];

      for(int i=1;i<dim;i++)
	std::cout
	  << ","
	  << coord[i];
      
      std::cout
	<< ")  tnumber: "
	<< tnumber
	<< "   tpos: "
	<< tpos
	<< "   value: "
	<< tpos
	<< "   value: "
	<< fval
	<< "\n";
    }
  
  return fval;
}



void DataIO_sparky::print()
{
  print_base();
  std::cout << "type:      " << get_type() << "\n";
  std::cout << "open_flg:  " << is_open() << "\n";
  std::cout << "swap_flg:  " << is_swap() << "\n";
  for (unsigned int i=0; i<get_dim(); i++)
    std::cout 
      << "tilesize["
      << i 
      << "]: "
      << tilesize[i]
      << "  "
      << "tilebufskip["
      << i 
      << "]: "
      << tilebufskip[i]
      << "  "
      << "tilecount["
      << i 
      << "]: "
      << tilecount[i]
      << "  "
      << "tileskip["
      << i 
      << "]: "
      << tileskip[i]
      << "  "
      << "\n";
  for (unsigned int i=0; i<get_dim(); i++)
    std::cout 
      << "sfreq["
      << i 
      << "]: "
      << sfreq[i]
      << "  "
      << "swidth["
      << i 
      << "]: "
      << swidth[i]
      << "  "
      << "ppmcenter["
      << i 
      << "]: "
      << ppmcenter[i]
      << "\n";

  std::cout << "tile_bsize:    " << tile_bsize << "\n";
  std::cout << "tile_initskip_bsize:    " << tile_initskip_bsize << "\n";
  if (current_tile!=0)
    std::cout << "current tile:    " << current_tile->get_tnumber() << "\n";
}



std::string DataIO_sparky::get_type()
{
  return ("sparky");
}


bool DataIO_sparky::is_swap()
{
  return (swap_flg);
}


bool DataIO_sparky::is_open()
{
  return (open_flg);
}


void DataIO_sparky::map_default(unsigned int idim, std::string& unit)  
{
  map_linear(idim, get_size(idim)/2, ppmcenter[idim], 
	     0, ppmcenter[idim]+0.5*swidth[idim]/sfreq[idim]);
}
