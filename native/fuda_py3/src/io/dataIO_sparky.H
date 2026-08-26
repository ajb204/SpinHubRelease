#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <cstdio>
#include <cmath>
#include "dataIO.H"

#ifndef _DATAIO_SPARKY_H
#define _DATAIO_SPARKY_H

// Byte size of initial sparky file header and of dimension headers. 
#define SPARKY_HEAD_SIZE ((size_t) 180)
#define SPARKY_DIMHEAD_SIZE ((size_t) 128)

// Predeclarations.
class DataIO_sparky_tile;


/************************************************************************/
/* sparky is a class for reading and writing sparky data files. */
class DataIO_sparky : public DataIO
{
  // Object members.
  FILE *stream; // This is the file stream.
  bool open_flg; // True when file is open.
  bool swap_flg; // Swap_flg is set if file is byte swapped.
  int debug; // Debug flag.
  char head[SPARKY_HEAD_SIZE]; // sparky file header.
  std::vector<char*> dimhead; // sparky dim headers.

  // Tile (submatrix) manipulation data structures.
  std::vector<unsigned int> tilesize; // dimensions of a tile.
  std::vector<unsigned int> tilebufskip; /* number of positions to
                                            skip in each dimension for
                                            a tile */
  std::vector<unsigned int> tilecount; // number of tiles in each dimension.
  std::vector<unsigned int> tileskip; /* number of tiles to skip in each
				         dimension. */
  unsigned int tile_bsize; // The size of a tile in bytes.
  unsigned int tile_initskip_bsize; /* the size of all headers in
				       bytes which precedes the
				       submatrix structure */
  // Tile buffering.
  DataIO_sparky_tile *current_tile; /* Ref. to the tile most recently
                                       addressed. Is identical to the
                                       first element in tilelist. */
  std::list<DataIO_sparky_tile*> tilelist; // the tile list buffer.
  unsigned int max_tilelist_size; // Maximum number of tiles in tilelist.
  
  std::vector<float> sfreq;
  std::vector<float> swidth;
  std::vector<float> ppmcenter;

  // Tile (private) methods.
  /* Calculate tile number (tiles are numbered from 0 and onwards) */
  unsigned int pt2tnumber(std::vector<unsigned int>& pt);
  /* Calculate position of point in tile. */
  unsigned int pt2tpos(std::vector<unsigned int>& pt);
  /* Calculate byte position for tile in file. */
  unsigned int tnumber2byteskip(unsigned int tnumber);
  
  
public:

  // Constructor. Opens an nmrpipe file.
  DataIO_sparky(std::string fname);

  // Destructor closes file and deallocates tile stack.
  ~DataIO_sparky();
    
  // Object methods.
  float read_pt_nooffset(std::vector<unsigned int> coord);
  void print();
  void set_debug(int dbg);
  std::string get_type();
  bool is_swap();
  bool is_open();
  void map_default(unsigned int dim, std::string& unit);  
};


class DataIO_sparky_tile
{
  unsigned int tnumber;
  float *buf;

public:

  DataIO_sparky_tile(FILE *stream,
		     unsigned int tile_number,
		     unsigned int byte_size,
		     unsigned int byte_skip);
  
  ~DataIO_sparky_tile();
  
  unsigned int get_tnumber();
  float *get_buf();
};





#endif
