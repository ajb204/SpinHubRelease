#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <cstdio>
#include <cmath>
#include "dataIO.H"

#ifndef _DATAIO_NMRPIPE_H
#define _DATAIO_NMRPIPE_H

/************************************************************************/
/* NMRPipe is a class for reading and writing NMRPipe data
   files. Currently, only real matrixes are supported and only read
   operations are implemented. */
class DataIO_nmrpipe : public DataIO
{
  // Object members.
  FILE *stream; // This is the file stream.
  bool open_flg; // True when file is open.
  bool swap_flg; // Swap_flg is set if file is byte swapped.
  int debug; // Debug flag.
  float *header; // NMRPipe file header pointer.
  std::vector<unsigned int> bsize; /* Accumulated sizes in
				      bytes. See constructor*/
  bool wdefined; // is set when wspace has well-defined contents.
  unsigned int wsize; // size of IO work space.
  std::vector<unsigned int> wcoord; // Coordinate of row populating wspace. 
  float *wspace; // Work space
  
public:

  // Constructor. Opens an nmrpipe file.
  DataIO_nmrpipe(std::string fname);

  // Destructor closes file and deallocates work buffer.
  ~DataIO_nmrpipe();
    
  // Object methods.
  float read_pt_nooffset(std::vector<unsigned int> coord);
  void print();
  void set_debug(int dbg);
  std::string get_type();
  bool is_swap();
  bool is_open();
  void map_default(unsigned int dim, std::string& unit);  
};

#endif
