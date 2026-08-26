#include "dataIO.H"
#include "dataIO_nmrpipe.H"

int main()
{
  // Open nmrpipe file.
  dataIO *io = new dataIO_nmrpipe("/kl5/smk/tmp/spec.pip");  
}
