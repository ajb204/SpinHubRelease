/***/
/* testsize.h: definitions for returning file and file system sizes, accommodating 64-bit systems.
/***/

#include "prec.h"

#define DEF_BLOCK_SIZE         1024

/***/
/* Returned by functions to extract file or file system size:
/***/

#define FILE_STATUS_FAIL         -1 /* Valid info could not be extracted.          */
#define FILE_STATUS_INT_OK        0 /* Valid size is returned correctly as an int. */
#define FILE_STATUS_INT_TOO_SMALL 1 /* Valid size can't be expressed as an NMR_INT */

/***/
/* Returned by functions to compare actual file sizes with predicted size,
/* or predicted size with available space on filesystem:
/***/

#define FILE_SIZE_FAIL           -1 /* Valid info could not be extracted.          */
#define FILE_SIZE_OK              0 /* File size is acceptable.                    */
#define FILE_SIZE_TOO_SMALL       1 /* File is smaller than prediction.            */ 
#define FILE_SIZE_TOO_LARGE       2 /* File is larger than prediction.             */

#define NMR_SINGLE_FILE           0 /* Size pertains to a single data file or stream. */
#define NMR_FILE_SERIES           1 /* Size pertains to a multi-file data set.        */

/* iTotalBytes = blockSize*blockCount + extraBytes */

struct FileSize
   {
    int     status;
    double  dTotalBytes;
    NMR_INT iTotalBytes;
    NMR_INT blockSize;
    NMR_INT blockCount;
    NMR_INT extraBytes;
   };

double getFileBytes(), getFreeBytes(), getNMRBytes();
int    testNMRin(), testNMRout(), fileSizeInit();
char   *GetNMRIntStr();

#define getNMRIntStr( I ) GetNMRIntStr( (NMR_INT)((NMR_INT)I) )

