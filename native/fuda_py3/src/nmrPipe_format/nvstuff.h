
#ifndef _PARAMH
#define _PARAMH

#ifndef TRUE
#define TRUE        1
#endif

#ifndef FALSE
#define FALSE       0
#endif
/*
#ifndef Boolean
#define Boolean     int
#endif
*/
#ifndef NIL
#define NIL         0
#endif

#ifndef NULL
#define NULL        0
#endif

#ifndef X
#define X 0
#endif

#ifndef Y
#define Y 1
#endif

#ifndef BLKREAD 
#define  BLKREAD 0
#endif

#ifndef BLKWRITE 
#define  BLKWRITE 1
#endif

#endif


#ifndef _NMRCONSTANTSH
#define _NMRCONSTANTSH

#ifndef NMR_NDIM
#define NMR_NDIM 	4
#endif

/* was 131072 4/2017 */

#ifndef NMR_NBLKS
#define NMR_NBLKS 16777216
#endif

#endif


#ifndef _NMRFILEH
#define _NMRFILEH

typedef struct {
	int      ndim;
	float   *blk[NMR_NBLKS];
	int      size[NMR_NDIM];
	int      vsize[NMR_NDIM];
	int      blksize[NMR_NDIM];
	int      nblks[NMR_NDIM];
	int      nvblks[NMR_NDIM];
	int	 offblk[NMR_NDIM];
	int	 offvblk[NMR_NDIM];
	int	 blkmask[NMR_NDIM];
	int	 offpt[NMR_NDIM];
	float	 sf[NMR_NDIM];
	float	 sw[NMR_NDIM];
	float	 refpt[NMR_NDIM];
	float	 ref[NMR_NDIM];
	float	 foldup[NMR_NDIM];
	float	 folddown[NMR_NDIM];
	int	 refunits[NMR_NDIM];
	int      blkelems;
        int      fheadersz;
        int      bheadersz;
	int	 open;
	char	 label[NMR_NDIM][16];
} MATRIXFILE; 

#endif

#ifndef _DATASETSH
#define _DATASETSH

typedef struct {
	unsigned long id;
	char          filename[128]; 
	char          dirname[128];
	MATRIXFILE    matfile;
	FILE         *fileptr;
	int           fi;
} DATASET; 
#endif

#ifndef _ERRORSH
#define _ERRORSH

#ifndef ERROR
#define ERROR	    1
#define NOERROR	    0
#endif

DATASET *Make_Dataset();

int DatasetDumpPar();
int DimDataset();
int HeaderWriteNV();

int SetBlockSize();
int WriteEmptyFile();
int nvCloseFile();
int nvSetup();
int nvWrite();
int nvWriteBuffers();
int open_file();
int zerobuff();

float *RWBlk();
float *getbuf();

#endif
