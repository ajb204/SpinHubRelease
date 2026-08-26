
#define SMX_MAXDIM       8

#define SMX_ERR_MAXDIM   1
#define SMX_ERR_NULLSIZE 2
#define SMX_ERR_MATLIM   3
#define SMX_ERR_SMXLIM   4
#define SMX_ERR_EDGE     5
#define SMX_ERR_MISMATCH 6

int smxInit();
int smx2matrix();
int matrix2smx();
int getSMXLoc();
int getMatrixLoc();
