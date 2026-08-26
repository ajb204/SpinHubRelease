/*
COPYRIGHT (C) 1994 MERCK AND CO., INC.
Whitehouse Station, N.J. U.S.A
All rights reserved.

IN NO EVENT SHALL MERCK AND COMPANY BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF MERCK AND CO.,INC.
HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

MERCK AND CO. SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
ON AN "AS IS" BASIS, AND MERCK AND CO. HAS NO OBLIGATION TO
PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prec.h"
#include "nvstuff.h"
#include "nvheader.h"
#include "nv_mf.h"

#define FPR (void)fprintf

#ifdef SUN
#define SIZE_T unsigned
#define VOID   char 
#else
#define SIZE_T size_t
#define VOID   void 
#endif

float          *sb, *bb;
int             ib, jb;
float          *RWBlk();

/* was 1024 4/2017 */

#define NMR_NBUFFER 4096

struct Buffer {
    float          *bufpt;
    DATASET        *data;
    int             fileblk;
    int             size;
}               buffer[NMR_NBUFFER];
FILEHEADER      fhead;

/****************************************************************************/

float *getbuff(size, data, blk)
    int             size, blk;
    DATASET        *data;
{
    static int      n = 0;
    float          *buff;

    /*
     * Is the datafile too big. 
     */
#ifndef S_SPLINT_S
    if (blk >= NMR_NBLKS) {
	FPR(stderr, "File too large !!! \n");
	exit(1);
    }
    if (blk < 0) {
	FPR(stderr, "Invalid Block # %d !!! \n", blk);
	exit(1);
    }
#endif

    if (n == NMR_NBUFFER)
	n = 0;

    if (buffer[n].bufpt != NULL) {
	(void)RWBlk(buffer[n].data, buffer[n].fileblk, BLKWRITE);
	if (size != buffer[n].size) {
	    FPR(stderr, "Freeing buffer %d \n", n);
	    free((VOID*)buffer[n].bufpt);
	    FPR(stderr, "Allocating buffer %d of size %d\n", n, size);
	    buffer[n].bufpt = (float *) calloc((SIZE_T)size, sizeof(float));
	}
	(buffer[n].data)->matfile.blk[buffer[n].fileblk] = NULL;
    } else {
	/* FPR(stderr,"Allocating buffer %d of size %d\n",n,size); */
	buffer[n].bufpt = (float *) calloc((SIZE_T)size, sizeof(float));
    }
    buffer[n].data = data;
    buffer[n].fileblk = blk;
    buffer[n].size = size;

    if (buffer[n].bufpt == NULL)
	FPR(stderr, "Error allocating memory for buffer\n");
    n++;
    return (buffer[n - 1].bufpt);
}

/****************************************************************************/

int zerobuff()
{
    int             i;
    for (i = 0; i < NMR_NBUFFER; i++) {
	buffer[i].bufpt = NULL;
	buffer[i].data = NULL;
	buffer[i].fileblk = -1;
    }

    return(0);
}

int nvCloseFile(data)
    DATASET        *data;
{
    (void) fclose(data->fileptr);
    return(0);
}

int nvWriteBuffers(data)
    DATASET        *data;
{
    int i;

/* delaglio */

    if (data->matfile.ndim == 1) return( 0 );

    for (i = 0; i < NMR_NBUFFER; i++) {
	if ((buffer[i].data == data) && (buffer[i].bufpt != NULL)) {
	    if (buffer[i].fileblk >= 0)
		(void)RWBlk(data, buffer[i].fileblk, BLKWRITE);
	}
    }
    return(0);
}



/********************************************************************************/

DATASET *Make_Dataset(dirname, name)
    char           *dirname, *name;
/*
 * Create data structure for a new NMRView dataset. 
 */
{
    int             i;
    DATASET        *data;
    /*
     * Allocate memory for the new dataset. 
     */
    if (!(data = (DATASET *) malloc(sizeof(DATASET))))
	return (NULL);

    /*
     * Set all variables in the structure. 
     */
#ifndef S_SPLINT_S
    data->id = (unsigned long) data;
#endif
    (void) strcpy(data->filename, name);
    (void) strcpy(data->dirname, dirname);
    data->matfile.open = 0;
    data->fileptr = NULL;
    data->fi = 0;
    for (i = 0; i < NMR_NDIM; i++) {
	data->matfile.refunits[i] = 3;
	data->matfile.sw[i] = 7000.0;
	data->matfile.sf[i] = 600.0;
	data->matfile.refpt[i] = -99000;
	data->matfile.ref[i] = 4.73;
	data->matfile.label[i][0] = '\0';
    }

    return (data);
}

int SetBlockSize(data)
    DATASET        *data;
{
    NMR_INT         npoints;
    int             i, ptspblk, blksize, blkspdim, blog;
    int             size[4], blksiz[4], nbdim[4];
    int             nblks, ndim, j, bmax, imax, bmin, imin,iSize2;
    npoints = 1;
    ndim = data->matfile.ndim;

    for (i = 0; i < ndim; i++) {
	size[i] = data->matfile.size[i];
        iSize2 = 1;
	for (j=0;j<20;j++) {
		if (size[i]==iSize2) {
			break;
		}
		if (size[i]<iSize2) {
			size[i]=iSize2;
			break;
		}
		iSize2 *= 2;
	}
#ifndef S_SPLINT_S
	if (j>=20) {
		FPR(stderr, "Dimension size too large !!! \n");
		exit(1);
	}
#endif

	/* FPR(stderr, "%d\n", size[i]); */
	npoints *= size[i];
    }
    nblks = npoints / 4096;
 /*   FPR(stderr, "%d %ld\n", nblks, (long)npoints);*/
    blkspdim = exp((1.0 / ndim) * log((double) nblks));
 /*   FPR(stderr, "blkspdim %d\n", blkspdim);*/
#ifdef S_SPLINT_S
    blog = 1.0 + 0.5;
#else
    blog = log((double) blkspdim) / log(2.0) + 0.5;
#endif
    blkspdim = exp(blog * log(2.0)) + 0.5;
    if (blkspdim < 1) blkspdim = 1; /* fd 3/9/2020 */
/*    FPR(stderr, "blkspdim %d\n", blkspdim);*/
    for (i = 0; i < ndim; i++) {
	nbdim[i] = blkspdim;
    }

    for (j = 0; j < 2; j++) {
	blksize = 1;
	bmax = 0;
	bmin = 1e6;
	for (i = 0; i < ndim; i++) {
	    /*
	     * FPR(stderr, "blksize %d %d\n", size[i] / nbdim[i], size[i]
	     * % nbdim[i]); 
	     */
	    blksize *= size[i] / nbdim[i];
	    if (nbdim[i] > bmax) {
		bmax = nbdim[i];
		imax = i;
	    }
	    if (nbdim[i] < bmin) {
		bmin = nbdim[i];
		imin = i;
	    }
	}
	/* FPR(stderr, "%d\n", blksize); */
	if (blksize > 4096)
	    nbdim[imax] = nbdim[imax] * 2;
	if (blksize < 4096)
	    nbdim[imax] = nbdim[imax] / 2;
	if (nbdim[imin] < 2) {
	    nbdim[imin] = 2;
	    nbdim[imax] = nbdim[imax] / 2;
	}
    }
    for (i = 0; i < ndim; i++) {
	data->matfile.blksize[i] = size[i] / nbdim[i];
	data->matfile.vsize[i] = size[i];
    }
    return(0);
}

/****************************************************************************/

float *RWBlk(data, i, mode)
    DATASET        *data;
    int             i, mode;
{
    int             nblks,vblk,blk, buffn,iBlock[NMR_NDIM],dim;
    long            offset;
    float          *buff;

    /*
     * Calculate which block to read. 
     */
    vblk = i;
    /* FPR(stderr,"%d mode for block %d \n",mode,i); */
    if (vblk < 0) {
	FPR(stderr, "Invalid block %d\n", vblk);
	return ((float *) NULL);
    }
    /*
     * Calculate offset of block in file. 
     */

    nblks = vblk;
    for (dim=(data->matfile.ndim-1);dim>=0;dim--) {
	iBlock[dim] = nblks /data->matfile.offvblk[dim];
	nblks  = nblks % data->matfile.offvblk[dim];
    }

    blk=0;
    for (dim=0;dim<data->matfile.ndim;dim++) {
	blk += iBlock[dim]*data->matfile.offblk[dim];
          /*FPR(stderr," iBlock %d offblk %d ",iBlock[dim],data->matfile.offblk[dim]); */
        /* FPR(stderr," %d ",iBlock[dim]); */
    }

/*     FPR(stderr," vblk %d blk %d %d \n",vblk,blk,mode); */
     /* FPR(stderr," blk %d \n",blk); */

    offset = (long)data->matfile.fheadersz + ((long)(blk + 1))*data->matfile.bheadersz + ((long)blk)*data->matfile.blkelems*4;

    /*
     * Which buffer. 
     */
    if (mode == BLKREAD)
	buff = getbuff(data->matfile.blkelems, data, i);
    else
	buff = data->matfile.blk[i];

    /*
     * Move file pointer to offset. 
     */
    if (fseek(data->fileptr, offset, 0))
       {
	FPR( stderr,
                 "Error Seeking block %d offset %ld blkelems %d\n",
                 blk, offset, data->matfile.blkelems );

	return ((float *) NULL);
       }


    /*
     * Read data into buffer 
     */
    if (mode == BLKREAD) {
	if (fread((void *)buff, (size_t)4, (size_t)data->matfile.blkelems, data->fileptr) != (size_t)data->matfile.blkelems) {
	    FPR(stderr, "Error reading data from block %d\n", blk);
	    return ((float *) NULL);
	}
	/*
	 * Set ptr to buffer 
	 */
	data->matfile.blk[i] = buff;
	/*
	 * Write data from buffer 
	 */
    } else {
	if (fwrite((void *)buff, (size_t)4, (size_t)data->matfile.blkelems, data->fileptr) != (size_t)data->matfile.blkelems) {
	    FPR(stderr, "Error writing data to block %d\n", blk);
	}
    }

    /*
     * Return ptr to buffer. 
     */
#ifdef S_SPLINT_S
    return( (float *)NULL );
#else
    return (data->matfile.blk[i]);
#endif
}

int nvSetup(data)
    DATASET        *data;
{
    int             filedim;
    filedim = data->matfile.ndim;
    if (filedim == 2) {
	M_2Setup(data);
    }
    if (filedim == 3) {
	M_3Setup(data);
    }
    if (filedim == 4) {
	M_4Setup(data);
    }
    return(0);
}

int nvWrite(data, mat, ix, iy, iz, iz2, islices)
    DATASET        *data;
    float          *mat;
    int             ix, iy, iz, iz2, islices;
{
    int             filedim, xsize, ky1, ky2;
    int             i, j, k, l;
    long            offset;
    filedim = data->matfile.ndim;
    xsize = data->matfile.size[0];
    ky1 = iy;
    ky2 = ky1 + islices;
    offset = data->matfile.fheadersz;

    /* FPR(stderr, "%d  %d  %d %d %d %d\n",filedim,ix,iy,iz,iz2,islices); */
    switch (filedim) {
/* delaglio */
    case 1:
       if (fseek(data->fileptr, offset, SEEK_SET))
          {
           FPR( stderr, "NV Error Seeking to Offset %ld\n", offset );
           return( 1 );
          }

        if (fwrite((void *)mat, sizeof(float), (size_t)xsize, data->fileptr) != (size_t)xsize)
           {
            FPR(stderr, "NV Error Writing %d 1D Data Points.\n", xsize );
            return( 1 );
           }

       break; 
    case 2:
	for (j = ky1; j < ky2; j++) {
	    for (i = 0; i < xsize; i++) {
		M_2F(data, i, j);
		*sb = *mat;
		mat++;
	    }
	}
	break;
    case 3:
	k = iz;
	for (j = ky1; j < ky2; j++) {
	    for (i = 0; i < xsize; i++) {
		M_3F(data, i, j, k);
		*sb = *mat;
		mat++;
	    }
	}
	break;
    case 4:
	l = iz2;
	k = iz;
	for (j = ky1; j < ky2; j++) {
	    for (i = 0; i < xsize; i++) {
		M_4F(data, i, j, k, l);
		*sb = *mat;
		mat++;
	    }
	}
	break;
    }
    return(0);
}
int open_file(data)
    DATASET        *data;
{
    char            fullname[255], parname[255];
    int             len, i;

    /*
     * Return if the dataset is already open. 
     */
    if (data->matfile.open == 1) {
	data->fi = 1;
	return (ERROR);
    }
    (void) strcpy(fullname, data->dirname);
    (void) strcat(fullname, "/");
    (void) strcat(fullname, data->filename);

#ifndef S_SPLINT_S
    if ((data->fileptr = fopen(fullname, "w+")) == NULL) {
	FPR(stderr, "Could not open file %s \n", fullname);
	data->fi = -1;
	return (ERROR);
    }
#endif

    data->matfile.open = 1;

    (void) DimDataset(data);
    (void) HeaderWriteNV(data);
    (void) DatasetDumpPar(data);
    (void) WriteEmptyFile(data);
    return (NOERROR);
}

int WriteEmptyFile(data)
    DATASET        *data;
{
    int             i, dim, ndim, delem, nblks_tot,size,bsize,nblks;
    float          *eblock;
    ndim = data->matfile.ndim;
    nblks_tot = 1;
    for (dim = 0; dim < ndim; dim++) {
	nblks_tot *= data->matfile.nblks[dim];
    }

    delem = data->matfile.blkelems;
    eblock = (float *) calloc((SIZE_T)delem, sizeof(float));
    for (i = 0; i < nblks_tot; i++) {
	if (fwrite((void *)eblock, sizeof(float), (size_t)delem, data->fileptr) != (size_t)delem) {
	    FPR(stderr, "Could not write block %d\n", i);
	}
    }
    free((VOID*)eblock);
    return(0);
}

int DimDataset(data)
    DATASET        *data;
{
    int             ndim;
    int             size[4];
    int             vsize[4];
    int             blksize[4];
    int             i, dim, nblks_tot;

    ndim = data->matfile.ndim;
    data->matfile.blkelems = 1;

    for (dim = 0; dim < ndim; dim++) {
	size[dim] = data->matfile.size[dim];
	vsize[dim] = data->matfile.vsize[dim];
	blksize[dim] = data->matfile.blksize[dim];

	data->matfile.nblks[dim] = size[dim] / blksize[dim];
        if (blksize[dim]*data->matfile.nblks[dim] < size[dim])
		data->matfile.nblks[dim] += 1;

	data->matfile.nvblks[dim] = vsize[dim] / blksize[dim];
	if (dim > 0) {
	    data->matfile.offblk[dim] = data->matfile.nblks[dim - 1] * data->matfile.offblk[dim - 1];
	    data->matfile.offvblk[dim] = data->matfile.nvblks[dim - 1] * data->matfile.offvblk[dim - 1];
	    data->matfile.offpt[dim] = data->matfile.blksize[dim - 1] * data->matfile.offpt[dim - 1];
	} else {
	    data->matfile.offblk[dim] = 1;
	    data->matfile.offvblk[dim] = 1;
	    data->matfile.offpt[dim] = 1;
	}
	data->matfile.blkmask[dim] = data->matfile.blksize[dim] - 1;
	data->matfile.blkelems = data->matfile.blkelems * data->matfile.blksize[dim];


	data->matfile.foldup[dim] = 0.0;
	data->matfile.folddown[dim] = 0.0;
    }
    for (dim = 1; dim < ndim; dim++) {
#ifdef S_SPLINT_S
	data->matfile.offpt[dim] = (int) (0.5 + (double)1.0);
#else
	data->matfile.offpt[dim] = (int) (0.5 + (log((double) data->matfile.offpt[dim]) / log(2.0)));
#endif
    }

    nblks_tot = data->matfile.nblks[0];
    for (dim = 1; dim < ndim; dim++)
	nblks_tot *= data->matfile.nblks[dim];

#ifndef S_SPLINT_S
    for (i = 0; i < nblks_tot; i++)
	data->matfile.blk[i] = NULL;
#endif

    return(0);
}

int DatasetDumpPar(data)
    DATASET        *data;
{
    int             i;
    FPR(stderr, "%d %d %d\n", data->matfile.fheadersz,
	    data->matfile.bheadersz,
	    data->matfile.ndim);
    for (i = 0; i < data->matfile.ndim; i++) {
	FPR(stderr, "%5d %5d %3d %3d %3d %4d",
		data->matfile.size[i],
		data->matfile.vsize[i],
		data->matfile.blksize[i],
		data->matfile.nblks[i],
		data->matfile.nvblks[i],
		data->matfile.offblk[i]);
	FPR(stderr, " %4d %3d",
		data->matfile.blkmask[i],
		data->matfile.offpt[i]);
	FPR(stderr, "%7.2f %8.2f %7.2f %7.2f %s \n",
		data->matfile.sf[i],
		data->matfile.sw[i],
		data->matfile.ref[i],
		data->matfile.refpt[i],
		data->matfile.label[i]);

    }
    return(0);
}

int HeaderWriteNV(data)
    DATASET        *data;
{
    int             delem, i;
    fhead.magic = 0x3418abcd;
    delem = (int)sizeof(fhead);
    data->matfile.fheadersz = delem;
    data->matfile.bheadersz = 0;
    fhead.fheadersz = data->matfile.fheadersz;
    fhead.bheadersz = data->matfile.bheadersz;
    fhead.blkelems = data->matfile.blkelems;
    fhead.ndim = data->matfile.ndim;
    FPR(stderr, "Header size %d\n", delem);
    for (i = 0; i < fhead.ndim; i++) {
	fhead.dim[i].sw = data->matfile.sw[i];
	fhead.dim[i].sf = data->matfile.sf[i];
	fhead.dim[i].ref = data->matfile.ref[i];
	fhead.dim[i].refpt = data->matfile.refpt[i];
	fhead.dim[i].refunits = data->matfile.refunits[i];
	fhead.dim[i].sw = data->matfile.sw[i];
	fhead.dim[i].size = data->matfile.size[i];
	fhead.dim[i].blksize = data->matfile.blksize[i];
	fhead.dim[i].size = data->matfile.size[i];
	fhead.dim[i].offblk = data->matfile.offvblk[i];
	fhead.dim[i].blkmask = data->matfile.blkmask[i];
	fhead.dim[i].offpt = data->matfile.offpt[i];
	fhead.dim[i].nblks = data->matfile.nvblks[i];
	(void) strcpy(fhead.dim[i].label, data->matfile.label[i]);
    }
    if (fwrite((void *)&fhead, (size_t)1, (size_t)delem, data->fileptr) != (size_t)delem) {
	FPR(stderr, "Could not write header\n");
    }
    return(0);
}
