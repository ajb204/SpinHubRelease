
int   msize0,msize1,msize2,msize3,offb1,offb2,offb3,bmask0,bmask1,bmask2,bmask3,offp1,offp2,offp3;
float **mblk;

#ifndef S_SPLINT_S

#define M_2Setup(data)	msize0 = (int)(0.5+(log((double)data->matfile.blksize[0])/log(2.0))); \
			msize1 = (int)(0.5+(log((double)data->matfile.blksize[1])/log(2.0))); \
			offb1 = (int)(0.5+(log((double)data->matfile.offvblk[1])/log(2.0))); \
			mblk = data->matfile.blk;\
			bmask0 = data->matfile.blkmask[0];\
			bmask1 = data->matfile.blkmask[1];\
			offp1 = data->matfile.offpt[1];

#define M_2F(data,I,J)  jb = (I>>msize0) + ((J>>msize1)<<offb1);\
		        sb = mblk[jb];\
			if (sb == (float *) NULL) sb=RWBlk(data,jb,BLKREAD);\
			sb += (I & bmask0) + ((J & bmask1)<<offp1);

#define M_2B(data,I,J)  jb = (I>>msize0) + ((J>>msize1)<<offb1);\
		        bb = mblk[jb];\
			if (bb == (float *) NULL) {\
			if ((bb = RWBlk(data,jb,BLKREAD)) ==(float *) NULL) {\
				printf("Error at %d %d \n",I,J);\
				return(1);\
				}}

#define M_2BD(I,J) sb = bb + (I & bmask0) + ((J & bmask1)<<offp1);


#define M_3Setup(data)	msize0 = (int)(0.5+(log((double)data->matfile.blksize[0])/log(2.0))); \
			msize1 = (int)(0.5+(log((double)data->matfile.blksize[1])/log(2.0))); \
			msize2 = (int)(0.5+(log((double)data->matfile.blksize[2])/log(2.0))); \
			offb1 = (int)(0.5+(log((double)data->matfile.offvblk[1])/log(2.0))); \
			offb2 = (int)(0.5+(log((double)data->matfile.offvblk[2])/log(2.0))); \
			mblk = data->matfile.blk;\
			bmask0 = data->matfile.blkmask[0];\
			bmask1 = data->matfile.blkmask[1];\
			bmask2 = data->matfile.blkmask[2];\
			offp1 = data->matfile.offpt[1];\
			offp2 = data->matfile.offpt[2];

#define M_3F(data,I,J,K) jb = (I>>msize0) + ((J>>msize1)<<offb1) + ((K>>msize2)<<offb2);\
		       sb = mblk[jb];\
			if (sb == (float *) NULL) sb=RWBlk(data,jb,BLKREAD);\
			sb += (I & bmask0) + ((J & bmask1)<<offp1)+((K & bmask2)<<offp2);

#define M_3B(data,I,J,K) jb = (I>>msize0) + ((J>>msize1)<<offb1) + ((K>>msize2)<<offb2);\
		       bb = mblk[jb];\
			if (bb == (float *) NULL) {\
			bb = RWBlk(data,jb,BLKREAD);\
			if (bb == (float *) NULL) {\
				printf("%d %d\n",I,J);\
			}\
			}

#define M_3BD(I,J,K) sb = bb + (I & bmask0) + ((J & bmask1)<<offp1)+((K & bmask2)<<offp2);


#define M_4Setup(data)	msize0 = (int)(0.5+(log((double)data->matfile.blksize[0])/log(2.0))); \
			msize1 = (int)(0.5+(log((double)data->matfile.blksize[1])/log(2.0))); \
			msize2 = (int)(0.5+(log((double)data->matfile.blksize[2])/log(2.0))); \
			msize3 = (int)(0.5+(log((double)data->matfile.blksize[3])/log(2.0))); \
			offb1 = (int)(0.5+(log((double)data->matfile.offvblk[1])/log(2.0))); \
			offb2 = (int)(0.5+(log((double)data->matfile.offvblk[2])/log(2.0))); \
			offb3 = (int)(0.5+(log((double)data->matfile.offvblk[3])/log(2.0))); \
			mblk = data->matfile.blk;\
			bmask0 = data->matfile.blkmask[0];\
			bmask1 = data->matfile.blkmask[1];\
			bmask2 = data->matfile.blkmask[2];\
			bmask3 = data->matfile.blkmask[3];\
			offp1 = data->matfile.offpt[1];\
			offp2 = data->matfile.offpt[2];\
			offp3 = data->matfile.offpt[3];

#define M_4F(data,I,J,K,L) jb = (I>>msize0) + ((J>>msize1)<<offb1) + ((K>>msize2)<<offb2) + ((L>>msize3)<<offb3);\
		       sb = mblk[jb];\
			if (sb == (float *) NULL) sb=RWBlk(data,jb,BLKREAD);\
			sb += (I & bmask0) + ((J & bmask1)<<offp1)+((K & bmask2)<<offp2) + ((L & bmask3) << offp3);

#else

#define M_2Setup(data)	msize0 = (int)(0.5+(double)1.0); \
			msize1 = (int)(0.5+(double)1.0); \
			offb1  = (int)(0.5+(double)1.0); \
			mblk = data->matfile.blk;\
			bmask0 = data->matfile.blkmask[0];\
			bmask1 = data->matfile.blkmask[1];\
			offp1 = data->matfile.offpt[1];

#define M_2F(data,I,J)  jb = (I>>msize0) + ((J>>msize1)<<offb1);\
		        sb = mblk[jb];\
			if (sb == (float *) NULL) sb=RWBlk(data,jb,BLKREAD);\
			sb += (I & bmask0) + ((J & bmask1)<<offp1);

#define M_2B(data,I,J)  jb = (I>>msize0) + ((J>>msize1)<<offb1);\
		        bb = mblk[jb];\
			if (bb == (float *) NULL) {\
			if ((bb = RWBlk(data,jb,BLKREAD)) ==(float *) NULL) {\
				printf("Error at %d %d \n",I,J);\
				return(1);\
				}}

#define M_2BD(I,J) sb = bb + (I & bmask0) + ((J & bmask1)<<offp1);


#define M_3Setup(data)	msize0 = (int)(0.5+(double)1.0); \
			msize1 = (int)(0.5+(double)1.0); \
			msize2 = (int)(0.5+(double)1.0); \
			offb1  = (int)(0.5+(double)1.0); \
			offb2  = (int)(0.5+(double)1.0); \
			mblk = data->matfile.blk;\
			bmask0 = data->matfile.blkmask[0];\
			bmask1 = data->matfile.blkmask[1];\
			bmask2 = data->matfile.blkmask[2];\
			offp1 = data->matfile.offpt[1];\
			offp2 = data->matfile.offpt[2];

#define M_3F(data,I,J,K) jb = (I>>msize0) + ((J>>msize1)<<offb1) + ((K>>msize2)<<offb2);\
		       sb = mblk[jb];\
			if (sb == (float *) NULL) sb=RWBlk(data,jb,BLKREAD);\
			sb += (I & bmask0) + ((J & bmask1)<<offp1)+((K & bmask2)<<offp2);

#define M_3B(data,I,J,K) jb = (I>>msize0) + ((J>>msize1)<<offb1) + ((K>>msize2)<<offb2);\
		       bb = mblk[jb];\
			if (bb == (float *) NULL) {\
			bb = RWBlk(data,jb,BLKREAD);\
			if (bb == (float *) NULL) {\
				printf("%d %d\n",I,J);\
			}\
			}

#define M_3BD(I,J,K) sb = bb + (I & bmask0) + ((J & bmask1)<<offp1)+((K & bmask2)<<offp2);


#define M_4Setup(data)	msize0 = (int)(0.5+(double)1.0); \
			msize1 = (int)(0.5+(double)1.0); \
			msize2 = (int)(0.5+(double)1.0); \
			msize3 = (int)(0.5+(double)1.0); \
			offb1  = (int)(0.5+(double)1.0); \
			offb2  = (int)(0.5+(double)1.0); \
			offb3  = (int)(0.5+(double)1.0); \
			mblk = data->matfile.blk;\
			bmask0 = data->matfile.blkmask[0];\
			bmask1 = data->matfile.blkmask[1];\
			bmask2 = data->matfile.blkmask[2];\
			bmask3 = data->matfile.blkmask[3];\
			offp1 = data->matfile.offpt[1];\
			offp2 = data->matfile.offpt[2];\
			offp3 = data->matfile.offpt[3];

#define M_4F(data,I,J,K,L) jb = (I>>msize0) + ((J>>msize1)<<offb1) + ((K>>msize2)<<offb2) + ((L>>msize3)<<offb3);\
		       sb = mblk[jb];\
			if (sb == (float *) NULL) sb=RWBlk(data,jb,BLKREAD);\
			sb += (I & bmask0) + ((J & bmask1)<<offp1)+((K & bmask2)<<offp2) + ((L & bmask3) << offp3);

#endif

