/* C source code is found in dgemm_example.c */

#define min(x,y) (((x) < (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mkl.h"

#include <boost/numeric/ublas/matrix.hpp>
//#include <mkl_boost_ublas_matrix_prod.hpp>

#include <iostream>

using namespace std;

namespace ublas = boost::numeric::ublas;

int main()
{
  double *A, *B, *C,*A1,*B1,*C1;
    int m, n, k, i, j;
    double alpha, beta;

    printf ("\n This example computes real matrix C=alpha*A*B+beta*C using \n"
            " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
            " alpha and beta are double precision scalars\n\n");

    //m = 2000, k = 200, n = 1000;
    m = 2, k = 2, n = 2;
    printf (" Initializing data for matrix multiplication C=A*B for matrix \n"
            " A(%ix%i) and matrix B(%ix%i)\n\n", m, k, k, n);
    alpha = 1.0; beta = 0.0;

    printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
            " performance \n\n");

    //A = (double *)mkl_malloc( m*k*sizeof( double ), 64 );
    //B = (double *)mkl_malloc( k*n*sizeof( double ), 64 );
    //C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );


    A1 = new double[m*k]; //(double *) malloc ( m*k*sizeof( double ) );
    B1 =  new double[k*n]; //(double *) malloc( k*n*sizeof( double ) );
    C1 =  new double[m*n]; //(double *) malloc( m*n*sizeof( double ) );

    ublas::matrix<double> A2(2,2);
    ublas::matrix<double> B2(2,2);
    ublas::matrix<double> C2(2,2);
    
    /*  
    if (A == NULL || B == NULL || C == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      mkl_free(A);
      mkl_free(B);
      mkl_free(C);
      return 1;
    }
*/
    printf (" Intializing matrix data \n\n");
    for (i = 0; i < (m*k); i++) {
      // A[i] = (double)(i+1);
	A1[i]=(double) (i+1);
    }

    for (i = 0; i < (k*n); i++) {
      //  B[i] = (double)(-i-1);
	B1[i] = (double)(-i-1);
    }

    for (i = 0; i < (m*n); i++) {
      // C[i] = 0.0;
	C1[i] = 0.0;
    }

    for(i=0;i< n;++i)
      {
	for(j=0;j< n;++j)
	  {
	    A2(i,j)=i*j+1;
	    B2(i,j)=-i*j-1;
	    C2(i,j)=0.0;

	    A1[i*n+j]=i*j+1;
	    B1[i*n+j]=-i*j-1;
	    C1[i*n+j]=0;

	  }
      }
    
    int num=1000000;

    clock_t timec1=clock();
    printf (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");
    /*
    for (int i=0;i<num;++i)
      {
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	//      m, n, k, alpha, A, k, B, n, beta, C, n);
      }
    */
    clock_t timec2=clock();
    /*
    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", C[j+i*n]);
      }
      printf ("\n");
    }
*/

    
    clock_t timec3=clock();    
    for (i=0;i<num;++i)
      {
		for (j=0;j<10;++j)
      {
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
        //        m, n, k, alpha, A, k, B, n, beta, C, n);
	//C[0]=A[0]*B[0] + A[1]*B[2];
       	//C[1]=A[0]*B[1] + A[1]*B[3];
	//C[2]=A[2]*B[0] + A[3]*B[2];
	//C[3]=A[2]*B[1] + A[3]*B[3];

		C1[0]=A1[0]*B1[0] + A1[1]*B1[2];
		C1[1]=A1[0]*B1[1] + A1[1]*B1[3];
		C1[2]=A1[2]*B1[0] + A1[3]*B1[2];
		C1[3]=A1[2]*B1[1] + A1[3]*B1[3];

      }	
      }
      
    clock_t timec4=clock();
    // 	cout << C[0]+C[1]+C[2]+C[3] << endl;
    
    /*printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", C[j+i*n]);
      }
      printf ("\n");
      }*/


    clock_t timec5=clock();    

    for (i=0;i<num;++i)
      {
	for (j=0;j<10;++j)
      {


	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
        //        m, n, k, alpha, A, k, B, n, beta, C, n);

	//C1[0]=A1[0]*B1[0] + A1[1]*B1[2];
	//C1[1]=A1[0]*B1[1] + A1[1]*B1[3];
	//C1[2]=A1[2]*B1[0] + A1[3]*B1[2];
	//C1[3]=A1[2]*B1[1] + A1[3]*B1[3];

	C2(0,0)=A2(0,0)*B2(0,0) + A2(0,1)*B2(1,0);
	C2(0,1)=A2(0,0)*B2(0,1) + A2(0,1)*B2(1,1);
	C2(1,0)=A2(1,0)*B2(0,0) + A2(1,1)*B2(1,0);
	C2(1,1)=A2(1,0)*B2(0,1) + A2(1,1)*B2(1,1);

	//double g=C1[0]+C1[1]+C1[2]+C1[3];
	//C[0]=A[0]*B[0] + A[1]*B[2];
	//C[1]=A[0]*B[1] + A[1]*B[3];
	//C[2]=A[2]*B[0] + A[3]*B[2];
	//C[3]=A[2]*B[1] + A[3]*B[3];
	
	
	//	t+=1;
      
      
	//cout << i << endl;
	//	C[0]=A[0]*B[0] + A[1]*B[2];
	//C[1]=A[0]*B[1] + A[1]*B[3];
	//C[2]=A[2]*B[0] + A[3]*B[2];
	//C[3]=A[2]*B[1] + A[3]*B[3];
	//cout << C1[0]+C1[1]+C1[2]+C1[3] << endl;
      
      }
      }

    clock_t timec6=clock();


    for (i=0;i<num;++i)
      {
	for (j=0;j<10;++j)
      {


	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		    m, n, k, alpha, &(A2.data()[0]), k, &(B2.data()[0]), n, beta, &(C2.data()[0]), n);

	//C1[0]=A1[0]*B1[0] + A1[1]*B1[2];
	//C1[1]=A1[0]*B1[1] + A1[1]*B1[3];
	//C1[2]=A1[2]*B1[0] + A1[3]*B1[2];
	//C1[3]=A1[2]*B1[1] + A1[3]*B1[3];

	//C2(0,0)=A2(0,0)*B2(0,0) + A2(0,1)*B2(1,0);
	//C2(0,1)=A2(0,0)*B2(0,1) + A2(0,1)*B2(1,1);
	//C2(1,0)=A2(1,0)*B2(0,0) + A2(1,1)*B2(1,0);
	//C2(1,1)=A2(1,0)*B2(0,1) + A2(1,1)*B2(1,1);

	//double g=C1[0]+C1[1]+C1[2]+C1[3];
	//C[0]=A[0]*B[0] + A[1]*B[2];
	//C[1]=A[0]*B[1] + A[1]*B[3];
	//C[2]=A[2]*B[0] + A[3]*B[2];
	//C[3]=A[2]*B[1] + A[3]*B[3];
	
	
	//	t+=1;
      
      
	//cout << i << endl;
	//	C[0]=A[0]*B[0] + A[1]*B[2];
	//C[1]=A[0]*B[1] + A[1]*B[3];
	//C[2]=A[2]*B[0] + A[3]*B[2];
	//C[3]=A[2]*B[1] + A[3]*B[3];
	//cout << C1[0]+C1[1]+C1[2]+C1[3] << endl;
      
      }
      }
    clock_t timec7=clock();

    cout << "time:  " << double(timec4-timec3)/CLOCKS_PER_SEC << endl;
    cout << "time:  " << double(timec6-timec5)/CLOCKS_PER_SEC << endl;
    cout << "time:  " << double(timec7-timec6)/CLOCKS_PER_SEC << endl;
    
    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", C1[j+i*n]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", C2(i,j));
      }
      printf ("\n");
    }

    
    //cout << t << " " << num <<  endl;
    /*
    
    //std::cout << std::endl << "n:   " << n << std::endl;
    cout << "time1:  " << double(timec2-timec1)/CLOCKS_PER_SEC << endl;
    cout << "time2:  " << double(timec4-timec3)/CLOCKS_PER_SEC << endl;
    cout << "time3:  " << double(timec6-timec4)/CLOCKS_PER_SEC << endl;
    //cout << "time1:  " << double(timec2-timec1) << endl;
    //cout << "time2:  " << double(timec4-timec3) << endl;
    //cout << "time3:  " << double(timec6-timec4) << endl;
    
    printf ("\n Computations completed.\n\n");

    printf (" Top left corner of matrix A: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(k,6); j++) {
        printf ("%12.0f", A[j+i*k]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix B: \n");
    for (i=0; i<min(k,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.0f", B[j+i*n]);
      }
      printf ("\n");
    }
    
    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", C[j+i*n]);
      }
      printf ("\n");
    }
    */
    //printf ("\n Deallocating memory \n\n");
    //mkl_free(A);
    //mkl_free(B);
    //mkl_free(C);

    //printf (" Example completed. \n\n");
    return 0;
}
