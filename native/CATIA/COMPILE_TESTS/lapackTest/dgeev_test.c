/* Calling DGEEEV using col-major order */
/*15th Sept 2020                        */
/* avoding the need to use boost lapack bindings */

#include <stdio.h>
#include <time.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <MatrixExponential.hpp>

//#ifdef __APPLE__
//#endif

#include <lapack.h>

namespace ublas = boost::numeric::ublas;
using ublas::prod; //matrix product

typedef double real_t;
typedef std::complex<real_t> cmplx_t;
typedef ublas::vector<cmplx_t> cvec_t;
typedef ublas::matrix<real_t, ublas::row_major> rrmat_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> rcmat_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> ccmat_t;

//print matrices and vectors to the screen
void PrintMatrix(char *tag,rrmat_t a)
{
  printf("%s\n",tag);
  for(int i=0;i<a.size1();i++)
    {
      for(int j=0;j<a.size2();j++)
	printf("%lf ",a(i,j));
      printf("\n");
    }
}
void PrintMatrix(char *tag,ccmat_t a)
{
  printf("%s\n",tag);
  for(int i=0;i<a.size1();i++)
    {
      for(int j=0;j<a.size2();j++)
	printf("%lf+i%lf  ",real(a(i,j)),imag(a(i,j)));
      printf("\n");
    }
}
void PrintMatrix(char *tag,rcmat_t a)
{
  printf("%s\n",tag);
  for(int i=0;i<a.size1();i++)
    {
      for(int j=0;j<a.size2();j++)
	printf("%lf+i%lf  ",real(a(i,j)),imag(a(i,j)));
      printf("\n");
    }
}
void PrintVector(char *tag,cvec_t a)
{
  printf("%s\n",tag);
  for(int i=0;i<a.size();i++)
    printf("%lf+i%lf  ",real(a(i)),imag(a(i)));
  printf("\n");
}


//**************************//
//Eigenvalue expm functions
void expm_eig(rrmat_t a,   rcmat_t &G,  std::vector<double> par)
{
  int n=a.size1();//set size of matrix
  double wr[n];   //declare workspace
  double wi[n];   //declare workspace
  double *vl;     //declare workspace (null pointer)
  double vr[n*n]; //declare matrix for real right eigenvectors
  int ipv[n];     //pivot for inverse
  cvec_t val(n);  //for vector multiplciation
  rcmat_t vr_mat(n,n); //complex right eigenvectors
  rcmat_t vri_mat(n,n); //complex right inverse eigenvectors
  rcmat_t tempDiagMat(n,n); //temp array
  double work[4*n];  //workspace for dgeev
  std::complex<double> workT[4*n]; //workspace for moving columns around for zgetri
  int lwork=n*4;    //workspace for dgeev/zgetri
  int info;         //output of calc.
  dgeev_("N","V", &n,&(a.data()[0]), &n,wr, wi, vl,&n,vr,&n,work,&lwork,&info); //EV
  //POSSIBLE PLACE WHERE WE COULD GET ERRORS: //looks risky as one is double and other is complex double...seems to work...
  std::copy(vr, vr+n*n, vr_mat.data().begin());  //copy up right eigenvectors 
  for(int i=0;i<n-1;++i) //update complex eigenvalues if we have conjugate pairs.
    if(wi[i]+wi[i+1]==0 && abs(wi[i])>0) //conjugate pair test
      //std::cout << i << " " << i+1 << std::endl; //print pairs.
      for(int j=0;j<n;++j) //update complex eigenvalues
	{
	  vr_mat(i,j)=cmplx_t (vr[i*n+j],vr[(j)+(i+1)*n]);
	  vr_mat(i+1,j)=cmplx_t (vr[i*n+j],-vr[(j)+(i+1)*n]);
	}
  std::copy(vr_mat.data().begin(),vr_mat.data().end(),vri_mat.data().begin());  //copy up right eigenvectors to inverse
  zgetrf_(&n,&n,&(vri_mat.data()[0]),&n,ipv,&info);   //take inverse in two stages (fast)...
  zgetri_(&n,&(vri_mat.data()[0]),&n,ipv,workT,&lwork,&info); //and complete inverse.
  for (unsigned int i = 0; i < par.size(); i++) //for each parameter in list...
    {
      for (unsigned int j = 0; j < n; ++j) //exp(val*par) . vr
	{
	  std::complex<double> pp=cmplx_t (par[i]*wr[j],par[i]*wi[j]);
	  //std::cout << "values: " << real(val(j)) << " "  << imag(val(j)) << " " <<par[i] << " " << real(pp) << " " << imag(pp) << std::endl;
	  for (unsigned int k = 0; k < n; ++k)
	    tempDiagMat(j,k) = exp(pp)*vr_mat(j,k);
	}
      G = prod(vri_mat, tempDiagMat); //vr^-1.exp(val*par) .vr
    }
}


void expm_eig(rcmat_t a,   rcmat_t &G,  std::vector<double> par)
{
  int n=a.size1(); //set size of matrix
  std::complex<double> w[n];  //declare workspace
  std::complex<double> *vl;   //declare workspace (null pointer)
  std::complex<double> work[2*n];   //declare workspace
  int lwork=2*n;
  double rwork[2*n];
  int ipv[n];
  int info;
  rcmat_t tempDiagMat(n,n);  
  rcmat_t vr_mat(n,n); //right eigenvectors for matrix multiplcaition
  rcmat_t vri_mat(n,n); //right eigenvector inverse for matrix multipication
  zgeev_("N","V", &n,&(a.data()[0]), &n,w, vl,&n,&(vr_mat.data()[0]),&n,work,&lwork,rwork,&info); //EV
  std::copy(vr_mat.data().begin(), vr_mat.data().end(), vri_mat.data().begin());  //copy up right eigenvectors
  zgetrf_(&n,&n,&(vri_mat.data()[0]),&n,ipv,&info);   //take inverse in two stages (fast)...
  zgetri_(&n,&(vri_mat.data()[0]),&n,ipv,work,&lwork,&info);
  for (int i = 0; i < par.size(); ++i) //for each parameter in list...
    {
      for (int j = 0; j < n; ++j) //exp(val*par) . vr
	{
	  std::complex<double> pp=par[i]*w[j];
	  for (int k = 0; k < n; ++k)
	    tempDiagMat(j, k) = exp(pp)*vr_mat(j,k);
	}
      G = prod(vri_mat, tempDiagMat); //vr^-1.exp(val*par) .vr
    }
}




//function to compare eig versus pade matrix exponentials
//both in accuracy and in speed for various sizes.
int main (int argc, const char * argv[])
{

  //double aaa[3][3]= {1,2,3,4,5,6,7,8,9}; //col major ((1,2,3),(4,5,6),(7,8,9))

  double test,diag; //test real matricies
  std::vector<double> par;
  par.push_back(1);
  
  int grid=10000;
  int Ngrid=15;
  FILE* fp;
  fp=fopen("test.out","w");
  for(int k=0;k<Ngrid;++k)
    {
      int n=k+2;//matrix size from 2 to 2+Ngrid

      //rrmat_t aa(n,n);
      rrmat_t aa1(n,n); //real test matricies
      rcmat_t aa2(n,n); //complex test matricies
      
      rcmat_t G1(n,n); //real expm_eig
      rcmat_t G2(n,n); //real expm_pade
      rcmat_t G3(n,n); //complex expm_eig
      rcmat_t G4(n,n); //complex expm_pade
      
      for(int i=0;i<n;++i) //set aa
	for(int j=0;j<n;++j)
	  {
	    aa1(i,j)=(j+i*n+1)*-1E-6;
	    aa2(i,j)=cmplx_t ((j+i*n+1)*-1E-6,(j+i*n+1)*-1E-6);
	    //std::cout << aa2(i,j) << std::endl;
	  }
      
      clock_t time1=clock();
      for(int i=0;i<grid;++i)
	expm_eig(aa1,G1,par); //real, eig
      
      clock_t time2=clock();
      for(int i=0;i<grid;++i)
	G2 = ublas::expm_pad(aa1); //real, pade
      clock_t time3=clock();
      
      for(int i=0;i<grid;++i)
	expm_eig(aa2,G3,par); //complex, eig
      clock_t time4=clock();
      
      for(int i=0;i<grid;++i)
	G4 = ublas::expm_pad(aa2); //complex, pade
      clock_t time5=clock();
      
      std::cout << std::endl << "n:   " << n << std::endl;
      std::cout << "Realexpm_pad: " << double(time3-time2)/CLOCKS_PER_SEC << std::endl;
      std::cout << "expm_eig: " << double(time2-time1)/CLOCKS_PER_SEC << std::endl;
      std::cout << "speedUp:  " << double( (time3-time2))/double((time2-time1)) << std::endl;
      std::cout << "Complexexpm_pad: " << double(time5-time4)/CLOCKS_PER_SEC << std::endl;
      std::cout << "expm_eig: " << double(time4-time3)/CLOCKS_PER_SEC << std::endl;
      std::cout << "speedUp:  " << double( (time5-time4))/double((time4-time3)) << std::endl;
      
      fprintf(fp,"%i\t%f\t%f\t%f\t%f\n",n,double(time2-time1)/CLOCKS_PER_SEC,double(time3-time2)/CLOCKS_PER_SEC,double(time4-time3)/CLOCKS_PER_SEC,double(time5-time4)/CLOCKS_PER_SEC);
      
      test=0;diag=0; //test real results
      for(int i=0;i<n;++i)
	{
	  diag+=(abs(G1(i,i))+abs(G2(i,i)))/2.;
	  for(int j=0;j<n;++j)
	    test+=abs(G1(i,j)-G2(i,j));
	  
	}
      std::cout << "real test : " << test << " diag: " << diag << " err: " << test/diag << std::endl;
      if(test/diag>1E-9)
	{
	  PrintMatrix("G1",G1);
	  PrintMatrix("G2",G2);
	  std::cout << "FAIL: the two matrix exponentials are very different" << std::endl;
	  exit(100);
	}
      
      test=0;diag=0; //test complex results
      for(int i=0;i<n;++i)
	{
	  diag+=(abs(G3(i,i))+abs(G4(i,i)))/2.;
	  for(int j=0;j<n;++j)
	    test+=abs(G3(i,j)-G4(i,j));
	  
	}
      std::cout << "complex test : " << test << " diag: " << diag << " err: " << test/diag << std::endl;
      if(test/diag>1E-6)
	{
	  PrintMatrix("G3",G3);
	  PrintMatrix("G4",G4);
	  std::cout << "FAIL: the two matrix exponentials are very different" << std::endl;
	  exit(100);
	}
    }
  fclose(fp);
  
  fp=fopen("gnu.gp","w");
  fprintf(fp,"set term post eps enh color solid\n");
  fprintf(fp,"set output 'plot.eps'\n");
  fprintf(fp,"set xlabel 'N'\n");
  fprintf(fp,"set title 'real/complex pade versus eig matrix exponentials'\n");
  fprintf(fp,"set ylabel 'time'\n");
  fprintf(fp,"set key left\n");
  fprintf(fp,"plot 'test.out' u 1:2 ti 'real eig','' u 1:3 ti 'real pade','' u 1:4 ti 'complex eig','' u 1:5 ti 'complex pade'\n");
  fprintf(fp,"set output 'ratio.eps'\n");
  fprintf(fp,"set xlabel 'N'\n");
  fprintf(fp,"set title 'real/complex pade versus eig matrix exponentials'\n");
  fprintf(fp,"set ylabel 'pade/eig'\n");
  fprintf(fp,"set key left\n");
  fprintf(fp,"plot 'test.out' u 1:($3/$2) ti 'real','' u 1:($5/$4) ti 'complex'\n");
  fclose(fp);
  system("gnuplot gnu.gp");
  
  return 0;
}
