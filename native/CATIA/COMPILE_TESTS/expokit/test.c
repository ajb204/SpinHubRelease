/* Calling DGEEEV using col-major order */
/*15th Sept 2020                        */
/* avoding the need to use boost lapack bindings */

#include <stdio.h>
#include <time.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <MatrixExponential.hpp> //real expm in c++
#include <lapack.h>              //expm via eigenvalues and lapack
#include <lapackexpo.h>         //expm using expokit



namespace ublas = boost::numeric::ublas;
using ublas::prod; //matrix product
using ublas::diagonal_matrix; //matrix product

typedef double real_t;
typedef std::complex<real_t> cmplx_t;
typedef ublas::vector<cmplx_t> cvec_t;
typedef ublas::vector<real_t> rvec_t;
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
void PrintVector(char *tag,rvec_t a)
{
  printf("%s\n",tag);
  for(int i=0;i<a.size();i++)
    printf("%lf  ",a(i));
  printf("\n");
}


//need to transpse a if we are doing matrix times vec and using rrmat.
//void prodNew2(ublas::matrix<double> &a,ublas::vector<double> &b,ublas::vector<double> &c)
void prodNew2(rrmat_t &a,rvec_t &b,rrmat_t &c)
{
  int n=a.size1(); //rows/cols of a, rows of b/c
  int k=1;         //cols of b and c
  double alpha=1.;  //alpha . a . b + beta.c. To be saved in C
  double beta=0;
  //a=trans(a);
  double cc[n];
  dgemv_("T",&n,&n,&alpha,&(a.data()[0]),&n,&(b.data()[0]),&k,&beta,&(c.data()[0]),&k);
}


//need to transpse a if we are doing matrix times vec and using rrmat.
void prodNew(ublas::matrix<double> &a,ublas::vector<double> &b,ublas::vector<double> &c)
{
  int n=a.size1(); //rows/cols of a, rows of b/c
  int k=1;         //cols of b and c
  double alpha=1.;  //alpha . a . b + beta.c. To be saved in C
  double beta=0.;
  //a=trans(a);
  dgemm_("T","N",&n,&k,&n,&alpha,&(a.data()[0]),&n,&(b.data()[0]),&n,&beta,&(c.data()[0]),&n);
}
//need to transpse a if we are doing matrix times vec and using rrmat.
void prodNew(rrmat_t &a,rvec_t &b,rrmat_t &c)
{
  int n=a.size1(); //rows/cols of a, rows of b/c
  int k=1;         //cols of b and c
  double alpha=1.;  //alpha . a . b + beta.c. To be saved in C
  double beta=0.;
  //a=trans(a);
  dgemm_("T","N",&n,&k,&n,&alpha,&(a.data()[0]),&n,&(b.data()[0]),&n,&beta,&(c.data()[0]),&n);
}

//no need to transpose a or b
void prodNew(rrmat_t &a,rrmat_t &b,rrmat_t &c)
{
  int n=a.size1(); //rows/cols of a, rows of b/c
  int k=b.size2();         //cols of b and c
  double alpha=1.;  //alpha . a . b + beta.c. To be saved in C
  double beta=0.;
  //a=trans(a);
  dgemm_("N","N",&n,&k,&n,&alpha,&(a.data()[0]),&n,&(b.data()[0]),&n,&beta,&(c.data()[0]),&n);
}
//no need to transpose a or b
void prodNew(rcmat_t &a,rcmat_t &b,rcmat_t &c)
{
  int n=a.size1(); //rows/cols of a, rows of b/c
  int k=b.size2();         //cols of b and c
  std::complex <double> alpha=1.;  //alpha . a . b + beta.c. To be saved in C
  std::complex<double> beta=0.;
  //a=trans(a);
  zgemm_("N","N",&n,&k,&n,&alpha,&(a.data()[0]),&n,&(b.data()[0]),&n,&beta,&(c.data()[0]),&n);
}





//**************************//
//Eigenvalue expm functions

void expm_eig(rrmat_t a,   rcmat_t &G,  std::vector<double> &par)
{
  int n=a.size1();//set size of matrix
  double wr[n];   //declare workspace
  double wi[n];   //declare workspace
  double *vl;     //declare workspace (null pointer)
  double vr[n*n]; //declare matrix for real right eigenvectors
  int ipv[n];     //pivot for inverse
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
      //std::vector<std::complex<double> > parss;
      //rcmat_t test(n,n);
      //for (unsigned int j = 0; j < n; ++j) //exp(val*par) . vr
	//parss[j]=cmplx_t (par[i]*wr[j],par[i]*wi[j]);
      //diagonal_matrix <std::complex<double> > test;
      //test=diagMatrix(parss.size(),parss.data());

      for (unsigned int j = 0; j < n; ++j) //exp(val*par) . vr
	{
	  std::complex<double> epp=exp(cmplx_t (par[i]*wr[j],par[i]*wi[j]));
	  //std::cout << "values: " << real(val(j)) << " "  << imag(val(j)) << " " <<par[i] << " " << real(pp) << " " << imag(pp) << std::endl;
	  for (unsigned int k = 0; k < n; ++k)
	    //tempDiagMat(j,k) = exp(pp)*vr_mat(j,k);
	    tempDiagMat(j,k) = epp*vr_mat(j,k);
	}
      //G = prod(vri_mat, tempDiagMat); //vr^-1.exp(val*par) .vr
      //G = prod(vri_mat, tempDiagMat); //vr^-1.exp(val*par) .vr
      prodNew(tempDiagMat,vri_mat,G);//reversing the direction: not sure why... but seems to work.
      //rcmat_t GG(n,n);
      //G=prod(vri_mat, tempDiagMat); //vr^-1.exp(val*par) .vr
      //PrintMatrix("GG",GG);
      //PrintMatrix("G",G);
      //exit(100);
    }
}


void expm_eig(rcmat_t a,   rcmat_t &G,  std::vector<double> &par)
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
	  std::complex<double> epp=exp(par[i]*w[j]);
	  for (int k = 0; k < n; ++k)
	    tempDiagMat(j, k) = epp*vr_mat(j,k);
	}
      prodNew(tempDiagMat,vri_mat,G);//reversing the direction: not sure why... but seems to work.
      //G = prod(vri_mat, tempDiagMat); //vr^-1.exp(val*par) .vr
    }
}





void expm_padeExpo(rrmat_t &a,   rcmat_t &G,  std::vector<double> &par)
{
  int n=a.size1();//set size of matrix
  int ideg=6;              //the degree of the diagonal pade to be used. 6 typically good.
  double work[4*n*n+ideg+1];  //workspace for dgeev
  int lwork=4*n*n+ideg+1;    //workspace for dgeev/zgetri
  int info;         //output of calc.
  double ipiv[n];  //workspace?
  int iexph;  //will store position of output in work
  int ns;     //number of scaling-squaring used
  for(int i=0;i<par.size();++i)
    {
      dgpadm_(&ideg,&n,&par[i],&(a.data()[0]), &n,work,&lwork,ipiv,&iexph,&ns,&info); //EV
      std::copy(work+iexph-1, work+iexph+n*n-1, G.data().begin());  //copy up result
    }
}


void expm_padeExpo(rcmat_t &a,   rcmat_t &G,  std::vector<double> &par)
{
  int n=a.size1();//set size of matrix
  int ideg=6;              //the degree of the diagonal pade to be used. 6 typically good.
  std::complex<double> work[4*n*n+ideg+1];  //workspace for dgeev
  int lwork=4*n*n+ideg+1;    //workspace for dgeev/zgetri
  int info;         //output of calc.
  std::complex<double> ipiv[n];  //workspace?
  int iexph;  //will store position of output in work
  int ns;     //number of scaling-squaring used
  for(int i=0;i<par.size();++i)
    {
      zgpadm_(&ideg,&n,&par[i],&(a.data()[0]), &n,work,&lwork,ipiv,&iexph,&ns,&info); //EV
      std::copy(work+iexph-1, work+iexph+n*n-1, G.data().begin());  //copy up result
    }
}



void expm_pad(rrmat_t &a,rcmat_t &G,std::vector<double> &par)
{
  for(int i=0;i<par.size();++i)
    {
      rrmat_t aa=a;
      aa*=par[i];
      G = ublas::expm_pad(aa); //real, pade
    }
}
void expm_pad(rcmat_t &a,rcmat_t &G,std::vector<double> &par)
{
  for(int i=0;i<par.size();++i)
    {
      rcmat_t aa=a;
      aa*=par[i];
      G = ublas::expm_pad(aa); //real, pade
    }
}



//function to compare eig versus pade matrix exponentials
//both in accuracy and in speed for various sizes.
int main (int argc, const char * argv[])
{

  //double aaa[3][3]= {1,2,3,4,5,6,7,8,9}; //col major ((1,2,3),(4,5,6),(7,8,9))

  double test,diag; //test real matricies

  int grid=100; //number of repeats
  int Ngrid=15;  //matrix sizes...
  int Pgrid=20;  //matrix sizes...


  FILE* fp;
  fp=fopen("test.out","w");fclose(fp);
  for(int p=0;p<Pgrid;++p)
    {
      int pars=p*5+5;
      std::vector<double> par;
      for(int i=0;i<pars;++i)
	par.push_back(-1E-6);

      
      for(int k=0;k<Ngrid;++k)
	{
	  int n=k+2;//matrix size from 2 to 2+Ngrid
	  //rrmat_t aa(n,n);
	  rrmat_t aa1a(n,n); //real test matricies *no par mult
	  rrmat_t aa1b(n,n); //real test matricies *no par mult

	  rcmat_t aa2a(n,n); //complex test matricies

	  rrmat_t C1(n,1); //real test matricies *no par mult
	  rrmat_t C2(n,1); //real test matricies *no par mult
	  rrmat_t C3(n,1); //real test matricies *no par mult
	  
	  rcmat_t G1(n,n); //real expm_eig
	  rcmat_t G1a(n,n); //real expm_eig
	  rcmat_t G2(n,n); //real expm_pade
	  
	  rcmat_t G3(n,n); //complex expm_eig
	  rcmat_t G3a(n,n); //complex expm_eig
	  rcmat_t G4(n,n); //complex expm_pade

	  rvec_t b(n);  //vector for multiplying
	  rrmat_t b2(n,1);
	  for(int i=0;i<n;++i) //set aa
	    {
	      b(i)=i*1.;
	      b2(i,0)=i*1;
	      for(int j=0;j<n;++j)
		{
		  aa1a(i,j)=(j+i*n+1);
		  aa2a(i,j)=cmplx_t ((j+i*n+1),(j+i*n+1));
		  //std::cout << aa2(i,j) << std::endl;
		}
	    }
	  

	  //PrintMatrix("a",aa1a);
	  //PrintMatrix("b",b);
	  


	  clock_t timec1=clock();
	  for(int i=0;i<grid*100;++i)
	    prodNew(aa1a,b,C1);
	  clock_t timec2=clock();
	  for(int i=0;i<grid*100;++i)
	    C2=prod(aa1a,b2);
	  clock_t timec3=clock();
	  for(int i=0;i<grid*100;++i)
	    prodNew2(aa1a,b,C3);
	  clock_t timec4=clock();


	  
	  std::cout << std::endl << "n:   " << n << std::endl;
	  std::cout << "prodNew:  " << double(timec2-timec1)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "prod:     " << double(timec3-timec2)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "speedUp:  " << double( (timec3-timec2))/double((timec2-timec1)) << std::endl;
	  std::cout << "prodNew2: " << double(timec4-timec3)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "speedUp:  " << double( (timec2-timec1))/double((timec4-timec3)) << std::endl;
	  std::cout << std::endl;

	  
	  //PrintMatrix("a",aa1a);
	  //PrintMatrix("b",b);

	  //PrintMatrix("c1",C1);
	  //PrintMatrix("c2",C2);

	  test=0;diag=0; //test complex results
	  for(int i=0;i<C1.size1();++i)
	    {
	      if(C1.size2()==1)
		diag+=(abs(C1(i,0))+abs(C2(i,0)))/2.;
	      else
		diag+=(abs(C1(i,i))+abs(C2(i,i)))/2.;
	      for(int j=0;j<C1.size2();++j)
		test+=abs(C1(i,j)-C2(i,j));
	      
	    }
	  std::cout << "complex test : " << test << " diag: " << diag << " err: " << test/diag << std::endl;
	  if(test/diag>1E-6)
	    {
	      PrintMatrix("C1",C1);
	      PrintMatrix("C2",C2);
	      std::cout << "FAIL: the two matrix exponentials are very different" << std::endl;
	      exit(100);
	    }
	  test=0;diag=0; //test complex results
	  for(int i=0;i<C1.size1();++i)
	    {
	      if(C1.size2()==1)
		diag+=(abs(C1(i,0))+abs(C3(i,0)))/2.;
	      else
		diag+=(abs(C1(i,i))+abs(C3(i,i)))/2.;
	      for(int j=0;j<C1.size2();++j)
		test+=abs(C1(i,j)-C3(i,j));
	      
	    }
	  std::cout << "complex test : " << test << " diag: " << diag << " err: " << test/diag << std::endl;
	  if(test/diag>1E-6)
	    {
	      PrintMatrix("C1",C1);
	      PrintMatrix("C3",C3);
	      std::cout << "FAIL: the two matrix exponentials are very different" << std::endl;
	      exit(100);
	    }

	  
	  

	  clock_t time1=clock();
	  for(int i=0;i<grid;++i)
	    expm_eig(aa1a,G1,par); //real, eig
	  clock_t time2=clock();
	  for(int i=0;i<grid;++i)
	    expm_pad(aa1a,G2,par);
	  clock_t time3=clock();
	  for(int i=0;i<grid;++i)
	    expm_padeExpo(aa1a,G1a,par); //real, eig
	  clock_t time3a=clock();
	  
	  
	  for(int i=0;i<grid;++i)
	    expm_eig(aa2a,G3,par); //complex, eig
	  clock_t time4=clock();
	  for(int i=0;i<grid;++i)
	    expm_pad(aa2a,G4,par);
	  //G4 = ublas::expm_pad(aa2); //complex, pade
	  clock_t time5=clock();
	  for(int i=0;i<grid;++i)
	    expm_padeExpo(aa2a,G3a,par); //complex, eig
	  clock_t time5a=clock();
	  
	  
	  std::cout << std::endl << "n:   " << n << std::endl;
	  std::cout << "Realexpm_pad:    " << double(time3-time2)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "expokit_pad:     " << double(time3a-time3)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "expm_eig:        " << double(time2-time1)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "speedUpEig:      " << double( (time3-time2))/double((time2-time1)) << std::endl;
	  std::cout << "speedUpExpokit:  " << double( (time3-time2))/double((time3a-time3)) << std::endl;
	  std::cout << "Expokit/Eig:     " << double( (time2-time1))/double((time3a-time3)) << std::endl;
	  
	  
	  std::cout << "Complexexpm_pad: " << double(time5-time4)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "expokit_pad:     " << double(time5a-time5)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "expm_eig:        " << double(time4-time3)/CLOCKS_PER_SEC << std::endl;
	  std::cout << "speedUpEig:      " << double( (time5-time4))/double((time4-time3)) << std::endl;
	  std::cout << "speedUpExpokit:  " << double( (time5-time4))/double((time5a-time5)) << std::endl;
	  std::cout << "Expokit/Eig:     " << double( (time4-time3))/double((time3a-time3)) << std::endl;
	  
	  
	  fp=fopen("test.out","a");
	  fprintf(fp,"%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\n",n,pars,double(time2-time1)/CLOCKS_PER_SEC,double(time3-time2)/CLOCKS_PER_SEC,double(time3a-time3)/CLOCKS_PER_SEC,double(time4-time3)/CLOCKS_PER_SEC,double(time5-time4)/CLOCKS_PER_SEC,double(time5a-time5)/CLOCKS_PER_SEC);
	  fclose(fp);
	  
	  //test g1/g2
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
	  
	  //test g1/g1a
	  test=0;diag=0; //test real results
	  for(int i=0;i<n;++i)
	    {
	      diag+=(abs(G1(i,i))+abs(G1a(i,i)))/2.;
	      for(int j=0;j<n;++j)
		test+=abs(G1(i,j)-G1a(i,j));
	      
	    }
	  std::cout << "real test : " << test << " diag: " << diag << " err: " << test/diag << std::endl;
	  if(test/diag>1E-9)
	    {
	      PrintMatrix("G1",G1);
	      PrintMatrix("G1a",G1a);
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
	  
	  
	  test=0;diag=0; //test complex results
	  for(int i=0;i<n;++i)
	    {
	      diag+=(abs(G3(i,i))+abs(G3a(i,i)))/2.;
	      for(int j=0;j<n;++j)
		test+=abs(G3(i,j)-G3a(i,j));
	      
	    }
	  std::cout << "complex test : " << test << " diag: " << diag << " err: " << test/diag << std::endl;
	  if(test/diag>1E-6)
	    {
	      PrintMatrix("G3",G3);
	      PrintMatrix("G3a",G3a);
	      std::cout << "FAIL: the two matrix exponentials are very different" << std::endl;
	      exit(100);
	    }

	}
      fp=fopen("test.out","a");
      fprintf(fp,"\n");
      fclose(fp);
    }
  


  /*
  fp=fopen("gnu.gp","w");
  fprintf(fp,"set term post eps enh color solid\n");
  fprintf(fp,"set output 'plot.eps'\n");
  fprintf(fp,"set xlabel 'N'\n");
  fprintf(fp,"set title 'real/complex pade versus eig matrix exponentials'\n");
  fprintf(fp,"set ylabel 'time'\n");
  fprintf(fp,"set key left\n");
  fprintf(fp,"plot 'test.out' u 1:2 ti 'real eig','' u 1:3 ti 'real pade','' u 1:4 ti 'real pad expkit','' u 1:5 ti 'complex eig','' u 1:6 ti 'complex complex pade','' u 1:7 ti 'complex pade expo'\n");
  fprintf(fp,"set output 'ratio.eps'\n");
  fprintf(fp,"set xlabel 'N'\n");
  fprintf(fp,"set title 'real/complex pade versus eig matrix exponentials'\n");
  fprintf(fp,"set ylabel 'New/OldPade'\n");
  fprintf(fp,"set key left\n");
  fprintf(fp,"plot 'test.out' u 1:($3/$2) ti 'real eig','' u 1:($3/$4) ti 'real pade expo','' u 1:($6/$5) ti 'complex eig','' u 1:($6/$7) ti 'complex pade expkit' \n");
  fclose(fp);
  system("gnuplot gnu.gp");
  */
  return 0;
}
