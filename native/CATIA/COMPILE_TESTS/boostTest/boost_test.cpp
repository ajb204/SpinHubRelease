#include <complex>
#include <boost/math/special_functions/sign.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

//#include <boost-numeric-bindings/boost/numeric/bindings/lapack/workspace.hpp>
//#include <boost-numeric-bindings/boost/numeric/bindings/lapack/geev.hpp>

//#include <boost/numeric/bindings/lapack/workspace.hpp>
//#include <boost/numeric/bindings/lapack/geev.hpp>

//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

//#include <boost/numeric/bindings/lapack/workspace.hpp>
//#include <boost/numeric/bindings/lapack/driver/geev.hpp>

#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/lapack/driver/geev.hpp>

//#include "include/boost/numeric/bindings/lapack/driver/geev.hpp"







namespace ublas = boost::numeric::ublas;

using ublas::prod;
typedef double real_t;
typedef std::complex<real_t> cmplx_t;
typedef ublas::vector<cmplx_t> cvec_t;
typedef ublas::matrix<real_t, ublas::row_major> rrmat_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> rcmat_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> ccmat_t;

//ccmat_t' (aka 'matrix<complex<double>, basic_column_major<> >'






#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/lapack/driver/geev.hpp>

namespace lapack = boost::numeric::bindings::lapack;
//namespace lapack = boost::numeric::lapack;


//#include <boost/numeric/bindings/lapack/computational/getri.hpp>
//#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
//#include <MatrixExponential.hpp>

//#include <Catia.h>
//#include <Dataset.h>
//#include <Abort.h>




//#include <iostream>
//#include <iomanip>
//#include <complex>
//#include <cmath>







int main()
{

  const unsigned int n = 2;

cvec_t eigValVect(n);
ccmat_t lEigVectMat(n, n);
ccmat_t rEigVectMat(n, n);
ccmat_t freeEvolMatColMaj(n,n);

 cmplx_t fff(3,4);
 cmplx_t ggg(5,6);
 

 freeEvolMatColMaj(0,0)=fff;
 freeEvolMatColMaj(1,1)=ggg;
 
 std::cout << real(freeEvolMatColMaj(0,0)) << std::endl;
  std::cout << imag(freeEvolMatColMaj(0,0)) << std::endl;
  std::cout << freeEvolMatColMaj(1,1) << std::endl;
 
  //lapack::geev('N', 'V',freeEvolMatColMaj, eigValVect, NULL,NULL,lapack::optimal_workspace());
  lapack::geev('N', 'V',freeEvolMatColMaj, eigValVect, lEigVectMat,rEigVectMat,lapack::optimal_workspace());
  // lapack::geev("N","V",freeEvolMatColMaj, eigValVect,NULL,NULL,lapack::optimal_workspace());


  
  //std::cout << lapack::poo(6) << std::endl;

  //identity_matrix<double> m (3);
  //std::cout << m << std::endl;

  
  // int  lapack::geev (A &a, W &w, V *vl, V *vr, optimal_workspace);
// lapack::geev (freeEvolMatColMaj,eigValVect, lEigVectMat, rEigVectMat, lapack::optimal_workspace());
 
 return 0;
}
