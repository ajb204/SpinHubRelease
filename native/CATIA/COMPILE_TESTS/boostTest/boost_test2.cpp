#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/bindings/lapack/syevd.hpp>

#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

typedef boost::numeric::ublas::matrix<int> iMatrix;
typedef boost::numeric::ublas::matrix<double> dMatrix;
typedef boost::numeric::ublas::vector<int> iVector;
typedef boost::numeric::ublas::vector<double> dVector;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

int main() {
    int n = 10;
    dMatrix jacobi(n,n); // then actually initialize it
    dVector eigenvals(n);


    //int error = lapack::gesvd('S','S', jacobi, eigenvals, eigenvects1, eigenvects2);
    int error = lapack::syevd('V','L', jacobi, eigenvals, lapack::optimal_workspace() );

    std::cout << eigenvals << std::endl;
    std::cout << jacobi << std::endl;
    return 0;
}
