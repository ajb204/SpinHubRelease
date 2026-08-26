#include<Catia.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <lapack.h>

namespace ublas = boost::numeric::ublas;

//types of data structure
typedef double real_t;
typedef std::complex<real_t> cmplx_t;
typedef ublas::vector<cmplx_t> cvec_t;
typedef ublas::matrix<real_t, ublas::row_major> rrmat_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> rcmat_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> ccmat_t;


double Catia::MajorPeakIntensity(double pb, double kex, double DeltaOmega, double IntA, double IntB) {
	// Calculate the intensity of the major peak
	// From Hansen & Led, JMR, 2003 _163_, p. 215-227
	std::complex<double> Intensity;
	std::complex<double> kA = kex * pb;
	std::complex<double> kB = kex * (1. - pb);
	std::complex<double> k2A = kA - sqrt(std::complex<double>(-1.)) * 0.;
	std::complex<double> k2B = kB - sqrt(std::complex<double>(-1.)) * DeltaOmega;

	std::complex<double> theta1 = 0.5 * (-(k2A + k2B) + sqrt((k2A - k2B) * (k2A - k2B) + 4. * kA * kB));
	std::complex<double> theta2 = 0.5 * (-(k2A + k2B) - sqrt((k2A - k2B) * (k2A - k2B) + 4. * kA * kB));

	// THETA1+THETA3
	std::complex<double> Intensity0 = (IntA * (kA - theta2 - k2A) + IntB * (kB + theta1 + k2A)) / (theta1 - theta2);

	// THETA2+THETA4
	std::complex<double> Intensity1 = (IntA * (-kA + theta1 + k2A) + IntB * (-kB - theta2 + k2A)) / (theta1 - theta2);
	if (fabs(real(Intensity1)) > fabs(real(Intensity0))) {
		Intensity = Intensity1;
	} else {
		Intensity = Intensity0;
	}

	return abs(Intensity);
}

double Catia::MajorPeakIntensity3States(double pb, double pc, double kex_ab, double kex_ac, double kex_bc, double deltaOmega_ab, double deltaOmega_ac, double IntA, double IntB,
		double IntC) {

	const double pa = 1.0 - pb - pc;

	const double kex_p_ab = kex_ab / (pa + pb);
	const double kex_p_ac = kex_ac / (pa + pc);
	const double kex_p_bc = kex_bc / (pb + pc);

	double k_ab = kex_p_ab * pb;
	double k_ba = kex_p_ab * pa;
	double k_ac = kex_p_ac * pc;
	double k_ca = kex_p_ac * pa;
	double k_bc = kex_p_bc * pc;
	double k_cb = kex_p_bc * pb;

	std::complex<double> I(0.0, 1.0);

	ublas::matrix<std::complex<double>, ublas::column_major> R(3, 3);
	//ccmat_t R(3,3);

	R(0, 0) =  +k_ab + k_ac;
	R(1, 0) =  -k_ab;
	R(2, 0) =  -k_ac;

	R(0, 1) =  -k_ba;
	R(1, 1) =  +k_ba + k_bc - I * deltaOmega_ab;
	R(2, 1) =  -k_bc;

	R(0, 2) =  -k_ca;
	R(1, 2) =  -k_cb;
	R(2, 2) =  +k_ca + k_cb - I * deltaOmega_ac;



	/*
	ublas::vector<std::complex<double> > EVal(3);
	ublas::matrix<std::complex<double>, ublas::column_major> EigVecL(1, 1);
	ublas::matrix<std::complex<double>, ublas::column_major> EigVecR(3, 3);

    lapack::geev('N','V', R, EVal, EigVecL, EigVecR, lapack::optimal_workspace());

	ublas::matrix<std::complex<double>, ublas::column_major > IEigVecR(EigVecR);

    ublas::vector<int> ipiv (3);   // pivot vector
    lapack::getrf (IEigVecR, ipiv);  // no lu_factor() alias for getrf() available
    lapack::getri (IEigVecR, ipiv);  // no lu_invert() alias for getrf() available

	unsigned int happyFaceIndex = 0;
//	double bestEVal = fabs(real(EVal[0]));
//	for (unsigned int i = 1; i < 3; ++i) {
//		double rEVal = fabs(real(EVal[i]));
//		if (rEVal < bestEVal) {
//			bestEVal = rEVal;
//			happyFaceIndex = i;
//		}
//	}





	ublas::matrix<std::complex<double> > M(3, 3);
	for (unsigned int i = 0; i < 3; ++i) {
		for (unsigned int j = 0; j < 3; ++j) {
			M(i, j) = EigVecR(i, happyFaceIndex) * IEigVecR(happyFaceIndex, j);
		}
	}
	*/

	ublas::matrix<std::complex<double> > M(3, 3); //array to save....
	{
	  int n=R.size1(); //set size of matrix


	  std::complex<double> w[n];  //declare workspace
	  std::complex<double> *vl;   //declare workspace (null pointer)
	  std::complex<double> vr[n*n];   //declare workspace
	  std::complex<double> work[4*n]; //workspace for moving columns around
	  int lwork=4*n;
	  double rwork=3*n-2; //max value from documentation
	  int info;
	  zgeev_("N","V", &n,&(R.data()[0]), &n,w, vl,&n,vr,&n,work,&lwork,&rwork,&info); //EV
	  cvec_t val(n);  //for vector multiplciation
	  for(int i=0;i<n;++i) //update complex eigenvalues
	    val(i)=w[i];//cmplx_t (wr[i],wi[i]); 
	  rcmat_t vr_mat(n,n); //for matrix multiplcaition
	  std::copy(vr, vr+n*n, vr_mat.data().begin());  //copy up right eigenvectors
	  int ipv[n];
	  zgetrf_(&n,&n,vr,&n,ipv,&info);   //take inverse in two stages (fast)...
	  zgetri_(&n,vr,&n,ipv,work,&lwork,&info);
	  rcmat_t vri_mat(n,n); //for matrix multipication
	  std::copy(vr, vr+n*n, vri_mat.data().begin()); //copy up inverse of right eigenvectors

	  unsigned int happyFaceIndex = 0;

	  for (unsigned int i = 0; i < n; ++i) //exp(val*par) . vr
	    for (unsigned int j = 0; j < n; ++j)
	      M(i,j)=vri_mat(i, happyFaceIndex) * vr_mat(happyFaceIndex, j);
	  //tempDiagMat(j, k) = exp(par[i] * val(j))*vr_mat(j,k);
	  //MG[i] = ublas::prod(vri_mat, tempDiagMat); //vr^-1.exp(val*par) .vr


	}
	
		  
	std::vector<std::complex<double> > theta(3);

	std::complex<double> Intensity;
	for (unsigned int i = 0; i < 3; ++i) {
		Intensity += M(i, 0) * IntA;
		Intensity += M(i, 1) * IntB;
		Intensity += M(i, 2) * IntC;
	}

	return abs(Intensity);

	return 0.0;
}
