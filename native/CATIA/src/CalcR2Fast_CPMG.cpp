/*
 * CalcR2Fast.cpp
 *
 *  Created on: May 3, 2010
 *      Author: guillaume
 */


#include <complex>
//#include <boost/math/special_functions/sign.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <lapack.h>

//#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>

namespace ublas = boost::numeric::ublas;
//#include <MatrixExponentialLapack.hpp>
using ublas::prod;

typedef double real_t;
typedef std::complex<real_t> cmplx_t;
typedef ublas::vector<real_t> rvec_t;
typedef ublas::vector<cmplx_t> cvec_t;
typedef ublas::matrix<real_t> rmat_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> cmat_t;




size_t GoodEigenValueIndex(const rvec_t &x);

void Catia::CalcR2Fast_CPMG(Dataset& dset, int GlobalAtom) {
	// GlobalAtom is the AtomNumber of the Global set of Atoms;
  
    int LocalAtom = -1;
	for (unsigned int i = 0; i < dset._localToGlobalAtomIndex.size(); i++) {
		if (GlobalAtom == dset._localToGlobalAtomIndex[i]) {
			LocalAtom = i;
			break;
		}
	}
	if (LocalAtom == -1) {
		std::cerr << " Atom :" << AtomNumber2AtomName(GlobalAtom);
		std::cerr << " is not found in the Dataset named :" << dset._id;
		std::cerr << std::endl;
		Abort(1);
	}
	//
	//Check that ncyc is not 0:
	for (unsigned int i = 0; i < dset.ncyc[LocalAtom].size(); i++) {
		if (dset.ncyc[LocalAtom][i] < 1) {
			std::cerr << " This function returns R2eff - i.e.,\n";
			std::cerr << " one must have ncyc>0 \n";
			std::cerr << " Function: CalcR2CW_CPMG();\n";
			Abort(1);
		}
	}

	
	///////////////////////////////////////////////////
	std::vector<double> taucp;
	for (unsigned int i = 0; i < dset.ncyc[LocalAtom].size(); i++) {
		taucp.push_back(2.0 / (4.0 * dset.ncyc[LocalAtom][i]));
	}
	
	ublas::matrix<std::complex<double> >* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<std::complex<double> >[taucp.size()];


	
	CalcMatrix(Gfree, dset, GlobalAtom, "free", taucp);

	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {

		const int n = Gfree[ncycC].size1();

		std::complex<double> eigValVect[n];
		std::complex<double> *vl;
		std::complex<double> vr[n*n];
		std::complex<double> work[4*n]; //workspace for moving columns around
		int lwork=4*n;
		double rwork=3*n-2; //max value from documentation
		int info;

		//cmat_t rEigVectMat(n, n);
		cmat_t freeEvolMatColMaj(prod(conj(Gfree[ncycC]), Gfree[ncycC]));
		
		//lapack::geev('N', 'V', freeEvolMatColMaj, eigValVect, lEigVectMat, rEigVectMat, lapack::optimal_workspace());
		//LAPACKE_cgeev(LAPACK_COL_MAJOR,'N','V', n,&(aa.data()[0]), n,eigValVect, vl,n,vr,n); //EV
		zgeev_("N","V", &n,&(freeEvolMatColMaj.data()[0]), &n,eigValVect, vl,&n,vr,&n,work,&lwork,&rwork,&info); //EV
		
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.0f", dset._sfrq);
		std::string fieldMarker(line);
		ClearBuf(line, sizeof(line));
		//
		std::string temperatureMarker("");
		if (_multipleTemperatures) {
			sprintf(line, "_%.1f", dset._temperature);
			temperatureMarker = line;
			ClearBuf(line, sizeof(line));
		}
		std::string marker(fieldMarker + temperatureMarker);

		const double R0 = ResolveParam(LocalParam[GlobalAtom], "R0" + marker);
//		std::cerr << R0 << " | ";
		double R2 = 1.0e+16;
		for (size_t i = 0; i < n; i++) {
			double tmp = -1.0 / (2.0 * taucp[ncycC]) * real(log(eigValVect[i]));
//			std::cerr << tmp << " ";
			if (tmp < R2 && tmp > R0 + 1.0e-10) {
				R2 = tmp;
			}
		}
//		std::cerr << "| " << R2 << std::endl;

		dset.R2_calc[LocalAtom][ncycC] = R2;

	}

	delete[] Gfree;

}


// To get the right eigenvalue, the best is probably to look at the intensities.
	//		cmat_t rEigVectMatInv(rEigVectMat);
	//		ublas::vector<int> ipiv(n); // pivot vector
	//		lapack::getrf(rEigVectMatInv, ipiv); // no lu_factor() alias for getrf() available
	//		lapack::getri(rEigVectMatInv, ipiv); // no lu_invert() alias for getrf() available
	//
	//		const double pb = Pb_3st(dset);
	//		const double pc = Pc_3st(dset);
	//		const double pa = 1.0 - pb - pc;
	//		rvec_t p(n);
	//		p(0) = pa; p(1) = pb; p(2) = pc;
	//		cvec_t coeff = ublas::prod(rEigVectMatInv, p);
	//
	//		size_t goodIndex = 0;
	//		double bigger = 0.0;
	//		for(size_t i = 0; i < n; ++i) {
	//			double tmp = real(sum(ublas::matrix_column<cmat_t>(rEigVectMat, i)) * coeff(i));
	//			if (tmp > bigger) {
	//				bigger = tmp;
	//				goodIndex = i;
	//			}
	//		}

