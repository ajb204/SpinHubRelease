/*
 This function returns the evolution Matrix, thus, if M is the relative
 components of the coherences in the vector M and

 dM/dt = -G M, then, M(t)=exp(-t*G)*M(0)
 and the Matrix returned is exp(-t*G)

 written by D.F Hansen

 */

#include <complex>
//#include <boost/math/special_functions/sign.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>



//#include <MatrixExponential.hpp>     		//G = ublas::expm_pad(freeEvolMat); 
#include <MatrixExponentialLapack.hpp>             //expm_eig   expm_padeExpo

#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>

namespace ublas = boost::numeric::ublas;
//using ublas::prod;

typedef double real_t;
typedef std::complex<real_t> cmplx_t;
typedef ublas::vector<cmplx_t> cvec_t;
typedef ublas::matrix<real_t, ublas::row_major> rrmat_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> rcmat_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> ccmat_t;

void Catia::CalcMatrix(ublas::matrix<double>& G, Dataset& dset, int Atom, std::string opt, double par) {
  //std::cout << opt << std::endl;
	if (opt == "Iinv") {

		if (dset.vpar("basis")[1] == "tratr_13" || dset.vpar("basis")[1] == "n_nh_13") {
			G = ublas::zero_matrix<double>(13, 13);
			for (unsigned int i = 0; i < 3; i++) {
				G(1 + i, 4 + i) = 1.0;
				G(4 + i, 1 + i) = 1.0;
				G(1 + i + 6, 4 + i + 6) = 1.0;
				G(4 + i + 6, 1 + i + 6) = 1.0;
			}
			G(0, 0) = 1.0;

		} else if (dset.vpar("basis")[1] == "iphaph_13" || dset.vpar("basis")[1] == "iphaph_13_deltar2" || dset.vpar("basis")[1] == "iphaph_13_dr2") {
			G = ublas::zero_matrix<double>(13, 13);
			for (unsigned int i = 0; i < 3; i++) {
				G(1 + i, 1 + i) = 1.0;
				G(1 + 6 + i, 1 + 6 + i) = 1.0;
				G(4 + i, 4 + i) = -1.0;
				G(4 + i + 6, 4 + i + 6) = -1.0;
			}
			G(0, 0) = 1.0;

		} else if (dset.vpar("basis")[1] == "iph_7" || dset.vpar("basis")[1] == "n_7") {
			G = ublas::identity_matrix<double>(7);

		} else if (dset.vpar("basis")[1] == "3st_iph_10") {
			G = ublas::identity_matrix<double>(10);

		} else if (dset.vpar("basis")[1] == "tr_7") {
			// This is not really rigid - but what else can we do.
			// the 180(H) actually converts all trosy to antitrosy,
			// and the whole basis is gone.
			G = ublas::identity_matrix<double>(7);

		} else if (dset.vpar("basis")[1] == "aph_7") {
			G = ublas::zero_matrix<double>(7, 7);
			G(0, 0) = 1.0;
			for (unsigned int i = 1; i < 7; i++) {
				G(i, i) = -1.0;
			}

		} else if (dset.vpar("basis")[1] == "c_ch3_25") {
			G = ublas::zero_matrix<double>(25, 25);
			G(0, 0) = 1.0;
			for (unsigned int i = 0; i < 3; i++) {
				for (unsigned int j = 0; j < 2; j++) {
					G(1 + i + 12 * j, 10 + i + 12 * j) = 1.0; // A <-> D
					G(10 + i + 12 * j, 1 + i + 12 * j) = 1.0; // D <-> A
					//
					G(4 + i + 12 * j, 7 + i + 12 * j) = 1.0; // B <-> C
					G(7 + i + 12 * j, 4 + i + 12 * j) = 1.0; // C <-> B
				}
			}

		} else {
			std::cerr << " Basis:" << dset.vpar("basis")[1] << " is undefined\n";
			std::cerr << " Function CalcMatrix();\n";
			std::cerr << std::endl;
			Abort(1);
		}
			    
	} else if (opt == "free") { //This is valid for all bases.
		rrmat_t freeEvolMat;
		//std::cout << "tootdles" << std::endl;
		//FreePrecess(freeEvolMat, dset, Atom);
		//freeEvolMat *= -par;
		//G = ublas::expm_pad(freeEvolMat);
		//expm_eig(freeEvolMat,G,par);
		expm_padeExpo(freeEvolMat,G,par);


	} else if (opt == "spinlock_x") { //This is valid for all bases.
	  	rrmat_t freeEvolMat;
		FreePrecess(freeEvolMat, dset, Atom);
		//freeEvolMat *= -par;
		//G = ublas::expm_pad(freeEvolMat);
		//expm_eig(freeEvolMat,G,par);
		expm_padeExpo(freeEvolMat,G,par);
	  
	} else if (opt == "90x" || opt == "90mx") {
	  	rrmat_t Gf;
		FreePrecess(Gf, dset, Atom);

		// calculate the pulsing field.
		double w1 = 2 * DFH_PI / (4 * par);
		if (opt == "90mx") {
			w1 = -w1;
		}

		rrmat_t P = ublas::zero_matrix<double>(Gf.size1(), Gf.size2());

		// basis depend
		if (dset.vpar("basis")[1] == "tratr_13" || dset.vpar("basis")[1] == "n_nh_13") {
			for (unsigned int k = 0; k < 2; k++) {
				//trosy
				P(2 + 6 * k, 3 + 6 * k) = w1;
				P(3 + 6 * k, 2 + 6 * k) = -w1;
				//anti trosy
				P(2 + 3 + 6 * k, 3 + 3 + 6 * k) = w1;
				P(3 + 3 + 6 * k, 2 + 3 + 6 * k) = -w1;
			}

		} else if (dset.vpar("basis")[1] == "iphaph_13" || dset.vpar("basis")[1] == "iphaph_13_deltar2" || dset.vpar("basis")[1] == "iphaph_13_dr2") {
			for (unsigned int k = 0; k < 2; k++) {
				//inphase
				P(2 + 6 * k, 3 + 6 * k) = w1;
				P(3 + 6 * k, 2 + 6 * k) = -w1;
				//antiphase
				P(2 + 3 + 6 * k, 3 + 3 + 6 * k) = w1;
				P(3 + 3 + 6 * k, 2 + 3 + 6 * k) = -w1;
			}

		} else if (dset.vpar("basis")[1] == "iph_7" || dset.vpar("basis")[1] == "n_7") {
			for (unsigned int k = 0; k < 2; k++) {
				//inphase
				P(2 + 3 * k, 3 + 3 * k) = w1;
				P(3 + 3 * k, 2 + 3 * k) = -w1;
			}

		} else if (dset.vpar("basis")[1] == "3st_iph_10") {
			for (unsigned int k = 0; k < 3; k++) {
				//inphase
				P(2 + 3 * k, 3 + 3 * k) = w1;
				P(3 + 3 * k, 2 + 3 * k) = -w1;
			}

		} else if (dset.vpar("basis")[1] == "aph_7") {
			for (unsigned int k = 0; k < 2; k++) {
				//inphase
				P(2 + 3 * k, 3 + 3 * k) = w1;
				P(3 + 3 * k, 2 + 3 * k) = -w1;
			}

		} else if (dset.vpar("basis")[1] == "tr_7") {
			for (unsigned int k = 0; k < 2; k++) {
				//inphase
				P(2 + 3 * k, 3 + 3 * k) = w1;
				P(3 + 3 * k, 2 + 3 * k) = -w1;
			}

		} else if (dset.vpar("basis")[1] == "c_ch3_25") {
			for (unsigned int k = 0; k < 2; k++) {
				//A
				P(2 + 12 * k, 3 + 12 * k) = w1;
				P(3 + 12 * k, 2 + 12 * k) = -w1;
				//B
				P(2 + 3 + 12 * k, 3 + 3 + 12 * k) = w1;
				P(3 + 3 + 12 * k, 2 + 3 + 12 * k) = -w1;
				//C
				P(2 + 6 + 12 * k, 3 + 6 + 12 * k) = w1;
				P(3 + 6 + 12 * k, 2 + 6 + 12 * k) = -w1;
				//D
				P(2 + 9 + 12 * k, 3 + 9 + 12 * k) = w1;
				P(3 + 9 + 12 * k, 2 + 9 + 12 * k) = -w1;
			}

		} else {
			std::cerr << " Basis:" << dset.vpar("basis")[1] << " is undefined\n";
			std::cerr << " Function CalcMatrix();\n";
			std::cerr << std::endl;
			Abort(1);
		}

		Gf += P;
		//Gf *= -par;
		//G = ublas::expm_pad(Gf);
		//expm_eig(Gf,G,par);
		expm_padeExpo(Gf,G,par);
		
		
	} else if (opt == "90y" || opt == "90my") {
	  	rrmat_t freeEvolMat;
		FreePrecess(freeEvolMat, dset, Atom);

		// Calculate the pulsing field.
		//double w1 = boost::math::sign(dset._gamma) * 2 * DFH_PI / (4 * par);
		double w1 = dset._gamma/abs(dset._gamma) * 2 * DFH_PI / (4 * par);

		//double w1=2*DFH_PI/(4*par);
		if (opt == "90my") {
			w1 = -w1;
		}

		rrmat_t P = ublas::zero_matrix<double>(freeEvolMat.size1(), freeEvolMat.size2());

		if (dset.vpar("basis")[1] == "tratr_13" || dset.vpar("basis")[1] == "n_nh_13") {
			for (unsigned int k = 0; k < 2; k++) {
				//trosy
				P(1 + 6 * k, 3 + 6 * k) = -w1;
				P(3 + 6 * k, 1 + 6 * k) = w1;
				//anti trosy
				P(1 + 3 + 6 * k, 3 + 3 + 6 * k) = -w1;
				P(3 + 3 + 6 * k, 1 + 3 + 6 * k) = w1;
			}

		} else if (dset.vpar("basis")[1] == "iphaph_13" || dset.vpar("basis")[1] == "iphaph_13_deltar2"|| dset.vpar("basis")[1] == "iphaph_13_dr2") {
			for (unsigned int k = 0; k < 2; k++) {
				//inphase
				P(1 + 6 * k, 3 + 6 * k) = -w1;
				P(3 + 6 * k, 1 + 6 * k) = w1;
				//antiphase
				P(1 + 3 + 6 * k, 3 + 3 + 6 * k) = -w1;
				P(3 + 3 + 6 * k, 1 + 3 + 6 * k) = w1;
			}

		} else if (dset.vpar("basis")[1] == "iph_7" || dset.vpar("basis")[1] == "n_7") {
			for (unsigned int k = 0; k < 2; k++) {
				//inphase
				P(1 + 3 * k, 3 + 3 * k) = -w1;
				P(3 + 3 * k, 1 + 3 * k) = w1;
			}

		} else if (dset.vpar("basis")[1] == "3st_iph_10") {
			for (unsigned int k = 0; k < 3; k++) {
				//inphase
				P(1 + 3 * k, 3 + 3 * k) = -w1;
				P(3 + 3 * k, 1 + 3 * k) = w1;
			}

		} else if (dset.vpar("basis")[1] == "aph_7") {
			for (unsigned int k = 0; k < 2; k++) {
				//inphase
				P(1 + 3 * k, 3 + 3 * k) = -w1;
				P(3 + 3 * k, 1 + 3 * k) = w1;
			}

		} else if (dset.vpar("basis")[1] == "tr_7") {
			for (unsigned int k = 0; k < 2; k++) {
				//inphase
				P(1 + 3 * k, 3 + 3 * k) = -w1;
				P(3 + 3 * k, 1 + 3 * k) = w1;
			}

		} else if (dset.vpar("basis")[1] == "c_ch3_25") {
			for (unsigned int k = 0; k < 2; k++) {
				//A
				P(1 + 12 * k, 3 + 12 * k) = -w1;
				P(3 + 12 * k, 1 + 12 * k) = w1;
				//B
				P(1 + 3 + 12 * k, 3 + 3 + 12 * k) = -w1;
				P(3 + 3 + 12 * k, 1 + 3 + 12 * k) = w1;
				//C
				P(1 + 6 + 12 * k, 3 + 6 + 12 * k) = -w1;
				P(3 + 6 + 12 * k, 1 + 6 + 12 * k) = w1;
				//D
				P(1 + 9 + 12 * k, 3 + 9 + 12 * k) = -w1;
				P(3 + 9 + 12 * k, 1 + 9 + 12 * k) = w1;
			}

		} else {
			std::cerr << " Basis:" << dset.vpar("basis")[1] << " is undefined\n";
			std::cerr << " Function CalcMatrix();\n";
			std::cerr << std::endl;
			Abort(1);
		}
		//
		freeEvolMat += P;
		//freeEvolMat *= -par;
		//G = ublas::expm_pad(freeEvolMat);
		//std::cout << "tddsaf" << std::endl;
		//expm_eig(freeEvolMat,G,par);
		expm_padeExpo(freeEvolMat,G,par);
		
	} else {
		std::cerr << " The option " << opt << " is unknown\n";
		std::cerr << " Function CalcMatrix();\n";
		std::cerr << std::endl;
		Abort(1);
	}
}



void Catia::CalcMatrix(rrmat_t* G, Dataset& dset, int Atom, std::string opt, std::vector<double> par) {
  //std::cout << opt << std::endl;
  if (opt == "free") {

    rrmat_t freeEvolMat;
		FreePrecess(freeEvolMat, dset, Atom);
		expm_eig(freeEvolMat,G,par); //update propagators

		//freeEvolMat *= -par;
		//G = ublas::expm_pad(freeEvolMat); 
		
		//std::cout << "tiddles" << std::endl;
		//NOTE: not sure why but padeexpot fails here//
		//AJB 16th Jan 2023
		//expm_padeExpo(freeEvolMat,G,par);

		/*
		for(int p=0;p<par.size();++p)
		  {
		    rrmat_t nom;
		    nom=freeEvolMat;
		    for(int i=0;i<freeEvolMat.size1();++i)
		      {
			for(int j=0;j<freeEvolMat.size1();++j)
			  {
			    nom(i,j)*=-par[p];
			  }
		      }
		    G[p]=real(ublas::expm_pad(nom));
		  }*/
		
		/*
		const unsigned int n = freeEvolMat.size1();
		cvec_t eigValVect(n);
		ccmat_t lEigVectMat(1, 1);
		ccmat_t rEigVectMat(n, n);
		ccmat_t freeEvolMatColMaj(freeEvolMat);
		//lapack::geev('N', 'V', freeEvolMatColMaj, eigValVect, lEigVectMat, rEigVectMat, lapack::optimal_workspace());
		ccmat_t rEigVectMatInv(rEigVectMat);
		ublas::vector<int> ipiv(n); // pivot vector
		//lapack::getrf(rEigVectMatInv, ipiv); // no lu_factor() alias for getrf() available
		//lapack::getri(rEigVectMatInv, ipiv); // no lu_invert() alias for getrf() available
		ublas::zero_matrix<cmplx_t, ublas::column_major> zeroMat(n);
		ccmat_t tempDiagMat(n, n);
		// we should allocate space on G 
		for (unsigned int i = 0; i < par.size(); i++) {
			tempDiagMat = zeroMat;
			for (unsigned int j = 0; j < n; ++j) {
				tempDiagMat(j, j) = exp(-par[i] * eigValVect(j));
			}
			G[i] = real(prod(rEigVectMat, ccmat_t(prod(tempDiagMat, rEigVectMatInv))));
		}
		*/
		
	} 
	else if (opt == "w1_free") {
	  	rrmat_t freeEvolMat;
		FreePrecess(freeEvolMat, dset, Atom);
		//ADD ON B1
		//expm_eig(freeEvolMat,G,par);
		expm_padeExpo(freeEvolMat,G,par);
		/*
		const unsigned int n = freeEvolMat.size1();
		cvec_t eigValVect(n);
		ccmat_t lEigVectMat(1, 1);
		ccmat_t rEigVectMat(n, n);
		ccmat_t freeEvolMatColMaj(freeEvolMat);
		lapack::geev('N', 'V', freeEvolMatColMaj, eigValVect, lEigVectMat, rEigVectMat, lapack::optimal_workspace());
		ccmat_t rEigVectMatInv(rEigVectMat);
		ublas::vector<int> ipiv(n); // pivot vector
		lapack::getrf(rEigVectMatInv, ipiv); // no lu_factor() alias for getrf() available
		lapack::getri(rEigVectMatInv, ipiv); // no lu_invert() alias for getrf() available
		ublas::zero_matrix<cmplx_t, ublas::column_major> zeroMat(n);
		ccmat_t tempDiagMat(n, n);
		// we should allocate space on G 
		for (unsigned int i = 0; i < par.size(); i++) {
			tempDiagMat = zeroMat;
			for (unsigned int j = 0; j < n; ++j) {
				tempDiagMat(j, j) = exp(-par[i] * eigValVect(j));
			}
			G[i] = real(prod(rEigVectMat, ccmat_t(prod(tempDiagMat, rEigVectMatInv))));
		}
		*/

	}
	else {
		std::cerr << " Can only use array of free precession, i.e., you cannot array " << opt << " !\n";
		std::cerr << " Function .CalcMatrix();\n";
		std::cerr << " Dataset " << dset._id << std::endl;
		Abort(1);
	}

}

void Catia::CalcMatrix(rcmat_t* G, Dataset& dset, int Atom, std::string opt, std::vector<double> par) {
  //std::cout << opt << std::endl;
	if (opt == "free") {
	  	rcmat_t freeEvolMat;
		//std::cout << "tddsaf" << std::endl;
		FreePrecess(freeEvolMat, dset, Atom);
		//expm_eig(freeEvolMat,G,par);
		expm_padeExpo(freeEvolMat,G,par);
		
		/*
		const unsigned int n = freeEvolMat.size1();
		cvec_t eigValVect(n);
		ccmat_t lEigVectMat(1, 1);
		ccmat_t rEigVectMat(n, n);
		ccmat_t freeEvolMatColMaj(freeEvolMat);
		//lapack::geev('N', 'V', freeEvolMatColMaj, eigValVect, lEigVectMat, rEigVectMat, lapack::optimal_workspace());
		ccmat_t rEigVectMatInv(rEigVectMat);
		ublas::vector<int> ipiv(n); // pivot vector
		//lapack::getrf(rEigVectMatInv, ipiv); // no lu_factor() alias for getrf() available
		//lapack::getri(rEigVectMatInv, ipiv); // no lu_invert() alias for getrf() available
		ublas::zero_matrix<cmplx_t, ublas::column_major> zeroMat(n);
		ccmat_t tempDiagMat(n, n);
		//we should allocate space on G 
		for (unsigned int i = 0; i < par.size(); i++) {
			tempDiagMat = zeroMat;
			for (unsigned int j = 0; j < n; ++j) {
				tempDiagMat(j, j) = exp(-par[i] * eigValVect(j));
			}
			G[i] = prod(rEigVectMat, ccmat_t(prod(tempDiagMat, rEigVectMatInv)));
		}
		*/

	} else {
		std::cerr << " Can only use array of free precession, i.e., you cannot array " << opt << " !\n";
		std::cerr << " Function .CalcMatrix();\n";
		std::cerr << " Dataset " << dset._id << std::endl;
		Abort(1);
	}

}
