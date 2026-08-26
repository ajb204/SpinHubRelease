#include <Catia.h>
#include <Dataset.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

/*
 This is ready to be read by the "Numerical Recipies" in C
 minimization routines

 D.F. Hansen, September 11 2007.

 Modified by Flemming on October 12 2007
 - now ready for the Flemming-version of the lm fitting function
 */

void Catia::CalcR2drv(VecDoub x, VecDoub& y, MatDoub& dyda, const std::vector<double*> a, std::vector<int*> ai, const std::vector<int> AtomNumber,
		const std::vector<std::string> ParamName) {
	// x: x-coordinate:
	//    This is now refering to combination of
	//    int()Dataset,int()atom,float(ncyc),
	//    Thus, we keep x as an integer and use the map:
	//    (int) x-> int[3]
	//
	// a: array of pointers to the parameters
	// y: R2_calc
	// dyda: ..
	// ma: number of parameters
	//
	//////////////////////////////////////////////////////////////////
	//
	int ndata(x.size());
	int ma(a.size());
	//
	std::map<std::string, double>::iterator itDS; //iterator:double<-string
	//
	// sort out which (Dataset,atom) we have to calculate
	//
	std::vector<std::vector<int> > R22Calc; // (Dataset,atom);
	std::vector<int> dset_atom(3);
	for (unsigned int i = 0; i < ndata; i++) {
		dset_atom = X2dset_atom[(int) floor(x[i] + 0.1)];
		if (i == 0) {
			R22Calc.push_back(dset_atom);
			continue;
		}
		//
		if ((dset_atom[0] != R22Calc[R22Calc.size() - 1][0]) || (dset_atom[1] != R22Calc[R22Calc.size() - 1][1])) {
			R22Calc.push_back(dset_atom);
		}
	}
	//
	// Calc R2
	//
	for (unsigned int i = 0; i < R22Calc.size(); i++) {
		CalcR2(Datasets[R22Calc[i][0]], R22Calc[i][1]);
	}
	//
	//  std::cerr<<" Calculated R2!"<<std::endl;
	for (unsigned int i = 0; i < ndata; i++) {
		dset_atom = X2dset_atom[(int) floor(x[i] + 0.1)];
		y[i] = Datasets[dset_atom[0]].R2_calc[Datasets[dset_atom[0]]._atomNameToLocalAtomNumber[AtomNumber2AtomName(dset_atom[1])]][dset_atom[2]];
		//std::cerr<<AtomNumber2AtomName(dset_atom[1])<<" "<<Datasets[dset_atom[0]].id<<" "<<y[i]<<std::endl;
	}
	//std::cerr<<" R2 uploaded to heap "<<std::endl;
	//
	// Calculate the derivative.
	//
	if (CalcDeriv != 0) {
		//
		std::map<std::string, double>::iterator itDS; //iterator:double<-string
		double* yfita;
		double* da;
		yfita = new double[ndata];
		da = new double[ma];
		std::vector<int> dset_atom(3);
		//
		// store the original parameters;
		for (unsigned int k = 0; k < ma; k++) {
			da[k] = *(a[k]);
		}
		//
		//
		for (unsigned int i = 0; i < ma; i++) {
			//
			// Do we know a priori that this derivative is 0 (zero)
			bool global = ((AtomNumber[i] == -1) ? true : false);
			int free = (*(ai[i]));
			signed int Atom = AtomNumber[i];
			std::string name = ParamName[i];
			double SD; // square deviation
			//
			//
			// we know have the following information achieved:
			//
			// free = 1 or 0 or -2 (Force fixed )
			// name (ParamName_AtomName)
			// Atom GlobalAtomNumber
			// global true or false;
			//
			if (free == 1) {
				//	cblas_scopy(ma, a, 1, da, 1);
				for (unsigned int k = 0; k < ma; k++) {
					*(a[k]) = da[k];
				}
				*(a[i]) = (fabs(da[i] * (1. + 1E-4)) > 1E-5 ? da[i] * (1. + 1E-4) : 1E-5);
				//
				// Calc R2 (only when needed though )
				//
				int ndataC = 0;
				for (unsigned int j = 0; j < R22Calc.size(); j++) {
					//
					int AtomNumberInDataSet = Datasets[R22Calc[j][0]]._atomNameToLocalAtomNumber[AtomNumber2AtomName(R22Calc[j][1])];
					//
					// Check that GlobalAtom (R22Calc[j][1]) depend of param i
					if (global || Atom == R22Calc[j][1]) {
						CalcR2(Datasets[R22Calc[j][0]], R22Calc[j][1]);
						//
						// Now calculate dyda and put on array dyda[]
						//
						for (int jj = 0; jj < Datasets[R22Calc[j][0]].ncyc[AtomNumberInDataSet].size(); jj++) {
							dyda[ndataC][i] = (Datasets[R22Calc[j][0]].R2_calc[AtomNumberInDataSet][jj] - y[ndataC]) / (*(a[i]) - da[i]);
							ndataC++;
						}
					} else {
						for (int jj = 0; jj < Datasets[R22Calc[j][0]].ncyc[AtomNumberInDataSet].size(); jj++) {
							dyda[ndataC][i] = 0.;
							ndataC++;
						}
					}
				}
				//START - check if we have some dummy parameter
				SD = 0.;
				for (unsigned int j = 0; j < ndata; j++) {
					SD += pow(dyda[j][i] / y[j], 2);
				}
				SD = sqrt(SD) / ndata;
				if (SD < _fixParamLimit) {
					if (global) {
						std::cerr << " The global parameter: " << name << " has changed status to fixed" << std::endl;
						*(ai[i]) = -2;
					} else {
						std::cerr << " The local parameter: " << name << " of atom " << AtomNumber2AtomName(Atom) << " has changed status to fixed" << std::endl;
						*(ai[i]) = -2;
					}
				}
				//END - check if we have a dummy parameter
			}
		}
		// Reset the parameters to the original value after calulating the derivatives
		for (unsigned int k = 0; k < ma; k++) {
			*(a[k]) = da[k];
		}
		// Clean up
		delete[] da;
		delete[] yfita;
	}
	//
	return;
}

