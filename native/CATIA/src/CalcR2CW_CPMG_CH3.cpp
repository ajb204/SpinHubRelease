/*
 This is the pulse sequence to measure CH3 cpmg dispersion of the four
 individual carbon lines of the CH3 group.

 based on HtoC_CH3_exchange_600_DC_lek.c

 Written by Flemming on September 24 2008
 */
#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;
#include <MatrixExponentialLapack.hpp>
using boost::to_lower_copy;

void Catia::CalcR2CW_CPMG_CH3(Dataset& dset, int GlobalAtom) {
	//
	const double kex = Kex(dset);
	const double pb = Pb(dset);
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
	//Check that nu(CPMG) is not 0 or less
	for (unsigned int i = 0; i < dset.ncyc[LocalAtom].size(); i++) {
		if (dset.ncyc[LocalAtom][i] < 1) {
			std::cerr << " This function returns R2eff - i.e.,\n";
			std::cerr << " one must have nu(CPMG)>0 \n";
			std::cerr << " Function: CalcR2CW_CPMG_CH3();\n";
			Abort(1);
		}
	}
	// We need the following matrices:
	// ProtonInversion, Xpulsing, -Xpulsing,
	// FreePrecessing for (time_T2/(4*ncyc),time_equil,
	///////////////////////////////////////////////////
	std::vector<double> taucp;
	double time_T2 = atof(dset.vpar("time_t2")[1].c_str());
	double pwx_cp = atof(dset.vpar("pwx_cp")[1].c_str());
	for (unsigned int i = 0; i < dset.ncyc[LocalAtom].size(); i++) {
		if (fabs(dset.ncyc[LocalAtom][i] * time_T2 - 0.50) < 1e-6) {
			taucp.push_back(time_T2 / 2.);
		} else {
			taucp.push_back(1 / (4 * dset.ncyc[LocalAtom][i]) - pwx_cp);
		}
	}
	taucp.push_back(atof(dset.vpar("time_equil")[1].c_str()));
	taucp.push_back(-2. * pwx_cp / DFH_PI);
	std::string select_flg = to_lower_copy(dset.vpar("select_flg")[1]);

	ublas::matrix<double> p90x, p90mx, p90y;
	ublas::matrix<double>* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<double>[taucp.size()];
	//
	CalcMatrix(p90x, dset, GlobalAtom, "90x", pwx_cp); //90(x) pulse on S spin
	CalcMatrix(p90mx, dset, GlobalAtom, "90mx", pwx_cp); //90(-x) pulse on S spin
	CalcMatrix(p90y, dset, GlobalAtom, "90y", pwx_cp); // ..
	CalcMatrix(Gfree, dset, GlobalAtom, "free", taucp); // Free precession for all delays in taucp
	//
	//Equilibrium vector
	//
	ublas::vector<double> Sigma(Gfree[0].size1());
	ublas::vector<double> SigmaInit(Gfree[0].size1());
	ublas::vector<double> Signal(Gfree[0].size1());
	//
	// Which line are we looking at:
	if (select_flg == "aaa") {
		CalcCoherence("CH3_A(z)", Sigma, dset);
	} else if (select_flg == "aab") {
		CalcCoherence("CH3_B(z)", Sigma, dset);
	} else if (select_flg == "abb") {
		CalcCoherence("CH3_C(z)", Sigma, dset);
	} else if (select_flg == "bbb") {
		CalcCoherence("CH3_D(z)", Sigma, dset);
	} else {
		std::cerr << " select_flg only takes the following values:\n";
		std::cerr << " aaa  : The |aaa><aaa| transistion\n";
		std::cerr << " aab  : The |aab><aab| + |aba><aba| + ... transistion\n";
		std::cerr << " abb  : The |abb><abb| transistion\n";
		std::cerr << " bbb  : |bbb><bbb|\n";
		std::cerr << " Function: CalcR2CW_CPMG_CH3();\n";
		std::cerr << " while evaluating Dataset " << dset._id << std::endl;
		Abort(1);
	}
	//
	// Store the initial sigma;
	for (unsigned int i = 0; i < Gfree[0].size1(); i++) {
		SigmaInit[i] = Sigma[i];
		Signal[i] = 0.;
	}
	//
	// ncyc=0 condition.
	for (unsigned int p = 0; p < 4; p++) {
		if (p == 0 || p == 1) {
			propagate(p90x, Sigma);
		} else {
			propagate(p90mx, Sigma);
		}
		if (p == 0 || p == 2) { //pi x
			propagate(p90x, Sigma);
			propagate(p90x, Sigma);
		} else {
			propagate(p90mx, Sigma);
			propagate(p90mx, Sigma);
		}
		propagate(p90x, Sigma);
		propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //time_equil
		//
		for (unsigned int j = 0; j < Gfree[0].size1(); j++) {
			if (p == 0 || p == 1) {
				Signal[j] += Sigma[j] / 4.;
			} else {
				Signal[j] -= Sigma[j] / 4.;
			}
			Sigma[j] = SigmaInit[j];
		}
	}
	double MajorInt, MinorInt;
	if (select_flg == "aaa") {
		MajorInt = Detect("CH3_A(z)", 0, dset, Signal);
		MinorInt = Detect("CH3_A(z)", 1, dset, Signal);
	} else if (select_flg == "aab") {
		MajorInt = Detect("CH3_B(z)", 0, dset, Signal);
		MinorInt = Detect("CH3_B(z)", 1, dset, Signal);
	} else if (select_flg == "abb") {
		MajorInt = Detect("CH3_C(z)", 0, dset, Signal);
		MinorInt = Detect("CH3_C(z)", 1, dset, Signal);
	} else if (select_flg == "bbb") {
		MajorInt = Detect("CH3_D(z)", 0, dset, Signal);
		MinorInt = Detect("CH3_D(z)", 1, dset, Signal);
	}
	//
	// Measure the intensity
	double InitMag = MajorPeakIntensity(pb, kex, (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), MajorInt, MinorInt);
	//
	//
	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {
		//
		int ncyc = (int) floor(time_T2 * dset.ncyc[LocalAtom][ncycC] + 0.1);
		//check if we have the ncyc==-1 condition.
		//std::cerr<< dset.ncyc[LocalAtom][ncycC]*time_T2 <<std::endl;
		if (fabs(dset.ncyc[LocalAtom][ncycC] * time_T2 - 0.50) < 1e-6) {
			ncyc = -1;
		}
		for (unsigned int s = 0; s < Gfree[0].size1(); s++) {
			Signal[s] = 0.;
		}
		//
		for (unsigned int p = 0; p < 4; p++) {
			//
			// Reset the intensity.
			for (unsigned int s = 0; s < Gfree[0].size1(); s++) {
				Sigma[s] = SigmaInit[s];
			}
			if (p == 0 || p == 1) {
				propagate(p90x, Sigma);
			} else {
				propagate(p90mx, Sigma);
			}
			if (ncyc == -1) {
				propagate(Gfree[ncycC], Sigma);
				propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
			} else {
				//
				//CPMG Train
				for (unsigned int ns = 0; ns < ncyc; ns++) {
					propagate(Gfree[ncycC], Sigma);
					if (ns == 0) {
						propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
					}
					propagate(p90y, Sigma);
					propagate(p90y, Sigma);
					propagate(Gfree[ncycC], Sigma);
				}
			}
			// 180 pulse to remove pulse inperfections
			if (p == 0 || p == 2) {
				propagate(p90x, Sigma);
				propagate(p90x, Sigma);
			} else {
				propagate(p90mx, Sigma);
				propagate(p90mx, Sigma);
			}
			if (ncyc == -1) {
				propagate(Gfree[ncycC], Sigma);
				propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
			} else {
				for (unsigned int ns = 0; ns < ncyc; ns++) {
					propagate(Gfree[ncycC], Sigma);
					propagate(p90y, Sigma);
					propagate(p90y, Sigma);
					propagate(Gfree[ncycC], Sigma);
					if (ns == ncyc - 1) {
						propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
					}
				}
			}
			propagate(p90x, Sigma);
			propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
			for (unsigned int j = 0; j < Gfree[0].size1(); j++) {
				if (p == 0 || p == 1) {
					Signal[j] += Sigma[j] / 4.;
				} else {
					Signal[j] -= Sigma[j] / 4.;
				}
				Sigma[j] = 0.;
			}
		}
		//
		//
		if (select_flg == "aaa") {
			MajorInt = Detect("CH3_A(z)", 0, dset, Signal);
			MinorInt = Detect("CH3_A(z)", 1, dset, Signal);
		} else if (select_flg == "aab") {
			MajorInt = Detect("CH3_B(z)", 0, dset, Signal);
			MinorInt = Detect("CH3_B(z)", 1, dset, Signal);
		} else if (select_flg == "abb") {
			MajorInt = Detect("CH3_C(z)", 0, dset, Signal);
			MinorInt = Detect("CH3_C(z)", 1, dset, Signal);
		} else if (select_flg == "bbb") {
			MajorInt = Detect("CH3_D(z)", 0, dset, Signal);
			MinorInt = Detect("CH3_D(z)", 1, dset, Signal);
		}
		//
		// Measure the intensity
		double Int = MajorPeakIntensity(pb, kex, (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), MajorInt, MinorInt);
		//std::cerr<<ncyc<<"\t"<<Int<<std::endl;
		dset.R2_calc[LocalAtom][ncycC] = -log(Int / InitMag) / time_T2;
	}
	delete[] Gfree;
}
