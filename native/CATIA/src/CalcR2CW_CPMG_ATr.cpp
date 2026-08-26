#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;
#include <MatrixExponentialLapack.hpp>

void Catia::CalcR2CW_CPMG_ATr(Dataset& dset, int GlobalAtom) {
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
			std::cerr << " one must have ncyc>0 \n";
			std::cerr << " Function: CalcR2CW_CPMG_ATr();\n";
			Abort(1);
		}
	}
	// We need the following matrices:
	// ProtonInversion, Xpulsing, -Xpulsing,
	// Ypulsing(protonDEC),
	// FreePrecessing for (time_T2/(4*ncyc),time_equil,
	///////////////////////////////////////////////////
	std::vector<double> taucp;
	double time_T2 = atof(dset.vpar("time_t2")[1].c_str());
	double pwx_cp = atof(dset.vpar("pwx_cp")[1].c_str());
	for (unsigned int i = 0; i < dset.ncyc[LocalAtom].size(); i++) {
		taucp.push_back(1 / (4 * dset.ncyc[LocalAtom][i]) - pwx_cp);
	}
	taucp.push_back(atof(dset.vpar("time_equil")[1].c_str()));

	ublas::matrix<double> p90x, p90mx, p90y, Iinv;
	ublas::matrix<double>* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<double>[taucp.size()];
	//
	CalcMatrix(Iinv, dset, GlobalAtom, "Iinv", 0.);
	CalcMatrix(p90x, dset, GlobalAtom, "90x", pwx_cp); //90(x) pulse on S spin
	CalcMatrix(p90mx, dset, GlobalAtom, "90mx", pwx_cp); //90(-x) pulse on S spin
	CalcMatrix(p90y, dset, GlobalAtom, "90y", pwx_cp); // ..
	CalcMatrix(Gfree, dset, GlobalAtom, "free", taucp); // Free precession for all delays in taucp
	//
	//Equilibrium vector
	//
	ublas::vector<double> Sigma(Gfree[0].size1());
	ublas::vector<double> SigmaInit(Gfree[0].size1());
	//
	// Do we have the initial conditions ?
	if (dset._haveInitIntensity) {
		if (pow(pow(dset.initialIntensities[LocalAtom][0], 2.) + pow(dset.initialIntensities[LocalAtom][1], 2.), 0.5) < 1E-15) {
			CalcCoherence("2IzSy", Sigma, dset);
		} else {
			//
			double TrInt = dset.initialIntensities[LocalAtom][0];
			double ATrInt = dset.initialIntensities[LocalAtom][1];
			TrInt = TrInt / (ATrInt + TrInt);
			ATrInt = 1 - TrInt;
			ublas::vector<double> SigmaTr(Gfree[0].size1());
			ublas::vector<double> SigmaATr(Gfree[0].size1());
			CalcCoherence("Tr(y)", SigmaTr, dset);
			CalcCoherence("ATr(y)", SigmaATr, dset);
			for (unsigned int i = 0; i < Gfree[0].size1(); i++) {
				Sigma[i] = TrInt * SigmaTr[i] + ATrInt * SigmaATr[i];
			}
		}
	} else {
		CalcCoherence("2IzSy", Sigma, dset);
	}
	//
	// Store the initial sigma;
	for (unsigned int i = 0; i < Gfree[0].size1(); i++) {
		SigmaInit[i] = Sigma[i];
	}
	//
	propagate(p90y, Sigma);
	propagate(p90x, Sigma);
	propagate(p90x, Sigma);
	propagate(p90y, Sigma);
	propagate(p90x, Sigma);
	propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma);
	//
	// Measure the trosy peak intensity
	double InitMag = MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), (-Detect("Sz", 0, dset, Sigma)
			+ Detect("2IzSz", 0, dset, Sigma)) / 2., (-Detect("Sz", 1, dset, Sigma) + Detect("2IzSz", 1, dset, Sigma)) / 2.);
	//
	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {
		int ncyc = (int) floor(time_T2 * dset.ncyc[LocalAtom][ncycC] + 0.5);
		//
		double Intensity[2][2];
		//
		for (unsigned int p = 0; p < 2; p++) {
			//
			// Reset the intensity.
			for (unsigned int s = 0; s < Gfree[0].size1(); s++) {
				Sigma[s] = SigmaInit[s];
			}
			// convert Trosy <-> AntiTrosy;
			propagate(Iinv, Sigma);
			//
			//CPMG Train
			for (unsigned int ns = 0; ns < ncyc; ns++) {
				propagate(Gfree[ncycC], Sigma);
				propagate(p90y, Sigma);
				propagate(p90y, Sigma);
				propagate(Gfree[ncycC], Sigma);
				//Maybe 2pwn/pi
			}
			// HVK Element - START
			if (p == 0) {
				propagate(p90x, Sigma);
				propagate(p90x, Sigma);
			} else {
				propagate(p90mx, Sigma);
				propagate(p90mx, Sigma);
			}
			// HVK Element - END
			for (unsigned int ns = 0; ns < ncyc; ns++) {
				propagate(Gfree[ncycC], Sigma);
				propagate(p90y, Sigma);
				propagate(p90y, Sigma);
				propagate(Gfree[ncycC], Sigma);
			}
			propagate(p90x, Sigma);
			// convert Trosy <-> AntiTrosy;
			propagate(Iinv, Sigma);

			propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
			//
			// store Trosy magnetization
			Intensity[p][0] = (Detect("2IzSz", 0, dset, Sigma) - Detect("Sz", 0, dset, Sigma)) / 2.;
			Intensity[p][1] = (Detect("2IzSz", 1, dset, Sigma) - Detect("Sz", 1, dset, Sigma)) / 2.;
		}
		//
		dset.R2_calc[LocalAtom][ncycC] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),
				(Intensity[0][0] + Intensity[1][0]) / 2., (Intensity[0][1] + Intensity[1][1]) / 2.) / InitMag) / time_T2;
		//    std::cerr<<dset.ncyc[LocalAtom][ncycC]<<"\t"<<dset.R2_calc[LocalAtom][ncycC]
		//	     <<"\t"<<dset.R2_exp[LocalAtom][ncycC]<<std::endl;

	}
	delete[] Gfree;
	//
	return;
}
