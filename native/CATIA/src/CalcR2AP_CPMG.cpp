#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>

#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;
#include <MatrixExponentialLapack.hpp>
/*
 This sequence is the seqfil="AP_CPMG"

 Essentially it is build up of
 2IzSy -> {CPMG}_y - 180_x  - {CPMG}_y -> 2IzSy

 */
void Catia::CalcR2AP_CPMG(Dataset& dset, int GlobalAtom) {
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
			std::cerr << " Function: CalcR2AP_CPMG();\n";
			Abort(1);
		}
	}
	// We need the following matrices:
	// Xpulsing,
	// Ypulsing,
	// FreePrecessing for (time_T2/(4*ncyc),
	///////////////////////////////////////////////////
	std::vector<double> taucp;
	double time_T2 = atof(dset.vpar("time_t2")[1].c_str());
	double pwx_cp = atof(dset.vpar("pwx_cp")[1].c_str());
	for (unsigned int i = 0; i < dset.ncyc[LocalAtom].size(); i++) {
		taucp.push_back(1 / (4 * dset .ncyc[LocalAtom][i]) - pwx_cp);
	}
	//
	ublas::matrix<double> p90x, p90y;
	ublas::matrix<double>* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<double>[taucp.size()];
	//
	CalcMatrix(p90x, dset, GlobalAtom, "90x", pwx_cp); //90(x) pulse on S spin
	CalcMatrix(p90y, dset, GlobalAtom, "90y", pwx_cp); // ..
	CalcMatrix(Gfree, dset, GlobalAtom, "free", taucp); // Free precession for all delays in taucp
	//
	//Equilibrium vector
	//
	ublas::vector<double> Sigma(Gfree[0].size1());
	//
	CalcCoherence("2IzSy", Sigma, dset); //the dset contains the basis information
	//
	propagate(p90x, Sigma);
	propagate(p90x, Sigma);
	//
	// Measure the peak intensity
	double InitMag = MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), Detect("2IzSy", 0, dset, Sigma),
			Detect("2IzSy", 1, dset, Sigma));
	//
	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {
		int ncyc = (int) floor(time_T2 * dset.ncyc[LocalAtom][ncycC] + 0.5);

		double Intensity[2];

		CalcCoherence("2IzSy", Sigma, dset); //the dset contains the basis information
		//
		//CPMG Train
		for (unsigned int ns = 0; ns < ncyc; ns++) {
			propagate(Gfree[ncycC], Sigma);
			propagate(p90y, Sigma);
			propagate(p90y, Sigma);
			propagate(Gfree[ncycC], Sigma);
			//Maybe 2pwn/pi
		}
		//180_x;
		propagate(p90x, Sigma);
		propagate(p90x, Sigma);
		//
		for (unsigned int ns = 0; ns < ncyc; ns++) {
			propagate(Gfree[ncycC], Sigma);
			propagate(p90y, Sigma);
			propagate(p90y, Sigma);
			propagate(Gfree[ncycC], Sigma);
		}
		//
		// store Trosy magnetization
		Intensity[0] = Detect("2IzSy", 0, dset, Sigma);
		Intensity[1] = Detect("2IzSy", 1, dset, Sigma);

		//
		dset.R2_calc[LocalAtom][ncycC] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),
				Intensity[0], Intensity[1]) / InitMag) / time_T2;

	}
	delete[] Gfree;
}
