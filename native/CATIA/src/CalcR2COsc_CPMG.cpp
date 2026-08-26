#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric::ublas;
#include <MatrixExponentialLapack.hpp>
/*
 This is the sequence of CO CPMG, with inversion of the
 side-chain CO in the middle of the experiment, to remove
 to a first order, to CO(bb) - CO(sc) couplings.

 sequence file: hnco_CO_CPMG_SCFilter_600_dfh.c

 Flemming, April 2008

 */

void Catia::CalcR2COsc_CPMG(Dataset& dset, int GlobalAtom) {
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
			std::cerr << " Function: CalcR2COsc_CPMG();\n";
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
	double taucc = atof(dset.vpar("taucc")[1].c_str()); // this is 1/2*JCC - approx 9ms
	double pwx_cp = atof(dset.vpar("pwx_cp")[1].c_str());
	for (unsigned int i = 0; i < dset.ncyc[LocalAtom].size(); i++) {
		taucp.push_back(1 / (4 * dset.ncyc[LocalAtom][i]) - pwx_cp);
	}
	taucp.push_back(atof(dset.vpar("time_equil")[1].c_str()));
	taucp.push_back(-2. * pwx_cp / DFH_PI);
	taucp.push_back(taucc);

	ublas::matrix<double> Hinv, p90x, p90mx, p90y, p90my;
	ublas::matrix<double>* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<double>[taucp.size()];

	CalcMatrix(Hinv, dset, GlobalAtom, "Iinv", 180.); //invert the I spin
	CalcMatrix(p90x, dset, GlobalAtom, "90x", pwx_cp);
	CalcMatrix(p90mx, dset, GlobalAtom, "90mx", pwx_cp);
	CalcMatrix(p90y, dset, GlobalAtom, "90y", pwx_cp);
	CalcMatrix(p90my, dset, GlobalAtom, "90my", pwx_cp);
	CalcMatrix(Gfree, dset, GlobalAtom, "free", taucp);
	//
	//Equilibrium vector
	//
	ublas::vector<double> Sigma(Gfree[0].size1());
	ublas::vector<double> Intensity(2);
	//
	// INITIAL CONDITION AND NCYC=0 PLANE
	//
	CalcCoherence("Sz", Sigma, dset); //the dset contains the basis information
	// evolve the taucc element.
	propagate(Gfree[dset.ncyc[LocalAtom].size() + 2], Sigma); //taucc delay
	//180(y);
	propagate(p90y, Sigma);
	propagate(Hinv, Sigma);
	propagate(p90y, Sigma);
	//
	propagate(Gfree[dset.ncyc[LocalAtom].size() + 2], Sigma); //taucc delay
	//
	double InitMag = MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), Detect("Sz", 0, dset, Sigma),
			Detect("Sz", 1, dset, Sigma));
	//
	// START THE REAL EXPERIMENT
	//
	//
	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {
		int ncyc = (int) floor(time_T2 * dset.ncyc[LocalAtom][ncycC] + 0.5);
		// phase cycle
		for (unsigned int p = 0; p < 2; p++) {
			CalcCoherence("Sz", Sigma, dset); //the dset contains the basis information
			//
			//CPMG Train
			propagate(p90y, Sigma);
			for (unsigned int ns = 0; ns < ncyc; ns++) {
				propagate(Gfree[ncycC], Sigma);
				if (ns == 0) {
					propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //the -2*pw/PI delay.
				}
				propagate(p90x, Sigma);
				propagate(p90x, Sigma);
				propagate(Gfree[ncycC], Sigma);
				if (ns == ncyc - 1) {
					propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //the -2*pw/PI delay.
				}
			}
			//
			//inversion element begin.
			propagate(p90y, Sigma);
			propagate(Gfree[dset.ncyc[LocalAtom].size() + 2], Sigma); //taucc delay
			//180(y);
			if (p == 0) {
				propagate(p90y, Sigma);
				propagate(Hinv, Sigma);
				propagate(p90y, Sigma);
			} else if (p == 1) {
				propagate(p90my, Sigma);
				propagate(Hinv, Sigma);
				propagate(p90my, Sigma);
			}
			//
			propagate(Gfree[dset.ncyc[LocalAtom].size() + 2], Sigma); //taucc delay
			propagate(p90my, Sigma);
			// inversion element end
			//
			for (unsigned int ns = 0; ns < ncyc; ns++) {
				propagate(Gfree[ncycC], Sigma);
				if (ns == 0) {
					propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //the -2*pwn_cp/PI delay.
				}
				propagate(p90x, Sigma);
				propagate(p90x, Sigma);
				propagate(Gfree[ncycC], Sigma);
				if (ns == ncyc - 1) {
					propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //the -2*pwc_cp/PI delay.
				}
			}
			//
			// We need the equilibration delay
			propagate(p90y, Sigma);
			//
			//Delay of time_equil
			propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); // time_equil
			Intensity[p] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), Detect("Sz", 0, dset,
					Sigma), Detect("Sz", 1, dset, Sigma)) / InitMag) / time_T2;
		}
		dset.R2_calc[LocalAtom][ncycC] = (Intensity[0] + Intensity[1]) / 2.;
	}
	delete[] Gfree;
}
