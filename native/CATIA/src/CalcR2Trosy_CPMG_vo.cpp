#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;
#include <MatrixExponentialLapack.hpp>

/*
 This is the function for the seqfil=Trosy_cpmg_vo
 The sequence is based on Vladislav Yu. Orekhov; May 2003
 (Single quantum proton part of NH_cpmg_vo.c)

 The sequence is:
 Trosy(y) -> { CPMG }_{phase(n)} -> Detect(Trosy(y))

 where the phase follows a xy super cycle:
 CPMGphase={x,y,x,y,y,x,y,x,-x,-y,-x,-y,-y,-x,-y,-x,.....}

 Writting by Flemming on September 26 2007

 Modified on September 27 2007 by Flemming
 - Added phase cycle: Instead of phase cycling the pulses
 and receiver, the CPMG block is phase cycled.
 Phase cycle removed :(

 */
void Catia::CalcR2Trosy_CPMG_vo(Dataset& dset, int GlobalAtom) {
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
			std::cerr << " Function: CalcR2CW_CPMG();\n";
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

	ublas::matrix<double> p90x, p90mx, p90y, p90my;
	ublas::matrix<double>* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<double>[taucp.size()];

	CalcMatrix(p90x, dset, GlobalAtom, "90x", pwx_cp); //90(x) pulse on S spin
	CalcMatrix(p90mx, dset, GlobalAtom, "90mx", pwx_cp); //90(-x) pulse on S spin
	CalcMatrix(p90y, dset, GlobalAtom, "90y", pwx_cp); // ..
	CalcMatrix(p90my, dset, GlobalAtom, "90my", pwx_cp);
	CalcMatrix(Gfree, dset, GlobalAtom, "free", taucp); // Free precession for all delays in taucp
	//
	//Equilibrium vector
	//
	ublas::vector<double> Sigma(Gfree[0].size1());
	//
	CalcCoherence("Tr(y)", Sigma, dset); //the dset contains the basis information
	//
	//
	propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma);
	propagate(p90y, Sigma);
	propagate(p90y, Sigma);
	propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma);
	//
	//
	// Measure the trosy peak intensity
	double InitMag = MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), (-Detect("Sy", 0, dset, Sigma)
			+ Detect("2IzSy", 0, dset, Sigma)) / 2., (-Detect("Sy", 1, dset, Sigma) + Detect("2IzSy", 1, dset, Sigma)) / 2.);
	//
	double Intensity[2][2];
	// We use the 16 xy step phase cycle for the CPMG train.

	const int CPMGphase[] = { 0, 1, 0, 1, 1, 0, 1, 0, 2, 3, 2, 3, 3, 2, 3, 2 };

	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {
		for (unsigned p = 0; p < 1; p++) {
			int ncyc = (int) floor(time_T2 * dset.ncyc[LocalAtom][ncycC] + 0.5);
			if (ncyc == 1 || ((ncyc % 2) != 0)) {
				std::cerr << " ncyc must be even!\n Function .CalcR2Trosy_CPMG_vo\n" << std::endl;
				std::cerr << " DataPoint #" << ncycC << " has ncyc=" << ncyc << "=" << time_T2 << "*" << dset.ncyc[LocalAtom][ncycC] << std::endl;
				Abort(1);
			}
			//Dump the trosy equilibrium from previous onto Sigma[]
			CalcCoherence("Tr(y)", Sigma, dset); //the dset contains the basis information
			//
			//CPMG Train
			for (unsigned int ns = 0; ns < 2 * ncyc; ns++) {
				propagate(Gfree[ncycC], Sigma);
				if (CPMGphase[(ns + 8 * p) % 16] == 0) {
					propagate(p90x, Sigma);
					propagate(p90x, Sigma);
				} else if (CPMGphase[(ns + 8 * p) % 16] == 1) {
					propagate(p90y, Sigma);
					propagate(p90y, Sigma);
				} else if (CPMGphase[(ns + 8 * p) % 16] == 2) {
					propagate(p90mx, Sigma);
					propagate(p90mx, Sigma);
				} else if (CPMGphase[(ns + 8 * p) % 16] == 3) {
					propagate(p90my, Sigma);
					propagate(p90my, Sigma);
				} else {
					std::cerr << " Wrong phase ?\n Function .CalcR2Trosy_CPMG_vo();\n" << std::endl;
					Abort(1);
				}
				propagate(Gfree[ncycC], Sigma);
			}
			propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
			propagate(p90y, Sigma);
			propagate(p90y, Sigma);
			propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
			//
			// store Trosy magnetizatio
			Intensity[p][0] = (Detect("2IzSy", 0, dset, Sigma) - Detect("Sy", 0, dset, Sigma)) / 2.;
			Intensity[p][1] = (Detect("2IzSy", 1, dset, Sigma) - Detect("Sy", 1, dset, Sigma)) / 2.;
		}
		//
		dset.R2_calc[LocalAtom][ncycC] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),
				Intensity[0][0], Intensity[0][1]) / InitMag) / time_T2;

		// std::cerr<<dset.id<<" ";
		//std::cerr<<"dset.R2_calc["<<LocalAtom<<"]["<<ncycC<<"] = "<<dset.R2_calc[LocalAtom][ncycC]<<std::endl;
	}
	delete[] Gfree;
	//  exit(10);
}
