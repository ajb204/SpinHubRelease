#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;
#include <MatrixExponentialLapack.hpp>

void Catia::CalcR2AntiTrosySel_CPMG(Dataset& dset, int GlobalAtom) {
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
	taucp.push_back(atof(dset.vpar("taub")[1].c_str()));

	ublas::matrix<double> Hinv, p90x, p90mx, p90y, p90my;
	ublas::matrix<double>* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<double>[taucp.size()];

	CalcMatrix(Hinv, dset, GlobalAtom, "Iinv", 180.); //invert the I spin
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
	double Intensity0[2][2];
	for (unsigned int p = 0; p < 2; p++) {
		CalcCoherence("ATr(y)", Sigma, dset); //the dset contains the basis information
		propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
		propagate(p90mx, Sigma);
		propagate(p90y, Sigma);
		propagate(Hinv, Sigma);
		propagate(p90y, Sigma);
		propagate(p90mx, Sigma);
		propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
		if (p == 0) {
			propagate(p90x, Sigma);
		} else {
			propagate(p90mx, Sigma);
		}
		propagate(Hinv, Sigma);
		propagate(p90y, Sigma);
		propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma);
		Intensity0[p][0] = (Detect("Sz", 0, dset, Sigma) + Detect("2IzSz", 0, dset, Sigma)) / 2.;
		Intensity0[p][1] = (Detect("Sz", 1, dset, Sigma) + Detect("2IzSz", 1, dset, Sigma)) / 2.;
	}
	//
	// Measure the trosy peak intensity
	double InitMag = MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),
			(Intensity0[0][0] + Intensity0[1][0]) / 2., (Intensity0[0][1] + Intensity0[1][1]) / 2.);
	//
	//
	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {
		int ncyc = (int) floor(time_T2 * dset.ncyc[LocalAtom][ncycC] + 0.5);
		if (fabs(float(ncyc) - time_T2 * dset.ncyc[LocalAtom][ncycC]) > 0.05) {
			std::cerr << "\n There is a mismatch between nu_cpmg(Hz) and time_T2 (sec)\n";
			std::cerr << " Observed for Dataset " << dset._id << std::endl;
			std::cerr << " Number of cycles calculated: " << time_T2 * dset.ncyc[LocalAtom][ncycC] << " which is not an integer\n";
			std::cerr << " Function: CalcR2AntiTrosy_CPMG();\n";
			std::cerr << std::endl;
			Abort(1);
		}

		double Intensity[2][2];

		for (unsigned int p = 0; p < 2; p++) {
			CalcCoherence("ATr(y)", Sigma, dset); //the dset contains the basis information
			//
			//CPMG Train
			for (unsigned int ns = 0; ns < ncyc; ns++) {
				propagate(Gfree[ncycC], Sigma);
				propagate(p90y, Sigma);
				propagate(p90y, Sigma);
				propagate(Gfree[ncycC], Sigma);
				//Maybe 2pwn/pi
			}
			// Palmer Element - START
			if (p == 0) {
				propagate(p90y, Sigma);
			} else {
				propagate(p90my, Sigma);
			}
			propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //taub
			//180(compositeX);
			propagate(p90mx, Sigma);
			propagate(p90y, Sigma);
			propagate(Hinv, Sigma);
			propagate(p90y, Sigma);
			propagate(p90mx, Sigma);
			propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //taub
			propagate(Hinv, Sigma);
			if (p == 0) {
				propagate(p90x, Sigma);
			} else {
				propagate(p90mx, Sigma);
			}
			// Palmer Element - END
			for (unsigned int ns = 0; ns < ncyc; ns++) {
				propagate(Gfree[ncycC], Sigma);
				propagate(p90x, Sigma);
				propagate(p90x, Sigma);
				propagate(Gfree[ncycC], Sigma);
			}
			propagate(p90y, Sigma);
			propagate(Hinv, Sigma);
			propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
			//
			// store Trosy magnetization
			Intensity[p][0] = (Detect("2IzSz", 0, dset, Sigma) - Detect("Sz", 0, dset, Sigma)) / 2.;
			Intensity[p][1] = (Detect("2IzSz", 1, dset, Sigma) - Detect("Sz", 1, dset, Sigma)) / 2.;
		}
		//
		dset.R2_calc[LocalAtom][ncycC] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),
				(Intensity[0][0] + Intensity[1][0]) / 2., (Intensity[0][1] + Intensity[1][1]) / 2.) / InitMag) / time_T2;

		// std::cerr<<dset.id<<" ";
		//std::cerr<<"dset.R2_calc["<<LocalAtom<<"]["<<ncycC<<"] = "<<dset.R2_calc[LocalAtom][ncycC]<<std::endl;
	}
	delete[] Gfree;
	//exit(10);
}
