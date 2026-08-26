#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;
#include <MatrixExponentialLapack.hpp>

/*
 Modified on May 13 2009 by Flemming
 - changed the phase cycle to include phases on the initial 90 and
 the 180 in the middle of the CPMG element.


 */

void Catia::CalcR2CW_CPMG(Dataset& dset, int GlobalAtom) {
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
		taucp.push_back(1 / (4 * dset .ncyc[LocalAtom][i]) - pwx_cp);
	}
	taucp.push_back(atof(dset.vpar("time_equil")[1].c_str()));
	taucp.push_back(-2 * pwx_cp / DFH_PI);

	ublas::matrix<double> Hinv, p90x, p90mx, p90y;
	ublas::matrix<double>* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<double>[taucp.size()];

	CalcMatrix(Hinv, dset, GlobalAtom, "Iinv", 180.); //invert the I spin
	CalcMatrix(p90x, dset, GlobalAtom, "90x", pwx_cp);
	CalcMatrix(p90mx, dset, GlobalAtom, "90mx", pwx_cp);
	CalcMatrix(p90y, dset, GlobalAtom, "90y", pwx_cp);
	CalcMatrix(Gfree, dset, GlobalAtom, "free", taucp);
	//
	//Equilibrium vector
	//
	ublas::vector<double> Sigma(Gfree[0].size1());
	//
	double Intensity[2];
	Intensity[0] = 0;
	Intensity[1] = 0;

	for (unsigned int p1 = 0; p1 < 2; p1++) {
		for (unsigned int p2 = 0; p2 < 2; p2++) {
			CalcCoherence("Sz", Sigma, dset); //the dset contains the basis information
			if (p1 == 0) {
				propagate(p90x, Sigma);
			} else {
				propagate(p90mx, Sigma);
			}
			if (p2 == 0) {
				propagate(p90x, Sigma);
				propagate(Hinv, Sigma);
				propagate(p90x, Sigma);
			} else {
				propagate(p90mx, Sigma);
				propagate(Hinv, Sigma);
				propagate(p90mx, Sigma);
			}
			propagate(p90mx, Sigma);
			propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
			//
			Intensity[0] += Detect("Sz", 0, dset, Sigma) * (p1 - 0.5) * 2.;
			Intensity[1] += Detect("Sz", 1, dset, Sigma) * (p1 - 0.5) * 2.;
		}
	}
	//
	double InitMag =
			MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), Intensity[0] / 4., Intensity[1] / 4.);
	//
	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {
		Intensity[0] = 0.;
		Intensity[1] = 0.;
		//
		int ncyc = (int) floor(time_T2 * dset.ncyc[LocalAtom][ncycC] + 0.5);
		if (fabs(float(ncyc) - time_T2 * dset.ncyc[LocalAtom][ncycC]) > 0.05) {
			std::cerr << "\n There is a mismatch between nu_cpmg(Hz) and time_T2 (sec)\n";
			std::cerr << " Observed for Dataset " << dset._id << std::endl;
			std::cerr << " Number of cycles calculated: " << time_T2 * dset.ncyc[LocalAtom][ncycC] << " which is not an integer\n";
			std::cerr << " Function: CalcR2CW_CPMG();\n";
			std::cerr << std::endl;
			Abort(1);
		}
		for (unsigned int p1 = 0; p1 < 2; p1++) {
			for (unsigned int p2 = 0; p2 < 2; p2++) {

				CalcCoherence("Sz", Sigma, dset); //the dset contains the basis information
				//
				// Initial excitation pulse
				if (p1 == 0) {
					propagate(p90x, Sigma);
				} else {
					propagate(p90mx, Sigma);
				}
				//
				//CPMG Train
				for (unsigned int ns = 0; ns < ncyc; ns++) {
					propagate(Gfree[ncycC], Sigma);
					if (ns == 0) {
						propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); // -2.*pw/PI delay
					}
					propagate(p90y, Sigma);
					propagate(p90y, Sigma);
					propagate(Gfree[ncycC], Sigma);
				}
				//180(x);
				if (p2 == 0) {
					propagate(p90x, Sigma);
					propagate(Hinv, Sigma);
					propagate(p90x, Sigma);
				} else {
					propagate(p90mx, Sigma);
					propagate(Hinv, Sigma);
					propagate(p90mx, Sigma);
				}
				for (unsigned int ns = 0; ns < ncyc; ns++) {
					propagate(Gfree[ncycC], Sigma);
					propagate(p90y, Sigma);
					propagate(p90y, Sigma);
					propagate(Gfree[ncycC], Sigma);
					if (ns == ncyc - 1) {
						propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); // -2.*pw/PI delay
					}
				}
				// We need the equilibration delay
				propagate(p90mx, Sigma);
				//
				propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
				//
				Intensity[0] += -Detect("Sz", 0, dset, Sigma) * (p1 - 0.5) * 2.;
				Intensity[1] += -Detect("Sz", 1, dset, Sigma) * (p1 - 0.5) * 2.;
			}
		}

		dset.R2_calc[LocalAtom][ncycC] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),
				Intensity[0] / 4., Intensity[1] / 4.) / (InitMag)) / time_T2;

	}
	delete[] Gfree;
}
