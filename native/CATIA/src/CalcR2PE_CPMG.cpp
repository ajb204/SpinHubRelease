#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

#include <MatrixExponentialLapack.hpp>

/*
 This sequence is the seqfil="PE_CPMG"

 Essentially it is build up of
 2IzSy -> {CPMG}_y - taub - 180 - taub - {CPMG}_x -> Sx

 */
void Catia::CalcR2PE_CPMG(Dataset& dset, int GlobalAtom) {
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
	taucp.push_back(atof(dset.vpar("taub")[1].c_str()));
	//
	ublas::matrix<double> Hinv, p90x, p90mx, p90y, p90my;
	ublas::matrix<double>* Gfree;
	//Declare the Gfree according to the taucp array;
	Gfree = new ublas::matrix<double>[taucp.size()];
	//
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
	//ublas::vector<double> Sogma(Gfree[0].size1());
	//
	CalcCoherence("2IzSy", Sigma, dset); //the dset contains the basis information
	//CalcCoherence("2IzSy", Sogma, dset); //the dset contains the basis information
	//
	//

	
	/*
	Sigma = ublas::prod(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
	Sigma = ublas::prod(p90y, Sigma);
	Sigma = ublas::prod(p90x, Sigma);
	Sigma = ublas::prod(Hinv, Sigma);
	Sigma = ublas::prod(p90x, Sigma);
	Sigma = ublas::prod(p90y, Sigma);
	Sigma = ublas::prod(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
	Sigma = ublas::prod(p90y, Sigma);
	Sigma = ublas::prod(Gfree[dset.ncyc[LocalAtom].size()], Sigma);
	*/
	
	propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
	propagate(p90y, Sigma);
	propagate(p90x, Sigma);
	propagate(Hinv, Sigma);
	propagate(p90x, Sigma);
	propagate(p90y, Sigma);
	propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma);
	propagate(p90y, Sigma);
	propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma);
	

	

	//
	//
	// Measure the trosy peak intensity
	double InitMag = MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), Detect("Sz", 0,
			dset, Sigma), Detect("Sz", 1, dset, Sigma));
	//
	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) {
		int ncyc = (int) floor(time_T2 * dset.ncyc[LocalAtom][ncycC] + 0.5);

		double Intensity[2];

		CalcCoherence("2IzSy", Sigma, dset); //the dset contains the basis information
		//
		//CPMG Train
		for (unsigned int ns = 0; ns < ncyc; ns++) {
		  /*
		    Sigma = ublas::prod(Gfree[ncycC], Sigma);
			Sigma = ublas::prod(p90y, Sigma);
			Sigma = ublas::prod(p90y, Sigma);
			Sigma = ublas::prod(Gfree[ncycC], Sigma);
			//Maybe 2pwn/pi
			*/	

		  propagate(Gfree[ncycC], Sigma);
		  propagate(p90y, Sigma);
		  propagate(p90y, Sigma);
		  propagate(Gfree[ncycC], Sigma);
		  //Maybe 2pwn/pi


		}
		/*
		Sigma = ublas::prod(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //taub
		//
		//180(compositeX);
		Sigma = ublas::prod(p90y, Sigma);
		Sigma = ublas::prod(p90x, Sigma);
		Sigma = ublas::prod(Hinv, Sigma);
		Sigma = ublas::prod(p90x, Sigma);
		Sigma = ublas::prod(p90y, Sigma);
		Sigma = ublas::prod(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //taub
		// Palmer Element - END
		*/
		


		propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //taub
		//180(compositeX);
		propagate(p90y, Sigma);
		propagate(p90x, Sigma);
		propagate(Hinv, Sigma);
		propagate(p90x, Sigma);
		propagate(p90y, Sigma);
		propagate(Gfree[dset.ncyc[LocalAtom].size() + 1], Sigma); //taub
		// Palmer Element - END

		for (unsigned int ns = 0; ns < ncyc; ns++) {
		  /*
		    Sigma = ublas::prod(Gfree[ncycC], Sigma);
		    Sigma = ublas::prod(p90x, Sigma);
		    Sigma = ublas::prod(p90x, Sigma);
		    Sigma = ublas::prod(Gfree[ncycC], Sigma);
		  */
		  
		  propagate(Gfree[ncycC], Sigma);
		  propagate(p90x, Sigma);
		  propagate(p90x, Sigma);
		  propagate(Gfree[ncycC], Sigma);
		  
		}
		/*
		Sigma = ublas::prod(p90y, Sigma);
		Sigma = ublas::prod(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
		*/
		
		propagate(p90y, Sigma);
		propagate(Gfree[dset.ncyc[LocalAtom].size()], Sigma); //Delay of time_equil
		
		//
		// store Trosy magnetization
		Intensity[0] = Detect("Sz", 0, dset, Sigma);
		Intensity[1] = Detect("Sz", 1, dset, Sigma);
		//
		dset.R2_calc[LocalAtom][ncycC] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),
				Intensity[0], Intensity[1]) / InitMag) / time_T2;

	}
	delete[] Gfree;
	//  exit(10);
}
