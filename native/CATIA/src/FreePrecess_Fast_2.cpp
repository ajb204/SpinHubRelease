/*
 * FreePrecessFast.cpp
 *
 *  Created on: May 3, 2010
 *      Author: guillaume
 */

#include <complex>
#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>

namespace ublas = boost::numeric::ublas;

void Catia::FreePrecess_Fast_2(ublas::matrix<std::complex<double> >& G, Dataset& dset, int Atom) {
	//
	G.resize(2, 2);
	G.clear();

	std::vector<std::string> basis=dset.vpar("basis");
	// Now fetch the parameters from the different Dataset/Atom parameters
	//
	const double kex = Kex(dset);
	const double pb = Pb(dset);
	//
	std::string Nucl = dset._nucleus;
	double gammaS = 0.;
	double gammaI = 0.;
	double delta_csa = 0;
	double r_is = 0.;
	double hbar = 6.626075e-34;
	if (Nucl == "N") {
		gammaS = _gammaN;
	} else if (Nucl == "C") {
		gammaS = _gammaC;
	} else if (Nucl == "H") {
		gammaS = _gammaH;
	} else if (Nucl == "D") {
		gammaS = _gammaD;
	} else if (Nucl == "F") {
		gammaS = _gammaF;
	} else {
		std::cerr << " Could not resolve the nucleus of type " << Nucl << "\n";
		std::cerr << " in Dataset " << dset._id << "\n";
		std::cerr << " Functions: FreePrecess_Iph_7()\n";
		Abort(1);
	};
	//Store it for later use
	dset._gamma = gammaS;
	//
	char line[MAX_STRING_LENGTH];
	sprintf(line, "_%.0f", dset._sfrq);
	std::string fieldMarker(line);
	ClearBuf(line, sizeof(line));
	//
	std::string temperatureMarker("");
	if (_multipleTemperatures) {
		sprintf(line, "_%.1f", dset._temperature);
		temperatureMarker = line;
		ClearBuf(line, sizeof(line));
	}
	std::string marker(fieldMarker + temperatureMarker);
	//

	//Check that the parameters are there.
	std::vector<std::string> RequiredParam;
	RequiredParam.push_back("R0" + marker);
	for (unsigned int i = 0; i < RequiredParam.size(); i++) {
		if (!(BResolveParam(LocalParam[Atom], RequiredParam[i]))) {
			std::cerr << " The parameters:" << RequiredParam[i] << " is required by the basisset";
			std::cerr << " " << basis[1] << "\n but is not provided for atom " << AtomNumber2AtomName(Atom) << "\n";
			std::cerr << " please provide in the LocalParameter set\n";
			std::cerr << " Function FreePrecess_Fast_2()\n" << std::endl;
			Abort(1);
		}
	}
	//
	double R0 = ResolveParam(LocalParam[Atom], "R0" + marker);
	double DeltaO = (gammaS / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, Atom);
	//
	// kex additions
	std::complex<double> I(0, 1);
	G(0, 0) =  kex * pb + R0;	G(0, 1) = -kex * (1 - pb)                  ;
	G(1, 0) = -kex * pb     ;	G(1, 1) =  kex * (1 - pb) + R0 + I * DeltaO;
	//
}
