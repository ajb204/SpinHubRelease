/*
 * FreePrecess_Fast_3.cpp
 *
 *  Created on: May 13, 2010
 *      Author: guillaume
 */

#include <complex>
#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>

namespace ublas = boost::numeric::ublas;

void Catia::FreePrecess_Fast_3(ublas::matrix<std::complex<double> >& G, Dataset& dset, int Atom) {
	//
	G.resize(3, 3);
	G.clear();

	std::vector<std::string> basis=dset.vpar("basis");
	// Now fetch the parameters from the different Dataset/Atom parameters
	//
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
	const double R0 = ResolveParam(LocalParam[Atom], "R0" + marker);
	const double DeltaO_ab = (gammaS / _gammaH) * dset._sfrq * 2.0 * DFH_PI * DeltaOmega_ab(dset, Atom);
	const double DeltaO_ac = (gammaS / _gammaH) * dset._sfrq * 2.0 * DFH_PI * DeltaOmega_ac(dset, Atom);
	//
	// kex additions
	const double kex_ab = Kex_ab_3st(dset);
	const double kex_ac = Kex_ac_3st(dset);
	const double kex_bc = Kex_bc_3st(dset);

	const double pb = Pb_3st(dset);
	const double pc = Pc_3st(dset);
	const double pa = 1.0 - pb - pc;

	const double kab = kex_ab * pb / (pa + pb);
	const double kba = kex_ab * pa / (pa + pb);

	const double kac = kex_ac * pc / (pa + pc);
	const double kca = kex_ac * pa / (pa + pc);

	const double kbc = kex_bc * pc / (pb + pc);
	const double kcb = kex_bc * pb / (pb + pc);

	const std::complex<double> I(0, 1);

	G(0, 0) =  kab + kac + R0;	G(0, 1) = -kba                           ;	G(0, 2) = -kca                           ;
	G(1, 0) = -kab           ;	G(1, 1) =  kba + kbc + R0 + I * DeltaO_ab;	G(1, 2) = -kcb                           ;
	G(2, 0) = -kac           ;	G(2, 1) = -kbc                           ;	G(2, 2) =  kca + kcb + R0 + I * DeltaO_ac;
	//
}
