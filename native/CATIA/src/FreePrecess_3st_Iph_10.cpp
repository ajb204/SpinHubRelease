/*
 * FreePrecess_3st_Iph_10.cpp
 *
 *  Created on: Jul 3, 2009
 *      Author: guillaume
 */

#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>
namespace ublas = boost::numeric::ublas;
void Catia::FreePrecess_3st_Iph_10(ublas::matrix<double>& G, Dataset& dset, int Atom) {
	/*
	 NameOfBasis: 3st_Iph_10

	 The basis is:
	 { E,
	 [ Sx ,Sy ,Sz] x { Site A, Site B, Site C }
	 }

	 Thus - a 10 element basis
	 */
	std::vector<std::string> basis = dset.vpar("basis");
	//
	G.resize(10, 10);
	G.clear();
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
		delta_csa = -172E-6;
		r_is = 1.02E-10;
	} else if (Nucl == "C") {
		gammaS = _gammaC;
		delta_csa = -25E-6; //Calpha
		r_is = 1.098e-10; // From Kowalewski, JMR, 2002, p171-177
	} else if (Nucl == "H") {
		gammaS = _gammaH;
		delta_csa = -10E-6;
		r_is = 1.02E-10;
	} else if (Nucl == "D") {
		gammaS = _gammaD;
		delta_csa = 0.;
		r_is = 1.0E-10;
	} else if (Nucl == "F") {
		gammaS = _gammaF;
		delta_csa = 0.;
		r_is = 1.0E-10;
	} else {
		std::cerr << " Could not resolve the nucleus of type " << Nucl << "\n";
		std::cerr << " in Dataset " << dset._id << "\n";
		std::cerr << " Functions: FreePrecess_3st_Iph_10()\n";
		Abort(1);
	}
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
	RequiredParam.push_back("R1iph" + marker);
	RequiredParam.push_back("R0iph" + marker);
	RequiredParam.push_back("Omega" + temperatureMarker);
	for (unsigned int i = 0; i < RequiredParam.size(); i++) {
		if (!(BResolveParam(LocalParam[Atom], RequiredParam[i]))) {
			std::cerr << " The parameters:" << RequiredParam[i] << " is required by the basisset";
			std::cerr << " " << basis[1] << "\n but is not provided for atom " << AtomNumber2AtomName(Atom) << "\n";
			std::cerr << " please provide in the LocalParameter set\n";
			std::cerr << " Function FreePrecess_3st_Iph_10()\n" << std::endl;
			Abort(1);
		}
	}
	double R1iph = ResolveParam(LocalParam[Atom], "R1iph" + marker);
	double R0iph = ResolveParam(LocalParam[Atom], "R0iph" + marker);
	double Omega = (gammaS / _gammaH) * dset._sfrq * 2 * DFH_PI * (ResolveParam(LocalParam[Atom], "Omega" + temperatureMarker) - dset._xcar);
	double DeltaO_ab = (gammaS / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega_ab(dset, Atom);
	double DeltaO_ac = (gammaS / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega_ac(dset, Atom);
	//
	gammaI = _gammaH;
	//
	// Approximate value for equlibrium value
	const double pb = Pb_3st(dset);
	const double pc = Pc_3st(dset);
	const double pa = 1 - pb - pc;
	const double p[] = { pa, pb, pc };
	const double omega[] = { Omega, Omega + DeltaO_ab, Omega + DeltaO_ac };
	for (unsigned int i = 0; i < 3; ++i) {
		// site i
		//
		const double offset = 3 * i;
		// Auto Relaxations
		G(1 + offset, 1 + offset) = R0iph;
		G(2 + offset, 2 + offset) = R0iph;
		G(3 + offset, 3 + offset) =  R1iph;
		G(3 + offset, 0) = -p[i] * (R1iph) * gammaS / (gammaI * 2.);
		//
		// Omega
		G(1 + offset, 2 + offset) =  omega[i];
		G(2 + offset, 1 + offset) = -omega[i];

	}

	// kex additions

	const double kex_ab = Kex_ab_3st(dset);
	const double kex_ac = Kex_ac_3st(dset);
	const double kex_bc = Kex_bc_3st(dset);

	const double kex_p_ab = kex_ab / (pa + pb);
	const double kex_p_ac = kex_ac / (pa + pc);
	const double kex_p_bc = kex_bc / (pb + pc);

	const double k_ab = kex_p_ab * pb;
	const double k_ba = kex_p_ab * pa;
	const double k_ac = kex_p_ac * pc;
	const double k_ca = kex_p_ac * pa;
	const double k_bc = kex_p_bc * pc;
	const double k_cb = kex_p_bc * pb;

	for (unsigned int i = 0; i < 3; i++) {
		G(i + 1, i + 1) = G(i + 1, i + 1) + k_ab;
		G(i + 4, i + 1) = G(i + 4, i + 1) - k_ab;
		//
		G(i + 4, i + 4) = G(i + 4, i + 4) + k_ba;
		G(i + 1, i + 4) = G(i + 1, i + 4) - k_ba;
		//
		G(i + 1, i + 1) = G(i + 1, i + 1) + k_ac;
		G(i + 7, i + 1) = G(i + 7, i + 1) - k_ac;
		//
		G(i + 7, i + 7) = G(i + 7, i + 7) + k_ca;
		G(i + 1, i + 7) = G(i + 1, i + 7) - k_ca;
		//
		G(i + 4, i + 4) = G(i + 4, i + 4) + k_bc;
		G(i + 7, i + 4) = G(i + 7, i + 4) - k_bc;
		//
		G(i + 7, i + 7) = G(i + 7, i + 7) + k_cb;
		G(i + 4, i + 7) = G(i + 4, i + 7) - k_cb;
		//
	}
}
