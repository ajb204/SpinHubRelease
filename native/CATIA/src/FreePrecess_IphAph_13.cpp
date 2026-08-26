//
// June 3, 2009:
//    Modified by G. Bouvignies:
//                  - to take into account the temperature dependance of some variables;
//                  - to use function for calculating kex, pb and delatO.
//

#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>
namespace ublas = boost::numeric::ublas;
using boost::to_lower_copy;

void Catia::FreePrecess_IphAph_13(ublas::matrix<double>& G, Dataset& dset, int Atom) {
	/*
	 NameOfBasis: IphAph_13

	 The basis is:
	 { E,
	 [ Ix ,Iy ,Iz,
	 2SxIz,2SyIz,2SzIz,
	 ] (x) { Site A, Site B }
	 }

	 Thus - a 13 element basis

	 Modified Guillaume and Flemming: April 2009:
	 Check if the relaxation of anti-phase is due to external spins
	 or caused by the directly attached spin.

	 */
	//
	std::vector<std::string> basis = dset.vpar("basis");
	//
	G.resize(13, 13);
	G.clear();
	// Fetch the parameters from the different Dataset/Atom parameters
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
		delta_csa = -172E-6;
		r_is = 1.02E-10;
	} else if (Nucl == "C") {
		gammaS = _gammaC;
		delta_csa = -25E-6;
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
		std::cerr << " Functions: FreePrecess_IphAph_13()\n";
		Abort(1);
	}
	//Store it for later use
	dset._gamma = gammaS;
	//
	std::string cn = to_lower_copy(dset.vpar("couplednucleus")[1]);
	if (cn.find('n') < cn.npos) {
		gammaI = _gammaN;
	}
	if (cn.find('h') < cn.npos) {
		gammaI = _gammaH;
	}
	if (cn.find('c') < cn.npos) {
		gammaI = _gammaC;
	}
	if (fabs(gammaI) < 1E-6) {
		std::cerr << " Cannot convert " << cn << " into a nucleus type\n";
		std::cerr << " for the coupled nucleus\n";
		std::cerr << " Function FreePrecess_IphAph_13();\n";
		Abort(1);
	}
	if (fabs(gammaI - gammaS) < 1E-6) {
		std::cerr << " You have defined a homo-nuclear spin-system\n";
		std::cerr << " program is not ready for this yet\n";
		std::cerr << " Function FreePrecess_IphAph();\n";
		Abort(1);
	}
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
	RequiredParam.push_back("R1aph" + marker);
	RequiredParam.push_back("R0iph" + marker);
	//  RequiredParam.push_back("R0aph"+marker); estimate from R1aph-R1iph
	RequiredParam.push_back("Omega" + temperatureMarker);
	//RequiredParam.push_back("JIS");
	RequiredParam.push_back("DeltaJ");
	for (unsigned int i = 0; i < RequiredParam.size(); i++) {
		if (!(BResolveParam(LocalParam[Atom], RequiredParam[i]))) {
			std::cerr << " The parameter:" << RequiredParam[i] << " is required by the basisset";
			std::cerr << " " << basis[1] << "\n but is not provided for atom " << AtomNumber2AtomName(Atom) << "\n";
			std::cerr << " please provide in the LocalParameter set\n";
			std::cerr << " Function FreePrecess_IphAph_13()\n" << std::endl;
			Abort(1);
		}
	}
	if (!(BResolveParam(LocalParam[Atom], "JIS") || BResolveParam(LocalParam[Atom], "JIS" + marker))) {
		std::cerr << " The parameter JIS or JIS" << marker << " is required by the basisset";
		std::cerr << " " << basis[1] << "\n but is not provided for atom " << AtomNumber2AtomName(Atom) << "\n";
		std::cerr << " please provide one of these parameters in the LocalParameter set\n";
		std::cerr << " Function FreePrecess_IphAph_13()\n" << std::endl;
		Abort(1);
	}
	double R1iph = ResolveParam(LocalParam[Atom], "R1iph" + marker);
	double R1aph = ResolveParam(LocalParam[Atom], "R1aph" + marker);
	double R0iph = ResolveParam(LocalParam[Atom], "R0iph" + marker);
	//
	//  double R0aph=R0iph+(R1aph-R1iph);
	double R0aph = 0;
	if (Nucl == "N" && cn.find('h') < cn.npos) {// The scenario that H relaxes due to external spins
		R0aph = R0iph + (R1aph - R1iph);
	} else if (Nucl == "H" && (cn.find('n') < cn.npos || cn.find('c') < cn.npos)) { // Nz relaxes due to the directly attached H
		R0aph = R0iph - (R1aph - R1iph);
	} else if (Nucl == "C" && cn.find('n') < cn.npos) {
		// COyNz relaxes due to CO(csa) and R1(N)
		//   approximation invoked:  R(CzNz) ~ 1.5*R1(N), Read in R1(N) as R1aph
		//                           R(Cz) ~ 0.6 - 1.1/s, Read in ~0.6 as R1iph (500/800 MHz, ~6.5 ns tauc)
		//                           R(Cy) ~ R(CO,CSA),   Read in R2(inf) estimate as R0iph
		//                           R(CyNz) ~ R(CO,CSA) + R1(N), R0aph = R0iph - R1aph
		R0aph = R0iph - R1aph;
	} else if (Nucl == "C" && cn.find('h') < cn.npos) {
		R0aph = R0iph - R1aph;
	} else {
		std::cerr << " The basis IphAph_13 only allows for 15N-1H, 1H-15N, 13CO-15N and methyl spin systems.";
		std::cerr << " Function .FreePrecess_IphAph_13()" << std::endl;
		Abort(1);
	}
	double Omega = (gammaS / _gammaH) * dset._sfrq * 2 * DFH_PI * (ResolveParam(LocalParam[Atom], "Omega" + temperatureMarker) - dset._xcar);

	double DeltaO = (gammaS / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, Atom);
	double JIS = 0.;
	if (BResolveParam(LocalParam[Atom], "JIS" + marker)) {
		JIS = ResolveParam(LocalParam[Atom], "JIS" + marker);
	} else {
		JIS = ResolveParam(LocalParam[Atom], "JIS");
	}
	double DeltaJ = ResolveParam(LocalParam[Atom], "DeltaJ");
	if (dset.Bvpar("deltajscaling")) {
		double djs = atof(dset.vpar("deltajscaling")[0].c_str());
		DeltaJ = DeltaJ * djs;
		sprintf(line, "Scaling used for DeltaJ, i.e., DeltaJ(%s)=%g*DeltaJ", dset._id.c_str(), djs);
		std::string note(line);
		ClearBuf(line, sizeof(line));
		LocalNotes[Atom]["DeltaJ@" + dset._id] = note;
	}
	double EtaZ = 0.;
	double EtaXY = 0.;
	bool estimate_etaZ = false;
	bool estimate_etaXY = false;
	for (unsigned int i = 0; i < basis.size(); i++) {
		if (to_lower_copy(basis[i]) == "estimate_etaz") {
			estimate_etaZ = true;
		}
		if (to_lower_copy(basis[i]) == "estimate_etaxy") {
			estimate_etaXY = true;
		}
	}
	if (estimate_etaZ || estimate_etaXY) {
		//We estimate from R1iph and assume
		//CSA of -172 for Nitrogen, 25 for Carbon, and -10ppm for H
		double B0 = dset._sfrq * 1E6 * 2. * DFH_PI / _gammaH;
		double cc = B0 * gammaS * delta_csa / sqrt(3.);
		double dd = (1.e-7) * hbar * gammaI * gammaS * pow(r_is, -3.) / (DFH_PI * 2.);
		double JwS = R1iph / (cc * cc + dd * dd * 0.75);
		if (estimate_etaZ) {
			EtaZ = -sqrt(3.) * cc * dd * JwS;
		}
		if (estimate_etaXY) {
			double J0 = 0.25 * (R0iph / (dd * dd / 8. + cc * cc / 6.) - 3. * JwS);
			EtaXY = -(sqrt(3.) / 6.) * cc * dd * (4. * J0 + 3. * JwS);
		}
	}
	// Approximate value for equlibrium value
	// site A
	//
	for (unsigned int i = 0; i < 2; i++) {
		// Auto Relaxations
		G(1 + i * 6, 1 + i * 6) = R0iph;
		G(2 + i * 6, 2 + i * 6) = R0iph;
		G(3 + i * 6, 3 + i * 6) = R1iph;
		G(3 + i * 6, 0) = -(1 - pb + i * (2 * pb - 1)) * R1iph * gammaS / (_gammaH * 2.);
		G(4 + i * 6, 4 + i * 6) = R0aph;
		G(5 + i * 6, 5 + i * 6) = R0aph;
		G(6 + i * 6, 6 + i * 6) = R1aph;
		//    G(6+i*6,0) =    -(1-pb+i*(2*pb-1))*R1aph*gammaS/(gammaI*2.);

		G(4 + i * 6, 1 + i * 6) = EtaXY;
		G(5 + i * 6, 2 + i * 6) = EtaXY;
		G(6 + i * 6, 3 + i * 6) = EtaZ;

		G(1 + i * 6, 4 + i * 6) = EtaXY;
		G(2 + i * 6, 5 + i * 6) = EtaXY;
		G(3 + i * 6, 6 + i * 6) = EtaZ;
		//
		// J couplings
		G(1 + i * 6, 5 + i * 6) = (JIS + i * DeltaJ) * DFH_PI;
		G(5 + i * 6, 1 + i * 6) = -(JIS + i * DeltaJ) * DFH_PI;
		G(2 + i * 6, 4 + i * 6) = -(JIS + i * DeltaJ) * DFH_PI;
		G(4 + i * 6, 2 + i * 6) = (JIS + i * DeltaJ) * DFH_PI;
		//
		// Omega;
		G(1 + i * 6, 2 + i * 6) = (Omega + i * DeltaO);
		G(2 + i * 6, 1 + i * 6) = -(Omega + i * DeltaO);
		G(4 + i * 6, 5 + i * 6) = (Omega + i * DeltaO);
		G(5 + i * 6, 4 + i * 6) = -(Omega + i * DeltaO);
	}
	//
	// kex additions
	for (unsigned int i = 0; i < 6; i++) {
		G(i + 1, i + 1) = G(i + 1, i + 1) + kex * pb;
		G(i + 7, i + 1) = G(i + 7, i + 1) - kex * pb;
		//
		G(i + 7, i + 7) = G(i + 7, i + 7) + kex * (1 - pb);
		G(i + 1, i + 7) = G(i + 1, i + 7) - kex * (1 - pb);
	}
}
