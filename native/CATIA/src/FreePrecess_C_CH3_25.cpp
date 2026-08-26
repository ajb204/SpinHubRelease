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

void Catia::FreePrecess_C_CH3_25(ublas::matrix<double>& G, Dataset& dset, int Atom) {
	/*
	 NameOfBasis: C_CH3_25

	 The basis is:
	 { E,
	 [ A(x) ,A(y) ,A(z),
	 B(x), B(y), B(z),
	 C(x), C(y), C(z),
	 D(x), D(y), D(z),
	 ] x { Site A, Site B }
	 }

	 |-proton-|
	 where      A(x) is the Cx|aaa><aaa| carbon transistion =
	 (1/8 + Hz/4 + HzHz/2 + HzHzHz )

	 B(x) is the Cx(|aab><aab| + |aba><aba| + |baa><baa| ) =
	 (3/8 + Hz/4 - HzHz/2 - 3HzHzHz )

	 C(x) is the Cx(|abb><abb| + |bab><bab| + |bba><bba| ) =
	 (3/8 - Hz/4 - HzHz/2 + 3HzHzHz )

	 D(x) is the Cx|bbb><bbb|
	 (1/8 - Hz/4 + HzHz/2 - HzHzHz )

	 Thus - a 25 element basis
	 */
	std::vector<std::string> basis = dset.vpar("basis");
	//
	G.resize(25, 25);
	G.clear();
	// Now fetch the parameters from the different Dataset/Atom parameters
	//
	const double kex = Kex(dset);
	double pb = Pb(dset);


	if (_multipleTemperatures) //This has gotten back to front somewhere between the rate types
	    pb=1-pb;


	//
	std::string Nucl = dset._nucleus;
	double gammaS = 0.;
	double gammaI = 0.;
	double delta_csa = 0;
	double r_is = 0.;
	double hbar = 6.626075e-34 / (2 * DFH_PI);
	double mu0 = 4. * DFH_PI * 1e-7;
	if (Nucl != "C") {
		std::cerr << " The basis set C_CH3_25 requires that the nucleus is of type C13\n";
		std::cerr << " in Dataset " << dset._id << "\n";
		std::cerr << " Function: FreePrecess_C_CH3_25()\n;" << std::endl;
		Abort(1);
	} else {
		gammaS = _gammaC;
		r_is = 1.091E-10; // Acta Cryst (1988) C44 p. 439-443
		// This value is not explicitly used - we use
		// P2()/r3 from Mittermaier&Kay, J Biomol NMR 2002, 23, p 35-45
		// P2(cos(theta))= -0.228
	}
	dset._gamma = gammaS;
	//
	std::string cn = to_lower_copy(dset.vpar("couplednucleus")[1]);
	if (cn.find('h') < cn.npos) {
		gammaI = _gammaH;
	} else {
		std::cerr << " The basis set C_CH3_25 requires that the coupled nucleus is\n";
		std::cerr << " of type H1\n";
		std::cerr << " in Dataset " << dset._id << "\n";
		std::cerr << " Function: FreePrecess_C_CH3_25()\n;" << std::endl;
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
	RequiredParam.push_back("R1iph" + marker); // R(Cz)
	RequiredParam.push_back("R1aph" + marker); // R(2CzHz)
	//RequiredParam.push_back("R1aph"); // R(2CzHz)
	RequiredParam.push_back("Rex" + temperatureMarker);
	RequiredParam.push_back("tauMe" + temperatureMarker);
	RequiredParam.push_back("S2tc" + temperatureMarker);
	RequiredParam.push_back("CSA");
	RequiredParam.push_back("Omega" + temperatureMarker);
	RequiredParam.push_back("DeltaJ");
	RequiredParam.push_back("DeltaS2tc" + temperatureMarker);
	for (unsigned int i = 0; i < RequiredParam.size(); i++) {
		if (!(BResolveParam(LocalParam[Atom], RequiredParam[i]))) {
			std::cerr << " The parameters:" << RequiredParam[i] << " is required by the basisset";
			std::cerr << " " << basis[1] << "\n but is not provided for atom " << AtomNumber2AtomName(Atom) << "\n";
			std::cerr << " please provide in the LocalParameter set\n";
			std::cerr << " Function FreePrecess_C_CH3_25()\n" << std::endl;
			Abort(1);
		}
	}
	if (!(BResolveParam(LocalParam[Atom], "JIS") || BResolveParam(LocalParam[Atom], "JIS" + marker))) {
		std::cerr << " The parameter JIS or JIS" << marker << " is required by the basisset";
		std::cerr << " " << basis[1] << "\n but is not provided for atom " << AtomNumber2AtomName(Atom) << "\n";
		std::cerr << " please provide one of these parameters in the LocalParameter set\n";
		std::cerr << " Function FreePrecess_C_CH3_25()\n" << std::endl;
		Abort(1);
	}
	double R1iph = ResolveParam(LocalParam[Atom], "R1iph" + marker);
	double R1aph = ResolveParam(LocalParam[Atom], "R1aph" + marker);
	//double R1aph=ResolveParam(LocalParam[Atom],"R1aph");
	double Rex = ResolveParam(LocalParam[Atom], "Rex"+ temperatureMarker);
	double tauMe = ResolveParam(LocalParam[Atom], "tauMe"+ temperatureMarker);
	double CSA = ResolveParam(LocalParam[Atom], "CSA");
	double Omega = ((gammaS) / _gammaH) * dset._sfrq * 2 * DFH_PI * (ResolveParam(LocalParam[Atom], "Omega" + temperatureMarker) - dset._xcar);
	double DeltaO = ((gammaS) / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, Atom);
	double JIS = 0.;
	if (BResolveParam(LocalParam[Atom], "JIS" + marker)) {
		JIS = ResolveParam(LocalParam[Atom], "JIS" + marker);
	} else {
		JIS = ResolveParam(LocalParam[Atom], "JIS");
	}
	double DeltaJ = ResolveParam(LocalParam[Atom], "DeltaJ");


	double DeltaS2tc = ResolveParam(LocalParam[Atom], "DeltaS2tc");
    double S2tc = ResolveParam(LocalParam[Atom], "S2tc");

    //double DeltaS2tc= ResolveParam(LocalParam[0], "DeltaS2tc"+ temperatureMarker);
    //double S2tc = ResolveParam(LocalParam[0], "S2tc"+ temperatureMarker);



	//
	double S2tcA = S2tc; //*S2tc;
	double S2tcB = (S2tc + DeltaS2tc); //*(S2tc+DeltaS2tc);
	//
	double B0 = dset._sfrq * 1E6 * 2. * DFH_PI / _gammaH;
	//
	// Dipolar-dipole relaxation.
	//double DD_CH=1.e-9*pow(mu0*gammaI*gammaS*hbar/(4.*DFH_PI*pow(r_is,3)),2);
	//DD_CH=DD_CH*0.840831; //using H-C-H angle of 107.8 Acta Cryst (1988) C44 p. 439-443
	// --> C-C-H angle is 111.097 and P2/P2 is 0.84

	//            /--- the 9 removes the previously assumed 109.5 degrees
	double DD_CH = 9 * 1.e-9 * pow(mu0 * gammaI * gammaS * hbar / (4. * DFH_PI * pow(1e-10, 3)), 2) * (-0.228) * (-0.228);
	//
	double DD = 1e-9 * pow(mu0 * gammaI * gammaI * hbar / (4. * DFH_PI * pow(1.76e-10, 3)), 2);
	//
	// Four lines goes as  S2tm [ 1./5.  1./45. 1./45 1./5 ]
	//
	// Dipole-csa cross relaxation (J. biomol, NMR, tugarinov, .. _30_, p 397-406 (2004).
	// The factor of '2' is added and converted to some-thing-wierd units
	double D_CSA = 1. * (-(4. / 15.) * (-gammaS * B0) * // omegaC
			gammaI * gammaS * hbar * (-0.228) * // P2(cos(theta))/r^3.
			CSA * // CSA is in ppm
			1e8); // (1e-6*1e30*1e-9*1e-7)= (ppm, AA->meter, ns->s, stupid(CGS unit) -> nice metric units)
	// four lines goes as S2tm [ -3, -1, 1, 3 ];

	//
	// CSA-CSA relaxation
	double CSA_C = 1e-9 * pow(CSA * 1e-6 * gammaS * B0, 2.);
	// four lines goes as 4./45 [ S2tm + 3/4 S2tm/{1+(W*tc)^2} ] ~ 4./45 [ S2tm ] (x) 4
	//
	bool estimate_etaZ = false;
	for (unsigned int i = 0; i < basis.size(); i++) {
		if (to_lower_copy(basis[i]) == "estimate_etaz") {
			estimate_etaZ = true;
			break;
		}
	}
	//
	// Set up the spin-flip matrix
	// methyl spin goes as (3./10)*(sin(2*beta)^2+sin(beta)^4)
	double RHz = (R1aph - R1iph) + 0.36264 * DD * tauMe; // assume same for site A and B
	ublas::matrix<double> SF;
	SF.resize(4, 4);
	//
	SF(0, 0) = 3. / 2.;
	SF(1, 0) = -3. / 2.;
	SF(2, 0) = 0.;
	SF(3, 0) = 0.;
	//
	SF(0, 1) = -1. / 2.;
	SF(1, 1) = 3. / 2.;
	SF(2, 1) = -1.0;
	SF(3, 1) = 0.;

	SF(0, 2) = 0.;
	SF(1, 2) = -1.;
	SF(2, 2) = 3. / 2.;
	SF(3, 2) = -1. / 2.;

	SF(0, 3) = 0.;
	SF(1, 3) = 0.;
	SF(2, 3) = -3. / 2.;
	SF(3, 3) = 3. / 2.;
	//
	//
	// Auto relaxations (A)
	G(1, 1) = S2tcA * (DD_CH / 5. + D_CSA * (-3) + CSA_C * 4. / 45.);
	G(2, 2) = S2tcA * (DD_CH / 5. + D_CSA * (-3) + CSA_C * 4. / 45.);
	G(3, 3) = G(3, 3) + R1iph;
	//
	// Auto relaxations (B)
	G(4, 4) = S2tcA * (DD_CH / 45. + D_CSA * (-1) + CSA_C * 4. / 45.);
	G(5, 5) = S2tcA * (DD_CH / 45. + D_CSA * (-1) + CSA_C * 4. / 45.);
	G(6, 6) = G(6, 6) + R1iph;
	//
	// Auto relaxations (C)
	G(7, 7) = S2tcA * (DD_CH / 45. + D_CSA * (1) + CSA_C * 4. / 45.);
	G(8, 8) = S2tcA * (DD_CH / 45. + D_CSA * (1) + CSA_C * 4. / 45.);
	G(9, 9) = G(9, 9) + R1iph;
	//
	// Auto relaxations (D)
	G(10, 10) = S2tcA * (DD_CH / 5. + D_CSA * (3) + CSA_C * 4. / 45.);
	G(11, 11) = S2tcA * (DD_CH / 5. + D_CSA * (3) + CSA_C * 4. / 45.);
	G(12, 12) = G(12, 12) + R1iph;
	//
	// Chemical shift and J coupling evolutions
	// (A)
	G(1, 2) = (-JIS * 3.) * DFH_PI + Omega;
	G(2, 1) = -(-JIS * 3.) * DFH_PI - Omega;
	// (B)
	G(4, 5) = (-JIS * 1.) * DFH_PI + Omega;
	G(5, 4) = -(-JIS * 1.) * DFH_PI - Omega;
	// (C)
	G(7, 8) = (JIS * 1.) * DFH_PI + Omega;
	G(8, 7) = -(JIS * 1.) * DFH_PI - Omega;
	// (D)
	G(10, 11) = (JIS * 3.) * DFH_PI + Omega;
	G(11, 10) = -(JIS * 3.) * DFH_PI - Omega;
	//
	//////////////////////// Minor state ////////////////////////////////////////
	// We assume that R1iph(major state) ~ R1iph(minor state)
	//
	// Auto relaxations (A)
	G(13, 13) = (S2tcB) * (DD_CH / 5. + D_CSA * (-3) + CSA_C * 4. / 45.);
	G(14, 14) = (S2tcB) * (DD_CH / 5. + D_CSA * (-3) + CSA_C * 4. / 45.);
	G(15, 15) = G(15, 15) + R1iph;
	//
	// Auto relaxations (B)
	G(16, 16) = (S2tcB) * (DD_CH / 45. + D_CSA * (-1) + CSA_C * 4. / 45.);
	G(17, 17) = (S2tcB) * (DD_CH / 45. + D_CSA * (-1) + CSA_C * 4. / 45.);
	G(18, 18) = G(18, 18) + R1iph;
	//
	// Auto relaxations (C)
	G(19, 19) = (S2tcB) * (DD_CH / 45. + D_CSA * (1) + CSA_C * 4. / 45.);
	G(20, 20) = (S2tcB) * (DD_CH / 45. + D_CSA * (1) + CSA_C * 4. / 45.);
	G(21, 21) = G(21, 21) + R1iph;
	//
	// Auto relaxations (D)
	G(22, 22) = (S2tcB) * (DD_CH / 5. + D_CSA * (3) + CSA_C * 4. / 45.);
	G(23, 23) = (S2tcB) * (DD_CH / 5. + D_CSA * (3) + CSA_C * 4. / 45.);
	G(24, 24) = G(24, 24) + R1iph;
	//
	// Chemical shift and J coupling evolutions
	// (A)
	G(13, 14) = (-(JIS + DeltaJ) * 3.) * DFH_PI + (Omega + DeltaO);
	G(14, 13) = -(-(JIS + DeltaJ) * 3.) * DFH_PI - (Omega + DeltaO);
	// (B)
	G(16, 17) = (-(JIS + DeltaJ) * 1.) * DFH_PI + (Omega + DeltaO);
	G(17, 16) = -(-(JIS + DeltaJ) * 1.) * DFH_PI - (Omega + DeltaO);
	// (C)
	G(19, 20) = ((JIS + DeltaJ) * 1.) * DFH_PI + (Omega + DeltaO);
	G(20, 19) = -((JIS + DeltaJ) * 1.) * DFH_PI - (Omega + DeltaO);
	// (D)
	G(22, 23) = ((JIS + DeltaJ) * 3.) * DFH_PI + (Omega + DeltaO);
	G(23, 22) = -((JIS + DeltaJ) * 3.) * DFH_PI - (Omega + DeltaO);
	//
	// Store the spin-flip matrix on the full matrix
	for (unsigned int c0 = 0; c0 < 4; c0++) {
		for (unsigned int c1 = 0; c1 < 4; c1++) {
			for (unsigned int i = 0; i < 3; i++) { // x,y,z
				//
				// Ground state
				G(c0 * 3 + i + 1, c1 * 3 + i + 1) += RHz * SF(c0, c1);
				//
				//scale the excited state with J(0)[excited_state]/J(0)[ground_state]
				G(c0 * 3 + i + 13, c1 * 3 + i + 13) += RHz * SF(c0, c1) * S2tcB / S2tcA;
			}
		}
	}
	//
	// Towards infinity
	G(3, 0) = -(1 - pb) * (G(3, 3)) * gammaS / (_gammaH * 2.);
	G(6, 0) = -3. * (1 - pb) * (G(6, 6)) * gammaS / (_gammaH * 2.);
	G(9, 0) = -3. * (1 - pb) * (G(9, 9)) * gammaS / (_gammaH * 2.);
	G(12, 0) = -(1 - pb) * (G(12, 12)) * gammaS / (_gammaH * 2.);
	//
	G(15, 0) = -(pb) * (G(15, 15)) * gammaS / (_gammaH * 2.);
	G(18, 0) = -3. * (pb) * (G(18, 18)) * gammaS / (_gammaH * 2.);
	G(21, 0) = -3. * (pb) * (G(21, 21)) * gammaS / (_gammaH * 2.);
	G(24, 0) = -(pb) * (G(24, 24)) * gammaS / (_gammaH * 2.);
	//
	// kex additions
	for (unsigned int i = 0; i < 12; i++) {
		G(i + 1, i + 1) = G(i + 1, i + 1) + kex * (pb);
		G(i + 13, i + 1) = G(i + 13, i + 1) - kex * (pb);
		//
		G(i + 13, i + 13) = G(i + 13, i + 13) + kex * (1-pb);
		G(i + 1, i + 13) = G(i + 1, i + 13) - kex * (1-pb);
	}
	//
	// Add the offset tauMe ( --> tau_me in nanoseconds)
	for (unsigned int i = 0; i < 4; i++) { //lines
		for (unsigned int j = 0; j < 2; j++) { //x,y
			if (i == 1 || i == 2) { //inner lines
				// contribution goes as: ((99/80)*DHH^2+(87/40)DCH^2)sin(beta)^4 +(87/40)DCH^2 sin(2*beta)^2
				//G(i*3+j+1, i*3+j+1, G(i*3+j+1, i*3+j+1) +tauMe*DD*3.0660 );
				//G(i*3+j+13,i*3+j+13,G(i*3+j+13,i*3+j+13)+tauMe*DD*3.0660 );
				G(i * 3 + j + 1, i * 3 + j + 1) = G(i * 3 + j + 1, i * 3 + j + 1) + tauMe * DD * 4.1523;
				G(i * 3 + j + 13, i * 3 + j + 13) = G(i * 3 + j + 13, i * 3 + j + 13) + tauMe * DD * 4.1423;
			} else {
				// contribution goes as: (9/8)DCH^2*sin(2*beta)^2 + {(27/16)DHH^2+ (9/8)DCH^2}sin(beta)^4
				//G(i*3+j+1, i*3+j+1, G(i*3+j+1, i*3+j+1) +tauMe*DD*1.7767 );
				//G(i*3+j+13,i*3+j+13,G(i*3+j+13,i*3+j+13)+tauMe*DD*1.7767 );
				G(i * 3 + j + 1, i * 3 + j + 1) = G(i * 3 + j + 1, i * 3 + j + 1) + tauMe * DD * 3.19517;
				G(i * 3 + j + 13, i * 3 + j + 13) = G(i * 3 + j + 13, i * 3 + j + 13) + tauMe * DD * 3.19517;
			}
		}
	}
	//
	// Add the exchange B0*B0*Rex
	for (unsigned int i = 0; i < 4; i++) { //lines
		for (unsigned int j = 0; j < 2; j++) { //x,y
			G(i * 3 + j + 1, i * 3 + j + 1) = G(i * 3 + j + 1, i * 3 + j + 1) + Rex * B0 * B0;
			G(i * 3 + j + 13, i * 3 + j + 13) = G(i * 3 + j + 13, i * 3 + j + 13) + 0.0;
		}
	}
	/*
	 std::cerr<<"-------------------------------------------"<<std::endl;
	 std::cerr<<" SF     \t"<< SF(0,0)*RHz<<std::endl;
	 std::cerr<<" DD-DD  \t"<< S2tcA*( DD_CH/5 )<<std::endl;
	 std::cerr<<" DD-CSA \t"<< S2tcA*( D_CSA*(-3) )<<std::endl;
	 std::cerr<<" CSA-CSA\t"<< S2tcA*( CSA_C*4./45. )<<std::endl;
	 std::cerr<<" tauMe  \t"<< tauMe*DD*3.19517 <<std::endl;
	 std::cerr<<" Rex    \t"<< Rex*B0*B0<<std::endl;
	 std::cerr<<" Sum    \t"
	 <<G(1,1)  -kex*pb<<"\t"
	 <<G(4,4)  -kex*pb<<"\t"
	 <<G(7,7)  -kex*pb<<"\t"
	 <<G(10,10)-kex*pb<<"\t"
	 <<std::endl;
	 std::cerr<<"-------------------------------------------"<<std::endl;
	 */
	//
}
