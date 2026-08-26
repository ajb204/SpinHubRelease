#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>

/*
 (1) The basis set is extracted from the dataset dset
 (2) The coherence with name 'coh' is stored in the vector *Eq

 Flemming September 2007
 */

namespace ublas = boost::numeric::ublas;

void Catia::CalcCoherence(std::string coh, ublas::vector<double>& Eq, Dataset dset) {
	std::vector<std::string> basis = dset.vpar("basis");
	double pa = 0.;
	double pb = 0.;
	double pc = 0.;

	//
	pb = Pb(dset);
	//
	if (basis[1] == "tratr_13" || basis[1] == "n_nh_13") {
		//check
		for (unsigned int i = 0; i < 13; i++) {
			Eq[i] = 0.;
		}
		if (coh == "Sz") {
			Eq[0] = 1.;
			Eq[3] = (1 - pb) / 2.;
			Eq[6] = (1 - pb) / 2.;
			Eq[9] = pb / 2.;
			Eq[12] = pb / 2.;
		} else if (coh == "Sy") {
			Eq[0] = 1.;
			Eq[2] = (1 - pb) / 2.;
			Eq[5] = (1 - pb) / 2.;
			Eq[8] = pb / 2.;
			Eq[11] = pb / 2.;
		} else if (coh == "2IzSz") {
			Eq[0] = 1.;
			Eq[3] = -(1 - pb) / 2.;
			Eq[6] = (1 - pb) / 2.;
			Eq[9] = -pb / 2.;
			Eq[12] = pb / 2.;
		} else if (coh == "2IzSy") {
			Eq[0] = 1.;
			Eq[2] = -(1 - pb) / 2.;
			Eq[5] = (1 - pb) / 2.;
			Eq[8] = -pb / 2.;
			Eq[11] = pb / 2.;
		} else if (coh == "Tr(y)") {
			Eq[0] = 1.;
			Eq[2] = 1 - pb;
			Eq[8] = pb;
		} else if (coh == "ATr(y)") {
			Eq[0] = 1.;
			Eq[5] = 1 - pb;
			Eq[11] = pb;
		} else if (coh == "Tr(z)") {
			Eq[0] = 1.;
			Eq[3] = 1 - pb;
			Eq[9] = pb;
		} else if (coh == "ATr(z)") {
			Eq[0] = 1.;
			Eq[6] = 1 - pb;
			Eq[12] = pb;
		} else {
			std::cerr << " Coherence " << coh << " is not available in the basis " << basis[1] << "" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "iphaph_13" || basis[1] == "iphaph_13_deltar2" || basis[1] == "iphaph_13_dr2" ) {
		//check
		for (unsigned int i = 0; i < 13; i++) {
			Eq[i] = 0.;
		}
		if (coh == "Sz") {
			Eq[0] = 1.;
			Eq[3] = (1 - pb);
			Eq[9] = pb;
		} else if (coh == "2IzSz") {
			Eq[0] = 1.;
			Eq[6] = (1 - pb);
			Eq[12] = pb;
		} else if (coh == "2IzSy") {
			Eq[0] = 1.;
			Eq[5] = (1 - pb);
			Eq[11] = pb;
		} else {
			std::cerr << " Coherence " << coh << " is not available in the basis " << basis[1] << "" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "iph_7" || basis[1] == "n_7") {
		//check
		for (unsigned int i = 0; i < 7; i++) {
			Eq[i] = 0.;
		}
		if (coh == "Sz") {
			Eq[0] = 1.;
			Eq[3] = (1 - pb);
			Eq[6] = pb;
		} else if (coh == "Sy") {
			Eq[0] = 1.;
			Eq[2] = (1 - pb);
			Eq[5] = pb;
		} else if (coh == "Sx") {
			Eq[0] = 1.;
			Eq[1] = (1 - pb);
			Eq[4] = pb;
		} else if (coh == "Tr(y)") {
			Eq[0] = 1.;
			Eq[2] = (1 - pb);
			Eq[5] = pb;
		} else {
			std::cerr << " Coherence " << coh << " is not available in the basis " << basis[1] << "" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "aph_7") {
		//check
		for (unsigned int i = 0; i < 7; i++) {
			Eq[i] = 0.;
		}
		if (coh == "2IzSz") {
			Eq[0] = 1.;
			Eq[3] = (1 - pb);
			Eq[6] = pb;
		} else if (coh == "2IzSy") {
			Eq[0] = 1.;
			Eq[2] = (1 - pb);
			Eq[5] = pb;
		} else if (coh == "Tr(y)") {
			Eq[0] = 1.;
			Eq[2] = (1 - pb);
			Eq[5] = pb;
		} else if (coh == "2IzSx") {
			Eq[0] = 1.;
			Eq[1] = (1 - pb);
			Eq[4] = pb;
		} else {
			std::cerr << " Coherence " << coh << " is not available in the basis " << basis[1] << "" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "tr_7") {
		//check
		for (unsigned int i = 0; i < 7; i++) {
			Eq[i] = 0.;
		}
		if (coh == "2IzSz") {
			Eq[0] = 1.;
			Eq[3] = -(1 - pb) / 2;
			Eq[6] = -pb / 2;
		} else if (coh == "2IzSy") {
			Eq[0] = 1.;
			Eq[2] = -(1 - pb) / 2;
			Eq[5] = -pb / 2;
		} else if (coh == "2IzSx") {
			Eq[0] = 1.;
			Eq[1] = -(1 - pb) / 2;
			Eq[4] = -pb / 2;
		} else if (coh == "Sz") {
			Eq[0] = 1.;
			Eq[3] = (1 - pb) / 2;
			Eq[6] = pb / 2;
		} else if (coh == "Sy") {
			Eq[0] = 1.;
			Eq[2] = (1 - pb) / 2;
			Eq[5] = pb / 2;
		} else if (coh == "Sx") {
			Eq[0] = 1.;
			Eq[1] = (1 - pb) / 2;
			Eq[4] = pb / 2;
		} else if (coh == "Tr(y)") {
			Eq[0] = 1.;
			Eq[2] = (1 - pb);
			Eq[5] = pb;
		} else {
			std::cerr << " Coherence " << coh << " is not available in the basis " << basis[1] << "" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "c_ch3_25") {
		//check
		for (unsigned int i = 0; i < 25; i++) {
			Eq[i] = 0.;
		}
		Eq[0] = 1.;
		if (coh == "CH3_A(x)") {
			Eq[1] = (1 - pb);
			Eq[1 + 12] = pb;
		} else if (coh == "CH3_A(y)") {
			Eq[2] = (1 - pb);
			Eq[2 + 12] = pb;
		} else if (coh == "CH3_A(z)") {
			Eq[3] = (1 - pb);
			Eq[3 + 12] = pb;
		} else if (coh == "CH3_B(x)") {
			Eq[4] = (1 - pb);
			Eq[4 + 12] = pb;
		} else if (coh == "CH3_B(y)") {
			Eq[5] = (1 - pb);
			Eq[5 + 12] = pb;
		} else if (coh == "CH3_B(z)") {
			Eq[6] = (1 - pb);
			Eq[6 + 12] = pb;
		} else if (coh == "CH3_C(x)") {
			Eq[7] = (1 - pb);
			Eq[7 + 12] = pb;
		} else if (coh == "CH3_C(y)") {
			Eq[8] = (1 - pb);
			Eq[8 + 12] = pb;
		} else if (coh == "CH3_C(z)") {
			Eq[9] = (1 - pb);
			Eq[9 + 12] = pb;
		} else if (coh == "CH3_D(x)") {
			Eq[10] = (1 - pb);
			Eq[10 + 12] = pb;
		} else if (coh == "CH3_D(y)") {
			Eq[11] = (1 - pb);
			Eq[11 + 12] = pb;
		} else if (coh == "CH3_D(z)") {
			Eq[12] = (1 - pb);
			Eq[12 + 12] = pb;
		} else {
			std::cerr << " The coherence " << coh << " has not been defined yet" << std::endl;
			std::cerr << " in the basis " << basis[1] << "" << std::endl;
			std::cerr << " Function .CalcCoherence();" << std::endl;
			Abort(1);
		}
	} else {
		std::cerr << " The basis " << basis[1] << " is not defined" << std::endl;
		std::cerr << " Function .CalcCoherence()";
		std::cerr << std::endl;
		Abort(1);
	}
}

void Catia::CalcCoherence_3st(std::string coh, ublas::vector<double>& Eq, Dataset dset) {
	std::vector<std::string> basis = dset.vpar("basis");
	double pa = 0.;
	double pb = 0.;
	double pc = 0.;

	bool hasKexP3st = (BResolveParam(GlobalParam, "kex_ab")
			&& BResolveParam(GlobalParam, "kex_ac")
			&& BResolveParam(GlobalParam, "kex_bc")
			&& BResolveParam(GlobalParam, "pb")
			&& BResolveParam(GlobalParam, "pc"));

	bool hasDeltaHDeltaS3st = (BResolveParam(GlobalParam, "deltaHb")
			&& BResolveParam(GlobalParam, "deltaHc")
			&& BResolveParam(GlobalParam, "deltaHab")
			&& BResolveParam(GlobalParam, "deltaHac")
			&& BResolveParam(GlobalParam, "deltaHbc")
			&& BResolveParam(GlobalParam, "deltaSb")
			&& BResolveParam(GlobalParam, "deltaSc")
			&& BResolveParam(GlobalParam, "deltaSab")
			&& BResolveParam(GlobalParam, "deltaSac")
			&& BResolveParam(GlobalParam, "deltaSbc"));

	if (hasKexP3st || hasDeltaHDeltaS3st) {
		pb = Pb_3st(dset);
		pc = Pc_3st(dset);
		pa = 1.0 - pb - pc;
	} else {
		std::cerr << " No unique global exchange parameters have been provided." << std::endl;
		std::cerr << " Please provide either {kex,pb} or {kab,kba} or {deltaHa,deltaSa,deltaHb,deltaSb}" << std::endl;
		std::cerr << " or {kex_ab,kex_ac,_kex_bc,pb,pc}" << std::endl;
		std::cerr << " in the global parameter set." << std::endl;
		std::cerr << " Function: .CalcCoherence()\n;" << std::endl;
		Abort(1);
	}
	//
	if (basis[1] == "3st_iph_10") {
		//check
		for (unsigned int i = 0; i < 10; i++) {
			Eq[i] = 0.;
		}
		if (coh == "Sz") {
			Eq[0] = 1.;
			Eq[3] = pa;
			Eq[6] = pb;
			Eq[9] = pc;
		} else if (coh == "Sy") {
			Eq[0] = 1.;
			Eq[2] = pa;
			Eq[5] = pb;
			Eq[8] = pc;
		} else if (coh == "Sx") {
			Eq[0] = 1.;
			Eq[1] = pa;
			Eq[4] = pb;
			Eq[7] = pc;
		} else if (coh == "Tr(y)") {
			Eq[0] = 1.;
			Eq[2] = pa;
			Eq[5] = pb;
			Eq[8] = pc;
		} else {
			std::cerr << " Coherence " << coh << " is not available in the basis " << basis[1] << "" << std::endl;
			Abort(1);
		}
	} else {
		std::cerr << " The basis " << basis[1] << " is not defined" << std::endl;
		std::cerr << " Function .CalcCoherence()";
		std::cerr << std::endl;
		Abort(1);
	}
}
