#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>

namespace ublas = boost::numeric::ublas;

void Catia::FreePrecess(ublas::matrix<double>& G, Dataset& dset, int Atom) {
	//
	// Atom is globalatom number
	//
	std::vector<std::string> basis = dset.vpar("basis");
	if (basis[1] == "tratr_13") {
		FreePrecess_TrATr_13(G, dset, Atom);
	} else if (basis[1] == "iphaph_13") {
		FreePrecess_IphAph_13(G, dset, Atom);
	} else if (basis[1] == "iphaph_13_deltar2") {
		FreePrecess_IphAph_13_deltaR2(G, dset, Atom);
	} else if (basis[1] == "iphaph_13_dr2") {
		FreePrecess_IphAph_13_deltaR2_simple(G, dset, Atom);
	} else if (basis[1] == "iph_7") {
		FreePrecess_Iph_7(G, dset, Atom);
	} else if (basis[1] == "aph_7") {
		FreePrecess_Aph_7(G, dset, Atom);
	} else if (basis[1] == "tr_7") {
		FreePrecess_Tr_7(G, dset, Atom);
	} else if (basis[1] == "c_ch3_25") {
		FreePrecess_C_CH3_25(G, dset, Atom);
	} else if (basis[1] == "n_nh_13") {
		FreePrecess_N_NH_13(G, dset, Atom);
	} else if (basis[1] == "n_7") {
		FreePrecess_N_7(G, dset, Atom);
	} else if (basis[1] == "3st_iph_10") {
		FreePrecess_3st_Iph_10(G, dset, Atom);
	} else {
		std::cerr << " The basis:" << basis[1] << " is not defined \n";
		std::cerr << " Function .FreePrecess();\n";
		Abort(1);
	}
	//G.print("%10.3f");
	//exit(10);

}

void Catia::FreePrecess(ublas::matrix<std::complex<double> >& G, Dataset& dset, int Atom) {
	//
	// Atom is globalatom number
	//
	std::vector<std::string> basis = dset.vpar("basis");
	if (basis[1] == "fast_2") {
		FreePrecess_Fast_2(G, dset, Atom);
	} else if (basis[1] == "fast_3") {
			FreePrecess_Fast_3(G, dset, Atom);
	} else {
		std::cerr << " The basis:" << basis[1] << " is not defined \n";
		std::cerr << " Function .FreePrecess();\n";
		Abort(1);
	}
	//G.print("%10.3f");
	//exit(10);
}
