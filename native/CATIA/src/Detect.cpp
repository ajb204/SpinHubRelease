#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
/*
 Detect a given coherence 'coherence' from the magnetization vector *Sigma.
 The site (\in 0,1) refers to the chemically exchanging sites (So far only 0 or 1)
 for site A or B, respectively.

 If the coherence does not exists, then a zero is returned, e.g., if one
 detects 2IzSz in the basis of inphase-coherences (E,Sx(A),Sy(A),Sz(A)..)
 a zero is returned.

 Flemming, September 2007

 Modified on May 13 2009 by Flemming
 Check for all basis sets, that the coherence we are asking for is actually there,
 thus we should either get a 'return' statement, or an 'abort'.

 */
double Catia::Detect(std::string coherence, int site, Dataset dset, ublas::vector<double>& Sigma) {
	std::vector<std::string> basis = dset.vpar("basis");
	if (basis[1] == "tratr_13" || basis[1] == "n_nh_13") {
		if (coherence == "Sz") {
			if (site == 0) {
				return Sigma[3] + Sigma[6];
			} else if (site == 1) {
				return Sigma[3 + 6] + Sigma[6 + 6];
			} else {
				std::cerr << " Site must be either 0 or 1\n";
				std::cerr << " Function .Detect();\n";
				std::cerr << std::endl;
				Abort(1);
			}
		} else if (coherence == "2IzSz") {
			if (site == 0) {
				return -Sigma[3] + Sigma[6];
			} else if (site == 1) {
				return -Sigma[3 + 6] + Sigma[6 + 6];
			} else {
				std::cerr << " Site must be either 0 or 1\n";
				std::cerr << " Function .Detect();\n";
				std::cerr << std::endl;
				Abort(1);
			}
		} else {
			std::cerr << " Coherence: " << coherence << " is not available in the basis:\n";
			std::cerr << " " << basis[1] << std::endl;
			std::cerr << " Function .Detect()\n" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "iphaph_13" || basis[1] == "iphaph_13_deltar2" || basis[1] == "iphaph_13_dr2" ) {
		if (coherence == "Sz") {
			if (site == 0) {
				return Sigma[3];
			} else if (site == 1) {
				return Sigma[3 + 6];
			} else {
				std::cerr << " Site must be either 0 or 1\n";
				std::cerr << " Function .Detect();\n";
				std::cerr << std::endl;
				Abort(1);
			}
		} else if (coherence == "2IzSz") {
			if (site == 0) {
				return Sigma[6];
			} else if (site == 1) {
				return Sigma[6 + 6];
			} else {
				std::cerr << " Site must be either 0 or 1\n";
				std::cerr << " Function .Detect();\n";
				std::cerr << std::endl;
				Abort(1);
			}
		} else if (coherence == "2IzSy") {
			if (site == 0) {
				return Sigma[5];
			} else if (site == 1) {
				return Sigma[5 + 6];
			} else {
				std::cerr << " Site must be either 0 or 1\n";
				std::cerr << " Function .Detect();\n";
				std::cerr << std::endl;
				Abort(1);
			}
		} else {
			std::cerr << " Coherence: " << coherence << " is not available in the basis:\n";
			std::cerr << " " << basis[1] << std::endl;
			std::cerr << " Function .Detect()\n" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "iph_7" || basis[1] == "n_7") {
		if (coherence == "Sz") {
			if (site == 0) {
				return Sigma[3];
			} else if (site == 1) {
				return Sigma[3 + 3];
			} else {
				return 0.;
			}
		} else if (coherence == "Sy") {
			if (site == 0) {
				return Sigma[2];
			} else if (site == 1) {
				return Sigma[2 + 3];
			} else {
				return 0.;
			}
		} else if (coherence == "Sx") {
			if (site == 0) {
				return Sigma[1];
			} else if (site == 1) {
				return Sigma[1 + 3];
			} else {
				return 0.;
			}
		} else if (coherence == "2IzSz" || coherence == "2IzSy" || coherence == "2IzSx") {
			return 0.;
		} else {
			std::cerr << " Coherence: " << coherence << " is not available in the basis:\n";
			std::cerr << " " << basis[1] << std::endl;
			std::cerr << " Function .Detect()\n" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "3st_iph_10") {
		if (coherence == "Sz") {
			if (site < 3) {
				return Sigma[3 + 3 * site];
			} else {
				return 0.;
			}
		} else if (coherence == "Sy") {
			if (site < 3) {
				return Sigma[2 + 3 * site];
			} else {
				return 0.;
			}
		} else if (coherence == "Sx") {
			if (site < 3) {
				return Sigma[1 + 3 * site];
			} else {
				return 0.;
			}
		} else if (coherence == "2IzSz" || coherence == "2IzSy" || coherence == "2IzSx") {
			return 0.;
		} else {
			std::cerr << " Coherence: " << coherence << " is not available in the basis:\n";
			std::cerr << " " << basis[1] << std::endl;
			std::cerr << " Function .Detect()\n" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "aph_7") {
		if (coherence == "2IzSz") {
			if (site == 0) {
				return Sigma[3];
			} else if (site == 1) {
				return Sigma[3 + 3];
			} else {
				std::cerr << " Site must be either 0 or 1\n";
				std::cerr << " Function .Detect();\n";
				std::cerr << std::endl;
				Abort(1);
			}
		} else if (coherence == "Sz") {
			return 0.;
		} else {
			std::cerr << " Coherence: " << coherence << " is not available in the basis:\n";
			std::cerr << " " << basis[1] << std::endl;
			std::cerr << " Function .Detect()\n" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "tr_7") {
		if (coherence == "2IzSz") {
			if (site == 0) {
				return -Sigma[3];
			} else if (site == 1) {
				return -Sigma[3 + 3];
			} else {
				return 0.;
			}
		} else if (coherence == "Sz") {
			if (site == 0) {
				return Sigma[3];
			} else if (site == 1) {
				return Sigma[3 + 3];
			} else {
				return 0.;
			}
		} else {
			std::cerr << " Coherence: " << coherence << " is not available in the basis:\n";
			std::cerr << " " << basis[1] << std::endl;
			std::cerr << " Function .Detect()\n" << std::endl;
			Abort(1);
		}
	} else if (basis[1] == "c_ch3_25") {
		if (site == 0 || site == 1) {
			if (coherence == "CH3_A(x)") {
				return Sigma[1 + site * 12];
			}
			if (coherence == "CH3_A(y)") {
				return Sigma[2 + site * 12];
			}
			if (coherence == "CH3_A(z)") {
				return Sigma[3 + site * 12];
			}

			if (coherence == "CH3_B(x)") {
				return Sigma[4 + site * 12];
			}
			if (coherence == "CH3_B(y)") {
				return Sigma[5 + site * 12];
			}
			if (coherence == "CH3_B(z)") {
				return Sigma[6 + site * 12];
			}

			if (coherence == "CH3_C(x)") {
				return Sigma[7 + site * 12];
			}
			if (coherence == "CH3_C(y)") {
				return Sigma[8 + site * 12];
			}
			if (coherence == "CH3_C(z)") {
				return Sigma[9 + site * 12];
			}

			if (coherence == "CH3_D(x)") {
				return Sigma[10 + site * 12];
			}
			if (coherence == "CH3_D(y)") {
				return Sigma[11 + site * 12];
			}
			if (coherence == "CH3_D(z)") {
				return Sigma[12 + site * 12];
			}
		} else {
			std::cerr << " The basis " << basis[1] << " only has two site exchange\n";
			std::cerr << " You tried to detect the signal of site #" << site << "\n";
			std::cerr << " Function: .Detect();\n" << std::endl;
			Abort(1);
		}
	} else {
		std::cerr << " The basis " << basis[1] << " is not defined\n";
		std::cerr << " Function .Detect()";
		std::cerr << std::endl;
		Abort(1);
	}
}
