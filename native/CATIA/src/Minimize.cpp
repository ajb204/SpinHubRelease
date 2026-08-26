#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>

using boost::to_lower_copy;
using boost::trim_copy;

typedef std::pair<std::string, int> pairSI_t;

/*
 For wrapper functions see the bottom

 June 2, 2009:
 Modified by G Bouvignies to implement regular expression as the way to select atoms

 */

double Catia::Minimize(std::vector<int> dset, std::vector<int> aset, std::string opt, bool& conv) {
  /*
	 Returns final X2.
	 minimizes by using the Datasets d and atoms a
	 */
	//Sort the atomset:
	signed int swap = -1;
	for (unsigned int i = 0; i < aset.size(); i++) {
		for (unsigned int j = 0; j < i; j++) {
			if (aset[i] < aset[j]) {
				swap = aset[i];
				aset[i] = aset[j];
				aset[j] = swap;
			}
		}
	}
	//Check that each atom is only defined once!
	for (unsigned int i = 1; i < aset.size(); i++) {
		if (aset[i - 1] == aset[i]) {
			std::cerr << "\n You have defined the atom: " << AtomNumber2AtomName(aset[i]) << " more than once in the vector\n";
			std::cerr << " of atoms to be included in the least-squares minimization\n";
			std::cerr << " Each atom cannot be included more than once\n\n";
			std::cerr << " Function .Minimize();" << std::endl;
			Abort(2);
		}
	}

	std::vector<std::string> ParamName; //StringName for each parameters
	std::vector<double> ParamVal; // ..
	std::vector<int> AtomNumber; // number of atom corresponding to the parameter
	//
	std::string minimizer = "lm";
	int maxiter = 500;
	bool print = true;
	double tol = 1.e-3;

	//Make Ready for NR3
	VecDoub x;
	VecDoub y;
	VecDoub sig;
	//  VecDoub a;
	//  VecBool ai;
	std::vector<double*> a;
	std::vector<int*> ai;
	std::vector<int> aa; //map parameter to atomnumber, (-1 is global)
	std::vector<std::string> an; // name(param)
	//
	int ma;
	//
	double ChiSq;
	//
	// solve the options;
	std::vector<std::string> temp;
	std::vector<std::string> temp1;
	opt = to_lower_copy(trim_copy(opt));
	temp = split(opt, ";");
	for (unsigned int i = 0; i < temp.size(); i++) {
		if (temp[i].find('=') < temp[i].npos) {
			temp1 = split(temp[i], "=");
			if (temp1[0].substr(0, 3) == "min") {
				minimizer = temp1[1];
			} else if (temp1[0].substr(0, 3) == "pri") {
				if (temp1[1].find('n') < temp1[1].npos) {
					print = false;
				} else {
					print = true;
				}
			} else if (temp1[0].substr(0, 4) == "maxi") {
				maxiter = atoi(temp1[1].c_str());
			} else if (temp1[0].substr(0, 3) == "tol") {
				tol = atof(temp1[1].c_str());
			} else {
				std::cerr << " The minimizer option: " << temp1[0] << " is unknown\n";
				std::cerr << " Function .Minimize();\n";
				std::cerr << std::endl;
				Abort(1);
			}
		} else {
			std::cerr << " The minimizer option-field: " << temp[i] << "\n";
			std::cerr << " is not of the type: param=val \n";
			std::cerr << " Function .Minimize();\n";
			std::cerr << std::endl;
			Abort(1);
		}
	}
	temp.clear();
	// Do some checking
	// (1) Do we have data for all atoms?
	//  std::vector<int> av;
	bool found = false;
	for (unsigned int i = 0; i < aset.size(); i++) {
		found = false;
		// is there a Dataset where we have data for atom aset[i]?
		for (unsigned int j = 0; j < dset.size(); j++) {
			if (Datasets[dset[j]]._atomNameToLocalAtomNumber.find(AtomNumber2AtomName(aset[i])) == Datasets[dset[j]]._atomNameToLocalAtomNumber.end()) {
				found = false;
			} else {
				found = true;
				break;
			}
		}
		if (!found) {
			std::cerr << " No experimental data is available for atom number " << aset[i] << " with name " << AtomNumber2AtomName(aset[i]) << "\n";
			std::cerr << std::endl;
			Abort(1);
		}
	}
	Atoms2Fit.clear();
	Atoms2Fit = aset;
	//av.clear();
	//
	// Management of the parameters:
	// Assign ParamName,AtomName,ParamVal to each of the parameters
	// first global parameters:
	std::map<std::string, double>::iterator itDS; //iterator:double<-string
	for (itDS = GlobalParam.begin(); itDS != GlobalParam.end(); ++itDS) {
		ParamName.push_back(itDS->first);
		AtomNumber.push_back(-1);
		ParamVal.push_back(itDS->second);
	}
	//Then local Parameters
	for (unsigned int i = 0; i < Atoms2Fit.size(); i++) {
		//atom a[i]
		for (itDS = LocalParam[Atoms2Fit[i]].begin(); itDS != LocalParam[Atoms2Fit[i]].end(); ++itDS) {
			ParamName.push_back(itDS->first);
			AtomNumber.push_back(Atoms2Fit[i]);
			ParamVal.push_back(itDS->second);
		}
	}
	// Allocate space for NR (and other) minimization functions
	ma = ParamVal.size();
	a.resize(ma);
	ai.resize(ma);
	aa.resize(ma);
	an.resize(ma);
	//
	// The initial parameters
	for (unsigned int i = 0; i < ma; i++) {
		//a[i]=ParamVal[i];
		if (AtomNumber[i] == -1) {
			a[i] = &GlobalParam[ParamName[i]];
			ai[i] = &GlobalParamF[ParamName[i]];
			aa[i] = -1;
			an[i] = ParamName[i];
		} else {
			a[i] = &LocalParam[AtomNumber[i]][ParamName[i]];
			ai[i] = &LocalParamF[AtomNumber[i]][ParamName[i]];
			aa[i] = AtomNumber[i];
			an[i] = ParamName[i];
		}
	}
	//
	// Now set up the X-coordinates.
	// (int)X \in {0,1,2,3,4....}
	// The full evolution matrix should be solved for each X only once,
	// thus, all datapoints in same Dataset and same atom (but with different ncyc)
	// uses the same evolution matrix (and the diagonalization of this, which is only
	// performed once).
	//
	X2dset_atom.clear();
	int ndata = 0;
	std::vector<int> dset_atom(3);
	for (unsigned int i = 0; i < dset.size(); i++) {
		for (unsigned int j = 0; j < Datasets[dset[i]]._localToGlobalAtomIndex.size(); j++) {
			//
			//is this atom in Atoms2Fit ?
			bool UseThisAtom = false;
			for (unsigned int k = 0; k < Atoms2Fit.size(); k++) {
				if (Datasets[dset[i]]._localToGlobalAtomIndex[j] == Atoms2Fit[k]) {
					UseThisAtom = true;
					break;
				}
			}
			if (!(UseThisAtom)) {
				continue;
			}
			for (unsigned int k = 0; k < Datasets[dset[i]].ncyc[j].size(); k++) {
				dset_atom[0] = dset[i];
				dset_atom[1] = Datasets[dset[i]]._localToGlobalAtomIndex[j];
				dset_atom[2] = k;
				X2dset_atom[ndata] = dset_atom;
				//
				ndata++;
			}
		}
	}
	//
	// We need to know number of data before allocating the memory
	x.resize(ndata);
	y.resize(ndata);
	sig.resize(ndata);
	//
	int dc = 0;
	for (unsigned int i = 0; i < dset.size(); i++) {
		for (unsigned int j = 0; j < Datasets[dset[i]]._localToGlobalAtomIndex.size(); j++) {
			bool UseThisAtom = false;
			for (unsigned int k = 0; k < Atoms2Fit.size(); k++) {
				if (Datasets[dset[i]]._localToGlobalAtomIndex[j] == Atoms2Fit[k]) {
					UseThisAtom = true;
					break;
				}
			}
			if (UseThisAtom) {
				for (unsigned int k = 0; k < Datasets[dset[i]].ncyc[j].size(); k++) {
					x[dc] = 1.0 * dc;
					y[dc] = Datasets[dset[i]].R2_exp[j][k];
					sig[dc] = Datasets[dset[i]].R2_esd[j][k];
					dc++;
				}
			}
		}
	}
	if (dc != ndata) {
		std::cerr << "dc=" << dc << " while ndata=" << ndata << std::endl;
		Abort(1);
	}
	//
	if (minimizer == "lm") {
		conv = false;
		//           x,y,sig,param,paramF,param(atom),object
		Fitmrq LMFIT(x, y, sig, a, ai, aa, an, (*this)); //passing this object on

		LMFIT.SetPrint(print);
		LMFIT.SetMaxIter(maxiter);
		LMFIT.SetTol(tol);
		LMFIT.fit();

		conv = LMFIT.converged();
		ChiSq = LMFIT.chisq;
		std::vector<int> tno;
		std::vector<std::string> tna;
		for (unsigned int i = 0; i < ParamName.size(); i++) {
			if (AtomNumber[i] < 0) { // the parameters is global
				GlobalParamE[ParamName[i]] = sqrt(LMFIT.covar[i][i] * LMFIT.RedChiSq);
				if (GlobalParamF[ParamName[i]]) {
					tno.push_back(i);
					tna.push_back(ParamName[i] + "_global");
				}
			} else {
				LocalParamE[AtomNumber[i]][ParamName[i]] = sqrt(LMFIT.covar[i][i] * LMFIT.RedChiSq);
				if (LocalParamF[AtomNumber[i]][ParamName[i]]) {
					tno.push_back(i);
					tna.push_back(ParamName[i] + "_" + AtomNumber2AtomName(AtomNumber[i]));
				}
			}
		}
		LastFitCovar.resize(tno.size(), tno.size());
		LastFitParamName.clear();
		for (unsigned int i = 0; i < tno.size(); i++) {
			LastFitParamName[tna[i]] = i;
			for (unsigned int j = 0; j < tno.size(); j++) {
				LastFitCovar(i, j) = LMFIT.covar[tno[i]][tno[j]];
			}
		}
	} else {
		std::cerr << " The minimizer: " << minimizer << " is not available\n";
		std::cerr << " Function .Minimize();\n";
		std::cerr << std::endl;
		Abort(1);
	}
	an.clear();
	a.clear();
	aa.clear();
	ai.clear();
	return (double) ChiSq;
}

//Wrapper functions
double Catia::Minimize(std::vector<int> dv, std::vector<std::string> asv, std::string min, bool& conv) {
	std::vector<int> av;
	for (unsigned int i = 0; i < asv.size(); i++) {
		if (AtomName2AtomNumber(asv[i]) == -1) {
			std::cerr << " You are trying to include the atom -->" << asv[i] << "\n";
			std::cerr << " in the minimization, but this atom is not defined\n";
			std::cerr << " Function .Minimize();\n";
			std::cerr << std::endl;
			Abort(1);
		} else {
			av.push_back(AtomName2AtomNumber(asv[i]));
		}
	}
	return Minimize(dv, av, min, conv);
}

double Catia::Minimize(std::string d, std::vector<int> av, std::string min, bool& conv) {
	d = to_lower_copy(d);

	std::vector<int> dv;
	if (d == "all") {
		for (unsigned int i = 0; i < Datasets.size(); i++) {
			dv.push_back(i);
		}
	} else {
		std::cerr << " The Dataset option -->" << d << " is unknown\n";
		std::cerr << " Function Minimize();\n";
		std::cerr << std::endl;
		Abort(1);
	}
	return Minimize(dv, av, min, conv);
}

double Catia::Minimize(std::vector<int> dv, std::string a, std::string min, bool& conv) {
	a = to_lower_copy(a);
	std::vector<int> av;
	if (a == "all" || a == "*") {
		a = ".+";
	}
	const boost::regex regularExpression(a);
	BOOST_FOREACH(pairSI_t atom, Atoms)
{	if ( regex_match(atom.first, regularExpression)) {
		av.push_back(atom.second);
	}
}
return Minimize(dv,av,min,conv);
}

double Catia::Minimize(std::string d, std::string a, std::string min, bool& conv) {
	d = to_lower_copy(d);
	//a=to_lower_copy(a);
	std::vector<int> dv, av;
	if (d == "all") {
		for (unsigned int i = 0; i < Datasets.size(); i++) {
			dv.push_back(i);
		}
	} else {
		std::cerr << " The Dataset option -->" << d << "<- is unknown\n";
		std::cerr << " Function Minimize();\n";
		std::cerr << std::endl;
		Abort(1);
	}
	if (a == "all" || a == "*") {
		a = ".+";
	}
	const boost::regex regularExpression(a);
	BOOST_FOREACH(pairSI_t atom, Atoms)
{	if ( regex_match(atom.first, regularExpression)) {
		av.push_back(atom.second);
	}
}
return Minimize(dv,av,min,conv);
}

