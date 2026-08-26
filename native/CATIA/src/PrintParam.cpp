#include <standard.h>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>

using std::endl;

using boost::format;
using boost::io::group;
using boost::to_lower_copy;
using boost::trim;

std::string Catia::PrintParam(std::string opt) {
	std::ostringstream OutStream;
	char line[MAX_STRING_LENGTH];
	//
	trim(opt);
	// containers
	std::vector<std::string> ParamName;
	std::vector<double> ParamVal;
	std::vector<int> AtomNumber;
	std::map<std::string, double>::iterator itDS; //iterator:double<-string
	std::map<std::string, int>::iterator itIS; //iterator:double<-string
	std::map<std::string, std::string>::iterator itSS;
	//
	// is the string a global or local parameter?
	if (BResolveParam(GlobalParam, opt)) {
		OutStream << "#";
		OutStream << PrintParam(true, -1, opt);
		return OutStream.str();
	}
	//
	// is the parameter a local parameter ?
	bool optIsLocal = false;
	bool optIsAtomName = false;
	for (itIS = Atoms.begin(); itIS != Atoms.end(); ++itIS) {
		if ((itIS->second) < LocalParam.size()) {
			if (BResolveParam(LocalParam[itIS->second], opt)) {
				sprintf(line, " %9s", (itIS->first).c_str());
				OutStream << line;
				ClearBuf(line, sizeof(line));
				OutStream << PrintParam(false, itIS->second, opt);
				optIsLocal = true;
			}
		}
	}
	//
	// is the option equal to an atomname ?
	if (AtomName2AtomNumber(opt) != -1) {
		if (LocalParam.size() > AtomName2AtomNumber(opt)) {
			optIsAtomName = true;
			OutStream << "#";
			OutStream << "# Atom: " << opt << "\n";
			for (itSS = LocalNotes[AtomName2AtomNumber(opt)].begin(); itSS != LocalNotes[AtomName2AtomNumber(opt)].end(); ++itSS) {
				OutStream << "# Notes regarding: " << itSS->first << ": " << itSS->second << std::endl;
			}
			for (unsigned int i = 0; i < 39; i++)
				OutStream << "#";
			OutStream << "\n";
			for (itDS = LocalParam[AtomName2AtomNumber(opt)].begin(); itDS != LocalParam[AtomName2AtomNumber(opt)].end(); ++itDS) {
				OutStream << "#" << PrintParam(false, AtomName2AtomNumber(opt), itDS->first);
			}
		}
	}
	if (optIsLocal || optIsAtomName)
		return OutStream.str();
	//
	opt = to_lower_copy(opt);
	if (opt.substr(0, 3) == "all" || opt.substr(0, 1) == "*") { //print all parameters
		//
		// first global parameters:
		OutStream << "# Global Parameters:\n";
		for (itDS = GlobalParam.begin(); itDS != GlobalParam.end(); ++itDS) {
			OutStream << "#";
			OutStream << PrintParam(true, -1, itDS->first);
		}
		//
		//Then local Parameters
		std::map<std::string, int>::iterator itIS; //iterator:int<-string
		for (itIS = Atoms.begin(); itIS != Atoms.end(); ++itIS) {
			OutStream << "#\n";
			OutStream << "# Atom: " << itIS->first << std::endl;
			for (unsigned int i = 0; i < 39; i++)
				OutStream << "#";
			OutStream << "\n";
			//
			for (itDS = LocalParam[itIS->second].begin(); itDS != LocalParam[itIS->second].end(); ++itDS) {
				OutStream << "#";
				OutStream << PrintParam(false, itIS->second, itDS->first);
			}
		}
	} else if (opt.substr(0, 6) == "global") {
		for (itDS = GlobalParam.begin(); itDS != GlobalParam.end(); ++itDS) {
			OutStream << PrintParam(true, -1, itDS->first);
		}
	} else if (opt.substr(0, 9) == "arrhenius") {
		OutStream << PrintArrheniusParam();
	}
	/*
	 else {
	 std::cerr<<" The option "<<opt<<" is not available\n";
	 std::cerr<<" and no Global nor Local parameters exists with that name\n";
	 std::cerr<<" Function .PrintParam();\n";
	 std::cerr<<std::endl;
	 Abort(1);
	 }
	 */
	return OutStream.str();
}

std::string Catia::PrintParam(bool global, int GlobalAtom, std::string ParamName) {
	std::ostringstream OutStream;
	char line[MAX_STRING_LENGTH];
	if (global) {
		sprintf(line, "%11s %13.6e", (ParamName).c_str(), GlobalParam[ParamName]);
		OutStream << line;
		ClearBuf(line, sizeof(line));
		if (fabs(GlobalParamE[ParamName] + 1) < 1E-6) {
			sprintf(line, "%13s", "NotFitted");
			OutStream << line;
			ClearBuf(line, sizeof(line));
		} else if (fabs(GlobalParamF[ParamName] + 2) < 1E-6) {
			sprintf(line, "%13s", "InternalFix");
			OutStream << line;
			ClearBuf(line, sizeof(line));
		} else {
			if (GlobalParamF[ParamName]) {
				sprintf(line, "%13.6e", GlobalParamE[ParamName]);
				OutStream << line;
				ClearBuf(line, sizeof(line));
			} else {
				sprintf(line, "%13s", "fixed");
				OutStream << line;
				ClearBuf(line, sizeof(line));
			}
		}
		OutStream << std::endl;
	} else {
		sprintf(line, "%11s %13.6e", (ParamName).c_str(), LocalParam[GlobalAtom][ParamName]);
		OutStream << line;
		ClearBuf(line, sizeof(line));
		if (fabs(LocalParamE[GlobalAtom][ParamName] + 1) < 1E-6) {
			sprintf(line, "%13s", "NotFitted");
			OutStream << line;
			ClearBuf(line, sizeof(line));
		} else if (fabs(LocalParamF[GlobalAtom][ParamName] + 2) < 1E-6) {
			sprintf(line, "%13s", "InternalFix");
			OutStream << line;
		} else {
			if (LocalParamF[GlobalAtom][ParamName]) {
				sprintf(line, "%13.6e", LocalParamE[GlobalAtom][ParamName]);
				OutStream << line;
				ClearBuf(line, sizeof(line));
			} else {
				sprintf(line, "%13s", "fixed");
				OutStream << line;
				ClearBuf(line, sizeof(line));
			}
		}
		OutStream << std::endl;
	}
	return OutStream.str();
}

std::string Catia::PrintArrheniusParam(void) const {

	bool hasDeltaHDeltaS2st = (BResolveParam(GlobalParam, "deltaHa")
			&& BResolveParam(GlobalParam, "deltaHb")
			&& BResolveParam(GlobalParam, "deltaSa")
			&& BResolveParam(GlobalParam, "deltaSb"));

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


	if(hasDeltaHDeltaS2st){
		return PrintArrheniusParam2st();
	} else if (hasDeltaHDeltaS3st){
		return PrintArrheniusParam3st();
	}

}

std::string Catia::PrintArrheniusParam3st(void) const {

	const double deltaHb  = ResolveParam(GlobalParam, "deltaHb");
	const double deltaHc  = ResolveParam(GlobalParam, "deltaHc");
	const double deltaHab = ResolveParam(GlobalParam, "deltaHab");
	const double deltaHac = ResolveParam(GlobalParam, "deltaHac");
	const double deltaHbc = ResolveParam(GlobalParam, "deltaHbc");
	const double deltaSb  = ResolveParam(GlobalParam, "deltaSb");
	const double deltaSc  = ResolveParam(GlobalParam, "deltaSc");
	const double deltaSab = ResolveParam(GlobalParam, "deltaSab");
	const double deltaSac = ResolveParam(GlobalParam, "deltaSac");
	const double deltaSbc = ResolveParam(GlobalParam, "deltaSbc");

	const double jouleToCal = 0.239005736;
	const double jouleToKCal = jouleToCal * 1e-3;

	std::ostringstream OutStream;

	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaH(B)" % (deltaHb * jouleToKCal) % "[kcal/mol]";
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaS(B)" % (deltaSb * jouleToCal)  % "[cal/mol/K]" << endl;
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaH(C)" % (deltaHc * jouleToKCal) % "[kcal/mol]";
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaS(C)" % (deltaSc * jouleToCal) % "[cal/mol/K]" << endl;
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaH(A<->B)" % (deltaHab * jouleToKCal) % "[kcal/mol]";
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaS(A<->B)" % (deltaSab * jouleToCal) % "[cal/mol/K]" << endl;
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaH(A<->C)" % (deltaHac * jouleToKCal) % "[kcal/mol]";
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaS(A<->C)" % (deltaSac * jouleToCal) % "[cal/mol/K]" << endl;
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaH(B<->C)" % (deltaHbc * jouleToKCal) % "[kcal/mol]";
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaS(B<->C)" % (deltaSbc * jouleToCal) % "[cal/mol/K]" << endl;

	std::vector<double> temperatureMemory;
	BOOST_FOREACH(const Dataset& dset, Datasets) {

		bool tempInTempMemory = false;
		BOOST_FOREACH(const double temperature, temperatureMemory) {
			if (temperature == dset._temperature) {
				tempInTempMemory = true;
			}
		}
		if(tempInTempMemory) {
			continue;
		} else {
			temperatureMemory.push_back(dset._temperature);
		}

		OutStream << "#######################################################" << endl;
		OutStream << format("# %-7s = %10.4f %-5s") % "Temp" % dset._temperature % "[C]";
		OutStream << format("# %-7s = %10.4f %-5s") % "kex(AB)" % Kex_ab_3st(dset) % "[/s]" << endl;
		OutStream << format("# %-7s = %10.4f %-5s") % "p(B)" % (Pb_3st(dset) * 100) % "[%]";
		OutStream << format("# %-7s = %10.4f %-5s") % "kex(AC)" % Kex_ac_3st(dset) % "[/s]" << endl;
		OutStream << format("# %-7s = %10.4f %-5s") % "p(C)" % (Pc_3st(dset) * 100) % "[%]";
		OutStream << format("# %-7s = %10.4f %-5s") % "kex(BC)" % Kex_bc_3st(dset) % "[/s]" << endl;
	}
	OutStream << "#######################################################" << endl;

	OutStream << format("%-10s = %+.6e") % "deltaHb" % deltaHb << endl;
	OutStream << format("%-10s = %+.6e") % "deltaSb" % deltaSb << endl;
	OutStream << format("%-10s = %+.6e") % "deltaHc" % deltaHc << endl;
	OutStream << format("%-10s = %+.6e") % "deltaSc" % deltaSc << endl;
	OutStream << format("%-10s = %+.6e") % "deltaHab" % deltaHab << endl;
	OutStream << format("%-10s = %+.6e") % "deltaSab" % deltaSab << endl;
	OutStream << format("%-10s = %+.6e") % "deltaHac" % deltaHac << endl;
	OutStream << format("%-10s = %+.6e") % "deltaSac" % deltaSac << endl;
	OutStream << format("%-10s = %+.6e") % "deltaHbc" % deltaHbc << endl;
	OutStream << format("%-10s = %+.6e") % "deltaSbc" % deltaSbc << endl;

	return OutStream.str();
}

std::string Catia::PrintArrheniusParam2st(void) const {

	const double deltaHa  = ResolveParam(GlobalParam, "deltaHa");
	const double deltaSa  = ResolveParam(GlobalParam, "deltaSa");
	const double deltaHb  = ResolveParam(GlobalParam, "deltaHb");
	const double deltaSb  = ResolveParam(GlobalParam, "deltaSb");

	const double jouleToCal = 0.239005736;
	const double jouleToKCal = jouleToCal * 1e-3;

	std::ostringstream OutStream;

	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaH(A->D)" % (deltaHa * jouleToKCal) % "[kcal/mol]";
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaS(A->D)" % (deltaSa * jouleToCal)  % "[cal/mol/K]" << endl;
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaH(B->D)" % (deltaHb * jouleToKCal) % "[kcal/mol]";
	OutStream << format("# %-15s = %10.4f %-15s") % "DeltaS(B->D)" % (deltaSb * jouleToCal)  % "[cal/mol/K]" << endl;

	std::vector<double> temperatureMemory;
	BOOST_FOREACH(const Dataset& dset, Datasets) {

		bool tempInTempMemory = false;
		BOOST_FOREACH(const double temperature, temperatureMemory) {
			if (temperature == dset._temperature) {
				tempInTempMemory = true;
			}
		}
		if(tempInTempMemory) {
			continue;
		} else {
			temperatureMemory.push_back(dset._temperature);
		}

		OutStream << "#######################################################" << endl;
		OutStream << format("# %-7s = %10.4f %-5s") % "Temp" % dset._temperature % "[C]" << endl;
		OutStream << format("# %-7s = %10.4f %-5s") % "p(B)" % (Pb(dset) * 100) % "[%]";
		OutStream << format("# %-7s = %10.4f %-5s") % "kex(AB)" % Kex(dset) % "[/s]" << endl;
		OutStream << format("# %-7s = %10.4f %-5s") % "k(A->B)" % (Pb(dset) * Kex(dset)) % "[/s]";
		OutStream << format("# %-7s = %10.4f %-5s") % "k(B->A)" % ((1.0 - Pb(dset)) * Kex(dset)) % "[/s]" << endl;
	}
	OutStream << "#######################################################" << endl;

	OutStream << format("%-10s = %+.6e") % "deltaHa" % deltaHa << endl;
	OutStream << format("%-10s = %+.6e") % "deltaSa" % deltaSa << endl;
	OutStream << format("%-10s = %+.6e") % "deltaHb" % deltaHb << endl;
	OutStream << format("%-10s = %+.6e") % "deltaSb" % deltaSb << endl;

	return OutStream.str();
}
