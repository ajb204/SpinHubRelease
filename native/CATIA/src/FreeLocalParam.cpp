/*
 D.F. Hansen,
 September 12 2007

 May 27 2009
 Modified by G. Bouvignies to manage regular expression for atom selection.

 April 16 2010
 Modified by G. Bouvignies to manage regular expression for parameter selection.
 */

#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <Catia.h>

using boost::to_lower_copy;

typedef std::pair<std::string, int> pairSI_t;
typedef std::pair<std::string, double> pairSD_t;

void Catia::FreeLocalParam(unsigned int atom, std::string name, int fix) {
  	if (atom < LocalParam.size()) {
		if (BResolveParam(LocalParam[atom], name)) {
			LocalParamF[atom][name] = fix;
		} else {
			std::cerr << " WARNING:\n";
			std::cerr << " The local parameter " << name << " of atom " << AtomNumber2AtomName(atom) << " is not defined\n";
			std::cerr << " Function .FreeLocalParam()\n";
			std::cerr << std::endl;
		}
		}
}


void Catia::FreeLocalParam(std::vector<int> atoms, std::string name, int fix) {
  	for (unsigned int i = 0; i < atoms.size(); i++) {
		FreeLocalParam(atoms[i], name, fix);
		}
}

void Catia::FreeLocalParam(std::string atomNameSelection, std::string paramNameSelection, int fix) {
  
	if (to_lower_copy(atomNameSelection) == "all") {
		atomNameSelection = ".+";
	}
	if (to_lower_copy(paramNameSelection) == "all") {
		paramNameSelection = ".+";
	}
	
	const boost::regex atomNameRegex(atomNameSelection);
	const boost::regex paramNameRegex(paramNameSelection);
	
	BOOST_FOREACH(pairSI_t atom, Atoms) {
	    if ( regex_match(atom.first, atomNameRegex) ) {
	    BOOST_FOREACH(pairSI_t param, LocalParam[atom.second]) {
	      if ( regex_match(param.first, paramNameRegex) ) {
			FreeLocalParam(atom.second, param.first,fix);
			      }
	    }
	  }
	}
	

}

void Catia::FreeLocalParam(std::vector<std::string> atomnames, std::string param, int fix) {
  	for (unsigned int i = 0; i < atomnames.size(); i++) {
		FreeLocalParam(atomnames[i], param, fix);
	}
}


