#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include <Catia.h>
#include <Abort.h>

using boost::to_lower_copy;

typedef std::pair<std::string, int> pairSI_t;
typedef std::pair<std::string, double> pairSD_t;

void Catia::SetLocalParam(std::string atomNameSelection, std::string paramNameSelection, double paramVal){

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
			BOOST_FOREACH(pairSD_t param, LocalParam[atom.second]) {
				if ( regex_match(param.first, paramNameRegex) ) {
					LocalParam[atom.second][param.first] = paramVal;
				}
			}
		}
	}

}
