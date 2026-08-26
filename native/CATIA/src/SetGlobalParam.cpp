#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include <Catia.h>
#include <Abort.h>

using boost::to_lower_copy;

typedef std::pair<std::string, double> pairSD_t;

void Catia::SetGlobalParam(std::string paramNameSelection, double paramVal) {

	if (to_lower_copy(paramNameSelection) == "all") {
		paramNameSelection = ".+";
	}

	const boost::regex paramNameRegex(paramNameSelection);

	BOOST_FOREACH(pairSD_t param, GlobalParam) {
		if ( regex_match(param.first, paramNameRegex) ) {
			GlobalParam[param.first] = paramVal;
		}
	}
}
