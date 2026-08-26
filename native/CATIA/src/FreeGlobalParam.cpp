#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <Abort.h>

using boost::to_lower_copy;

typedef std::pair<std::string, int> pairSI_t;

void Catia::FreeGlobalParam(std::string paramNameSelection, int fix) {

	if (to_lower_copy(paramNameSelection) == "all") {
		paramNameSelection = ".+";
	}

	const boost::regex paramNameRegex(paramNameSelection);

	BOOST_FOREACH(pairSI_t param, GlobalParam) {
		if ( regex_match(param.first, paramNameRegex) ) {
			GlobalParamF[param.first] = fix;
		}
	}

}
