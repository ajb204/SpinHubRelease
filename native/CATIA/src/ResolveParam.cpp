//
// June 3 2009:
//     Modified by G. Bouvignies to pass by const reference.
//         It speeds up the new code that makes a lot of call of ResolveParam.
//
//

#include <Catia.h>
#include <Abort.h>

double Catia::ResolveParam(const std::map<std::string, double>& map, const std::string& name) const {

	std::map<std::string, double>::const_iterator it(map.find(name));

	if (it == map.end()) {
		std::cerr << " The parameter " << name << " could not be located\n";
		std::cerr << " Function .ResolveParam();\n";
		Abort(1);
	} else {
		return it->second;
	}

}

bool Catia::BResolveParam(const std::map<std::string, double>& map, const std::string& name) const {

	if (map.find(name) == map.end()) {
		return false;
	} else {
		return true;
	}
}
