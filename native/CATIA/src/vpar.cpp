#include <Dataset.h>
#include <Abort.h>
/*
 This function returns a value of a specific parameter of the Dataset
 specified.
 */
std::vector<std::string> Dataset::vpar(const std::string par) const {
	std::vector<std::string> out;
	bool happyface = false;
	for (unsigned int i = 0; i < _procpar.size(); i++) {
		out = _procpar[i];
		if (out[0] == par) {
			happyface = true;
			break;
		}
	}
	if (!(happyface)) {
		std::cerr << " You requested the value of the parameter: " << par
				<< "\n";
		std::cerr
				<< " but this parameters is not present in the Dataset with ID\n";
		std::cerr << " " << _id << std::endl;
		std::cerr << " Function vpar();\n";
		std::cerr << std::endl;
		Abort(1);
	}
	return out;
}
