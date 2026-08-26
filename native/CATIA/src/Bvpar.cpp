#include <Dataset.h>

/*
 Thus bool function returns 'true' of the parameter par exists in
 the dataset dset

 Flemming, Aug. 07 2008
 */

bool Dataset::Bvpar(const std::string& par) {
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
		return false;
	} else {
		return true;
	}
}
