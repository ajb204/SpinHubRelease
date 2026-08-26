#include <Catia.h>

std::string Catia::AtomNumber2AtomName(int i) const {
	std::string Name;
	std::map<std::string, int>::const_iterator it;
	for (it = Atoms.begin(); it != Atoms.end(); ++it) {
		if (it->second == i) {
			Name = it->first;
			break;
		}
	}
	return Name;
}
