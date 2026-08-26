#include <Catia.h>
/*
 This function returns a value of a specific parameter of the dataset
 specified.
 */

using namespace boost::numeric::ublas;

ublas::matrix<double> Catia::CorrelationM(std::vector<std::string> n) {
	ublas::matrix<double> out(n.size(), n.size());
	out.clear();
	double cv = 0.;
	int in, jn;
	// dumb the matrik
	for (unsigned int i = 0; i < n.size(); i++) {
		for (unsigned int j = 0; j < n.size(); j++) {
			in = LastFitParamName[n[i]];
			jn = LastFitParamName[n[j]];
			cv = LastFitCovar(in, jn) / sqrt(LastFitCovar(in, in) * LastFitCovar(jn, jn));
			out(i, j) = cv;
		}
	}
	return out;
}
