#include <Catia.h>
#include <Dataset.h>
/*
 Return the Enorm = Sqrt(X2) for a given nucleus in a given Dataset
 */

double Catia::Enorm(Dataset& dset, int GlobalAtomNo) {
	double enorm = 0.;
	CalcR2(dset, GlobalAtomNo);
	int LocalAtomNo = dset._atomNameToLocalAtomNumber[AtomNumber2AtomName(GlobalAtomNo)];
	for (unsigned int i = 0; i < dset.R2_calc[LocalAtomNo].size(); i++) {
		enorm += pow((dset.R2_calc[LocalAtomNo][i] - dset.R2_exp[LocalAtomNo][i]) / dset.R2_esd[LocalAtomNo][i], 2);
	}
	return pow(enorm, 0.5);
}
