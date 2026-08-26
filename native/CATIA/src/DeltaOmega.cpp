/*
 *  DeltaOmega.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 24/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>

double Catia::DeltaOmega(const Dataset& dset, int globalAtomIndex) const {
	/*
	 Calculates the delta omega for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	if (dset.DeltaOmegaTempDepIs("standard")) {
		return DeltaOmegaStandard(dset, globalAtomIndex);
	} else if (dset.DeltaOmegaTempDepIs("linear")) {
		return DeltaOmegaLinear(dset, globalAtomIndex);
	} else if (dset.DeltaOmegaTempDepIs("harmonic")) {
		return DeltaOmegaHarmonic(dset, globalAtomIndex);
	} else {
		std::cerr << " No sequence compiled for deltaOmegaTempDep = " << dset.DeltaOmegaTempDep() << "\n";
		std::cerr << " Function: deltaOmega();\n";
		Abort(1);
	}
	return 0.0;
}

double Catia::DeltaOmega_ab(const Dataset& dset, int globalAtomIndex) const {
	/*
	 Calculates the delta omega for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	if (dset.DeltaOmegaTempDepIs("standard")) {
		return DeltaOmegaStandard_ab(dset, globalAtomIndex);
	} else if (dset.DeltaOmegaTempDepIs("linear")) {
		return DeltaOmegaLinear_ab(dset, globalAtomIndex);
	} else {
		std::cerr << " No sequence compiled for deltaOmegaType = " << dset.DeltaOmegaTempDep() << "\n";
		std::cerr << " Function: deltaOmega();\n";
		Abort(1);
	}
	return 0.0;
}

double Catia::DeltaOmega_ac(const Dataset& dset, int globalAtomIndex) const {
	/*
	 Calculates the delta omega for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	if (dset.DeltaOmegaTempDepIs("standard")) {
		return DeltaOmegaStandard_ac(dset, globalAtomIndex);
	} else if (dset.DeltaOmegaTempDepIs("linear")) {
		return DeltaOmegaLinear_ac(dset, globalAtomIndex);
	} else {
		std::cerr << " No sequence compiled for deltaOmegaType = " << dset.DeltaOmegaTempDep() << "\n";
		std::cerr << " Function: deltaOmega();\n";
		Abort(1);
	}
	return 0.0;
}
