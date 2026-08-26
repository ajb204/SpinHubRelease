/*
 *  Kex.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 25/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>

double Catia::Kex(const Dataset& dset) const {
	/*
	 Calculates the kex for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */
  double kexx=0.0;
	if (_rateType == "standard") {
		kexx= KexStandard(dset);
	} else if (_rateType == "arrhenius") {
		kexx= KexArrhenius(dset);
	} else if (_rateType == "kab_kba") {
		kexx= KexKabKba();
	} else {
		std::cerr << " Nothing compiled for rateType = " << _rateType << "\n";
		std::cerr << " Function: Kex();\n";
		Abort(1);
	}
       return fabs(kexx);
}

double Catia::Kex_ab_3st(const Dataset& dset) const {
	/*
	 Calculates the kex for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	if (_rateType == "standard") {
		return KexStandard_3st_ab(dset);
	} else if (_rateType == "arrhenius") {
		return KexArrhenius_3st_ab(dset);
	} else {
		std::cerr << " Nothing compiled for rateType = " << _rateType << "\n";
		std::cerr << " Function: Kex();\n";
		Abort(1);
	}
	return 0.0;
}

double Catia::Kex_ac_3st(const Dataset& dset) const {
	/*
	 Calculates the kex for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	if (_rateType == "standard") {
		return KexStandard_3st_ac(dset);
	} else if (_rateType == "arrhenius") {
		return KexArrhenius_3st_ac(dset);
	} else {
		std::cerr << " Nothing compiled for rateType = " << _rateType << "\n";
		std::cerr << " Function: Kex();\n";
		Abort(1);
	}
	return 0.0;
}

double Catia::Kex_bc_3st(const Dataset& dset) const {
	/*
	 Calculates the kex for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	if (_rateType == "standard") {
		return KexStandard_3st_bc(dset);
	} else if (_rateType == "arrhenius") {
		return KexArrhenius_3st_bc(dset);
	} else {
		std::cerr << " Nothing compiled for rateType = " << _rateType << "\n";
		std::cerr << " Function: Kex();\n";
		Abort(1);
	}
	return 0.0;
}

