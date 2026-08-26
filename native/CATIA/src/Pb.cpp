/*
 *  Pb.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 25/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>

double Catia::Pb(const Dataset& dset) const {
	/*
	 Calculates the pb for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	double pb = 0.0;

	if (_rateType == "standard") {
		pb =  PbStandard(dset);
	} else if (_rateType == "arrhenius") {
		pb =  PbArrhenius(dset);
	} else if (_rateType == "kab_kba") {
		pb =  PbKabKba();
	} else {
		std::cerr << " Nothing compiled for rateType = " << _rateType << "\n";
		std::cerr << " Function: Pb();\n";
		Abort(1);
	}

	return fabs(pb);
}

double Catia::Pb_3st(const Dataset& dset) const {
	/*
	 Calculates the pb for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	double pb = 0.0;

	if (_rateType == "standard") {
		pb =  PbStandard(dset);
	} else if (_rateType == "pbarrhenius") {
		pb = PbArrhenius_3st(dset);
	} else {
		std::cerr << " Nothing compiled for rateType = " << _rateType << "\n";
		std::cerr << " Function: Pb();\n";
		Abort(1);
	}

	return pb;
}

double Catia::Pc_3st(const Dataset& dset) const {
	/*
	 Calculates the pb for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */

	double pc = 0.0;

	if (_rateType == "standard") {
		pc = PcStandard(dset);
	} else if (_rateType == "pbarrhenius") {
		pc = PcArrhenius_3st(dset);
	} else {
		std::cerr << " Nothing compiled for rateType = " << _rateType << "\n";
		std::cerr << " Function: Pb();\n";
		Abort(1);
	}

	return pc;
}
