/*
 *  DeltaOmegaStandard.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 24/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>

double Catia::DeltaOmegaStandard(const Dataset& dset, int globalAtomIndex) const {
	//
	std::string deltaOStr("DeltaO");
	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		deltaOStr += line;
	}
	//
	return ResolveParam(LocalParam[globalAtomIndex], deltaOStr);
}

double Catia::DeltaOmegaStandard_ab(const Dataset& dset, int globalAtomIndex) const {
	//
	std::string deltaOStr("DeltaO_ab");
	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		deltaOStr += line;
	}
	return ResolveParam(LocalParam[globalAtomIndex], deltaOStr);
}

double Catia::DeltaOmegaStandard_ac(const Dataset& dset, int globalAtomIndex) const {
	//
	std::string deltaOStr("DeltaO_ac");
	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		deltaOStr += line;
	}
	//
	return ResolveParam(LocalParam[globalAtomIndex], deltaOStr);
}
