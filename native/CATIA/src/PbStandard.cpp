/*
 *  PbStandard.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 24/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>

double Catia::PbStandard(const Dataset& dset) const {

	std::string pStr("pb");

	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		pStr += line;
	}

	return ResolveParam(GlobalParam, pStr);

}

double Catia::PcStandard(const Dataset& dset) const {

	std::string pStr("pc");

	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		pStr += line;
	}

	return ResolveParam(GlobalParam, pStr);
}
