/*
 *  KexStandard.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 24/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>

double Catia::KexStandard(const Dataset& dset) const {

	std::string kexStr("kex");

	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		kexStr += line;
	}

	return ResolveParam(GlobalParam, kexStr);

}

double Catia::KexStandard_3st_ab(const Dataset& dset) const {

	std::string kexStr("kex_ab");

	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		kexStr += line;
	}

	return fabs(ResolveParam(GlobalParam, kexStr));

}

double Catia::KexStandard_3st_ac(const Dataset& dset) const {

	std::string kexStr("kex_ac");

	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		kexStr += line;
	}

	return fabs(ResolveParam(GlobalParam, kexStr));

}

double Catia::KexStandard_3st_bc(const Dataset& dset) const {

	std::string kexStr("kex_bc");

	if (_multipleTemperatures) {
		char line[MAX_STRING_LENGTH];
		sprintf(line, "_%.1f", dset._temperature);
		kexStr += line;
	}

	return fabs(ResolveParam(GlobalParam, kexStr));

}
