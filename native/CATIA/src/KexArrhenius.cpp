/*
 *  KexStandard.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 25/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>
#include <Dataset.h>

#define PLANCKS_H (6.62606896e-34) /* kg m^2 / s */
#define BOLTZMANN_KB (1.3806504e-23) /* kg m^2 / K s^2 */
#define MOLAR_GAS_R (8.314472e0) /* kg m^2 / K mol s^2 */

double Catia::KexArrhenius(const Dataset& dset) const {
	const double deltaSb = ResolveParam(GlobalParam, "deltaSb");
	const double deltaHb = ResolveParam(GlobalParam, "deltaHb");
	const double deltaSab = ResolveParam(GlobalParam, "deltaSab");
	const double deltaHab = ResolveParam(GlobalParam, "deltaHab");
	const double T = dset._temperature + 273.15;

	const double A = 3000. * T; // 1.6e-7 * BOLTZMANN_KB * T / PLANCKS_H;;
	const double RT = MOLAR_GAS_R * T;

	const double kab = A * exp(-(deltaHab - T * deltaSab) / RT);
	const double kba = A * exp(-((deltaHab - deltaHb) - T * (deltaSab - deltaSb)) / RT);

	const double kex = kab + kba;

	return kex;
}

double Catia::KexArrhenius_3st_ab(const Dataset& dset) const {
	const double deltaSb = ResolveParam(GlobalParam, "deltaSb");
	const double deltaSab = ResolveParam(GlobalParam, "deltaSab");

	const double deltaHb = ResolveParam(GlobalParam, "deltaHb");
	const double deltaHab = ResolveParam(GlobalParam, "deltaHab");

	const double T = dset._temperature + 273.15;

	const double A = 3000. * T; // 1.6e-7 * BOLTZMANN_KB * T / PLANCKS_H;;
	const double RT = MOLAR_GAS_R * T;

	const double kab = A * exp(-(deltaHab - T * deltaSab) / RT);
	const double kba = A * exp(-((deltaHab - deltaHb) - T * (deltaSab - deltaSb)) / RT);

	const double kex_ab = kab + kba;

	return kex_ab;
}

double Catia::KexArrhenius_3st_ac(const Dataset& dset) const {
	const double deltaSc = ResolveParam(GlobalParam, "deltaSc");
	const double deltaSac = ResolveParam(GlobalParam, "deltaSac");

	const double deltaHc = ResolveParam(GlobalParam, "deltaHc");
	const double deltaHac = ResolveParam(GlobalParam, "deltaHac");

	const double T = dset._temperature + 273.15;

	const double A = 3000. * T; // 1.6e-7 * BOLTZMANN_KB * T / PLANCKS_H;;
	const double RT = MOLAR_GAS_R * T;

	const double kac = A * exp(-(deltaHac - T * deltaSac) / RT);
	const double kca = A * exp(-((deltaHac - deltaHc) - T * (deltaSac - deltaSc)) / RT);

	const double kex_ac = kac + kca;

	return kex_ac;
}

double Catia::KexArrhenius_3st_bc(const Dataset& dset) const {
	const double deltaSb = ResolveParam(GlobalParam, "deltaSb");
	const double deltaSc = ResolveParam(GlobalParam, "deltaSc");
	const double deltaSbc = ResolveParam(GlobalParam, "deltaSbc");

	const double deltaHb = ResolveParam(GlobalParam, "deltaHb");
	const double deltaHc = ResolveParam(GlobalParam, "deltaHc");
	const double deltaHbc = ResolveParam(GlobalParam, "deltaHbc");

	const double T = dset._temperature + 273.15;

	const double A = 3000. * T; // 1.6e-7 * BOLTZMANN_KB * T / PLANCKS_H;;
	const double RT = MOLAR_GAS_R * T;

	const double kbc = A * exp(-((deltaHbc - deltaHb) - T * (deltaSbc - deltaSb)) / RT);
	const double kcb = A * exp(-((deltaHbc - deltaHc) - T * (deltaSbc - deltaSc)) / RT);

	const double kex_bc = kbc + kcb;

	return kex_bc;
}
