/*
 * SetDeltaOmegaType.cpp
 *
 *  Created on: Feb 17, 2010
 *      Author: guillaume
 */

#include <boost/foreach.hpp>

#include <Dataset.h>
#include <Catia.h>

void Catia::SetDeltaOmegaType(const std::string deltaOmegaTempDep) {

	BOOST_FOREACH(Dataset& dset, Datasets) {
		dset.SetDeltaOmegaTempDep(deltaOmegaTempDep);
	}

}
