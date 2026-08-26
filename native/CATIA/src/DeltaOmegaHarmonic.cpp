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

double Catia::DeltaOmegaHarmonic(const Dataset& dset,int globalAtomIndex) const {
    double deltaO_a = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_a");
    double deltaO_b = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_b");
    double deltaO_c = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_c");
    double temperature = dset._temperature;
    return deltaO_a * temperature * temperature + deltaO_b * temperature + deltaO_c;
}
