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

double Catia::DeltaOmegaLinear(const Dataset& dset, int globalAtomIndex) const {
    double deltaO_a = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_a");
    double deltaO_b = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_b");
    return deltaO_a * dset._temperature + deltaO_b;
}

double Catia::DeltaOmegaLinear_ab(const Dataset& dset, int globalAtomIndex) const {
    double deltaO_ab_a = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_ab_a");
    double deltaO_ab_b = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_ab_b");
    return deltaO_ab_a * dset._temperature + deltaO_ab_b;
}

double Catia::DeltaOmegaLinear_ac(const Dataset& dset, int globalAtomIndex) const {
    double deltaO_ac_a = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_ac_a");
    double deltaO_ac_b = ResolveParam(LocalParam[globalAtomIndex],"DeltaO_ac_b");
    return deltaO_ac_a * dset._temperature + deltaO_ac_b;
}
