/*
 *  PbKabKba.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 24/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>

double Catia::PbKabKba() const {
    double kab = ResolveParam(GlobalParam,"kab");
    double kba = ResolveParam(GlobalParam,"kba");
    return kab / (kab + kba);
}
