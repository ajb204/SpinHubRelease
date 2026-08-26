/*
 *  KexKabKba.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 24/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>

double Catia::KexKabKba() const {
    return ResolveParam(GlobalParam,"kab") + ResolveParam(GlobalParam,"kba");
}
