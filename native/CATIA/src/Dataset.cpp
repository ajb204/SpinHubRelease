/*
 * Dataset.cpp
 *
 *  Created on: Jun 15, 2009
 *      Author: guillaume
 */

#include <Dataset.h>

Dataset::Dataset()
: _deltaOmegaTempDep("standard")
, _inputFileName("")
, _id("")
, _sfrq(0.0)
, _xcar(0.0)
, _procpar()
, _temperature(0.0)
, _nucleus("")
, _gamma(0.0)
, _localToGlobalAtomIndex()
, _atomNameToLocalAtomNumber()
, _localAtomNumberToAtomName()
, _initialized(false)
, _haveInitIntensity(false)
, _intensityFileFormat()
, _dataDirectory("")
, ncyc()
, R2_exp()
, R2_esd()
, R2_calc()
, initialIntensities() {
	_minError[0] = 0.0;
	_minError[1] = 0.0;
}

Dataset::~Dataset() {}

