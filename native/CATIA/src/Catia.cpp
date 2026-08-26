/*
 CPMG And Trosy Intelligent Analysis (Catia)

 D. Flemming Hansen
 flemming@pound.med.utoronto.ca,
 August 2007


 */

#include <Catia.h>
#include <Dataset.h>

Catia::Catia()
: LastFitCovar()
, GlobalParam()
, GlobalParamF()
, GlobalParamE()
, LocalParam()
, LocalParamF()
, LocalParamE()
, LocalNotes()
, X2dset_atom()
, Datasets()
, Atoms()
, Atoms2Fit()
, CalcDeriv(1)
, _gammaH(2.67522128E8)
, _gammaN(-2.71261804E7)
, _gammaC(6.728284E7)
, _gammaD(4.10662791E7)
, _gammaF(2.518148E8)
, _fixParamLimit(1E-9)
, _stopFitting(false)
, _multipleTemperatures(false)
, _rateType("standard")
{}

Catia::~Catia() {}
