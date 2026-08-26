#include <Catia.h>
/*
  Set the threshold for, when a parameter is considered as a 
  'dummy' parameters.
*/
void Catia::SetInternalFixThres(double thres) {
  _fixParamLimit=thres;
  return;
};
